import sys,os
import copy
from gunfolds.tools import graphkit as gk
import warnings
from scipy import linalg, optimize
import numpy as np
import scipy
import pandas as pd
from libpgm.pgmlearner import PGMLearner
from libpgm.graphskeleton import GraphSkeleton
from libpgm.dyndiscbayesiannetwork import DynDiscBayesianNetwork
from libpgm.sampleaggregator import SampleAggregator
from gunfolds.tools import conversions as conv


def DiscretizeDataQuantiles(data,numvalues):
    #Discretize data my quantiles
    #numvalues is list of integers representing number of possible values for each variable
    data = np.asarray(np.r_[data[:, :-1], data[:, 1:]])
    numvalues=numvalues*2
    for i in range(len(data)):
        data[i]=pd.qcut(data[i],numvalues[i],labels = False)
    return data

def DiscretizeDataRanges(data,numvalues):
    data = np.asarray(np.r_[data[:, :-1], data[:, 1:]])
    numvalues=numvalues*2
    for i in range(len(data)):
        mx=max(data[i])
        mn=min(data[i])
        rng=float(mx-mn)/(numvalues[i])
        bins=np.array([mn+rng*j for j in range(1,numvalues[i]+1)])
        bins[-1]+=1
        data[i]=np.digitize(data[i],bins)
    return data

def MakeDataDictForBNPackage(data):
    #Formats data into Data Dictionary for BN Package
    ddata=[]
    for i in range(len(data[0])):
        a={}
        for j in range(len(data)):
            a[str(j+1)]=str(data[j][i])
        ddata.append(a)
    return ddata

def getBNparams(graph,ddata,n):
    # Gets Disc. BN parameters given a graph skeleton
    #skeleton should include t-1 and t nodes for each variable
    nodes=range(1,(n*2)+1)
    nodes=map(str,nodes)
    edges=gk.edgelist(graph)
    for i in range(len(edges)):
        edges[i]=list([edges[i][0],str(n+int(edges[i][1]))])
    skel=GraphSkeleton()
    skel.V=nodes
    skel.E=edges
    learner=PGMLearner()
    result=learner.discrete_mle_estimateparams(skel,ddata)
    return result
    
def GenProbTable(network):
    #Generates Generic probability table for a Bayesian Network.
    #Used for initialization of Dyn. Disc. BN
    agg= SampleAggregator()
    r=agg.aggregate(network.randomsample(5000))
    return r

def CreateDynDiscBN(vertices,edges,initialVdata,bnV_data):
    #Creates adn returns Dynamic Discrete BN
    d=DynDiscBayesianNetwork()
    d.V=vertices
    d.E=edges
    d.initial_Vdata=initialVdata
    d.twotbn_Vdata=bnV_data
    return d
    
def sampleBN(d,ss):
    #Samples BN and returns data in our generic format
    dictdata=d.randomsample(ss)
    data=[]
    count=0
    for i in dictdata:
        if count==0:
            for j in range(1,len(i)+1):
                data.append([i[str(j)]])
            count+=1
        else:
            for j in range(1,len(i)+1):
                data[j-1].append(i[str(j)])
    return np.array(data)
    
def alterinputsforBNtoDynBN(network,r,n,numvalues):
    #Takes in Bayesian network which includes both t-1 and t nodes
    #Outputs data needed for creating Dynamic Discrete Bayesian Network
    initialVdata=copy.deepcopy(network.Vdata)
    bnV_data=copy.deepcopy(network.Vdata)
    for i in range(1,n+1):
        del initialVdata[str(i+n)]
        del bnV_data[str(i+n)]
        initialVdata[str(i)]['cprob']=[r[str(i)][str(float(j))] for j in range(numvalues[i-1])]
        initialVdata[str(i)]['parents']=[]
        bnV_data[str(i)]['parents']=['past_'+str(int(j)-n) for j in bnV_data[str(i)]['parents']]
    return initialVdata,bnV_data

def makediscrete(graph,data,numvalues,ss):
    n=len(data)
    data=DiscretizeDataQuantiles(data,numvalues)
    ddata=MakeDataDictForBNPackage(data)
    result=getBNparams(graph,ddata,n)
    r=GenProbTable(result)
    initialVdata,bnV_data=alterinputsforBNtoDynBN(result,r,n,numvalues)
    trueEdges=gk.edgelist(graph)
    vertices=map(str,range(1,n+1))
    d=CreateDynDiscBN(vertices,trueEdges,initialVdata,bnV_data)
    data=sampleBN(d,ss)
    return data




def INITdata(graph,numvalues):
    initVdata={}
    V=[]
    E=[]
    for i in range(1,len(graph)+1):
        V.append(str(i))
        for j in graph[str(i)].keys():
            E.append([str(i),str(j)])
        initVdata[str(i)]={'children':None,'cprob':[.5,.5],
        'numoutcomes':numvalues[i-1],'ord':i-1,'parents':None,'vals':[str(q) for q in range(numvalues[i-1])]}
    return V,E,initVdata

import itertools
import random

def normalize(l):
    print l
    sm=sum(l)
    if sm==0:
        l=[1,1,1]
    sm=sum(l)
    for i in range(len(l)):
        l[i]=l[i]/sm
    return l

def gettfVdata(num_graph,numvalues,A):
    graph=conv.ian2g(num_graph)
    num_gt=gk.gtranspose(num_graph)
    gt=conv.ian2g(num_gt)
    tfVdata={}
    Cond={}
    for i in range(1,len(gt)+1):
        for j in gt[str(i)]:
            list1=range(numvalues[i-1])
            list2=range(numvalues[int(j)-1])
            maps=[zip(x,list2) for x in itertools.permutations(list1,len(list2))]
            mapp=random.choice(maps)
            Cond[str(j)]=dict.fromkeys(list(itertools.product(list1,list2)))
            for k in Cond[str(j)]:
                if k in mapp:
                    Cond[str(j)][k]=A
                else:
                    Cond[str(j)][k]=(1-A)/(numvalues[i-1]-1)
        lists=[]
        for j in gt[str(i)]:
            lists.append(range(numvalues[int(j)-1]))
        possible=list(itertools.product(*lists))
        conddict=dict.fromkeys([k for k in possible],[1 for v in range(numvalues[int(i)-1])])
        #conddictFinal=dict.fromkeys([str(i) for i in possible],1)
        #j is parent
        #k is (0,1) where 0 is parent value and 1 is node value
        key=sorted(gt[str(i)].keys())
        for j in range(len(key)):
            for k in Cond[key[j]]:
                for l in conddict:
                    if l[j]==k[0]:
                        a=copy.copy(conddict[l])
                        a[k[1]]=a[k[1]]*Cond[key[j]][k]
                        conddict[l]=a
        print(conddict)
        conddictfinal={}
        for j in conddict:
            a=conddict[j]
            q=[str(k) for k in j]
            conddictfinal[str(q)]=normalize(a)
        print("THIS IS CONDDICTFINAL:")
        print(conddictfinal)
        tfVdata[str(i)]={'children':graph[str(i)].keys(), 'parents':['past_'+p for p in sorted(gt[str(i)].keys())],
        'cprob':conddictfinal,'numoutcomes':numvalues[i-1],'ord':i-1,'vals':[str(v) for v in range(numvalues[i-1])]}
    return tfVdata       
        

def ConstructDynBN(num_graph,numvalues,A,ss):
    graph=conv.ian2g(num_graph)
    print(graph)
    V,E,initVdata=INITdata(graph,numvalues)
    tfVdata=gettfVdata(num_graph,numvalues,A)
    d=DynDiscBayesianNetwork()
    skel=GraphSkeleton()
    skel.V=V
    skel.E=E
    d.V=skel.V
    d.E=skel.E
    d.initial_Vdata=initVdata
    d.twotbn_Vdata=tfVdata
    print(d.V)
    print(d.E)
    print(d.initial_Vdata)
    print(d.twotbn_Vdata)
    data=sampleBN(d,ss)
    return data