# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 14:39:13 2016

@author: johncook
"""
from rpy2 import rinterface
from rpy2.robjects import r
import numpy as np
import itertools
from gunfolds.tools import graphkit as gk
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects

def initemptygraph(n):
    g={}
    for i in range(n):
        g[str(i+1)]={}
    return g

def no_string_superclique(n):
    g = {}
    for i in range(n):
        g[i + 1] = {j + 1: set([(0, 1), (2, 0)])
                         for j in range(n) if j != i}
        g[i + 1][i + 1] = set([(0, 1)])
    return g

def emptyGraph(n):
    '''purpose: empty graph formation for ges algorithm
       input: number of nodes
       output: the graph
       note: 0 = no edge between (symmetric), 
             1 = directed edge from 1st to 2nd index (not symmetric), 
            -1 = undirected edge between (symmetric)
    '''
    g = no_string_superclique(n)
    for i in g:
        for j in g[i]:
            g[i][j] = 0
    return g

def IntNoFlipsbedgelist(g,shift): # bidirected edge list with no flips
    l = []
    for n in g:
        l.extend([(int(n)+shift,int(e)+shift) for e in g[n] if (2,0) in g[n][e] if int(n)<int(e)])
    return l
    
def Intedgelist(g,shift): # directed
    '''
    return a list of tuples for edges of g
    '''
    l = []
    for n in g:
        print n
        l.extend([(int(n),int(e)+shift) for e in g[n] if (0,1) in g[n][e]])
    return l

def DiscreteMLL(data):
    #initialize
    r("""setwd('~/local/R_libs')""")
    basepy=importr('base')
    print r("""getwd()""")
    r("""library('rje',lib.loc='~/local/R_libs')""")
    r("""library('prodlim',lib.loc='~/local/R_libs')""")
    r("""library('ADMGs',lib.loc='~/local/R_libs')""")
    r("""library('rlist',lib.loc='~/local/R_libs')""")
    ADMG = importr("ADMGs",lib_loc="~/local/R_libs")
    rlist=importr("rlist", lib_loc="~/local/R_libs")
    r('dyn.load("~/MyEnv/foo.so")')
    r("""source('~/MyEnv/R-package.R')""")
    r_fitMEG=robjects.globalenv['fitMEG']
    def makeRlist(pylist):
        listR=basepy.list()
        for i,j in pylist:
            listR=rlist.list_append(listR,basepy.c(i,j))
        return listR
            
    def runScore(edgelist,bedgelist=basepy.list()):
        if bedgelist:
            bedgelist=makeRlist(bedgelist)
        edgelist=makeRlist(edgelist)
        gr=ADMG.makeGraph(n*2,basepy.list(),edgelist,bedgelist)
	robjects.globalenv['CPT']=CPT
	robjects.globalenv['gr']=gr
	r('CPT[CPT==0]=.000001')
        out=r('fitMEG(CPT,gr)')
        score=out[2][0]
        return score
            
    def getscore(currentgraph, x,y):
        y=y+n
        dedges=[(x,y)]
        return runScore(dedges)
        
    def Comparescore(direction,currentgraph, x, y):
        y=y+n
        dedges=Intedgelist(currentgraph,n)
        score1=runScore(dedges)
        if direction=="forward":
            dedges.append((x,y))
        else:
            dedges.remove((x,y))
        score2=runScore(dedges)
        return score2-score1

    def ComparescoreBiEdges(direction,currentgraph,x,y):
        y=y+n
        x=x+n
        dedges=Intedgelist(currentgraph,n)
        biedges=IntNoFlipsbedgelist(currentgraph,n)
        score1=runScore(dedges,biedges)
        if direction=="forward":
            biedges.append((x,y))
        else:
            biedges.remove((x,y))
        score2=runScore(dedges,biedges)
        return score2-score1

    def replacezeroes(ct):
        replace=r('function(x) { x[x==0]=.000001; return(x)}')
        return replace(ct)
    
    def makeCT(x):
        nr, nc = x.shape
        xvec = robjects.FloatVector(x.reshape(x.size))
        xr = r.matrix(xvec, nrow=nc, ncol=nr)
        df=basepy.data_frame(xr)
        ct=basepy.table(df)
        ct=replacezeroes(ct)
        return ct
    
    def ForwardDirected():
        allxycomb = [p for p in itertools.product(fr.keys(),repeat=2)]
        CandidateSet=[]
        for x,y in allxycomb:
            score=getscore(g, x, y)
            if score<0:
                CandidateSet.append((score,(x,y)))
        while CandidateSet:
            CandidateSet.sort(reverse=True)
            edgeToAdd=CandidateSet.pop()[1]
            g[str(edgeToAdd[0])][str(edgeToAdd[1])]=set([(0,1)])
            fr[edgeToAdd[1]][edgeToAdd[0]]=1
            reeval_node=edgeToAdd[1]
            newCandidateSet=[]
            for a,e in CandidateSet:
                if e[1]==reeval_node:
                    CandidateSet.remove((a,e))
                    score=Comparescore('forward',g, e[0], e[1])
                    if score<0:
                        newCandidateSet.append((score,(x,y)))
            CandidateSet=newCandidateSet+CandidateSet
        return g,fr
    
    def BackwardDirected():
        edges=gk.edgelist(g)
        CandidateSet=[]
        for e in edges:
            score=Comparescore('backward', g, int(e[0]),int(e[1]))
            if score<0:
                CandidateSet.append((score,(int(e[0]),int(e[1]))))
        while CandidateSet:
            CandidateSet.sort(reverse=True)
            edgetoDel=CandidateSet.pop()[1]
            fr[edgetoDel[1]][edgetoDel[0]]=0
            g[str(edgetoDel[0])][str(edgetoDel[1])]=set()
            reeval_node=edgetoDel[1]
            newCandidateSet=[]
            for a,e in CandidateSet:
                if e[1] ==reeval_node:
                    CandidateSet.remove((a,e))
                    score=Comparescore('backward',g, e[0],e[1])
                    if score<0:
                        newCandidateSet.append((score,(e[0],e[1])))
            CandidateSet = newCandidateSet+CandidateSet
        return g,fr
        
    def ForwardBi():
        fr=emptyGraph(n)
        Bixycomb = [p for p in itertools.combinations(fr.keys(),2)]
        CandidateSet=[]
        newCandidateSet=[]
        for x,y in Bixycomb:
            score=ComparescoreBiEdges('forward',g, x, y)
            if score<0:
                CandidateSet.append((score,(x,y)))
        while CandidateSet:
            CandidateSet.sort(reverse=True)
            edgeToAdd=CandidateSet.pop()[1]
            try:
                g[str(edgeToAdd[0])][str(edgeToAdd[1])].add((2,0))
            except:
                g[str(edgeToAdd[0])][str(edgeToAdd[1])]=set([(2,0)])
            try:
                g[str(edgeToAdd[1])][str(edgeToAdd[0])].add((2,0))
            except:
                g[str(edgeToAdd[1])][str(edgeToAdd[0])]=set([(2,0)])
            fr[edgeToAdd[1]][edgeToAdd[0]]=1
            fr[edgeToAdd[0]][edgeToAdd[1]]=1
            reeval_nodes=edgeToAdd
            newCandidateSet=[]
            for a,e in CandidateSet:
                if e[1] in reeval_nodes or e[0] in reeval_nodes:
                    CandidateSet.remove((a,e))
                    score=ComparescoreBiEdges('forward',g, e[0], e[1])
                    if score<0:
                        newCandidateSet.append((score,(e[0],e[1])))
            CandidateSet=newCandidateSet+CandidateSet
        return g,fr
        
    def BackwardBi():
        edges=gk.bedgelist(g)
        CandidateSet=[]
        for e in edges:
            if e[0]<e[1]:
                score=ComparescoreBiEdges('backward', g, int(e[0]),int(e[1]))
                if score<0:
                    CandidateSet.append((score,(int(e[0]),int(e[1]))))
        while CandidateSet:
            CandidateSet.sort(reverse=True)
            edgeToDel=CandidateSet.pop()[1]
            fr[edgeToDel[1]][edgeToDel[0]]=0
            fr[edgeToDel[0]][edgeToDel[1]]=0
            g[str(edgeToDel[0])][str(edgeToDel[1])].remove((2,0))
            g[str(edgeToDel[1])][str(edgeToDel[0])].remove((2,0))
            reeval_node=edgeToDel
            newCandidateSet=[]
            for a,e in CandidateSet:
                if e[1]==reeval_node[1] or e[0]==reeval_node[0]:
                    CandidateSet.remove((a,e))
                    score=ComparescoreBiEdges('backward',g, e[0], e[1])
                    if score<0:
                        newCandidateSet.append((score,(e[0],e[1])))
            CandidateSet = newCandidateSet+CandidateSet
        return g,fr
    
    import rpy2.robjects.numpy2ri
    rinterface.initr()
    rpy2.robjects.numpy2ri.activate()
    n=data.shape[0]
    data = np.asarray(np.r_[data[:, :-1], data[:, 1:]])
    CPT=makeCT(data)
    g = initemptygraph(n)
    fr=emptyGraph(n)
    g,fr=ForwardDirected()
    g,fr=BackwardDirected()
    g,fr=ForwardBi()
    g,fr=BackwardBi()
    return g
        
        
        
    
