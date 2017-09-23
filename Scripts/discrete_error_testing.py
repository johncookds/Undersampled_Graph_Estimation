import os, sys
from gunfolds.tools import traversal, bfutils
from gunfolds.tools import graphkit as gk
from gunfolds.tools import linear_model as lm
from gunfolds.tools import pc
from multiprocessing import Pool,Process, Queue, cpu_count, current_process
from functools import partial
from gunfolds.tools import conversions as conv
from gunfolds.tools import zickle as zkl
import time, socket
import scipy
import numpy as np
from gunfolds.tools import clingo as cg
from gunfolds.tools.clingo import eqclass
import timeout_decorator
from timeout_decorator import TimeoutError
from itertools import product
import pyfgs
import grangercausality as gc
import discrete_data_generation as dg
import gobnilp as GN
import DiscretePC as DPC
import DiscreteMLLscoringAlgo as MLL

#GRAPH/NETWORK PARAMETERS

ALL_MST=[.9]
ALL_DIST=['normal']
ALL_NOISE_STD=[.01]
ALL_NSAMPLES=[1000]
ALL_BURNIN=[100]
ALL_URATE=[2]
ALL_ALG=['DPC']
GOBNILPCommand='~/Desktop/MRN/GOBNILP/bin/gobnilp -f=dat -g=untitled ~/Desktop/MRN/gunfolds-master/scripts/data.txt'
ALL_CARDINALITY=[2]
densities = {6: [0.15],
    8: [0.15],
    10:[0.15],
    15:[0.15],
    20:[0.15],
    25:[0.15],
    30:[0.15],
    35:[0.15],
    40:[0.15],
    45:[0.15],
    50:[0.15],
    55:[0.15],
    60:[0.15],
    65:[0.15],
    70:[0.15]}

ALL_PARAMS=[ALL_MST,ALL_DIST,ALL_NOISE_STD,ALL_NSAMPLES,ALL_BURNIN,ALL_URATE,ALL_ALG,ALL_CARDINALITY]
#TESTING PARAMETERS
PARALLEL=False
estimate_gu=True
estimate_g1=True
TIMEOUT=3600 # seconds
MSLTIMEOUT = TIMEOUT
SATTIMEOUT = TIMEOUT
POSTFIX='sat4u_clingo'
UMAX = 1
INPNUM = 1 # number of randomized starts per graph
CAPSIZE= 1000 # stop traversing after growing equivalence class tothis size
REPEATS = 1
if socket.gethostname().split('.')[0] == 'leibnitz':
    PNUM=40
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'mars':
    PNUM=12
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'hooke':
    PNUM=21
    PNUM=max((1,PNUM/INPNUM))
else:
    # Setting the number  of parallel running processes  to the number
    # of cores minus 7% for breathing room
    PNUM=cpu_count()-int(0.07*cpu_count())
    PNUM=max((1,PNUM/INPNUM))
print 'processes: ',PNUM, INPNUM


#@timeout_decorator.timeout(MSLTIMEOUT, use_signals=True)
def msl_caller(g2):
    s = set()
    startTime = int(round(time.time() * 1000))
    s = traversal.v2g22g1(g2, capsize=CAPSIZE)
    endTime = int(round(time.time() * 1000))
    msl_time = endTime-startTime
    return s, msl_time

#@timeout_decorator.timeout(SATTIMEOUT, use_signals=False)
def sat_caller(g2, fold):
    startTime = int(round(time.time() * 1000))
    c = eqclass(g2, timeout=SATTIMEOUT+10, capsize=CAPSIZE)
    endTime = int(round(time.time() * 1000))
    sat_time = endTime-startTime
    c = {x[0] for x in c}
    return c, sat_time

def gu_estimate(data,alg,card=2):
    starttime=int(round(time.time() * 1000))
    if alg=='MLL':
        g2=MLL.DiscreteMLL(data)
    elif alg=='DPC':
        g2=DPC.dpc(data,card,.05)
    elif alg=='GOBNILP':
        g2=GN.gobnilp(data)
    endtime=int(round(time.time() * 1000))
    return g2,endtime-starttime

def fan_wrapper(fold='None',n=10,dens=.2,sst=.9,card=2):
    scipy.random.seed()
    msl_time = None
    sat_time = None
    s = None
    c = None
    k = bfutils.dens2edgenum(dens, n=nodes)
    output={}
    while True:
        try:
            output['gt'] = gk.ringmore(n,k)
            gdens = traversal.density(output['gt'])
            gu = bfutils.increment(output['gt'])
            if estimate_gu:
                data=dg.ConstructDynBN(output['gt'],[card]*n,sst,BURNIN+NSAMPLES*URATE)
                data=data[:,BURNIN:]
                gu_est,gu_time=gu_estimate(data[:,::URATE],ALG)
            else:
                gu_est=None
                gu_time=None
            if estimate_g1:
                try:
                    s, msl_time = msl_caller(gu_est)
                except TimeoutError:
                    s = None
                    msl_time = None
                try:
                    c, sat_time = sat_caller(gu_est, fold)
                except TimeoutError:
                    c = None
                    sat_time = None
                if msl_time is not None:
                    print "msl: {:2}: {:8} : {:4}  {:10} seconds".\
                        format(fold, round(sst,3), len(s),
                           round(msl_time/1000.,3))
                if sat_time is not None:
                    print "sat: {:2}: {:8} : {:4}  {:10} seconds".\
                        format(fold, round(sst,3), len(c),
                           round(sat_time/1000.,3))
            else:
                s=None
                c=None
                msl_time=None
                sat_time=None
            output['gu']={'est':gu_est,'truth':gu,'ms':gu_time}
            output['MSL']= {'eq':s,'ms':msl_time}
            output['SAT']= {'eq':c,'ms':sat_time}
        except MemoryError:
            print 'memory error... retrying'
            continue
        break
    return output

for nodes in np.sort(densities.keys())[1:3]:
    for PARAMS in list(product(*ALL_PARAMS)):
        MST,DIST,NOISE_STD,NSAMPLES,BURNIN,URATE,ALG,CARDINALITY=PARAMS
        print nodes, ': ----'
        print ''
        z = {}
        if PARALLEL: pool=Pool(processes=PNUM)
        for dens in densities[nodes]:
            print "{:2}: {:8} : {:10}  {:10}".format('id', 'density', 'eq class', 'time')
            z[dens]={}
            if PARALLEL:
                eqclasses = pool.map(partial(fan_wrapper, n=nodes, dens=dens,sst=MST,card=CARDINALITY), range(REPEATS))
                z[dens]['results'] = eqclasses
            else:
                z[dens]['results'] =[]
                for i in range(REPEATS):
                    z[dens]['results'].append(fan_wrapper(n=nodes,dens=dens,sst=MST,card=CARDINALITY))
            z[dens]['params']={'mst':MST,'dist':DIST,'noise_std':NOISE_STD,'nsamples':NSAMPLES,'burnin':BURNIN,'urate':URATE,'alg':ALG}
            zkl.save(z[dens],
                 socket.gethostname().split('.')[0]+\
                 '_nodes_'+str(nodes)+'_density_'+str(dens)+'_'+POSTFIX+'_.zkl')
            print ''
            print("Saved")
        if PARALLEL:
            pool.close()
            pool.join()
        zkl.save(z,socket.gethostname().split('.')[0]+'_nodes_'+str(nodes)+'_'+POSTFIX+'_.zkl')
        for dens in densities[nodes]:
            os.remove(socket.gethostname().split('.')[0]+\
                  '_nodes_'+str(nodes)+'_density_'+str(dens)+'_'+POSTFIX+'_.zkl')
        print("done")
