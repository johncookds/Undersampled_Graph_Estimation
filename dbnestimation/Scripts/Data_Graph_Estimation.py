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


TIMEOUT=3600 # seconds
MSLTIMEOUT = TIMEOUT
SATTIMEOUT = TIMEOUT
POSTFIX='sat4u_clingo'
UMAX = 1
INPNUM = 1 # number of randomized starts per graph
CAPSIZE= 1000 # stop traversing after growing equivalence class to this size
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

def gu_estimate(data,alg):
    starttime=int(round(time.time() * 1000))
    if alg=='pc':
        g2 = pc.dpc(data, pval=0.05)
    elif alg=='svar':
        g2 = lm.data2graph(data)
    elif alg=='gc':
        g2 = gc.gc(data,pval=.05)
    elif alg=='pyfgs':
        g2= pyfgs.pyfgs(data)
    endtime=int(round(time.time() * 1000))
    return g2,endtime-starttime

def Graph_Estimation(data,ALG,estimate_gu=True,estimate_g1=False,fold='None'):
    scipy.random.seed()
    msl_time = None
    sat_time = None
    s = None
    c = None
    output={}
    if estimate_gu:
        gu_est,gu_time=gu_estimate(data,ALG)
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
    return output


