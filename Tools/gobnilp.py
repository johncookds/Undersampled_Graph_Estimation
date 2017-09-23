import signal
import pprint
import time, socket
import numpy as np
import scipy
import functools, itertools
import progressbar as pb
import sys,os
from scipy.misc import comb
sys.path.append('/Users/johncook/Desktop/MRN/gunfolds-master/scripts')
import copy
from multiprocessing import Pool,Process, Queue, cpu_count, current_process
from progressbar import ProgressBar, Percentage, \
    Bar, RotatingMarker, ETA, FileTransferSpeed
import linear_model as lm
from gunfolds.tools import traversal as trv
import modbfutils as bfu
from gunfolds.tools import graphkit as gk
from gunfolds.tools import zickle as zkl
import pylab as plt
import modunknownrate as ur
import pandas as pd
import re
import DiscretePC as DPC
import DiscreteMLLscoringAlgo as MLL
import string
def gobnilp(d):
    d = np.asarray(np.r_[d[:, :-1], d[:, 1:]])
    d=pd.DataFrame(d)
    d.index=d.index+1
    d=d.T
    d.loc[-1]=[varst]*n*2
    d.index=d.index+1
    d=d.sort()
    d.columns=list(string.ascii_uppercase[:n*2])
    d.to_csv('data.txt',sep='\t',index=False)
    os.system(bashCommand)
    f=open('Output.txt')
    edges=[]
    bedges=[]
    for line in f:
        let=letters(line)
        if let=='BNscoreis':
            break
        count=0
        let=[string.uppercase.index(k) for k in let]
        for i in let:
            new=i+1
            if count==0:
                first=i+1
                if first <=n:
                    state="lt"
                else:
                    state="gt"
                count=1
            elif state=="lt":
                if new >n:
                    edges.append([first,new-n])
                if new <n:
                    bedges.append([first,new])
                    bedges.append([new,first])
            elif state=="gt":
                if new >n:
                    bedges.append([first-n,new-n])
                    bedges.append([new-n,first-n])
                if new <=n:
                    edges.append([new, first-n])
    g2=gk.emptyG(n)
    used=[]
    for i,j in edges:
        if [i,j] in bedges or [j,i] in bedges:
            g2[str(i)][str(j)]=set([(0,1),(2,0)])
            used.append([i,j])
        else:
           g2[str(i)][str(j)]=set([(0,1)])
    for i,j in bedges:
        if [i,j] not in used:
            g2[str(i)][str(j)]=set([(2,0)])