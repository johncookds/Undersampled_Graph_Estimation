from pycausal import search
from pycausal import pycausal as pc
from gunfolds.tools import conversions as conv
import re
import numpy as np
import pandas as pd

def pyfgs(d):
    '''
    fast greedy equivalence search from pycausal
    https://github.com/bd2kccd/py-causal
    Interpretation reasoning is provided in paper linked in repo README
    '''
    pc.start_vm()
    n = d.shape[0]
    d = np.asarray(np.r_[d[:, :-1], d[:, 1:]])
    d=pd.DataFrame(d)
    d=d.T
    d.columns=[str(x+1)+'*' for x in range(n)]+[str(x+1) for x in range(n)]
    fgs=search.fges(d,penaltydiscount = 2,
    faithfulnessAssumed = True, verbose = True)
    edg=fgs.edges
    print(edg)
    g2={}
    for i in np.arange(1,n+1):
        g2[str(i)]={}
    for i in edg:
        if re.match(r'[0-9]+\s.*[0-9]+\Z',i):
            y=re.findall(r'[0-9]+',i)
            if y[1] in g2[y[0]].keys():
                g2[y[0]][y[1]].add((2,0))
            else:
                g2[y[0]][y[1]]=set([(2,0)])
            if y[0] in g2[y[1]].keys():
                g2[y[1]][y[0]].add((2,0))
            else:
                g2[y[1]][y[0]]=set([(2,0)])
        elif re.match(r'[0-9]+\*.*[0-9]+\Z',i):
            y=re.findall(r'[0-9]+',i)
            if y[1] in g2[y[0]].keys():
                g2[y[0]][y[1]].add((0,1))
            else:
                g2[y[0]][y[1]]=set([(0,1)])
        elif re.match(r'[0-9]+\s.*[0-9]+\*',i):
            y=re.findall(r'[0-9]+',i)
            if y[0] in g2[y[1]].keys():
                g2[y[1]][y[0]].add((0,1))
            else:
                g2[y[1]][y[0]]=set([(0,1)])
    return conv.dict_format_converter(g2)