import string
import numpy as np
import pandas as pd
from gunfolds.tools import graphkit as gk
from gunfolds.tools import conversions as conv

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
    os.system()
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
    num_g2=gk.emptyG(n)
    g2 = conv.ian2g(num_g2)
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
    return conv.dict_format_converter(g2)
