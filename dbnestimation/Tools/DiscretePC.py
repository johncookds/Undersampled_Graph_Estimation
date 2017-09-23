from gunfolds.tools import bfutils as bfu
from gunfolds.tools import graphkit as gk
from gunfolds.tools import linear_model as lm
from gunfolds.tools import traversal as trv
from gunfolds.tools import zickle as zkl
import itertools
from itertools import combinations,permutations
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore, norm
import statsmodels.api as sm
from scipy.stats import chi2_contingency
import pandas as pd
import itertools
from scipy.stats import chisqprob
from gunfolds.tools import conversions as conv
import copy

def makeDF(data):
    data=pd.DataFrame(data)
    return data

def dpc(data, varst,pval=0.1):
    n = data.shape[0]
#    if n<200:
#        pval=.05
#    if n <1000:
#        pval=.1
#    elif n<2000:
#        pval=.1
    # stack the data: first n rows is t-1 slice, the next n are slice t
    data = np.asarray(np.r_[data[:, :-1], data[:, 1:]])

    def cindependent(y, x, counter, parents=[], pval=pval):
        for S in [j for j in combinations(parents, counter)]:
            print S
            if ChiSquaredTest(x, y, condset=list(S)):
                return True
        return False

    def bindependent(y, x, parents=[], pval=pval):
         print "done"
         return ChiSquaredTest(x, y, condset=parents, shift=n)

    def dir_prune(elist, mask, g):
        for e in mask:
            sett=copy.deepcopy(g[e[0]][e[1]])
            sett.remove((0,1))
            g[e[0]][e[1]]=sett
            elist.remove(e)

    def bi_prune(mask, g):
        for e in mask:
            sett=copy.deepcopy(g[e[0]][e[1]])
            sett.remove((2,0))
            g[e[0]][e[1]]=sett
            g[e[1]][e[0]]=sett

    def chisq_of_df_cols(df, c1, c2):
        groupsizes = df.groupby([c1, c2]).size()
        ctsum = groupsizes.unstack(c1)
    # fillna(0) is necessary to remove any NAs which will cause exceptions
        return(chi2_contingency(ctsum.fillna(0))) 
        
    def ChiSquaredTest(x,y,condset,shift=0):
        if condset:
            X = data[[shift+int(x)-1]+[n+int(y)-1]+condset,:].T
            df=makeDF(X)
            condnum=df.shape[1]-2
            for i in range(condnum):
                if i==0:
                    v=pd.Series.unique(df[i+2])
                else:
                    v=np.vstack([v,pd.Series.unique(df[i+2])])
            if condnum==1:
                condvalues=[v]
            else:
                condvalues=list(itertools.product(*v))
            chis=0
            dofs=0
            for i in condvalues:
                count=1
                for j in range(condnum):
                    if count==1:
                        newdf=df[df[j+2]==i[j]]
                        count=2
                    else:
                        newdf=newdf[newdf[j+2]==i[j]]
                try:
                    chi2, p, dof, ex = chisq_of_df_cols(newdf,0,1)
                except:
                    chi2=0
                    dof=(varst-1)**2
                chis+=chi2
                dofs+=dof
            val=chisqprob(chis,dofs)
        else:
            X = data[[shift+int(x)-1]+[n+int(y)-1],:].T
            df=makeDF(X)
            chi2, p, dof, ex = chisq_of_df_cols(df,0,1)
            val=chisqprob(chi2,dof)
        return val > pval

    def stringify(array):
        d=[]
        for i in array:
            d.append((str(i[0]),str(i[1])))
        return d
        
    num_g  = gk.superclique(n)
    el = gk.edgelist(num_g)
    el = stringify(el)
    print(el)
    num_gtr = gk.gtranspose(num_g)
    gtr=conv.ian2g(num_gtr)
    g = conv.ian2g(num_g)
    for counter in range(n):
        to_remove = []
        for e in el:
            ppp = [int(k)-1 for k in gtr[e[1]] if k != e[0]]
            if counter <= len(ppp):
                if cindependent(e[1], e[0], counter, parents=ppp, pval=pval):
                    to_remove.append(e)
                    gtr[e[1]].pop(e[0], None)
        dir_prune(el, to_remove, g)
    print(g)
    bel = [map(lambda k: str(k+1), x) for x in combinations(range(n), 2)]
    bi_list=[]
    for e in bel:
        ppp = list(set(gtr[e[0]].keys()) | set(gtr[e[1]].keys()))
        ppp = map(lambda x: int(x)-1, ppp)
        if bindependent(e[0], e[1], parents=ppp, pval=pval):
            bi_list.append(e)
    bi_prune(bi_list, g)
    g=conv.dict_format_converter(g)
    gk.clean_leaf_nodes(g)
    return g


# Local Variables:
# mode: python
# python-indent-offset: 4
# End:
