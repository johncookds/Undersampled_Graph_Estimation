# Undersampled_Graph_Estimation

Algorithms for graph estimation and data generation built on top of [gunfolds](https://gitlab.com/undersampling/gunfolds).

Built for Causal UAI 2017 paper [Learning Dynamic Structure from Undersampled Data](https://www.cs.purdue.edu/homes/eb/causal-uai17/papers/7.pdf).

Python 2.7

Packages/software used in the making of the paper and repo:  
gunfolds 
pycausal  
libpgm   
GOBNILP   

# Algorithm status:

(9/23/17)

## discrete_data_generation:
Working with latest gunfolds and latest libpgm

## continuous algorithms:

SVAR: working (uses gunfolds.tools version)

PC: working (uses gunfolds.tools version)

GC: working

pyfgs: working (not running in parallel due to java vm issues)

## discrete algorithms:

DPC: working

GOBNILP: working, requires [GOBNILP](https://www.cs.york.ac.uk/aig/sw/gobnilp/)

MLL: not working, excluded from paper due to poor performance

# To Do:
  1. explore new functions for fast greedy equivalence search now in latest pycausal release

Acknowledgment
========
This work was, in part, supported by  NSF IIS-1318759 grant.
