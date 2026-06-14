import numpy as np
from graph_base import *


'''
Calculate flip accuracy for a graph.
path: path to graph_topology.json
return: flip acc
'''
def cal_flipacc(path, verbose = False):
    graph = Graph(parse_one_log(path))
    nodes = graph.nodes
    flip = [n.is_inv for n in nodes]
    gt = [n.is_inv ^ n.correctness for n in nodes]
    flip = np.array(flip)
    gt = np.array(gt)
    flip_acc = max(sum(flip == gt)/len(flip), 1-sum(flip == gt)/len(flip))
    if verbose:
        print("flip acc: " + str(flip_acc))
    return flip_acc


'''
x: binary assignment vector of length n
A: n*n same-assignment weight matrix
B: n*n different-assignment weight matrix
return: weight_sum
'''
def cal_loss(x,A,B):
    n = len(x)
    # x = np.array(x, dtype=int)
    assert A.shape == (n,n)
    assert B.shape == (n,n)
    obj = 0
    for i in range(n):
        for j in range(n):
            obj += A[i,j]*(1 - (x[i]-x[j])*(x[i]-x[j])) + B[i,j]*(x[i]- x[j])*(x[i] - x[j])
    return obj
