import numpy as np
from graph_base import *


'''
计算一个graph的flip acc
path: graph_topology.json的路径
return: flip acc
'''
def cal_flipacc(path, verbose = False):
    graph = Graph(parse_one_log(path))
    nodes = graph.nodes
    flip = [n.is_inv for n in nodes]
    # gt 即 flip 亦或 correctness
    gt = [n.is_inv ^ n.correctness for n in nodes]
    flip = np.array(flip)
    gt = np.array(gt)
    flip_acc = max(sum(flip == gt)/len(flip), 1-sum(flip == gt)/len(flip))
    if verbose:
        print("flip acc: " + str(flip_acc))
    return flip_acc


'''
计算一个指派的weight_sum
x: 一个指派, 一个n维的向量,取值为0或1
A: 一个n*n的矩阵, A[i,j]表示x[i]和x[j]在指派相同的时候的权重
B: 一个n*n的矩阵, B[i,j]表示x[i]和x[j]在指派不同的时候的权重
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