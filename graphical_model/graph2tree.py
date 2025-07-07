import numpy as np
from graph_base import BiEdge,BidGraph,swmincut,Graph,ncut

    
class ClusterGraph:
    def __init__(self):
        self.V = [] # 
        self.E = [] # BiEdge, 其中u,v是在V中的idx, w为两个簇之间所有边的权重之和
        
    def add_node(self, cluster):
        '''
        node: NodeSet
        '''
        assert cluster!=set()
        id = len(self.V)
        for c in self.V:
            assert c.isdisjoint(cluster) 
        self.V.append(cluster)
        return id

    def add_edge(self, edge):
        self.E.append(edge)

    # 将ClusterGraph转化为BidGraph 其中每个Cluster被看成一个节点
    def toBiGraph(self):
        G = BidGraph()
        for i in range(len(self.V)):
            G.V.append(i)
        for edg in self.E:
            G.E.append(BiEdge(edg.u,edg.v,edg.w,edg.invw))
        return G

    '''
    转变为graph_base.graph
    其中每个簇的center是簇中心
    每个簇之间存在边,权重为两个簇之间所有边的权重之和
    '''
    def toGraph(self,g):
        class simpleNode:
            def __init__(self,id,center):
                self.id = id
                self.center = center
        class simpleEdge:
            def __init__(self,start,end,w):
                self.start = start
                self.end = end
                self.w = w
        new_g = Graph(None)
        
        
        for i in range(len(self.V)):
            id = i
            center = np.zeros(3)
            for node in self.V[i]:
                center = g.nodes[node].center
            # center /= len(self.V[i])
            new_g.nodes.append(simpleNode(id,center))
        for edg in self.E:
            new_g.edges.append(simpleEdge(edg.u,edg.v,edg.w))
        return new_g
        
    def __st_cut_on_longest_path(self):
        G = self.toBiGraph()
    
    def __stoer_wagner(self):
        assert len(self.V) > 1
        G = self.toBiGraph()
        # cut,f = ncut(G)
        cut,f = swmincut(G)
        S = [] 
        T = [] 
        for i in range(len(cut)):
            if cut[i] == 0:
                S.append(i)
            else:
                T.append(i)
        cut_edge = [] # BiEdge
        for edg in self.E:
            u = edg.u
            v = edg.v
            if u in S and v in T:
                cut_edge.append(edg)
            elif u in T and v in S:
                cut_edge.append(edg)
        if len(S)==0 or len(T) == 0:
            assert False
        return cut_edge,S,T,f

    # cut后,将A,B中的boundary簇合并
    def cut_and_comb(self):
        assert len(self.V) > 1
        cut_edge,A,B,_ = self.__stoer_wagner()
        boundaryA = set()
        boundaryB = set()
        is_boundary = np.zeros(len(self.V))
        
        for edg in cut_edge:
            if edg.u in A and edg.v in B:
                boundaryA = boundaryA.union(self.V[edg.u])
                boundaryB = boundaryB.union(self.V[edg.v])
                is_boundary[edg.u] = True
                is_boundary[edg.v] = True
            elif edg.u in B and edg.v in A:
                boundaryA = boundaryA.union(self.V[edg.v])
                boundaryB =boundaryB.union(self.V[edg.u])
                is_boundary[edg.u] = True
                is_boundary[edg.v] = True
            else:
                assert False
        resA = []
        resB = []
        for i in range(len(self.V)):
            if i in A:
                if not is_boundary[i]:
                    resA.append(self.V[i])
            elif i in B:
                if not is_boundary[i]:
                    resB.append(self.V[i])
            else:
                assert False                
        if boundaryA != set():
            resA.append(boundaryA)
        if boundaryB != set():
            resB.append(boundaryB)
        return resA,resB        

'''
G:BidGraph
A,B:两个簇,例如A = {1,2},B = {3,4}
返回A,B两个簇之间的边列表,边类型为BiEdge,其中u,v是在G.V中的idx
1,2,3,4是G.V中的idx
'''
def getEdgeBetweenCluster(G,A,B):
    edge_set = []
    for edg in G.E:
        u = edg.u
        v = edg.v
        if u in A and v in B or u in B and v in A:
            edge_set.append(edg)
    return edge_set
            
   
def getClusterGraph(G,V):
    # '''
    # G:BidGraph
    # V:Cluster list,每个Cluster是一个集合,其中每个元素都是G.V中的元素
    # '''
    edge_set = {}
    cluster_graph = ClusterGraph()
    node2id = np.array([-1 for _ in range(len(G.V))])
    
    for i in range(len(V)):
        id = cluster_graph.add_node(V[i])
        for node_id in V[i]:
            assert node_id in G.V
            node2id[node_id] = id 
    for edg in G.E:
        u = node2id[edg.u]
        v = node2id[edg.v]
        if u != -1 and v != -1 and u != v:
            if u > v:
                u,v = v,u
            if (u,v) not in edge_set:
                edge_set[(u,v)] = 0
            edge_set[(u,v)] = 1 / (len(V[u]) * len(V[v]))
            # edge_set[(u,v)] += edg.w
    for key in edge_set:
        u,v = key
        cluster_graph.add_edge(BiEdge(u,v,edge_set[key],edge_set[key]))
    return cluster_graph
    
    

'''
queue = [getClusterGraph(G,G.V)]
res = []
while queue is not empty:
    G = queue.pop()
    cut_edge,S,T = stoer_wagner(G)
    boundary = {cut_edge.u,cut_edge.v}
    A = S - boundary + comb(boundary ∩ S),B = T - boundary + comb(boundary ∩ T)
    if len(A) > 1:
        queue.append(G,A) # 图
    else:
        res.append(A) # 点集
    
    if len(B) > 1:
        queue.append(G,B)
    else:
        res.append(B)
return res
'''
def graph2tree(G):
    '''
    G:BidGraph
    '''
    res = []
    queue = []
    cid = G.get_clusters()
    cid_set = set(cid)
    for c in cid_set:
        sub_V = [{G.V[i]} for i in range(len(G.V)) if cid[i] == c]
        if len(sub_V) > 1:
            queue.append(getClusterGraph(G,sub_V))
        else:
            res.append(sub_V)
    print("联通分量数：",len(queue) + len(res))
    while len(queue) > 0:
        X = queue.pop()
        tg = X.toBiGraph()
        assert tg.get_cluster_num() == 1
        A,B = X.cut_and_comb()
        clusters = sorted(X.V,key=lambda x:len(x),reverse=True)
        if len(A) > 1:
            queue.append(getClusterGraph(G,A))
        else:
            res.append(A)
        if len(B) > 1:
            queue.append(getClusterGraph(G,B))
        else:
            res.append(B)
    return res  
    

        