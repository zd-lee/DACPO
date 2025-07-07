import os
import json
import numpy as np
import open3d as o3d
from draw_topology import add_topology, get_arrow ,get_sphere

class GraphVertex:
    def __init__(self,json_j):
        self.id = json_j['id']
        self.inv_time = json_j['log']['inv_times']
        self._size = json_j['metric']['total_count']
        self.loss = json_j['metric']['avg_nd_loss']
        self.is_inv = self.inv_time%2 == 1
        self.correctness = (json_j['metric']['avg_nd_loss'] == json_j['metric']['avg_loss'])
        self.center = json_j['center']
        self.center = np.array(self.center)

          
class Edge:
    def __init__(self,json_j):
        self.start = json_j['start']
        self.end = json_j['end']
        self.correctness = json_j['correctness']
        self.weight = json_j['weight']
        self.inv_weight =  json_j['inv_weight']
        self.confidence = json_j['confidence']
        self.node_loss = 0 #需要由外部计算
        
class Graph:
    def __init__(self,json_j):
        if json_j == None:
            self.nodes = []
            self.edges = []
            return
        
        nodes_counts = len(json_j['graph_topology'])
        self.nodes = []
        self.edges = []
        self.flip_acc = json_j['graph_metric']['flip_acc']
        self.edge_acc = json_j['graph_metric']['edge_acc']
        
        for i in range(nodes_counts):
            node_unit = json_j['graph_topology'][i]
            self.nodes.append(GraphVertex(node_unit['vertex']))
            
        
        for i in range(nodes_counts):
            node_unit = json_j['graph_topology'][i]
            edgesi = node_unit['edges']
            for i in range(len(edgesi)):
                e = Edge(edgesi[i])
                sloss = self.nodes[e.start].loss
                eloss = self.nodes[e.end].loss
                e.node_loss = max(sloss,eloss)
                self.edges.append(e)        

    def to_matrix(self,bool_inv = False):
        # 检验nodes的id是否是连续的
        set_id = set([i.id for i in self.nodes])
        for i in range(len(set_id)):
            assert i in set_id
        n = len(self.nodes)
        A = np.zeros((n,n))
        for e in self.edges:
            if bool_inv:
                A[e.start][e.end] = e.inv_weight
                A[e.end][e.start] = e.inv_weight
            else:
                A[e.start][e.end] = e.weight
                A[e.end][e.start] = e.weight
        return A
    
 
    def draw_topology(self,nodelabel=[],edgelabel=[]):
        '''
        用open3d画出图的拓扑结构
        其中 每个节点,在center处用get_sphere画一个球
        每个边,根据其起点和终点的位置,用get_arrow
        '''
        if len(nodelabel) == 0:
            nodelabel = np.zeros(len(self.nodes))
        if len(edgelabel) == 0: 
            edgelabel = np.zeros(len(self.edges))
        
        mesh = ([], [])
        colors = []
        unque_label = list(set(nodelabel).union(set(edgelabel)))
        label2color = {}
        for i in range(len(unque_label)):
            label2color[unque_label[i]] = np.random.rand(3)
        
    
        for i in range(len(self.nodes)):
            center = self.nodes[i].center
            assert i == self.nodes[i].id
            sp = get_sphere(center)
            add_topology(mesh,get_sphere(center))
            colors += [label2color[nodelabel[i]] for _ in range(len(sp[0]))]
        for i in range(len(self.edges)):
            start = self.nodes[self.edges[i].start].center
            end = self.nodes[self.edges[i].end].center
            arrow = get_arrow(start,end)
            add_topology(mesh,arrow)
            colors += [label2color[edgelabel[i]] for _ in range(len(arrow[0]))]    
        o3dmesh = o3d.geometry.TriangleMesh()
        o3dmesh.vertices = o3d.utility.Vector3dVector(mesh[0])
        o3dmesh.triangles = o3d.utility.Vector3iVector(mesh[1])
        o3dmesh.vertex_colors = o3d.utility.Vector3dVector(colors)
        o3d.visualization.draw_geometries([o3dmesh])
        return o3dmesh

def parse_one_log(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    return data 

class BiEdge:
    def __init__(self,u,v,w,invw):
        self.u = u
        self.v = v
        self.w = w
        self.invw = invw
    
    def __iter__(self):
        return iter([self.u,self.v,self.w,self.invw])

class BidGraph:
    def __init__(self):
        self.V = [] # id list
        self.E = [] # BiEdge list
    
    def checkVaild(self):
        # 检查E中的顶点是否都在V中
        for edg in self.E:
            u = edg.u
            v = edg.v
            assert u in self.V
            assert v in self.V

    def induced_subgraph(self,V):
        G = BidGraph()
        G.V = V
        for edg in self.E:
            u = edg.u
            v = edg.v
            if u in V and v in V:
                G.E.append(edg)
        return G
    
    def get_cut_edges(self,A,B):
        edgelist = []
        for edg in self.E:
            u = edg.u
            v = edg.v
            if u in A and v in B:
                edgelist.append(edg)
            if u in B and v in A:
                edgelist.append(edg)
        return edgelist
    
    def from_matrix(self,A):
        n = A.shape[0]
        self.V = [i for i in range(n)]
        for i in range(n):
            for j in range(i+1,n):
                if A[i,j] != 0:
                    assert A[i,j] == A[j,i] # 要求是无向图 即对称矩阵
                    self.E.append(BiEdge(i,j,A[i,j],A[i,j])) # 注意不要把invw与反向边混淆
        return self
    
    def to_matrix(self):
        # 检验nodes的id是否是连续的
        set_id = set(self.V)
        for i in range(len(set_id)):
            if not i in set_id:
                print("Error! the graph's vertex id is not continuous")
                assert i in set_id
        
        n = len(self.V)
        A = np.zeros((n,n))
        for edg in self.E:
            A[edg.u][edg.v] = edg.w
            A[edg.v][edg.u] = edg.w
        return A

    def get_clusters(self):
        '''
        返回图的连通分量
        return: 每个顶点所在的连通分量的id
        '''
        root = dict()
        for i in self.V:
            root[i] = i
        
        def find(x):
            if root[x] != x:
                root[x] = find(root[x])
            return root[x]
        def union(x,y):
            root[find(x)] = find(y)
        for edg in self.E:
            union(edg.u,edg.v)
        for i in self.V:
            find(i)
        res = [root[i] for i in self.V]    
        return res
    
    def get_cluster_num(self):
        return len(set(self.get_clusters()))  
    
    def has_cycle(self):
        '''
        判断是否有环
        return: bool
        '''
        mat = self.to_matrix()
        n = len(self.V)
        visited = np.zeros(n)
        def dfs(u,pre):
            visited[u] = 1
            for v in range(n):
                if mat[u,v] != 0 and visited[v] == 0:
                    if dfs(v,u):
                        return True
                elif mat[u,v] != 0 and visited[v] == 1 and v != pre:
                    return True
            return False
        for i in range(n):
            if visited[i] == 0 and dfs(i,-1):
                return True
        return False
       
    def get_longest_path(self):
        '''
        得到最长的最短路径
        return: pathlist, 例如[[0,1],[1,2],[2,3]]
        '''
        n = len(self.V)
        assert n<1000 #O(N^3)
        dp = np.full((n,n),1e9,dtype=int)
        for edg in self.E:
            u = edg.u
            v = edg.v
            dp[u,v] = 1
            dp[v,u] = 1
        for i in range(n):
            dp[i,i] = 0
        
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    dp[i,j] = min(dp[i,j],dp[i,k]+dp[k,j])

        max_length = 1e9
        s = -1
        t = -1
        for i in range(n):
            for j in range(n):
                if dp[i,j]<1e9:
                    if s==-1 and t==-1:
                        s,t = i,j
                    elif dp[i,j]>dp[s,t]:
                        s,t = i,j
                    
        # 从s到t的路径
        mat = np.zeros((n,n))
        for edg in self.E:
            u = edg.u
            v = edg.v
            mat[u,v] = 1
            mat[v,u] = 1
        dist = np.full(n,1e9)
        dist[s] = 0
        pre = np.zeros(n,dtype=int)
        for i in range(n):
            if i != s:
                pre[i] = s
        pq = [s]
        while len(pq) > 0:
            u = pq.pop(0)
            for v in range(n):
                if mat[u,v] != 0 and dist[v] > dist[u]+1:
                    dist[v] = dist[u]+1
                    pre[v] = u
                    pq.append(v)
        assert dist[t] == dp[s,t]
        # 从t到s的路径
        path = []
        u = t
        while u != s:
            path.append(u)
            u = pre[u]
        path.append(s)
        path.reverse()
        pathlist = []
        for i in range(len(path)-1):
            pathlist.append([path[i],path[i+1]])
        return pathlist


import maxflow 
def mincut_LR(G,L,R):
    '''
    G:BidGraph类型
    给定图G,两个点集L和R,求一个最小割,
    使得L中的点都在割的S中(lable = 0),R中的点都在割的T中(lable = 1)
    '''
    LR = L + R
    for i in LR:
        assert i in G.V
    for i in L:
        assert i not in R
    for i in R:
        assert i not in L
    
    n = len(G.V)
    mG = maxflow.Graph[float](n, n)
    nodes = mG.add_nodes(n)
    uv2node_id = {}
    for i in range(n):
        uv2node_id[G.V[i]] = nodes[i]
     
    for edg in G.E:
        a = uv2node_id[edg.u]
        b = uv2node_id[edg.v]
        mG.add_edge(a,b, edg.w, edg.w)# 无向图
    for i in L:
        mG.add_tedge(nodes[uv2node_id[i]], 1e9, 0)
    for i in R:
        mG.add_tedge(nodes[uv2node_id[i]], 0, 1e9)
    f = mG.maxflow()
    cut = mG.get_grid_segments(nodes)
    return cut, f

def mincut_S(G,S):
    propable_T = [i for i in G.V if i not in S]
    assert len(propable_T) > 0
    min_flow = 1e9
    for T in propable_T:
        _,flow = mincut_LR(G,S,[T])
        if flow < min_flow:
            min_flow = flow
            min_T = T
    return mincut_LR(G,S,[min_T])

# 粗糙版的stoer_wagner算法 嘿嘿
def swmincut(G):
    '''
    G:BidGraph类型
    return:cut, flow
    '''
    assert G.get_cluster_num() == 1
    if len(G.V) == 1:
        assert False
    S = [G.V[0]]
    cut,f =  mincut_S(G,S)
    return cut,f


from sklearn.cluster import KMeans
def ncut(G):
    x = G.V
    similarity = G.to_matrix()
    n = len(x)
    assert n>1
    assert similarity.shape[0] == n
    assert similarity.shape[1] == n
    D = np.diag(np.sum(similarity,axis=1))
    L = D - similarity
    D_inv = np.linalg.inv(D)
    L_sym = D_inv @ L
    eig_val,eig_vec = np.linalg.eig(L_sym)
    eig_val = eig_val.real
    eig_vec = eig_vec.real
    idx = np.argsort(eig_val)
    eig_val = eig_val[idx]
    eig_vec = eig_vec[:,idx]
    k = 2
    eig_vec = eig_vec[:,1:k]
    eig_vec = eig_vec / np.linalg.norm(eig_vec,axis=1).reshape(-1,1)
    kmeans = KMeans(n_clusters=k, random_state=0).fit(eig_vec)
    cluster = kmeans.labels_
    cluster = list(cluster)
    return cluster,0
    
    