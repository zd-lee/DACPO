#include <algorithm>
#include <bgl_api.h>
#include <queue>
#include <graph_arg.h>
#include <socket_api.h>

using namespace graph_arg;

void print_metric(nlohmann::json j) {
    std::cout << "weight_sum: " << j["weight_sum"] << "\t";
    std::cout << "edge_acc: " << j["edge_acc"] << "\t";
    std::cout << "flip_acc: " << j["flip_acc"] << std::endl;
}

graph_arg::FlipableNode::FlipableNode(nlohmann::json j):id(j["id"].get<int>())
{
    inv_time = 0;
    REAL loss = j["metric"]["avg_loss"].get<double>(),inv_loss = j["metric"]["avg_nd_loss"].get<double>();
    same_to_gt = (loss == inv_loss);
    _size = j["metric"]["total_count"];
}

graph_arg::FlipableNode::FlipableNode(nid id, bool if_gt) :id(id) {
    inv_time = 0;
}

bool graph_arg::FlipableNode::is_inv() const
{
	return inv_time % 2 == 1;
}

bool graph_arg::FlipableNode::is_same() const
{
    return same_to_gt != is_inv();
}

int graph_arg::FlipableNode::flip()
{
    inv_time++;
    return inv_time;
}

graph_arg::ConfidenceEdge::ConfidenceEdge(const FlipableNode* start, const FlipableNode* end, REAL weight, REAL inv_weight, REAL confidence):
     start(start), end(end), weight(weight), inv_weight(inv_weight), confidence(confidence)
 {
 }

void graph_arg::align_to_root::visit(graph_arg::Tree* root) {
    std::vector<std::vector<graph_arg::FlipableNode*>> g_n; // 存储每个子树的顶点列表
    std::vector<graph_arg::FlipableNode*> tmp_g_n; // 临时图的顶点列表
    std::vector<std::vector<ConfidenceEdge*>>  tmp_g_e; // 临时图的边矩阵
    
    g_n.push_back({root->_root});
    for(int i = 0;i<root->_childrens.size();i++){
        if(root->_childrens[i] == nullptr)continue;
        g_n.push_back(root->_childrens[i]->get_all_nodes());
    }
    tmp_g_n.resize(g_n.size());
    tmp_g_e.resize(g_n.size());
    for(int i = 0;i<g_n.size();i++){
        tmp_g_e[i].resize(g_n.size());
    }

    for(int i = 0;i<tmp_g_n.size();i++){
        tmp_g_n[i] = new graph_arg::FlipableNode(i);
    }
    for(int i = 0;i<tmp_g_n.size();i++){
        for(int j = 0;j<tmp_g_n.size();j++){
            if (j == i)continue;
            weight_pair w = graph->cal_weight(g_n[i],g_n[j]);
            tmp_g_e[i][j] = new graph_arg::ConfidenceEdge(tmp_g_n[i],tmp_g_n[j],w.first,w.second);
        }
    }
    graph_arg::FlipableGraph* tgraph = new graph_arg::FlipableGraph(tmp_g_n,tmp_g_e);
    BruteForceFlip bf;
    std::vector<bool> flip_res = bf.flip(tgraph);
    for(int i = 0;i<flip_res.size();i++){
        if(flip_res[i]){
            for(int j = 0;j<g_n[i].size();j++){
                g_n[i][j]->flip();
            }
        }
    }
    for(int i = 0;i<tmp_g_n.size();i++){
        delete tmp_g_n[i];
    }
    for(int i = 0;i<tmp_g_e.size();i++){
        for(int j = 0;j<tmp_g_e[i].size();j++){
            delete tmp_g_e[i][j];
        }
    }
    delete tgraph;
}

std::vector<bool> BruteForceFlip::flip(graph_arg::FlipableGraph* g){
    assert(g->_nodes.size()<33);
    std::vector<bool> res(g->_nodes.size(),false);
    REAL w = 1e10;
    std::vector<bool> ori_flip(g->_nodes.size(),false);
    for(int i = 0;i<g->_nodes.size();i++){
        ori_flip[i] = g->_nodes[i]->is_inv();
    }

    for(long long i = 0;i<(1<<g->_nodes.size());i++){
        for(int j = 0;j<g->_nodes.size();j++){
            if(i&(1<<j)){
                g->_nodes[j]->flip();
            }
        }
        REAL tmp = g->cal_current_weight_sum();
        if(tmp<w){
            w = tmp;
            for(int j = 0;j<g->_nodes.size();j++){
                res[j] = g->_nodes[j]->is_inv();
            }
        }
        //回退翻转
        for(int j = 0;j<g->_nodes.size();j++){
            if(i&(1<<j)){
                g->_nodes[j]->flip();
            }
        }

        // 检查是否有错误
        for(int j = 0;j<g->_nodes.size();j++){
            assert(g->_nodes[j]->is_inv() == ori_flip[j]);
        }
    }
    for(int i = 0;i<res.size();i++){
        g->_nodes[i]->flip();
    }
    return res;
}

nlohmann::json graph_arg::BruteForceFlip::get_config()
{
    nlohmann::json j;
    j["name"] = "BruteForceFlip";
    return j;
}

nlohmann::json graph_arg::BruteForceFlip::get_log()
{
    return nlohmann::json();
}

graph_arg::Tree::Tree(FlipableNode* root):_root(root)
{
}

void graph_arg::Tree::add_child(Tree* t){
    _childrens.push_back(t);
}

void graph_arg::Tree::PostTraverse(TreeVistor *vistor)
{
    for(int i = 0;i<_childrens.size();i++){
        if(_childrens[i] == nullptr)continue;
        _childrens[i]->PostTraverse(vistor);
    }
    vistor->visit(this);
}

std::vector<FlipableNode *> graph_arg::Tree::get_all_nodes()
{
    std::vector<FlipableNode*> res;
    res.push_back(_root);
    for(int i = 0; i < _childrens.size(); i++){
        std::vector<FlipableNode*> tmp = _childrens[i]->get_all_nodes();
        res.insert(res.end(), tmp.begin(), tmp.end());
    }
    return res;
}

void graph_arg::FlipableGraph::_align_flip(const std::vector<bool> &flip1, std::vector<bool> &flip2)
{
    assert(flip1.size() == flip2.size());
    assert(flip1.size() == _nodes.size());
    int dist = 0, inv_dist = 0;
    for(int i = 0; i < flip1.size(); i++){
        if(flip1[i] != flip2[i]){
            dist += _nodes[i]->_size;
        }
        else{
            inv_dist += _nodes[i]->_size;
        }
    }
    if(dist > inv_dist){
        for(int i = 0; i < flip2.size(); i++)flip2[i] = !flip2[i];
    }
}

void graph_arg::FlipableGraph::reset_flip_status()
{
    for (int i = 0; i < _nodes.size(); i++) {
        if (_nodes[i]->is_inv())_nodes[i]->flip();
    }
}

void graph_arg::FlipableGraph::apply_flip(const std::vector<bool> &node_status)
{
    assert(node_status.size() == _nodes.size());
    for(int i = 0; i < _nodes.size(); i++){
        if(node_status[i] == _nodes[i]->is_inv()){// 注意这里的逻辑
            _nodes[i]->flip();
        }
    }
}

graph_arg::FlipableGraph::FlipableGraph(std::vector<graph_arg::FlipableNode *> nodes, std::vector<std::vector<graph_arg::ConfidenceEdge *>> edges)
{
    _nodes = nodes;
    _edges = std::vector<std::vector<graph_arg::ConfidenceEdge*>>(_nodes.size(),std::vector<graph_arg::ConfidenceEdge*>(_nodes.size(),nullptr));
    for(int i = 0;i<_nodes.size();i++){
        for(int j = 0;j<_nodes.size();j++){
            if(edges[i][j] == nullptr)continue;
            _edges[i][j] = edges[i][j];
        }
    }
}

graph_arg::FlipableGraph::FlipableGraph(nlohmann::json j)
{
    using namespace nlohmann;
    json _items = j["graph_topology"];
    assert(_items.size()>0);

    _nodes.clear();
    _nodes.resize(_items.size());
    _edges.resize(_items.size(),std::vector<graph_arg::ConfidenceEdge*>(_items.size(),nullptr));
   
    for (int i = 0; i < _items.size(); i++)
    {
        json v = _items[i]["vertex"];
        _nodes[i] = new graph_arg::FlipableNode(v);
    }
    for(int i = 0;i<_items.size();i++){
        json edges_i = _items[i]["edges"];
        for(int j = 0;j<edges_i.size();j++){
            int s = edges_i[j]["start"].get<int>();
            int t = edges_i[j]["end"].get<int>();
            weight_pair w = {edges_i[j]["weight"].get<double>(),edges_i[j]["inv_weight"].get<double>() };
            REAL conf = edges_i[j]["confidence"].get<double>();
            _edges[s][t] = new graph_arg::ConfidenceEdge(_nodes[s],_nodes[t],w.first,w.second,conf);
        }
    }
    // boolize_weight();
    // multiply_conf();
}

inline int find(std::vector<int>& father,int x){
    if(father[x] == x)return x;
    return father[x] = find(father,father[x]);
}

Tree* Kruskal(std::vector<graph_arg::FlipableNode*> node, std::vector<graph_arg::ConfidenceEdge*> edges){
    
    // shuffle 
    std::random_shuffle(edges.begin(),edges.end());
    std::vector<graph_arg::ConfidenceEdge*> res_edges;
    std::map<int,int> nid2idx;
    for(int i = 0;i<node.size();i++){
        nid2idx[node[i]->id] = i;
    }
    // 先得到一个无向的树
    std::vector<std::vector<int>> tmptree(node.size(),std::vector<int>());
    std::vector<int> father(node.size(),-1);
    for(int i = 0;i<node.size();i++){
        father[i] = i;
    }
    
    
    int cnt = 0;
    for(int i = 0;i<edges.size();i++){
        int s = nid2idx[edges[i]->start->id], t = nid2idx[edges[i]->end->id];
        assert(s != t);
        if(find(father,s) == find(father,t))continue;
        tmptree[s].push_back(t);
        tmptree[t].push_back(s);
        father[find(father,s)] = find(father,t);
        cnt++;
    }    
    assert(cnt == node.size() - 1);
    // 从树中随机选择一个点作为根节点,并生成一个树
    int root = rand()%node.size();
    std::vector<std::vector<int>> tree(node.size(),std::vector<int>());
    std::vector<bool> visited = std::vector<bool>(node.size(),false);
    std::queue<int> q;
    q.push(root);
    visited[root] = true;
    while(!q.empty()){
        int cur = q.front();
        q.pop();
        for(int i = 0;i<tmptree[cur].size();i++){
            if(visited[tmptree[cur][i]])continue;
            tree[cur].push_back(tmptree[cur][i]);
            visited[tmptree[cur][i]] = true;
            q.push(tmptree[cur][i]);
        }
    }
    // 检查visited是否全为true
    for(int i = 0;i<visited.size();i++){
        assert(visited[i]);
    }
    // 生成树
    std::vector<Tree*> res_tree(node.size(),nullptr);
    for(int i = 0;i<node.size();i++){
        res_tree[i] = new Tree(node[i]);
    }
    for(int i = 0;i<node.size();i++){
        for(int j = 0;j<tree[i].size();j++){
            res_tree[i]->add_child(res_tree[tree[i][j]]);
        }
    }
    assert(res_tree[root]->get_all_nodes().size() == node.size());
    return res_tree[root];
}

/**
 * @brief 
 * 对给定的点集和边界生成一个最小生成树
 * @param add_times 最多增加的边的倍数 
 * @return Tree* 
 */
void AddEdgeByConf(const std::vector<FlipableNode*>& node, std::vector<graph_arg::ConfidenceEdge*>& edges,int add_times = 100){
    // 根据边的置信度,增加边
    REAL min_conf = 1e10, max_conf = 0;
    for(int i = 0;i<edges.size();i++){
        min_conf = (std::min)(min_conf,edges[i]->confidence);
        max_conf = (std::max)(max_conf,edges[i]->confidence);
    }
    assert(add_times > 0);
    double step = (max_conf - min_conf)/add_times;
    int esz = edges.size();
    for(int i = 0;i<esz;i++){        
        edges.push_back(edges[i]);
        if(step == 0)continue;
        int add_num = (edges[i]->confidence - min_conf)/step;
        for(int j = 0;j<add_num;j++){
            edges.push_back(edges[i]);
        }
    }
}

/**
 * @brief 
 * 将图中的每个连通分量生成一个随机生成树; **当且仅当两个顶点之间有双边时,两个顶点才会被连接**
 * @return std::vector<Tree *> 
 */
std::vector<Tree *> graph_arg::FlipableGraph::get_random_forest()
{
    std::vector<Tree*> res_forest;
    using namespace bgl_api;
    UGRAPH g(_nodes.size(),std::vector<int>(_nodes.size(),0));
    for(int i = 0;i<_nodes.size();i++){
        for(int j = 0;j<_nodes.size();j++){
            if(_edges[i][j] != NULL && _edges[j][i] != NULL){
                g[i][j] = 1;
            }
        }
    }
    std::vector<std::vector<int>> cpmponent_idx = bgl_api::bgl_connected_components(g);
    
    for(int i = 0; i< cpmponent_idx.size();i++){
        std::vector<FlipableNode*> tmp_n;
        std::vector<ConfidenceEdge*> tmp_e;
        std::vector<int> use_nodes = cpmponent_idx[i];
        for(int s = 0;s<use_nodes.size();s++){
            tmp_n.push_back(_nodes[use_nodes[s]]);
            for(int t = 0;t<use_nodes.size();t++){
                if(g[use_nodes[s]][use_nodes[t]] == 0)continue;
                tmp_e.push_back(_edges[use_nodes[s]][use_nodes[t]]);
            }
        }
        AddEdgeByConf(tmp_n,tmp_e);
        res_forest.push_back(Kruskal(tmp_n,tmp_e));
    }
    return res_forest;
}

nlohmann::json graph_arg::FlipableGraph::get_metric()
{
    nlohmann::json j;
    j["weight_sum"] = cal_current_weight_sum();
    j["edge_acc"] = cal_edge_acc(_edges);
    j["flip_acc"] = cal_flip_acc(_nodes);
    return j;
}

void graph_arg::FlipableGraph::boolize_weight()
{
    for (int i = 0; i < _edges.size(); i++) {
        for (int j = 0; j < _edges.size(); j++) {
            if (_edges[i][j] == NULL)continue;
            _edges[i][j]->boolize_weight();
        }
    }
}

// 将权重乘以置信度 以削弱坏边
void graph_arg::FlipableGraph::multiply_conf()
{
    for (int i = 0; i < _edges.size(); i++) {
        for (int j = 0; j < _edges.size(); j++) {
            if (_edges[i][j] == NULL)continue;
            _edges[i][j]->multiply_conf();
        }
    }
}

weight_pair graph_arg::FlipableGraph::cal_weight(const std::vector<FlipableNode *> &_node1, const std::vector<FlipableNode *> &_node2) const
{
    weight_pair res{0,0};
    for(int i = 0;i<_node1.size();i++){
        for(int j = 0;j<_node2.size();j++){
            int s = _node1[i]->id, t = _node2[j]->id;
            if(s == t){
                assert(false);
            }
            if(_edges[s][t] == NULL)continue;
            weight_pair tt = _edges[s][t]->get_weight();
            weight_pair tmpw = _edges[s][t]->get_weight(_node1[i]->is_inv() == _node2[j]->is_inv());
            assert(tt == tmpw);
            res.first += tmpw.first;
            res.second += tmpw.second;
        }
    }
    return res;
}

REAL graph_arg::FlipableGraph::cal_current_weight_sum() const 
{
    REAL res = 0;
    for(int i = 0;i<_nodes.size();i++){
        for(int j = 0;j<_nodes.size();j++){
            if(_edges[i][j] == NULL)continue;
            weight_pair tt = _edges[i][j]->get_weight();
            weight_pair tmpw = _edges[i][j]->get_weight(_nodes[i]->is_inv() == _nodes[j]->is_inv());
            assert(tt == tmpw);
            res += tmpw.first;
        }
    }
    return res;
}
graph_arg::align_to_root::align_to_root(FlipableGraph* graph):graph(graph)
{
    
}

std::vector<bool> graph_arg::ForestFlip::flip(FlipableGraph *g)
{
    std::vector<Tree*> trees = g->get_random_forest();
    std::vector<bool> res(g->_nodes.size(),false);
    for(int i = 0;i<trees.size();i++){
        TreeVistor *vistor = new align_to_root(g);
        trees[i]->PostTraverse(vistor);
    }
    
    for(int i = 0;i<g->_nodes.size();i++){
        res[i] = g->_nodes[i]->is_inv();
    }
    return res;
}

nlohmann::json graph_arg::ForestFlip::get_config()
{
    nlohmann::json j;
    j["name"] = "ForestFilp";    
    return j;
}

nlohmann::json graph_arg::ForestFlip::get_log()
{
    return nlohmann::json();
}

weight_pair graph_arg::ConfidenceEdge::get_weight(int if_same) const
{
    if (if_same == -1) {
        if_same = start->is_inv() == end->is_inv();
    }

    if(if_same){
        return {weight,inv_weight};
    }else{
        return {inv_weight,weight};
    }
}

bool graph_arg::ConfidenceEdge::is_good_edge()
{
    bool same = start->is_same() == end->is_same();
    auto w = get_weight();
    return same == (w.first < w.second);
}

void graph_arg::ConfidenceEdge::boolize_weight()
{
    if (weight > inv_weight) {
        weight = 1;
        inv_weight = 0;
    }
    else {
        weight = 0;
        inv_weight = 1;
    }
}

void graph_arg::ConfidenceEdge::multiply_conf(){
    assert(confidence >= 0.5);
    weight *= confidence - 0.5;
    inv_weight *= confidence - 0.5;
}


double graph_arg::cal_flip_acc(std::vector<FlipableNode *> _node)
{
    int cnt = 0;
    for(int i = 0;i<_node.size();i++){
        if(_node[i]->is_same()){
            cnt++;
        }
    }
    cnt = (std::max)(cnt,(int)_node.size() - cnt);
    return (double)cnt/_node.size();

}

double graph_arg::cal_edge_acc(std::vector<std::vector<ConfidenceEdge *>> _edges)
{
    int cnt = 0, all = 0;
    for(int i = 0;i<_edges.size();i++){
        for(int j = 0;j<_edges[i].size();j++){
            if(_edges[i][j] == NULL)continue;
            if(_edges[i][j]->is_good_edge()){
                cnt++;
            }
            all++;
        }
    }
    if (all == 0){
        printf("no edge!!\n");
        return 1;
    }
    return (double)cnt/ (double)all;
}

graph_arg::VoteFlip::VoteFlip(FlipGraph *base_alg,int times):_base_alg(base_alg),_times(times)
{
}



std::vector<bool> graph_arg::VoteFlip::flip(FlipableGraph *g)
{
    std::vector<std::vector<bool>> flip_res;
    for(int i = 0;i<_times;i++){
        auto tmpf = _base_alg->flip(g);
        printf("flip %d, metric: ",i);
        print_metric(g->get_metric());
        if(i!=0)g->_align_flip(flip_res[0],tmpf);
        _log_j.push_back(g->get_metric());
        flip_res.push_back(tmpf);
        g->reset_flip_status();
    }
    std::vector<bool> res(g->_nodes.size(),false);
    
    for(int i = 0;i<g->_nodes.size();i++){
        int cnt = 0;    
        for(int j = 0;j<flip_res.size();j++){
            if(flip_res[j][i]){
                cnt++;
            }
        }
        res[i] = cnt > _times/2;
    }
    g->apply_flip(res);
    printf("vote flip res:");
    print_metric(g->get_metric());
    return res;
}

nlohmann::json graph_arg::VoteFlip::get_config()
{
    nlohmann::json j;
    j["name"] = "VoteFlip";
    j["base_alg"] = _base_alg->get_config();
    j["times"] = _times;
    return j;
}

nlohmann::json graph_arg::VoteFlip::get_log()
{
    std::sort(_log_j.begin(),_log_j.end(),[](nlohmann::json a,nlohmann::json b){
        return a["weight_sum"].get<double>() < b["weight_sum"].get<double>();
    });
    return _log_j;
}


graph_arg::ChoseBestFlip::ChoseBestFlip(FlipGraph *base_alg, int times)
{
    _base_alg = base_alg;
    _times = times;
}

std::vector<bool> graph_arg::ChoseBestFlip::flip(FlipableGraph *g)
{
    std::vector<bool> flip_res;
    double min_w = 1e10;

    for(int i = 0;i<_times;i++){
        auto temp = _base_alg->flip(g);
        printf("flip %d, metric: ",i);
        print_metric(g->get_metric());
        double w = g->cal_current_weight_sum();
        if(w < min_w){
            min_w = w;
            flip_res = temp;
        }
        _log_j.push_back(g->get_metric());
        g->reset_flip_status();
    }
    g->apply_flip(flip_res);
    printf("ChoseBest flip res:");
    print_metric(g->get_metric());
    return flip_res;
}

nlohmann::json graph_arg::ChoseBestFlip::get_config()
{
    nlohmann::json j;
    j["name"] = "ChoseBestFlip";
    j["base_alg"] = _base_alg->get_config();
    j["times"] = _times;
    return j;
}

nlohmann::json graph_arg::ChoseBestFlip::get_log()
{
    // 根据"weight_sum"排序
    std::sort(_log_j.begin(),_log_j.end(),[](nlohmann::json a,nlohmann::json b){
        return a["weight_sum"].get<double>() < b["weight_sum"].get<double>();
    });
    return _log_j;
}

graph_arg::SearchBestFlip::SearchBestFlip(FlipGraph * base_alg, int times)
{
    _base_alg = base_alg;
    _times = times;
}

std::vector<bool> graph_arg::SearchBestFlip::flip(FlipableGraph *g)
{
    std::vector<bool> flip_res;
    double min_w = 1e10;
    for(int i = 0;i<_times;i++){
        auto temp = _base_alg->flip(g);
        printf("flip %d, metric: ",i);
        print_metric(g->get_metric());
        double w = g->cal_current_weight_sum();
        if(w < min_w){
            min_w = w;
            flip_res = temp;
        }
        _log_j.push_back(g->get_metric());
    }
    g->apply_flip(flip_res);
    printf("SearchBest flip res:");
    print_metric(g->get_metric());
    return flip_res;
}

nlohmann::json graph_arg::SearchBestFlip::get_config()
{
    nlohmann::json j;
    j["name"] = "SearchBestFlip";
    j["base_alg"] = _base_alg->get_config();
    j["times"] = _times;
    return j;
}

nlohmann::json graph_arg::SearchBestFlip::get_log()
{
    // 根据"weight_sum"排序
    std::sort(_log_j.begin(),_log_j.end(),[](nlohmann::json a,nlohmann::json b){
        return a["weight_sum"].get<double>() < b["weight_sum"].get<double>();
    });
    return _log_j;
}

std::vector<bool> graph_arg::OptimFlip::flip(FlipableGraph *g)
{
    // 根据grount truth的信息,对图进行翻转
    std::vector<bool> res(g->_nodes.size(),false);
    for(int i = 0;i<g->_nodes.size();i++){
        res[i] = !g->_nodes[i]->same_to_gt;
    }
    g->apply_flip(res);
    return res;
}

nlohmann::json graph_arg::OptimFlip::get_config()
{
    nlohmann::json j;
    j["name"] = "OptimFlip";
    return j;
}
nlohmann::json graph_arg::OptimFlip::get_log()
{
    return nlohmann::json();
}

graph_arg::MIQPFlip::MIQPFlip(nlohmann::json j)
{
    assert(j["name"] == "MIQPFlip");
    _ip = j["ip"].get<std::string>();
    _port = j["port"].get<int>();

}

std::vector<bool> graph_arg::MIQPFlip::flip(FlipableGraph* g)
{
    // 需要保证发送的单位数据是8字节的
    assert(sizeof(double) == 8);
    // 距离矩阵平铺到一维数组中后拼接 得到data
    std::vector<double> w(g->_nodes.size()*g->_nodes.size(),0);
    std::vector<double> inv_w(g->_nodes.size()*g->_nodes.size(),0);
    std::vector<double> data(w.size() + inv_w.size());
    

    g->reset_flip_status();
    for(int i = 0;i<g->_nodes.size();i++){
        for(int j = 0;j<g->_nodes.size();j++){
            if(g->_edges[i][j] == NULL)continue;
            weight_pair tmpw = g->_edges[i][j]->get_weight();
            weight_pair tmpw2 = g->_edges[i][j]->get_weight(g->_nodes[i]->is_inv() == g->_nodes[j]->is_inv());
            assert (tmpw == tmpw2);
            w[i*g->_nodes.size()+j] = tmpw.first;
            inv_w[i*g->_nodes.size()+j] = tmpw.second;
        }
    }

    for(int i = 0;i<w.size();i++){
        data[i] = w[i];
        data[i+w.size()] = inv_w[i];
    }

    lzd_tools::Socket* socket = new lzd_tools::Socket(_ip,_port);
    assert(socket->connectToServer());

    // 发送请求头
    nlohmann::json request;
    request["function_name"] = "MIQP";
    request["data_size"] = data.size();
    if (!socket->SendJson(request.dump())) {
        assert(false);
        std::cerr << "Failed to send request" << std::endl;
        return std::vector<bool>(g->_nodes.size());
    }
    // 接受确认
    std::string response_buff;
    nlohmann::json response;

    if (!socket->ReceiveJson(response_buff)) {
        assert(false);
        std::cerr << "Failed to receive response" << std::endl;
        return std::vector<bool>(g->_nodes.size());
    }
    response = nlohmann::json::parse(response_buff);
    if (response["status"] != "OK") {
        assert(false);
        std::cerr << "Error: " << response["error"] << std::endl;
        return std::vector<bool>(g->_nodes.size());
    }

    // 发送数据
    socket->SendDoubleArray(data);
    std::vector<int32_t> res(g->_nodes.size());
    socket->ReceiveIntArray(res, g->_nodes.size());
    delete socket;
    std::vector<bool> flip_res(g->_nodes.size(),false);
    for(int i = 0;i<res.size();i++){
        flip_res[i] = res[i] == 1;
    }
    g->apply_flip(flip_res);
    printf("MIQP res:");
    print_metric(g->get_metric());
    return flip_res;
}

nlohmann::json graph_arg::MIQPFlip::get_config()
{
    return nlohmann::json();
}

nlohmann::json graph_arg::MIQPFlip::get_log()
{
    return nlohmann::json();
}

FlipGraph* graph_arg::FlipGraph::get_flip_arg_from_json(nlohmann::json config_j)
{
    if(config_j["name"] == "BruteForceFlip"){
        return new BruteForceFlip();
    }
    if(config_j["name"] == "ForestFlip"){
        return new ForestFlip();
    }
    if(config_j["name"] == "VoteFlip"){
        return new VoteFlip(get_flip_arg_from_json(config_j["base_alg"]),config_j["times"].get<int>());
    }
    if(config_j["name"] == "ChoseBestFlip"){
        return new ChoseBestFlip(get_flip_arg_from_json(config_j["base_alg"]),config_j["times"].get<int>());
    }
    if(config_j["name"] == "SearchBestFlip"){
        return new SearchBestFlip(get_flip_arg_from_json(config_j["base_alg"]),config_j["times"].get<int>());
    }
    if(config_j["name"] == "OptimFlip"){
        printf("WARNING::USE OptimFlip, make sure you really want to use it\n");
        return new OptimFlip();
    }
    if(config_j["name"] == "MIQPFlip"){
        return new MIQPFlip(config_j);
    }
}