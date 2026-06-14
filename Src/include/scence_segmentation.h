#include <tools.h>
#include <socket_api.h>
#include <kdtree_api.h>

/**
 * @brief 
 * */
template<typename REAL, unsigned int DIM>
class Segmenter{
public:
    virtual nlohmann::json get_config() = 0;
    virtual std::vector<int> segment(const ORIENTED_POINTS& ops) = 0;
    
    virtual std::vector<int> nest_segment(const ORIENTED_POINTS& ops){
        return relabel(segment(ops));
    }

    std::vector<int> relabel(std::vector<int> labels){
        std::map<int, int> label_map;
        int label = 0;
        std::vector<int> res(labels.size());
        for(int i = 0;i<labels.size();i++){
            if(label_map.find(labels[i]) == label_map.end()){
                label_map[labels[i]] = label++;
            }
            res[i] = label_map[labels[i]];
        }
        return res;
    }
};

/**
 * @brief 
 */
class CombSegmenter : public Segmenter{
    std::vector<Segmenter*> _segmenters;
    double _mininum_rate;

public:
    CombSegmenter(std::vector<Segmenter*> segmenters, double mininum_rate){
        _segmenters = segmenters;
        _mininum_rate = mininum_rate;
    }
    nlohmann::json get_config(){
        nlohmann::json j;
        j["name"] = "CombSegmenter";
        j["mininum_rate"] = _mininum_rate;
        for(int i = 0;i<_segmenters.size();i++){
            j["segmenter" + std::to_string(i)] = _segmenters[i]->get_config();
        }
        return j;
    }

    std::vector<int> segment(const ORIENTED_POINTS& ops){
        std::vector<int> res(ops.size());
        for(auto segmenter:_segmenters){
            std::vector<int> L1_labels = segmenter->nest_segment(ops);
            int max_label = *std::max_element(L1_labels.begin(), L1_labels.end());
            std::vector<ORIENTED_POINTS> opss(max_label + 1);
            for(int i = 0;i<ops.size();i++){
                opss[L1_labels[i]].push_back(ops[i]);
            }
            for(int i = 0;i<=max_label;i++){
                if(opss[i].size() < _mininum_rate * ops.size()){
                    continue;
                }
                std::vector<int> L2_labels = segmenter->nest_segment(opss[i]);
                for(int j = 0;j<L2_labels.size();j++){
                    res[j] = i * (max_label + 1) + L2_labels[j];
                }
            }
        }
        return res;
    }

};


template<typename REAL, unsigned int DIM>
class DisjointSetCluster : public Segmenter<REAL, DIM>{
    double _eps;
    int _max_neighbers;
    
public:
    DisjointSetCluster(nlohmann::json j){
        _eps = j["eps"];
        _max_neighbers = j["max_neighbers"];
    }

    nlohmann::json get_config(){
        nlohmann::json j;
        j["name"] = "DisjointSetCluster";
        j["eps"] = _eps;
        j["max_neighbers"] = _max_neighbers;
        return j;
    }

    std::vector<int> segment(const ORIENTED_POINTS& ops){
        BasicKDTree<REAL, DIM>* kdt = new ktdKDTree<REAL, DIM>();
        kdt->build(opoints);
        std::vector<int> res;
        std::vector<int> roots(ops.size());

        auto find = [&roots](int x){
            if(roots[x] == x){
                return x;
            }else{
                return roots[x] = find(roots[x]);
            }
        };

        auto union_set = [&roots](int x, int y){
            int rx = find(x);
            int ry = find(y);
            if(rx != ry){
                roots[rx] = ry;
            }
        };

        for(int i = 0;i<ops.size();i++){
            roots[i] = i;
        }

        for(int i = 0;i<ops.size();i++){
            std::vector<int> neighbers = kdt->radius_search_with_klimit(ops[i].first, _eps, _max_neighbers);
            if(find(i) != i){
                continue;
            }
            for(auto j:neighbers){
                union_set(i, j);
            }
        }
        
        return roots;
    }
};



template<typename REAL, unsigned int DIM>
class SpectralSegmenter : public Segmenter<REAL, DIM>{
    float _mininum_rate;
    int _k_neighbers;

    // socket args
    std::string _ip_addr;
    int _port;

public:
    SpectralSegmenter(nlohmann::json j){
        _mininum_rate = j["mininum_rate"];
        _k_neighbers = j["k_neighbers"];
    }

    nlohmann::json get_config(){
        nlohmann::json j;
        j["name"] = "SpectralSegmenter";
        j["mininum_rate"] = _mininum_rate;
        j["k_neighbers"] = _k_neighbers;
        return j;
    }

    std::vector<int> segment(const ORIENTED_POINTS opoints){
        std::vector<Point<REAL, DIM>> points;
        for(auto op:opoints){
            points.push_back(op.first);
        }
        lzd_tools::Socket socket(_ip_addr, _port);
        assert(socket.connectToServer());
        assert(socket.SendJson(get_config().dump()));
        std::string resp;
        assert(socket.ReceiveJson(resp));
        nlohmann::json response = nlohmann::json::parse(resp);

        std::vector<int> res;
        std::vector<REAL> data;
        for(auto p:points){
            for(int i = 0;i<DIM;i++){
                data.push_back(p[i]);
            }
        }
        assert(socket.SendDoubleArray(data));
        assert(socket.ReceiveIntArray(res));
        return res;
    }



};

