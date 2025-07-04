#include <kdtree.h>
#include <tools>

template<typename REAL, int DIM>
class BasicKDTree{
    virtual int build(const std::vector<Point<REAL, DIM>>& points) = 0;
    virtual int build(const ORIENTED_POINTS& ops) = 0;
    virtual std::vector<int> knn_search(const Point<REAL,DIM>& q, int k) = 0;
    virtual std::vector<int> radius_search(const Point<REAL, DIM>& query, REAL radius) = 0;
    virtual std::vector<int> radius_search_with_klimit(const Point<REAL, DIM>& query, REAL radius, int k) = 0;
};

template<typename REAL, int DIM>
class ktdKDTree : public BasicKDTree<REAL, DIM>{
    kdt::KDTree<REAL, DIM> _kdt;
    kdt::KDTreePoint tokdtPoint(const Point<REAL,DIM> p){
        array<REAL, 3> _p_{p[0], p[1], p[2]};
        return kdt::KDTreePoint(_p_);
    }
    bool _is_built = false;

public:
    int build(const ORIENTED_POINTS& ops){
        vector<kdt::KDTreePoint> vertices;
        vertices.reserve(ops.size());
        for (size_t i = 0; i < ops.size(); ++i)
        {
            array<REAL, 3> _p_{ops[i].first[0], ops[i].first[1], ops[i].first[2]};
            vertices.push_back(kdt::KDTreePoint(_p_));
        }
        _kdt.build(vertices);
        _is_built = true;
    }

    
    int build(const std::vector<Point<REAL, DIM>>& points){
        vector<kdt::KDTreePoint> vertices;
        vertices.reserve(points.size());
        for (size_t i = 0; i < points.size(); ++i)
        {
            array<REAL, 3> _p_{points[i][0], points[i][1], points[i][2]};
            vertices.push_back(kdt::KDTreePoint(_p_));
        }
        _kdt.build(vertices);
        _is_built = true;
    }


    std::vector<int> knn_search(const Point<REAL,DIM>& q, int k){
        assert (_is_built);
        retutn _kdt.knnSearch(tokdtPoint(q),k);
    }

    std::vector<int> radius_search(const Point<REAL, DIM>& query, REAL radius){
        assert (_is_built);
        return _kdt.radiusSearch(tokdtPoint(query),radius);
    }

    std::vector<int> radius_search_with_klimit(const Point<REAL, DIM>& query, REAL radius, int k){
        assert (_is_built);
        std::vector<int> res =  _kdt.radiusSearch(tokdtPoint(query),radius);
        if(res.size() > k){
            res.resize(k);
        }
        return res;
    }



}