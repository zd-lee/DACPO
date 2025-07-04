// FieldPropagation.cpp
#include <limits>
#include <random>
#include <cmath>
#include <iostream>
#include <queue>
#include <algorithm>
#include <set>
#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/io/ply_io.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/normal_3d.h>
#include <pcl/search/kdtree.h>
#include "dipole.h"
#include <graph_arg.h>

namespace FIELD_PROPAGATION {

// Field implementation
Eigen::MatrixXd MaskField::mask_xie_intersection(const Eigen::MatrixXd& source,
                                        const Eigen::MatrixXd& target,
                                        const Eigen::MatrixXi& knn_mask) {
    int T = target.rows();
    int K = knn_mask.cols();
    Eigen::MatrixXd intersection(T, K);
    for (int i = 0; i < T; ++i) {
		auto ref_normal = mask_xie_field(source, target.row(i), knn_mask.row(i)); // K*3
		for (int j = 0; j < K; ++j) {
			intersection(i, j) = ref_normal.row(j).dot(target.row(i).tail<3>());
		}
    }
    return intersection;
}

Eigen::MatrixXd MaskField::mask_xie_field(const Eigen::MatrixXd& source,
                                          const Eigen::RowVectorXd& target,
                                          const Eigen::VectorXi& knn_mask) {
    int K = knn_mask.size();
    Eigen::MatrixXd result(K, 3);
    for (int j = 0; j < K; ++j) {
        Eigen::RowVectorXd diff = source.row(knn_mask[j]).head<3>() - target.head<3>();
        double norm = diff.norm();
        assert(norm != 0);
        Eigen::RowVectorXd unit_vec = diff / norm;
        Eigen::RowVectorXd normal_s = source.row(knn_mask[j]).tail<3>();
        Eigen::RowVectorXd ref_normal = normal_s - C * (normal_s.dot(unit_vec)) * unit_vec;
        ref_normal /= std::pow(norm, 3);
        result.row(j) = ref_normal;
    }
    return result;
}

// KnnGraph implementation
KnnGraph::KnnGraph(const Eigen::MatrixXd& pts, int k, double radius)
    : _pts(pts), _k(k), _radius(radius) {
    _k = std::min(k,(int)pts.rows());
    assert(_k > 1);
    int N = pts.rows();
    _graph_idx = std::vector<std::vector<int>>(N,std::vector<int>());
    _graph_weight = std::vector<std::vector<double>>(N,std::vector<double>());
    // build kd-tree
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    for (int i = 0; i < N; i++) {
        pcl::PointXYZ p;
        p.x = pts(i, 0);
        p.y = pts(i, 1);
        p.z = pts(i, 2);
        cloud->push_back(p);
    }
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
    tree->setInputCloud(cloud);
    // search knn
    std::vector<int> pointIdxNKNSearch(_k);
    std::vector<float> pointNKNSquaredDistance(_k);
    for (int i = 0; i < N; i++) {
        if (tree->nearestKSearch(cloud->points[i], _k, pointIdxNKNSearch, pointNKNSquaredDistance)) {
            for (int j = 1; j < _k; j++) // exclude itself
            {
                if (pointNKNSquaredDistance[j] <= _radius) {
                    _graph_idx[i].push_back(pointIdxNKNSearch[j]);
                    _graph_weight[i].push_back(pointNKNSquaredDistance[j]);
                }
            }
        }
    }
}

KnnGraph::KnnGraph(const Eigen::MatrixXd& pts, int k) : _pts(pts), _k(k)
{
    // TODO 可以优化一下
    _k = std::min(k, (int)pts.rows());
    assert(_k > 1); // 不会有大小为1的块吧...
    _radius = std::numeric_limits<double>::max();
    int N = pts.rows();
    _graph_idx = std::vector<std::vector<int>>(N, std::vector<int>(_k-1));
    // build kd-tree
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    for (int i = 0; i < N; i++) {
        pcl::PointXYZ p;
        p.x = pts(i, 0);
        p.y = pts(i, 1);
        p.z = pts(i, 2);
        cloud->push_back(p);
    }
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
    tree->setInputCloud(cloud);
    // search knn
    std::vector<int> pointIdxNKNSearch(_k);
    std::vector<float> pointNKNSquaredDistance(_k);
    for (int i = 0; i < N; i++) {
        if (tree->nearestKSearch(cloud->points[i], _k, pointIdxNKNSearch, pointNKNSquaredDistance)) {
            for (int j = 1; j < _k; j++) // exclude itself
            {
                _graph_idx[i][j-1] = pointIdxNKNSearch[j];
            }
        }
    }
}


std::vector<int> KnnGraph::get_BFS_route(int startpoint) {
    std::vector<bool> visited(_pts.rows(), false);
    std::vector<int> route;
    std::queue<int> q;
    q.push(startpoint);
    visited[startpoint] = true;
    while (!q.empty()) {
        int current = q.front();
        q.pop();
        route.push_back(current);
        for (int neighbor : _graph_idx[current]) {
            if (!visited[neighbor]) {
                q.push(neighbor);
                visited[neighbor] = true;
            }
        }
        if (q.size() == 0) {
            for (int i = 0; i < _pts.rows(); i++) {
                if (!visited[i]) {
                    q.push(i);
                    visited[i] = true;
                    break;
                }
            }
        }
    }
	assert(route.size() == _pts.rows());
    return route;
}

int KnnGraph::get_shortest_nei()
{
    int k = 1e9;
    for (int i = 0; i < _pts.rows(); i++)
    {
        if (_graph_idx[i].size() < k)
        {
            k = _graph_idx[i].size();
        }
    }
    assert(k != 1e9);
    return k;
}

std::vector<bool> MaskFieldPropagation::_xie_propagation_points_in_order(const Eigen::MatrixXd& pts, const std::vector<int>& order, 
    const Eigen::MatrixXd& pre_caled_field_intersaction_mat, const Eigen::MatrixXi& KnnMask) {
    assert(KnnMask.cols() == pre_caled_field_intersaction_mat.cols());
    assert(order.size() == pts.rows());
    assert(order.size() == pre_caled_field_intersaction_mat.rows());
    int N = pts.rows();
    std::vector<float> weight(N,0); // 0 for not visited, -1 for flip, 1 for not flip
    std::vector<bool> res(N);
    for (int i = 0; i < N; ++i) {
        int idx = order[i];
        double force_sum = 0;
        for (int j = 0 ;j < KnnMask.cols(); ++j) {
            force_sum += pre_caled_field_intersaction_mat(idx, j) * weight[KnnMask(idx, j)];
        }
        weight[idx] = force_sum > 0 ? 1 : -1;
    }

    if(_diffuse){
        for (int i = 0; i < N; ++i) {
            double force_sum = 0;
            for (int j = 0; j < KnnMask.cols(); ++j) {
                force_sum += pre_caled_field_intersaction_mat(i, j) * weight[KnnMask(i, j)];
            }
            weight[i] = force_sum > 0 ? 1 : -1;
        }
    }


    for (int i = 0; i < N; ++i) {
        res[i] = weight[i] > 0 ? 1 : 0;
    }
    return res;
}

static std::vector<int> get_random(int num,int min, int max){
    // 随机选择times个点作为起始点(使用固定的随机种子)
    std::set<int> start_points;
    std::mt19937 gen(0);
    std::uniform_int_distribution<int> dis(min, max);
    start_points.insert(0);
    while (start_points.size() < num){
        start_points.insert(dis(gen));
    }
    return std::vector<int>(start_points.begin(), start_points.end());
}
static void align_flip_status(std::vector<std::vector<bool>>& flip_status){
    int T = flip_status.size();
    int N = flip_status[0].size();
    std::vector<graph_arg::FlipableNode*> nodes(T);
    std::vector<std::vector<graph_arg::ConfidenceEdge*>> edges(T, std::vector<graph_arg::ConfidenceEdge*>(T, nullptr));
    for (int i = 0; i < T; i++) {
        nodes[i] = new graph_arg::FlipableNode(i);
    }
    
    for (int i = 0; i < T; i++) {
        for(int j = 0; j < T; j++){
            int distance_ij = 0;
            if(i == j){
                edges[i][j] = nullptr;
                continue;
            }
            if(j < i){
                distance_ij = edges[j][i]->weight;
            }
            else {
                for (int k = 0; k < N; k++) {
                    if (flip_status[i][k] != flip_status[j][k]) {
                        distance_ij++;
                    }
                }
            }
            edges[i][j] = new graph_arg::ConfidenceEdge(nodes[i], nodes[j], distance_ij, N - distance_ij);    
        }
    }

    //graph_arg::MIQPFlip bf("192.168.192.115", 11111);
    graph_arg::BruteForceFlip bf;
    graph_arg::FlipableGraph* tgraph = new graph_arg::FlipableGraph(nodes, edges);
    std::vector<bool> flip_res = bf.flip(tgraph);
    for(int i = 0; i < flip_res.size(); i++){
        if(flip_res[i]){
            for(int j = 0; j < N; j++){
                flip_status[i][j] = !flip_status[i][j];
            }
        }
    }
    return;
}

std::vector<bool> MaskFieldPropagation::xie_propagation_points_on_bfs_tree(const Eigen::MatrixXd &pts)
{
    int N = pts.rows();
    std::vector<bool> res(N, false);
    // 建立KnnGraph
    KnnGraph bfs_graph(pts, knngraph_k, knngraph_radius);
    KnnGraph neibor_graph(pts, mask_size);
    int K = neibor_graph.get_shortest_nei();
    Eigen::MatrixXi knn_mask(N, K);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < K; j++) {
			knn_mask(i, j) = neibor_graph._graph_idx[i][j];
		}
	}
    Eigen::MatrixXd pre_caled_field_intersaction_mat = _field.mask_xie_intersection(pts, pts, knn_mask);
    std::vector<std::vector<bool>> flip_status(_times, std::vector<bool>(N, 0));
    std::vector<std::vector<int>> orders(_times, std::vector<int>(N, 0));
    std::vector<int> start_points = get_random(_times, 0, N - 1);        

#pragma omp parallel for
    for (int i = 0; i < _times; ++i) {
        orders[i] = bfs_graph.get_BFS_route(start_points[i]);
    }

#pragma omp parallel for
    for (int i = 0; i < _times; ++i) {
        flip_status[i] = _xie_propagation_points_in_order(pts, orders[i], pre_caled_field_intersaction_mat, knn_mask);
    }

    align_flip_status(flip_status);

#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        int cnt = 0;
        for (int j = 0; j < _times; ++j) {
            if (flip_status[j][i]) {
                cnt++;
            }
        }
        res[i] = cnt > _times / 2;
    }
    return res;
}
} // namespace FieldPropagation
