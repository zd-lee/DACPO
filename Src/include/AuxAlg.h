#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <math.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/io/ply_io.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/normal_3d.h>
#include <pcl/search/kdtree.h>
#include <IPSR.h>

/**
*
* @brief
* 
*/


namespace aux_arg {
	template <typename REAL, unsigned int DIM>
	pcl::PointCloud<pcl::PointXYZ>::Ptr pn2pclxyz(const POINTS_NORMALS& pn);

	template<typename REAL, unsigned int DIM>
	class LabelGetter {

		pcl::search::KdTree<pcl::PointXYZ>::Ptr tree;
        std::vector<int> _label;
	
	public:
		LabelGetter(ORIENTED_POINTS op,std::vector<int> label){
			
			
			label = lzd_tools::squeeze_array(label,0);
			assert(label.size() == op.size());
			int max_label = *std::max_element(label.begin(), label.end());
			printf("label num = %d\n", max_label + 1);

			// 	_op_spilt_by_label[label[i]].push_back(op[i]);
			// 	_ops.insert(_ops.end(), _op_spilt_by_dbscan.begin(), _op_spilt_by_dbscan.end());
			// 		_ops.push_back(others);
			// label.clear();
			// pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
			// 		label.push_back(i);
			// 		pcl::PointXYZ p;
			// 			p.data[k] = _ops[i][j].first[k];
			// 		cloud->push_back(p);
			pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = pn2pclxyz<REAL, DIM>(op);
			_label = label;
			tree = pcl::search::KdTree<pcl::PointXYZ>::Ptr(new pcl::search::KdTree<pcl::PointXYZ>);
			tree->setInputCloud(cloud);
		}



		int get_label(Point<REAL, DIM> p){
			pcl::PointXYZ searchPoint;
			for(int i = 0; i < DIM; i++){
				searchPoint.data[i] = p[i];
			}
			std::vector<int> pointIdxNKNSearch(1);
			std::vector<float> pointNKNSquaredDistance(1);
			if(tree->nearestKSearch(searchPoint, 1, pointIdxNKNSearch, pointNKNSquaredDistance)){
				return _label[pointIdxNKNSearch[0]];
			}
			assert(0);
			return -1;
		}
        
    };

	template <typename REAL>
	using VV = std::vector<std::vector<REAL>>;

	template <typename REAL>
	std::vector<double> _boost_cal_node_betweeness(VV<REAL> graph) {
		return std::vector<double>(graph.size(), 0.0);
	}

	template <typename REAL>
	VV<double> _boost_cal_edge_betweenness(VV<REAL> graph) {
		VV<double> res(graph.size(), std::vector<double>(graph.size(), -1.0));
		for (int i = 0; i < graph.size(); i++) {
			for (int j = 0; j < graph[i].size(); j++) {
				if (graph[i][j] >= 0) {
					res[i][j] = 0.0;
				}
			}
		}
		return res;
	}

	template <typename REAL, unsigned int DIM>
	pcl::PointCloud<pcl::PointXYZ>::Ptr pn2pclxyz(const POINTS_NORMALS& pn) {
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->resize(pn.size());
		for (int i = 0; i < pn.size(); i++) {
			for (int j = 0; j < 3; j++) {
				cloud->points[i].data[j] = pn[i].first[j];
			}
		}
		return cloud;
	}

	/**
	 * @brief 
	 */
	template <typename REAL, unsigned int DIM>
	VV<REAL> pointcloud2graph(const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
	REAL radius = 0.01,int K = 10) {
		VV<REAL> graph;
		graph.resize(cloud->size());
		pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
		tree->setInputCloud(cloud);
		std::vector<int> pointIdxNKNSearch(K);
		std::vector<float> pointNKNSquaredDistance(K);
		for (int i = 0; i < cloud->size(); i++) {
			if (tree->nearestKSearch(cloud->points[i], K, pointIdxNKNSearch, pointNKNSquaredDistance)) {
				for (int j = 0; j < pointIdxNKNSearch.size(); j++) {
					if (pointIdxNKNSearch[j] != i && pointNKNSquaredDistance[j]<=radius) {
						graph[i].push_back(pointIdxNKNSearch[j]);
					}
				}
			}
		}
		return graph;
	}

	template <typename REAL, unsigned int DIM>
	int __BFSgrowing(const VV<REAL>& graph, std::vector<int>& arrived, float rate = 0.1) {
		assert(graph.size() > 0);
		assert(graph.size() == arrived.size());
		int arrived_num = 0, max_arrived_num = 0;
		std::queue<int> q;
		for(int i = 0; i < arrived.size(); i++){
			if(arrived[i] == 1){
				arrived_num++;
				q.push(i);
			}
			if (arrived[i] != -1) {
				max_arrived_num++;
			}
		}
		max_arrived_num = max_arrived_num * rate;
		

		while(!q.empty() && arrived_num < max_arrived_num){
			int i = q.front();
			q.pop();
			for(int j = 0; j < graph[i].size(); j++){
				if(graph[i][j] > 0 && arrived[graph[i][j]] == 0){
					arrived_num++;
					arrived[graph[i][j]] = 1;
					q.push(graph[i][j]);
					break;
				}
			}
			
		}
		if(arrived_num < max_arrived_num){
			printf("BFSgrowing: arrived_num < max_arrived_num\n");
			return -1;
		}
		return 0;
	}

	/**
	 * @brief 
	 * @param graph 
	 */
	template <typename REAL, unsigned int DIM>
	int BFS_growing(const POINTS_NORMALS& pn, std::vector<int>& arrived, float rate = 0.1) {
		VV<REAL> graph = pointcloud2graph<REAL, DIM>(pn2pclxyz<REAL, DIM>(pn));
		return __BFSgrowing<REAL, DIM>(graph, arrived, rate);
	}

	template<typename REAL>
	std::vector<int> _dijkstra_shortest_path(VV<REAL> graph, int s, int e){
		std::vector<int> path;
		std::vector<REAL> dist(graph.size(), -1);
		std::vector<int> prev(graph.size(), -1);
		std::vector<int> visited(graph.size(), 0);
		dist[s] = 0;
		for(int i = 0; i < graph.size(); i++){
			REAL min_dist = -1;
			int min_idx = -1;
			for(int j = 0; j < graph.size(); j++){
				if(visited[j] == 0 && dist[j] != -1 && (min_dist == -1 || dist[j] < min_dist)){
					min_dist = dist[j];
					min_idx = j;
				}
			}
			if(min_idx == -1){
				break;
			}
			visited[min_idx] = 1;
			for(int j = 0; j < graph[min_idx].size(); j++){
				if(graph[min_idx][j] != -1 && (dist[j] == -1 || dist[j] > dist[min_idx] + graph[min_idx][j])){
					dist[j] = dist[min_idx] + graph[min_idx][j];
					prev[j] = min_idx;
				}
			}
		}
		if(dist[e] == -1){
			return path;
		}
		int cur = e;
		while(cur != s){
			path.push_back(cur);
			cur = prev[cur];
		}
		path.push_back(s);
		std::reverse(path.begin(), path.end());
		return path;
	}	

};
