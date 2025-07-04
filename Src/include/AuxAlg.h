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
//#include <boost/graph/graph_traits.hpp>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/betweenness_centrality.hpp>
#include <IPSR.h>

/**
*
* @brief
* 辅助算法
* 
* * BFS_growing 扩大某个区域
* * pn2pclpc 将POINTS_NORMALS类型转化为pcl::PointCloud<pcl::PointXYZ>::Ptr
*/


namespace aux_arg {
	// 存储一个kdtree和_label。输入一个点是，返回其在空间中的领域内的点的label
	template<typename REAL, int DIM>
	class LabelGetter {

        // 存储一个kd树
		pcl::search::KdTree<pcl::PointXYZ>::Ptr tree;
        std::vector<int> _label;
	
	public:
		// @file_name 不带有后缀的文件名
		LabelGetter(ORIENTED_POINTS op,std::vector<int> label){
			
			
			label = lzd_tools::squeeze_array(label,0);
			assert(label.size() == op.size());
			int max_label = *std::max_element(label.begin(), label.end());
			printf("label num = %d\n", max_label + 1);

			// std::vector<ORIENTED_POINTS> _ops;
			// std::vector<ORIENTED_POINTS> _op_spilt_by_label(max_label + 1);
			// for(int i = 0; i < label.size(); i++){
			// 	_op_spilt_by_label[label[i]].push_back(op[i]);
			// }

			// // 对每个label的op进行dbscan再分类
			// for(int i = 0; i < _op_spilt_by_label.size(); i++){
			// 	if(_op_spilt_by_label[i].size() == 0){
			// 		continue;
			// 	}
			// 	std::vector<ORIENTED_POINTS> _op_spilt_by_dbscan;
			// 	ORIENTED_POINTS others;
			// 	o3d_dbscan_extraction(_op_spilt_by_label[i], _op_spilt_by_dbscan,others, eps, min_points);
			// 	_ops.insert(_ops.end(), _op_spilt_by_dbscan.begin(), _op_spilt_by_dbscan.end());
			// 	if (others.size() > 0) {
			// 		_ops.push_back(others);
			// 	}
			// }

			// label.clear();
			// pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
			// for(int i = 0; i < _ops.size(); i++){
			// 	for(int j = 0; j < _ops[i].size(); j++){
			// 		label.push_back(i);
			// 		pcl::PointXYZ p;
			// 		for(int k = 0; k < DIM; k++){
			// 			p.data[k] = _ops[i][j].first[k];
			// 		}
			// 		cloud->push_back(p);
			// 	}
			// }
			// _label = label;
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

	// 将POINT_NORMALS复制到一个pcl::PointCloud<pcl::PointXYZ>::Ptr中
	template <typename REAL, int DIM>
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
	 * 将点云转化为邻接矩阵。将每个点K近邻的点的距离作为边的权重。
	 * @param cloud pcl点云 
	 * @param graph 邻接矩阵
	 * @param radius 限制最大可连接距离
	 * @param K 限制最大可连接点数
	 */
	template <typename REAL, int DIM>
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

	template <typename REAL, int DIM>
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
			// 排除-1的情况
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
	 * BFS增长，从arrived中的点开始，
	 * 使用BFS算法，将所有与arrived中的点相连的点加入到arrived中,直到arrived中的点数达到rate*graph.size()
	 * @param graph 
	 * @param arrived 0表示未到达，1表示已到达，-1表示不可达
	 * @param rate 退出条件，已到达/（总点数-不可达点数）>= rate
	 */
	template <typename REAL, int DIM>
	int BFS_growing(const POINTS_NORMALS& pn, std::vector<int>& arrived, float rate = 0.1) {
		VV<REAL> graph = pointcloud2graph<REAL, DIM>(pn2pclxyz<REAL, DIM>(pn));
		return __BFSgrowing<REAL, DIM>(graph, arrived, rate);
	}

	// 找到图中，s到e的最短路径，并返回途径的点
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

	/**
	 * @brief 
	 * 计算图中每个边的betweenness
	 * 图不能太大；图中不可以出现负权边
	 * @param graph graph[i][j] = x (x>=0)表示边i->j的权重；x = -1表示这个边不存在
	 * @return graph[i][j] = x (x>=0)表示边i->j的betweenness=x；x = -1表示这个边不存在 
	 */
	// template <typename REAL>
	// VV<REAL> cal_betweenness(VV<REAL> graph){
	// 	VV<REAL> betweenness(graph.size(), std::vector<REAL>(graph.size(), 0));
	// 	for(int i = 0; i < graph.size(); i++){
	// 		for(int j = 0; j < graph.size(); j++){
	// 			if(graph[i][j] != -1){
	// 				std::vector<int> path = _dijkstra_shortest_path<REAL>(graph, i, j);
	// 				for(int k = 0; k < path.size() - 1; k++){
	// 					betweenness[path[k]][path[k+1]] += 1;
	// 				}
	// 			}
	// 		}
	// 	}
	// 	return betweenness;
	// }

	// 使用boost的brandes_betweeness_centrality算法计算node_betweenness
	//template <typename REAL>
	//std::vector<double> _boost_cal_node_betweeness(VV<REAL> graph){
	//	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> Graph;
	//	const int _num_vertices = graph.size();
	//	typedef std::pair<int, int> Edge;
	//	int num_edges = 0;
	//	// Create a graph object
	//	Graph g(_num_vertices);
	//	// Add edges to the graph object
	//	for (int i = 0; i < _num_vertices; i++) {
	//		for (int j = 0; j < _num_vertices; j++) {
	//			if (graph[i][j] != -1) {
	//				boost::add_edge(i, j, g);
	//				num_edges++;
	//			}
	//		}
	//	}

	//	boost::shared_array_property_map<double, boost::property_map<Graph, boost::vertex_index_t>::const_type>
 // 			centrality_map(num_vertices(g), get(boost::vertex_index, g));
	//	boost::brandes_betweenness_centrality(g, centrality_map);

	//	std::vector<double> res;
	//	boost::graph_traits<Graph>::vertex_iterator vi, vi_end;
	//	for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
	//		res.push_back(centrality_map[*vi]);
	//	}
	//	return res;
	//}

	// 使用boost的brandes_betweenness_centrality算法计算edge_betweenness
	// edge_betweeness = 
	//template <typename REAL>
	//VV<double> _boost_cal_edge_betweenness(VV<REAL> graph){
	//	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> Graph;
	//	const int _num_vertices = graph.size();
	//	typedef std::pair<int, int> Edge;
	//	int num_edges = 0;
	//	// Create a graph object
	//	Graph g(_num_vertices);
	//	// Add edges to the graph object
	//	for (int i = 0; i < _num_vertices; i++) {
	//		for (int j = 0; j < _num_vertices; j++) {
	//			if (graph[i][j] != -1) {
	//				boost::add_edge(i, j, g);
	//				num_edges++;
	//			}
	//		}
	//	}
	//	
	//	// typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
	//	// typedef boost::property_map<Graph,boost::vertex_index_t>::type IndexMap;
	//	// IndexMap index = boost::get(boost::vertex_index, g);
	//	// std::cout << "vertices(g) = ";
	//	// typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
	//	// std::pair<vertex_iter, vertex_iter> vp;
	//	// for (vp = vertices(g); vp.first != vp.second; ++vp.first)
	//	// 	std::cout << index[*vp.first] <<  " ";
	//	// std::cout << std::endl;
	//	// std::cout << "edges(g) = ";
 //  		// boost::graph_traits<Graph>::edge_iterator ei, ei_end;
 //   	// for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
 //       // 	std::cout << "(" << index[source(*ei, g)]<< "," << index[target(*ei, g)] << ") ";
 //   	// std::cout << std::endl;
	//	// compute the edge betweenness centrality
	//	// boost::shared_array_property_map<double, 
	//	
	//	using ECMap = std::map<Graph::edge_descriptor, double>;
	//	using ECEntry = ECMap::value_type;
	//	ECMap ecm;
	//	boost::brandes_betweenness_centrality(
	//		g, boost::edge_centrality_map(boost::make_assoc_property_map(ecm)));

	//	VV<double> res(graph.size());
	//	for(int i = 0; i < graph.size(); i++){
	//		res[i].resize(graph.size(), -1);
	//	}
	//	for (ECMap::iterator it = ecm.begin(); it != ecm.end(); ++it) {
	//		int s = boost::source(it->first, g);
	//		int t = boost::target(it->first, g);
	//		res[s][t] = it->second;
	//		assert(res[s][t] > 0);
	//	}
	//	return res;
	//}

};
