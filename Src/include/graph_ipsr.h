#pragma once
#define _UNSERTAIN //对于一些不确定的地方，使用assert进行检查
#include <IPSR.h>
#include <ipsr_controller.h>
#include "3dparty\nlohmann\json.hpp"
#include <AuxAlg.h>
#include <dipole.h>


extern nlohmann::json global_var::config_j;// 控制参数
extern nlohmann::json global_var::metric_j;// 效果、日志等

namespace GRAPH_IPSR {

/*********************************************GRAPH**********************************************************/
	// ipsr_handle 立方体
	template<typename T>
	using T_CUBE = std::vector<std::vector<std::vector<T>>>;

	template<typename REAL, int DIM>
	using IPSR_CUBE = T_CUBE<IPSR_HANDLE_P>;

	template<typename REAL, int DIM>
	using OP_CUBE = T_CUBE<POINTS_NORMALS>;

	template<typename T>
	void cube_resize(T_CUBE<T>& cube, int x_num, int y_num, int z_num){
		cube.resize(x_num);
		for (int i = 0; i < x_num;i++) {
			cube[i].resize(y_num);
			for (int j = 0; j < y_num; j++) {
				cube[i][j].resize(z_num);
			}
		}
	}
	
	// 顶点。每个顶点指向一个ipsr对象。
	// 后续可能会加入置信度等属性
	template<typename REAL, int DIM>
	class graph_vertex {
	private:
		std::vector<std::vector<float>> _diff_log;// _op中每个在_points_normals中的点的diff的log
		std::vector<int> _op_idx;//_op中的点在_handle->_points_normals中的索引
		bool init_flag = false;
		int inv_times = 0;//翻转的次数

		std::pair<REAL, REAL> _loss_according2gt;// first为inv_times == 0时的loss; second为inv_times == 1 时的loss

		lzd_tools::thread_safe_bool __is_changed;// 在_cal_weigth后是否再发生过改变

		void init() {
			assert(init_flag == false);
			init_flag = true;
			assert(_handle != NULL);
			// _handle->ixform 应该是一个只有平移或者放缩的变换 不应该有旋转
			std::vector<std::pair<int, int>> checkidx = {
				{0,1},{0,2},{1,2}
			};
			for(auto idx:checkidx){
				assert(_handle->ixform(idx.first,idx.second) == 0);
				assert(_handle->ixform(idx.second,idx.first) == 0);
			}

			_handle->update_target(_op, _handle->ixform.inverse(), _op_idx);
			_diff_log.resize(_op.size(), std::vector<float>());
			__is_changed.set(true);
		}

		// 求单个点的未归一化的置信度，方式为取最后五次var的总和sum以及均值avg，若sum不为零，则置信度为1/avg
		// 对于var很小的点，其置信度为1/0.01
		float _cal_confidence_by_var(std::vector<float> var){
			//assert(var.size());
			float sum = 0.0f;
			for(int i = var.size() - 1;i>=0&&i>= var.size()-5;i--){
				sum += var[i];
			}
			//if(sum == 0) assert(0);
			float avg = sum / std::min(5,(int)var.size());
			avg = max(1e-9f,avg);//防止除0,或者avg过小，导致某个点的置信度过大 FIX_ME diff的大小与点云大小有关
			avg = 1.0 / avg;
			// 取对数
			avg = log(avg);
			avg = max(0.0f,avg);
			return avg;
		}

		friend void comb_vertex(graph_vertex<REAL, DIM>* v1, graph_vertex<REAL, DIM>* v2, graph_vertex<REAL, DIM>* res);
		
	public:
		
		//enum GraphVertexOPType {
		//	SELF,//自身部分
		//	// 此处不需要区分与不同vertex的重合部分 
		//	OVERLAP //与其他Vertex重合部分。
		//};
		IPSR_HANDLE_P _handle;
		std::vector<int> _op_type;//_op中每个点的类型 0表示来自grid，1表示来自邻居grid
		POINTS_NORMALS _op;// 顶点对应的有向点云。注意这个点云与_handle->_points_normals中的点不一致，后者经过变换以及下采样
		std::vector<float> _op_confidence;//_op中的点的置信度
		Point<REAL, DIM> center_op;

		// 节点相较于最开始时的翻转状态
		bool is_inved(){
			return inv_times % 2;
		}	

		// 与gt相比的翻转状态,返回true表示与gt相同，1表示与gt相反
		int inved_status_according2gt(){
			auto loss_pair = _cal_loss();
			return loss_pair.first <= loss_pair.second;
		}

		int iter(IpsrController<REAL, DIM>* controller = NULL) {
			if (!init_flag) {
				init();
			}
			int flag;
			if(controller != NULL){
				flag = controller->iter(_handle);
			}else{
				flag = this->_handle->iter();
			}

			if(flag<0)return flag;
			
			for(int i = 0;i<_op.size();i++){
				_diff_log[i].push_back(_handle->_points_normals_diff[_op_idx[i]]);
			}
			return flag;
		}

		// 计算碎片中,每个op的置信度
		void cal_conf() {
			float max_confidence = 0;
			_op_confidence.clear();
			for (int i = 0; i < _op.size(); i++) {
				_op_confidence.push_back(_cal_confidence_by_var(_diff_log[i]));
				if (_op_confidence[i] > max_confidence) max_confidence = _op_confidence[i];
			}
			// 归一化置信度
			for (int i = 0; i < _op.size(); i++) {
				_op_confidence[i] /= max_confidence;
			}
		}

		// 计算碎片整体置信度
		REAL get_vertex_conf() {
			if (_handle->avg_max_var_log.size() > 0) {
				REAL w = _handle->avg_max_var_log.back();
				assert(w >= 0);
				w = max(w,1e-3); // 如果diff=0,近似认为其有微量的更新(防止除零)
				return (1 / w);
			}
			else return 0;
		}


		//更新 _handle->_points_normals => _op
		void update2op() {
			if (!init_flag) {
				init();
				init_flag = true;
			}
			assert(_handle->_points_normals.size() != 0);
			for(int i = 0;i<_op.size();i++){
				_op[i].second = _handle->_points_normals[_op_idx[i]].second;
			}
			__is_changed.set(true);
		}

		void iupdate2op(){
			assert(init_flag);
			assert(_handle->_points_normals.size() != 0);
			for(int i = 0;i<_op.size();i++){
				_handle->_points_normals[_op_idx[i]].second = _op[i].second;
			}
		}


		// 返回op_type为0的点的个数
		int cnt() {
			return std::count(_op_type.begin(), _op_type.end(), 0);
		}

		/**
		 * @brief 
		 * 计算小块的loss (需要有gt)
		 * 由于这个函数可能在不同阶段被调用，所以每次调用都需要重新计算
		 * @return std::pair<REAL,REAL> first表示小块当前与gt的loss;inv_loss表示如果翻转一次,所形成的新小块的loss的值 
		 */
		std::pair<REAL,REAL> _cal_loss(){
			if (!__is_changed.get()) {
				REAL ori_loss = _loss_according2gt.first, inv_loss = _loss_according2gt.second;
				if (is_inved()) {
					return std::pair<REAL, REAL>(inv_loss, ori_loss);
				}
				else {
					return _loss_according2gt;
				}
			}

			REAL loss = 0, inv_loss = 0;
			POINTS_NORMALS gt = get_self_gt();
			POINTS_NORMALS op = get_self_op();
			for(int i = 0;i<op.size();i++){
				REAL tloss = calculate_angle(gt[i].second,op[i].second);
				REAL inv_tloss = M_PI - tloss;
				loss += tloss;
				inv_loss += inv_tloss;
			}
			if (is_inved()) {
				_loss_according2gt = std::pair<REAL, REAL>(inv_loss, loss);
			}
			else {
				_loss_according2gt = std::pair<REAL, REAL>(loss, inv_loss);
			}
			__is_changed.set(false);
			return _cal_loss();
		}

		void inv_op(){
			for (auto& op : _op) {
				op.second *= -1;
			}
			inv_times++;
		}

		graph_vertex(){
			_handle = NULL;
			init_flag = false;
		}
		// graph_vertex(IPSR_HANDLE_P hp, std::vector<int> op_type, POINTS_NORMALS ori_op):_handle(hp),_op_type(op_type),_op(ori_op){
		// 	init();
		// }

		// 计算_op中_op_type为0的点的metric
		lzd_tools::PointNormalMetric<REAL,DIM> cal_metric(){
			POINTS_NORMALS gt = get_self_gt();
			POINTS_NORMALS op = get_self_op();
			return lzd_tools::PointNormalMetric<REAL,DIM>(gt,op);
		}

		nlohmann::json get_log(){
			nlohmann::json j;
			j["inv_times"] = inv_times;
			j["ipsr_log"] = _handle->get_log();
			return j;
		}

		nlohmann::json get_metric(){
			return cal_metric().to_json();
		}

		nlohmann::json to_json(){
			nlohmann::json j;
			j["id"] = -1; //需要在graph 的to_json中重新赋值
			Point<REAL, DIM> center;
			for(int i = 0;i<_op.size();i++){
				center += _op[i].first;
			}
			center /= _op.size();
			j["center"] = {center[0],center[1],center[2]};
			j["log"] = get_log();
			j["metric"] = get_metric();
			return j;
		}

		// 返回_op中_op_type为0的点	
		POINTS_NORMALS get_self_op(){
			POINTS_NORMALS res;
			for(int i = 0;i<_op.size();i++){
				if(_op_type[i] == 0){
					res.push_back(_op[i]);
				}
			}
			return res;
		}

		POINTS_NORMALS get_all_gt() {
			assert(init_flag);
			POINTS_NORMALS res;
			for (int i = 0; i < _op.size(); i++) {
				auto t_op = _handle->_gt_points_normals[_op_idx[i]];
				t_op.first = _handle->ixform * t_op.first;
				res.push_back(t_op);
			}
			return res;
		}

		// 返回_op中_op_type为0的点对应的gt
		POINTS_NORMALS get_self_gt(){
			POINTS_NORMALS res;
			for(int i = 0;i<_op.size();i++){
				if(_op_type[i] == 0){
					auto t_op = _handle->_gt_points_normals[_op_idx[i]];
					t_op.first = _handle->ixform * t_op.first;
					res.push_back(t_op);
				}
			}
			return res;
		}

		MESH get_self_mesh() {
			MESH tm;
			_handle->get_mesh(tm);
			for (int i = 0; i < tm.first.size(); i++) {
				tm.first[i] = _handle->ixform * tm.first[i];
			}
			auto op = get_self_op();
			std::vector<REAL> ttt;
			return shrink_boundary(tm, op, ttt, 20);
		}

	};
	

	// 有向边。
	// start、end都是指向VERTEX的指针
	// 后续可能会增加overlap等属性
	template<typename REAL, int DIM>
	class graph_edge {
		REAL weight;// 权重
		REAL inv_weight;// 翻转一个节点后的权重。TODO
		int cal_times = 0;// 计算权重的次数 确保权重只计算一次
		REAL length;// 边的长度,即两个顶点的center_op的距离
		
		REAL _confidence = 1; // 权重的置信度
		REAL _connection_strength = 1.0;// 边的强度。越高说明两个块之间的连接越强。使用边界点的数量作为边的连接强度

	public:
		static int _overlap_treshold; // 大于该值的overlap部分才会被使用

		graph_vertex<REAL, DIM>* start;
		graph_vertex<REAL, DIM>* end;
		// overlap部分在两个顶点的点云中的索引，overlap部分并不对称，只记录了start中来自end的点。
		std::vector<std::pair<int, int>> overlap_idxs;

		bool is_vaild() {
			return (int)overlap_idxs.size() >= graph_edge::_overlap_treshold && _connection_strength > 0; // 注意connetcion_strength必须大于零，否则所有顶点都会被认为是相连的
																									  // 但允许overlap_idxs等于零，因为可以不扩张
		}

		/**
		 * @brief
		 * 判断这个边是否是"误边"
		 * 这个边是"误边"有两种情况:
		 * 	1. weight<inv_weight, 但是start和end的翻转状态不一致
		 * 	2. weight>=inv_weight, 但是start和end的翻转状态一致
		 * 返回:true表示此边是一个"误边"
		 */
		bool is_wrong_edge() {
			int ss = start->inved_status_according2gt(), se = end->inved_status_according2gt();
			auto w = get_weight();
			if (w.first < w.second) {
				if (ss != se) return true;
			}
			else {
				if (ss == se) return true;
			}
			return false;
		}

		// 获得权重的置信度
		REAL get_brief() {
			return _confidence;
			//return std::sqrt(start->get_vertex_conf() * end->get_vertex_conf());
			//return overlap_idxs.size() * (max(weight, inv_weight) - min(weight, inv_weight));
			//return 1;
		}

		REAL get_connection_strength() { return _connection_strength; }

		void set_connection_strength(REAL strength) { _connection_strength = strength; }

		/**
		 * @brief 根据当前状态或者指定的状态，返回边的权重
		 * @param status 用于指定两个顶点的翻转状态是否一致，1表示一致,0表示不一致，其他值表示不指定。如果不指定，则根据start和end当前的翻转状态来决定
		 * @return std::pair<REAL, REAL> 
		 */
		std::pair<REAL, REAL> get_weight(int status = -1) {			
			assert(cal_times == 1);
			if(status != -1){
				if(status == 0){
					return std::pair<REAL, REAL>(inv_weight, weight);
				}else if (status == 1){
					return std::pair<REAL, REAL>(weight, inv_weight);
				}
			}

			if (start->is_inved() != end->is_inved()) {
				return std::pair<REAL, REAL>(inv_weight, weight);
			}
			else {
				return std::pair<REAL, REAL>(weight, inv_weight);
			}
		}

		/**
		 * @brief 设置边的权重weight
		 * @return REAL
		 */
		void set_weight(std::pair<REAL, REAL> weight, REAL confidence = 1) {
			assert(cal_times == 0);
			this->weight = weight.first;
			this->inv_weight = weight.second;
			_confidence = confidence;
			cal_times++;
		}

		graph_edge(graph_vertex<REAL, DIM>* s, graph_vertex<REAL, DIM>* e) :start(s), end(e) {
			weight = -1e9;
			inv_weight = -1e9;
			cal_times = 0;
			length = Distance(s->center_op, e->center_op);
			overlap_idxs = std::vector<std::pair<int, int>>();
		}

		graph_edge(graph_vertex<REAL, DIM>* s, graph_vertex<REAL, DIM>* e, std::vector<std::pair<int, int>> overlap_idxs) :start(s), end(e), overlap_idxs(overlap_idxs) {
			weight = -1e9;
			inv_weight = -1e9;
			cal_times = 0;
			length = Distance(s->center_op, e->center_op);
		}

		graph_edge() = delete;

		nlohmann::json to_json() {
			nlohmann::json j;
			j["start"] = -1;
			j["end"] = -1;
			auto w = get_weight();
			j["weight"] = w.first;
			j["inv_weight"] = w.second;
			j["overlap_size"] = overlap_idxs.size();
			j["correctness"] = is_wrong_edge() ? "incorrect" : "correct";
			j["confidence"] = get_brief();
			j["connection_strength"] = get_connection_strength();
			return j;
		}
	};

	template<typename REAL, int DIM>
	using VERTEX = graph_vertex<REAL, DIM>;

	template<typename REAL, int DIM>
	using VERTEX_P = VERTEX<REAL, DIM>*;

	template<typename REAL, int DIM>
	using EDGE = graph_edge<REAL, DIM>;

	template<typename REAL, int DIM>
	using EDGE_P = EDGE<REAL, DIM>*;

	template<typename REAL, int DIM>
	class edge_weight_calculator {
	protected:
		nlohmann::json caculate_res(const EDGE_P<REAL, DIM> edge){
			std::pair<REAL, REAL> weight;
			REAL conf;
			cal_weight(edge,weight,conf);
			nlohmann::json j;
			j["weight"] = weight;
			j["conf"] = conf;
			bool is_wrong = (edge->start->inved_status_according2gt() == edge->end->inved_status_according2gt()) == (weight.first > weight.second);
			j["correctness"] = is_wrong ? "incorrect" : "correct";
			j["gt_filp_status"] = { edge->start->inved_status_according2gt(),edge->end->inved_status_according2gt() };
			return j;
		}

	public:
		virtual void cal_weight(const EDGE_P<REAL, DIM> edge,std::pair<REAL,REAL>& weight, REAL& conf) = 0;
		virtual nlohmann::json get_config() = 0;
		// 将判断可视化; 建议在判断错误的case中使用
		virtual void visualization(const EDGE_P<REAL, DIM> edge,std::string floder) = 0;
	};

	template<typename REAL, int DIM>
	class linked_graph;

	/**
	 * @brief 
	 * weight = overlap部分中 start中的点与end中的点的平均角度差 * 两点的置信度的乘积
	 */
	template<typename REAL, int DIM>
	class overlap_diff_calculator : public edge_weight_calculator<REAL,DIM> {
		bool use_confidence;
	public:
		overlap_diff_calculator(bool use_confidence = true):use_confidence(use_confidence){}

		void cal_weight(const EDGE_P<REAL, DIM> edge,std::pair<REAL,REAL>& weight,REAL& conf){
			assert(edge->overlap_idxs.size() > 0); // 如果使用overlap来计算weight，则要求两个块之间必须在空间上接壤

			REAL diff = 0, inv_diff = 0;
			assert(edge->start->is_inved() == false);
			assert(edge->end->is_inved() == false);
			if(use_confidence){
				edge->start->cal_conf();
				edge->end->cal_conf();
			}
			for(auto p:edge->overlap_idxs){
				int s = p.first, e = p.second;
				REAL t = calculate_angle(edge->start->_op[s].second,edge->end->_op[e].second);
				REAL inv_t = M_PI - t;
				if(use_confidence){
					t     = t     * edge->start->_op_confidence[s] * edge->end->_op_confidence[e];
					inv_t = inv_t * edge->start->_op_confidence[s] * edge->end->_op_confidence[e];
				}
				diff     += t;
				inv_diff += inv_t;
			}
			diff /= edge->overlap_idxs.size();
			inv_diff /= edge->overlap_idxs.size();
			weight = std::pair<REAL,REAL>(diff/(diff+inv_diff), inv_diff/(diff + inv_diff));
			conf = max(diff, inv_diff) / (diff + inv_diff);
		}

		nlohmann::json get_config(){
			nlohmann::json j;
			j["name"] = "overlap_diff_calculator";
			j["use_op_confidence"] = use_confidence?"true":"false";
			return j;
		}

		void visualization(const EDGE_P<REAL, DIM> edge,std::string floder){
			rmkdir(floder);
			const int overlap_flag = 10;
			std::string op_path;
			POINTS_NORMALS sop = edge->start->_op;
			POINTS_NORMALS eop = edge->end->_op;
			POINTS_NORMALS op;
			std::vector<int> sop_type(edge->start->_op_type);
			std::vector<int> eop_type(edge->end->_op_type);
			std::vector<int> op_type;
			std::vector<std::pair<int,int>> overlap_idxs = edge->overlap_idxs;
			for(int i = 0;i<overlap_idxs.size();i++){
				int s = overlap_idxs[i].first, e = overlap_idxs[i].second;
				sop_type[s] = overlap_flag;
				eop_type[e] = overlap_flag;
			}
			for(int i = 0;i<sop.size();i++){
				if(sop_type[i] == overlap_flag || sop_type[i] == 0){
					op.push_back(sop[i]);
					op_type.push_back(sop_type[i]);
				}
			}
			op_path = floder + "/start.ply";
			lzd_tools::op2ply(op, op_path, XForm<REAL, DIM + 1>().Identity(), std::make_pair(op_type, lzd_tools::get_regular_colormap(0,overlap_flag)));
			op.clear();
			op_type.clear();
			op_path = floder + "/end.ply";			
			for(int i = 0;i<eop.size();i++){
				if(eop_type[i] == overlap_flag || eop_type[i] == 0){
					op.push_back(eop[i]);
					op_type.push_back(eop_type[i]);
				}
			}
			lzd_tools::op2ply(op,op_path,XForm<REAL,DIM+1>().Identity(), std::make_pair(op_type, lzd_tools::get_regular_colormap(0, overlap_flag)));
			std::string log_path = floder + "/log.json";
			nlohmann::json j;
			j["config"] = get_config();
			j["cal_res"] = caculate_res(edge);
			std::ofstream out(log_path);
			out << j.dump(4);
			out.close();
		}
	};

	// TODO 统一两种计算方式
	// template<typename REAL, int DIM>
	// class ConsistencyCaculator:public edge_weight_calculator<REAL,DIM>{
	// protected:
	// 	std::map<VERTEX_P<REAL, DIM>, int> _vertex2idx;
	// 	std::vector<std::vector<std::pair<REAL, REAL>>> _weight;
	// 	REAL _max_weight,_min_weight;
	// 	// 将consistency转为权重
	// 	std::pair<REAL, REAL> __toweight(REAL ori_consistency, REAL inv_consistency) {
	// 		return std::pair<REAL, REAL>(1 - ori_consistency / (ori_consistency + inv_consistency), 1 - inv_consistency / (ori_consistency + inv_consistency));
	// 	}
	// 	// 根据两个weight给出置信度
	// 	REAL _weight2conf(std::pair<REAL, REAL> weight) {
	// 		return max(weight.first, weight.second) / (weight.first + weight.second);
	// 	}
	// 	REAL _normalize_weight(REAL weight) {
	// 		REAL w =  (weight - _min_weight) / (_max_weight - _min_weight);
	// 		assert(w >= 0 && w <= 1);
	// 		return w;
	// 	}
	// 	virtual _cal_edge_weight(const EDGE_P<REAL, DIM> edge, std::pair<REAL, REAL>& weight) = 0;
	// public:
	// }

	/**
	 * @brief 
	 * weight为不同视角下，可见面片的朝向一致性的均值 
	 */
	template<typename REAL, int DIM>
	class consistency_calculator : public edge_weight_calculator<REAL, DIM> {
		MeshConsistency<REAL,DIM>* _oc;
		IPSR_Factory<REAL, DIM>* _factory;
		std::vector<MESH> _meshes;
		std::map<VERTEX_P<REAL, DIM>, int> _vertex2idx;
		std::vector<std::vector<std::pair<REAL, REAL>>> _weight;
		REAL _max_weight,_min_weight;
		// 将consistency转为权重
		std::pair<REAL, REAL> __toweight(REAL ori_consistency, REAL inv_consistency) {
			if (ori_consistency == inv_consistency) {
				assert(ori_consistency + inv_consistency == 0);
				return std::pair<REAL, REAL>(0.5, 0.5);
			}
			return std::pair<REAL, REAL>(1 - ori_consistency / (ori_consistency + inv_consistency), 1 - inv_consistency / (ori_consistency + inv_consistency));
		}

		std::vector<int> _get_mesh(const EDGE_P<REAL, DIM> edge, MESH& ori_mesh, MESH& inv_end_mesh){
			assert(ori_mesh.first.size() == 0 && inv_end_mesh.first.size() == 0);
			assert(ori_mesh.second.size() == 0 && inv_end_mesh.second.size() == 0);
			int s = _vertex2idx[edge->start], e = _vertex2idx[edge->end];
			MESH end_mesh = _meshes[e];
			lzd_tools::add_topology(ori_mesh,_meshes[s]);
			lzd_tools::add_topology(ori_mesh,end_mesh);
			lzd_tools::add_topology(inv_end_mesh,_meshes[s]);
			lzd_tools::revert_mesh(end_mesh);
			lzd_tools::add_topology(inv_end_mesh, end_mesh);
			std::vector<int> tri_type(ori_mesh.second.size(), 0);
			for (int i = 0; i < ori_mesh.second.size(); i++) {
				if (i < _meshes[s].second.size()) {
					tri_type[i] = 1;
				}
				else {
					tri_type[i] = 2;
				}
			}
			return tri_type;
		}

		// 根据两个weight给出置信度
		REAL _weight2conf(std::pair<REAL, REAL> weight) {
			return max(weight.first, weight.second) / (weight.first + weight.second);
		}

		void _cal_edge_weight(const EDGE_P<REAL, DIM> edge, std::pair<REAL, REAL>& weight) {
			MESH ori_mesh, inv_end_mesh;
			assert(edge->start->is_inved() == false);
			assert(edge->end->is_inved() == false);
			auto tri_type = _get_mesh(edge,ori_mesh, inv_end_mesh);
			std::vector<REAL> ori_consistencies;
			std::vector<REAL> inv_consistencies;
			REAL ori_consistency = _oc->cal_consistency(ori_mesh, ori_consistencies,"",tri_type);
			REAL inv_consistency = _oc->cal_consistency(inv_end_mesh, inv_consistencies,"", tri_type);
			weight = __toweight(ori_consistency, inv_consistency);
		}

		REAL _normalize_weight(REAL weight) {
			REAL w =  (weight - _min_weight) / (_max_weight - _min_weight);
			assert(w >= 0 && w <= 1);
			return w;
		}

		consistency_calculator():_oc(),_factory(NULL){
			printf("consistency_calculator construction called by son class...\n");
		}

	public:
		consistency_calculator(IPSR_Factory<REAL, DIM>* factory, const linked_graph<REAL, DIM>* lg, MeshConsistency<REAL,DIM>* _oc,int subdivision_times):_factory(factory),_oc(_oc){
			printf("consistency_calculator construction...\n");
			int num_nodes = lg->_nodes.size();
			_meshes.resize(num_nodes);
			lzd_tools::thread_safe_int finished(0);
#pragma omp parallel for
			for (int i = 0; i < num_nodes; i++) {
				if (subdivision_times > 0)_meshes[i] = subdivision_mesh(lg->_nodes[i]->get_self_mesh(),subdivision_times);
				else _meshes[i] = lg->_nodes[i]->get_self_mesh();
				_vertex2idx[lg->_nodes[i]] = i;
				++finished;
				// 打印进度 finish/total
				printf("\r%d/%d", finished.get(), (int)lg->_nodes.size());
			}

#pragma omp parallel for
			for (int i = 0; i < num_nodes; i++) {
				std::vector<REAL> temp;
				auto op = lg->_nodes[i]->get_self_op();
				//assert(_meshes[i].first.size() != 0 || op.size() < 10);
				_meshes[i] = clean_mesh(_meshes[i], op, temp,20 + 20 * 3 * subdivision_times);
			}

			printf("\npre_cal_weight...\n");
			_weight = std::vector<std::vector<std::pair<REAL, REAL>>>(num_nodes, std::vector<std::pair<REAL, REAL>>(num_nodes, std::pair<REAL, REAL>(-1, -1)));
			auto all_edges = lg->_edges;
			finished.set(0);
			int total_edges = 0;
			for (int i = 0; i < all_edges.size(); i++) {
				total_edges += all_edges[i].size();
			}
			_max_weight = 0;
			_min_weight = 1;
#pragma omp parallel for
			for (int i = 0; i < all_edges.size(); i++) {
				for (int j = 0; j < all_edges[i].size(); j++) {
					int s = _vertex2idx[all_edges[i][j]->start], e = _vertex2idx[all_edges[i][j]->end];
					if(s == e) continue;
					auto semi_weight = _weight[s][e];
					// 如果其逆边已经计算过，则直接使用
					if(semi_weight != std::pair<REAL,REAL>(-1,-1)){
						continue;
					}

					// 否则计算权重，并将其一并赋值给逆边
					_cal_edge_weight(all_edges[i][j], _weight[s][e]);
					_max_weight = max(_max_weight, max(_weight[s][e].first, _weight[s][e].second));
					_min_weight = min(_min_weight, min(_weight[s][e].first, _weight[s][e].second));
					_weight[e][s] = _weight[s][e];
					++finished;
					++finished;
					printf("\r%d/%d", finished.get(), total_edges);
				}
			}
			
			printf("\npre_cal_weight done!\n");
			if (finished.get() > total_edges) {
				printf("%d edge is not semi-symmetric\n", finished.get() - total_edges);
			}
			printf("consistency_calculator construction done!\n");
		}

		nlohmann::json get_config() {
			nlohmann::json j;
			j["name"] = "consistency_calculator";
			j["consistency_calculator_config"] = _oc->get_config();
			return j;
		}

		void cal_weight(const EDGE_P<REAL, DIM> edge, std::pair<REAL, REAL>& weight, REAL& conf) {
			// 翻转操作应该在cal_weight之后进行
			assert(edge->start->is_inved() == false);
			assert(edge->end->is_inved() == false);
			std::pair<REAL, REAL> unnormalized_weight = _weight[_vertex2idx[edge->start]][_vertex2idx[edge->end]];
			weight = std::pair<REAL, REAL>(_normalize_weight(unnormalized_weight.first), _normalize_weight(unnormalized_weight.second));
			conf = _weight2conf(weight);
		}

		void visualization(EDGE_P<REAL, DIM> edge, std::string floder) {
			rmkdir(floder);
			std::string ori_mesh_path = floder + "/ori_";
			std::string inv_mesh_path = floder + "/inv_";
			std::string log_path = floder + "/log.json";
			
			MESH ori_mesh, inv_end_mesh;
			auto tri_type = _get_mesh(edge,ori_mesh, inv_end_mesh);

			std::vector<REAL> ori_consistencies;
			std::vector<REAL> inv_consistencies;
			REAL ori_consistency = _oc->cal_consistency(ori_mesh,     ori_consistencies, ori_mesh_path,tri_type);
			REAL inv_consistency = _oc->cal_consistency(inv_end_mesh, inv_consistencies, inv_mesh_path,tri_type);
			

			nlohmann::json j;
			j["ori_consistencies"] = ori_consistencies;
			j["inv_consistencies"] = inv_consistencies;
			j["ori_avg_consistency"] = ori_consistency;
			j["inv_avg_consistency"] = inv_consistency;
	
			nlohmann::json log_j;
			log_j["config"] = get_config();
			log_j["cal_res"] = caculate_res(edge);
			log_j["consistency"] = j;	
			std::ofstream out(log_path);
			out << log_j.dump(4);
			out.close();
		}
	};

	/**
	 * @brief 
	 * 相比于consistency_calculator,其每个边单独生成一个mesh
	 * @tparam REAL 
	 * @tparam DIM 
	 */
	template<typename REAL, int DIM>
	class consistency_calculator_plus:public edge_weight_calculator<REAL,DIM>{
		MaxAB_Mesh<REAL, DIM> _oc;
		UnorderedPairMap<std::pair<MESH,MESH>> _mesh_graph; // 用于存储每个边对应的ori_mesh和inv_end_mesh. 由于不考虑overlap，所以边是无序的
		UnorderedPairMap<std::pair<REAL,REAL>> _weight;
		UnorderedPairMap<std::pair<std::vector<int>, std::vector<int>>> _tri_type;
		std::map<VERTEX_P<REAL, DIM>, int> _vertex2idx;
	
		REAL _max_weight, _min_weight;
		// 将consistency转为权重
		std::pair<REAL, REAL> __toweight(REAL ori_consistency, REAL inv_consistency) {
			if (ori_consistency == inv_consistency) {
				assert(ori_consistency + inv_consistency == 0);
				return std::pair<REAL, REAL>(0.5, 0.5);
			}
			assert(ori_consistency + inv_consistency != 0);
			return std::pair<REAL, REAL>(1 - ori_consistency / (ori_consistency + inv_consistency), 1 - inv_consistency / (ori_consistency + inv_consistency));
		}

		// 根据两个weight给出置信度
		REAL _weight2conf(std::pair<REAL, REAL> weight) {
			return max(weight.first, weight.second) / (weight.first + weight.second);
		}

		REAL _normalize_weight(REAL weight) {
			REAL w = (weight - _min_weight) / (_max_weight - _min_weight);
			assert(w >= 0 && w <= 1);
			return w;
		}

	 	std::pair<std::vector<int>,std::vector<int>> _get_mesh(const EDGE_P<REAL, DIM> edge, MESH& ori_mesh, MESH& inv_end_mesh) {
			std::pair<int, int> e = std::make_pair(_vertex2idx[edge->start], _vertex2idx[edge->end]);
			ori_mesh = _mesh_graph[e].first;
			inv_end_mesh = _mesh_graph[e].second;
			return _tri_type[e];
		}

		
		void _cal_edge_weight(const EDGE_P<REAL, DIM> edge, std::pair<REAL, REAL>& weight) {
			MESH ori_mesh, inv_end_mesh;
			assert(edge->start->is_inved() == false);
			assert(edge->end->is_inved() == false);
			auto tri_type = _get_mesh(edge,ori_mesh, inv_end_mesh);
			std::vector<REAL> ori_consistencies;
			std::vector<REAL> inv_consistencies;
			REAL ori_consistency = _oc.cal_consistency(ori_mesh, ori_consistencies,"",tri_type.first);
			REAL inv_consistency = _oc.cal_consistency(inv_end_mesh, inv_consistencies,"", tri_type.second);
			weight = __toweight(ori_consistency, inv_consistency);
		}	

	public:

		void cal_weight(const EDGE_P<REAL, DIM> edge, std::pair<REAL, REAL>& weight, REAL& conf) {
			// 翻转操作应该在cal_weight之后进行
			assert(edge->start->is_inved() == false);
			assert(edge->end->is_inved() == false);
			std::pair<REAL, REAL> unnormalized_weight = _weight[std::pair<int,int>(_vertex2idx[edge->start],_vertex2idx[edge->end])];
			weight = std::pair<REAL, REAL>(_normalize_weight(unnormalized_weight.first), _normalize_weight(unnormalized_weight.second));
			conf = _weight2conf(weight);
		}

		consistency_calculator_plus(IPSR_Factory<REAL, DIM>* factory, const linked_graph<REAL, DIM>* lg) : _oc() {
			printf("consistency_calculator_plus construction...\n");
			int num_nodes = lg->_nodes.size();
			lzd_tools::thread_safe_int finished(0);

			for (int i = 0; i < lg->_nodes.size(); i++) {
				_vertex2idx[lg->_nodes[i]] = i;
			}

			auto all_edges = lg->_edges;
			int total_edges = 0;
			for (int i = 0; i < lg->_edges.size(); i++) {
				total_edges += lg->_edges[i].size();
			}
			
			printf("\n construct mesh_graph...\n");
#pragma omp parallel for
			for (int i = 0; i < lg->_nodes.size(); i++) {
				for (int j = 0; j < lg->_edges[i].size(); j++) {
					printf("\r%d/%d", finished.get(), (int)total_edges);
					auto edge = lg->_edges[i][j];					
					std::pair<int,int> e = std::make_pair(_vertex2idx[edge->start], _vertex2idx[edge->end]);
					if (_mesh_graph.has(e)) {
						continue; // 跳过逆边
					}
					assert(_vertex2idx[edge->start] == i);
				 
				    POINTS_NORMALS ori_op, inv_end_op, ori_a, ori_b;
					ori_a = edge->start->get_self_op();
					ori_b  = edge->end->get_self_op();
					ori_op.insert(ori_op.end(), ori_a.begin(), ori_a.end());
					inv_end_op.insert(inv_end_op.end(), ori_a.begin(), ori_a.end());
					ori_op.insert(ori_op.end(), ori_b.begin(), ori_b.end());
					lzd_tools::revert_op(ori_b);
					inv_end_op.insert(inv_end_op.end(), ori_b.begin(), ori_b.end());
					
					MESH ori_mesh = clean_mesh(factory->get_mesh_from_points_normals(ori_op), ori_op, std::vector<REAL>());
					MESH inv_end_mesh = clean_mesh(factory->get_mesh_from_points_normals(inv_end_op), ori_op, std::vector<REAL>());
					_mesh_graph[e] = std::pair<MESH, MESH>(ori_mesh, inv_end_mesh);

										
					std::vector<int> op_type(ori_a.size(), 1);
					op_type.insert(op_type.end(), ori_b.size(), 2);
					// 由ori_op建立kdtree, tri_type = 最近的点的类型
					std::vector<int> _ori_tri_type(ori_mesh.second.size(), 0);
					std::vector<int> _inv_tri_type(inv_end_mesh.second.size(), 0);

					auto pcl_op = op2pcl(ori_op);
					pcl::KdTreeFLANN<pcl::PointNormal> kdtree;
					kdtree.setInputCloud(pcl_op);
					std::vector<int> pointIdxNKNSearch(1);
					std::vector<float> tp(1);
					for (int i = 0; i < ori_mesh.second.size(); i++) {
						Point<REAL,DIM> c = ori_mesh.first[ori_mesh.second[i][0]] + ori_mesh.first[ori_mesh.second[i][1]] + ori_mesh.first[ori_mesh.second[i][2]];
						c /= 3;
						pcl::PointNormal searchPoint{ (float)c[0],(float)c[1],(float)c[2] }; // TODO 把这个nc float类型的kdtree换掉	
						assert(kdtree.nearestKSearch(searchPoint, 1, pointIdxNKNSearch,tp) > 0);
						_ori_tri_type[i] = op_type[pointIdxNKNSearch[0]];
					
					}
					for (int i = 0; i < inv_end_mesh.second.size(); i++) {
						Point<REAL, DIM> c = inv_end_mesh.first[inv_end_mesh.second[i][0]] + inv_end_mesh.first[inv_end_mesh.second[i][1]] + inv_end_mesh.first[inv_end_mesh.second[i][2]];
						c /= 3;
						pcl::PointNormal searchPoint{ (float)c[0],(float)c[1],(float)c[2] };
						assert(kdtree.nearestKSearch(searchPoint, 1, pointIdxNKNSearch, tp) > 0);
						_inv_tri_type[i] = op_type[pointIdxNKNSearch[0]];
					}
					_tri_type[e] = std::pair<std::vector<int>, std::vector<int>>(_ori_tri_type, _inv_tri_type);
					++finished;
				}
			}
			printf("\nconstruct mesh_graph done!\n");

			printf("\npre_cal_weight...\n");
			finished.set(0);
#pragma omp parallel for
			for (int i = 0; i < all_edges.size(); i++) {
				for (int j = 0; j < all_edges[i].size(); j++) {
					printf("\r%d/%d", finished.get(), (int)total_edges);
					int s = _vertex2idx[all_edges[i][j]->start], e = _vertex2idx[all_edges[i][j]->end];
					std::pair<REAL, REAL> weight;
					_cal_edge_weight(all_edges[i][j], weight);
					_weight[std::make_pair(s, e)] = weight;
					_max_weight = max(_max_weight, max(weight.first, weight.second));
					_min_weight = min(_min_weight, min(weight.first, weight.second));
					++finished;
				}
			}
			printf("\npre_cal_weight done!\n");
			printf("consistency_calculator_plus construction done!\n");

		}
		
		nlohmann::json get_config() {
			nlohmann::json j;
			j["name"] = "consistency_calculator_plus";
			j["consistency_calculator_config"] = _oc.get_config();
			return j;
		}

		void visualization(EDGE_P<REAL, DIM> edge, std::string floder) {
			rmkdir(floder);
			std::string ori_mesh_path = floder + "/ori_";
			std::string inv_mesh_path = floder + "/inv_";
			std::string log_path = floder + "/log.json";

			MESH ori_mesh, inv_end_mesh;
			auto tri_type = _get_mesh(edge, ori_mesh, inv_end_mesh);

			std::vector<REAL> ori_consistencies;
			std::vector<REAL> inv_consistencies;
			REAL ori_consistency = _oc.cal_consistency(ori_mesh, ori_consistencies, ori_mesh_path, tri_type.first);
			REAL inv_consistency = _oc.cal_consistency(inv_end_mesh, inv_consistencies, inv_mesh_path, tri_type.second);


			nlohmann::json j;
			j["ori_consistencies"] = ori_consistencies;
			j["inv_consistencies"] = inv_consistencies;
			j["ori_avg_consistency"] = ori_consistency;
			j["inv_avg_consistency"] = inv_consistency;

			nlohmann::json log_j;
			log_j["config"] = get_config();
			log_j["cal_res"] = caculate_res(edge);
			log_j["consistency"] = j;
			std::ofstream out(log_path);
			out << log_j.dump(4);
			out.close();
		}
	};



	template<typename REAL, int DIM>
	class point_consistency_calculator:public edge_weight_calculator<REAL,DIM>{
		PointOrientationConsistency<REAL,DIM> _poc;
		
		std::map<VERTEX_P<REAL, DIM>, int> _vertex2idx;
		std::vector<std::vector<std::pair<REAL, REAL>>> _weight;
		REAL _max_weight,_min_weight;

		void _get_op(const EDGE_P<REAL, DIM> edge, ORIENTED_POINTS& ori_op, ORIENTED_POINTS& inv_op) {
			assert(ori_op.size() == 0 && inv_op.size() == 0);
			int s = _vertex2idx[edge->start], e = _vertex2idx[edge->end];
			auto start_op = edge->start->_op;
			auto end_op = edge->end->_op;
			for (int i = 0; i < start_op.size(); i++) {
				ori_op.push_back(start_op[i]);
				inv_op.push_back(start_op[i]);
			}
			for (int i = 0; i < end_op.size(); i++) {
				auto op = end_op[i];
				ori_op.push_back(op);
				op.second *= -1;
				inv_op.push_back(op);
			}
		}

		void _cal_edge_weight(const EDGE_P<REAL, DIM> edge, std::pair<REAL, REAL>& weight) {
			assert(edge->start->is_inved() == false);
			assert(edge->end->is_inved() == false);
			ORIENTED_POINTS ori_op, inv_op;
			_get_op(edge,ori_op, inv_op);
			std::vector<REAL> ori_consistencies;
			std::vector<REAL> inv_consistencies;
			REAL ori_consistency = _poc.cal_consistency(ori_op, ori_consistencies);
			REAL inv_consistency = _poc.cal_consistency(inv_op, inv_consistencies);
			weight =  std::pair<REAL, REAL>(1 - ori_consistency / (ori_consistency + inv_consistency), 1 - inv_consistency / (ori_consistency + inv_consistency));
		}

		REAL _normalize_weight(REAL weight) {
			REAL w =  (weight - _min_weight) / (_max_weight - _min_weight);
			assert(w >= 0 && w <= 1);
			return w;
		}

		// 根据两个weight给出置信度
		REAL _weight2conf(std::pair<REAL, REAL> weight) {
			return max(weight.first, weight.second) / (weight.first + weight.second);
		}

	public:
		point_consistency_calculator(double radius, const linked_graph<REAL, DIM>* lg){
			printf("point_consistency_calculator construction...\n");
			_poc = PointOrientationConsistency<REAL,DIM>(radius);
			int num_nodes = lg->_nodes.size();
			lzd_tools::thread_safe_int finished(0);
#pragma omp parallel for
			for (int i = 0; i < num_nodes; i++) {
				_vertex2idx[lg->_nodes[i]] = i;
				++finished;
				// 打印进度 finish/total
				printf("\r%d/%d", finished.get(), (int)lg->_nodes.size());
			}

			printf("\npre_cal_weight...\n");
			_weight = std::vector<std::vector<std::pair<REAL, REAL>>>(num_nodes, std::vector<std::pair<REAL, REAL>>(num_nodes, std::pair<REAL, REAL>(-1, -1)));
			auto all_edges = lg->_edges;
			finished.set(0);
			int total_edges = 0;
			for (int i = 0; i < all_edges.size(); i++) {
				total_edges += all_edges[i].size();
			}
			_max_weight = 0;
			_min_weight = 1;
#pragma omp parallel for
			for (int i = 0; i < all_edges.size(); i++) {
				for (int j = 0; j < all_edges[i].size(); j++) {
					int s = _vertex2idx[all_edges[i][j]->start], e = _vertex2idx[all_edges[i][j]->end];
					_cal_edge_weight(all_edges[i][j], _weight[s][e]);
					_max_weight = max(_max_weight, max(_weight[s][e].first, _weight[s][e].second));
					_min_weight = min(_min_weight, min(_weight[s][e].first, _weight[s][e].second));
					++finished;
					printf("\r%d/%d", finished.get(), total_edges);
				}
			}

			printf("\npre_cal_weight done!\n");
			printf("point_consistency_calculator construction done!\n");


		}
		nlohmann::json get_config(){
			nlohmann::json j;
			j["name"] = "point_consistency_calculator";
			j["CalculatorConfig"] = _poc.get_config();
			return j;
		}

		void cal_weight(const EDGE_P<REAL, DIM> edge,std::pair<REAL,REAL>& weight,REAL& conf){
			assert(edge->start->is_inved() == false);
			assert(edge->end->is_inved() == false);
			ORIENTED_POINTS ori_op, inv_end_op;
			_get_op(edge,ori_op, inv_end_op);
			std::pair<REAL, REAL> unnormalized_weight = _weight[_vertex2idx[edge->start]][_vertex2idx[edge->end]];
			weight = std::pair<REAL, REAL>(_normalize_weight(unnormalized_weight.first), _normalize_weight(unnormalized_weight.second));
			conf = _weight2conf(weight);
		}

		void visualization(const EDGE_P<REAL, DIM> edge,std::string floder){
			rmkdir(floder);
			std::string ori_op_path = floder + "/ori_op.ply";
			std::string inv_op_path = floder + "/inv_op.ply";
			std::string log_path = floder + "/log.json";
			ORIENTED_POINTS ori_op, inv_op;
			_get_op(edge,ori_op, inv_op);
			std::vector<REAL> ori_consistencies;
			std::vector<REAL> inv_consistencies;
			REAL ori_consistency = _poc.cal_consistency(ori_op, ori_consistencies, ori_op_path);
			REAL inv_consistency = _poc.cal_consistency(inv_op, inv_consistencies, inv_op_path);
			nlohmann::json j;
			j["ori_consistencies"] = ori_consistencies;
			j["inv_consistencies"] = inv_consistencies;
			j["ori_avg_consistency"] = ori_consistency;
			j["inv_avg_consistency"] = inv_consistency;
			nlohmann::json log_j;
			log_j["config"] = get_config();
			log_j["cal_res"] = caculate_res(edge);
			log_j["consistency"] = j;	
			std::ofstream out(log_path);
			out << log_j.dump(4);
			out.close();
		}
	};

	template<typename REAL, int DIM>
	class dipole_calculator:public edge_weight_calculator<REAL,DIM>{
	private:
		REAL _eps;
		REAL _max_weight; // 不需要记录_min_weight
		const linked_graph<REAL,DIM>* _lg;
		std::vector<std::vector<std::pair<REAL, REAL>>> _weight;
		std::map<VERTEX_P<REAL, DIM>, int> _vertex2idx;
		
		REAL _normalize_weight(REAL weight) {
			REAL w =  0.5 + weight  / (2 * _max_weight);
			assert(w >= 0 && w <= 1);
			return w;
		}

		void _cal_edge_weight(const EDGE_P<REAL, DIM> edge, std::pair<REAL, REAL>& weight) {
			assert(edge->start->is_inved() == false);
			assert(edge->end->is_inved() == false);
			REAL w = cal_field_intersaction(edge->start->_op, edge->end->_op, _eps);
			weight = std::pair<REAL, REAL>(w, -1 * w);
		}

		void _get_op(const EDGE_P<REAL, DIM> edge, ORIENTED_POINTS& ori_op, ORIENTED_POINTS& inv_op) {
			assert(ori_op.size() == 0 && inv_op.size() == 0);
			int s = _vertex2idx[edge->start], e = _vertex2idx[edge->end];
			auto start_op = edge->start->_op;
			auto end_op = edge->end->_op;
			for (int i = 0; i < start_op.size(); i++) {
				ori_op.push_back(start_op[i]);
				inv_op.push_back(start_op[i]);
			}
			for (int i = 0; i < end_op.size(); i++) {
				auto op = end_op[i];
				ori_op.push_back(op);
				op.second *= -1;
				inv_op.push_back(op);
			}
		}


	public:
		dipole_calculator(REAL eps, const linked_graph<REAL,DIM>* lg):_eps(eps),_lg(lg){
			printf("dipole_calculator construction...\n");
			int num_nodes = lg->_nodes.size();
			// 计算每个边的权重
			_weight = std::vector<std::vector<std::pair<REAL, REAL>>>(num_nodes, std::vector<std::pair<REAL, REAL>>(num_nodes, std::pair<REAL, REAL>(-1, -1)));

			lzd_tools::thread_safe_int finished(0);

			for (int i = 0; i < num_nodes; i++) {
				_vertex2idx[lg->_nodes[i]] = i;
			}

			printf("\npre_cal_weight...\n");
			auto all_edges = lg->_edges;
			finished.set(0);
			int total_edges = 0;
			for (int i = 0; i < all_edges.size(); i++) {
				total_edges += all_edges[i].size();
			}
			_max_weight = 1;
#pragma omp parallel for
			for (int i = 0; i < all_edges.size(); i++) {
				for (int j = 0; j < all_edges[i].size(); j++) {
					int s = _vertex2idx[all_edges[i][j]->start], e = _vertex2idx[all_edges[i][j]->end];
					_cal_edge_weight(all_edges[i][j], _weight[s][e]);
					_max_weight = max(_max_weight, max(_weight[s][e].first, _weight[s][e].second));
					++finished;
					printf("\r%d/%d", finished.get(), total_edges);
				}
			}
			printf("\npre_cal_weight done!\n");
			printf("dipole_calculator construction done!\n");
		}

		nlohmann::json get_config(){
			nlohmann::json j;
			j["name"] = "dipole_calculator";
			j["eps"] = _eps;
			return j;
		}

		void cal_weight(const EDGE_P<REAL, DIM> edge,std::pair<REAL,REAL>& weight,REAL& conf){
			assert(edge->start->is_inved() == false);
			assert(edge->end->is_inved() == false);
			std::pair<REAL, REAL> unnormalized_weight = _weight[_vertex2idx[edge->start]][_vertex2idx[edge->end]];
			weight = std::pair<REAL, REAL>(_normalize_weight(unnormalized_weight.first), _normalize_weight(unnormalized_weight.second));
			conf = max(weight.first, weight.second) / (weight.first + weight.second);
		}

		void visualization(const EDGE_P<REAL,DIM> edge,std::string floder){
			rmkdir(floder);
			std::string log_path = floder + "/log.json";
			nlohmann::json log_j;
			log_j["config"] = get_config();
			log_j["interaction"] = _weight[_vertex2idx[edge->start]][_vertex2idx[edge->end]].first;
			log_j["start_idx"] = _vertex2idx[edge->start];
			log_j["end_idx"] = _vertex2idx[edge->end];
			std::ofstream out(log_path);
			out << log_j.dump(4);
			out.close();
		}
	};

	/**
	 * @brief 
	 * 计算两种权重边和conf
	 * 定义一个权重边的conf_rate为其在所有边中的conf的排名比例(越大表示conf越大)
	 * 如果p_conf_rate与a_conf_rate差距大于gap,则使用aid的权重
	 * 注意使用的是排名比例，而不是conf的值,这样可以保证不同的conf计算方法可以比较;
	 * 但还是尽量保证conf、weight的值域接近，否则后续计算weight_sum时可能会出错
	 */
	template<typename REAL, int DIM>
	class CoCalculator:public edge_weight_calculator<REAL,DIM>{
		edge_weight_calculator<REAL,DIM>* _primary_calculator;
		edge_weight_calculator<REAL,DIM>* _aid_calculator;
		REAL _gap;
		std::vector<REAL> primary_edge_conf;
		std::vector<REAL> aid_edge_conf;

		// float __get_ration(REAL conf, std::vector<REAL>& confs){
		// 	auto it = std::lower_bound(confs.begin(),confs.end(),conf);
		// 	return (float)(it - confs.begin()) / (float)confs.size();
		// }

		// 判断是否使用primary的权重
		// 如果p_conf在primary_edge_conf中的位置比a_conf在aid_edge_conf中的位置靠前的比例大于_gap，则不使用primary的权重
		bool _use_primary(REAL p_conf,REAL a_conf){
			float p_conf_rate = (float)(std::lower_bound(primary_edge_conf.begin(),primary_edge_conf.end(),p_conf) - primary_edge_conf.begin()) / (float)primary_edge_conf.size();
			float a_conf_rate = (float)(std::lower_bound(aid_edge_conf.begin(),aid_edge_conf.end(),a_conf) - aid_edge_conf.begin()) / (float)aid_edge_conf.size();
			return a_conf_rate - p_conf_rate < _gap;
		}
		
	public:
		CoCalculator(edge_weight_calculator<REAL,DIM>* primary_calculator,edge_weight_calculator<REAL,DIM>* aid_calculator,REAL gap,const linked_graph<REAL,DIM>* lg):
		_primary_calculator(primary_calculator),_aid_calculator(aid_calculator),_gap(gap){
			assert(_gap >= 0 && _gap <= 1);
			auto all_edges = lg->_edges;
			for (int i = 0; i < all_edges.size(); i++) {
				for(int j = 0;j<all_edges[i].size();j++){
					std::pair<REAL,REAL> primary_weight,aid_weight;
					REAL tmp_conf = 0;
					_primary_calculator->cal_weight(all_edges[i][j],primary_weight,tmp_conf);
					primary_edge_conf.push_back(tmp_conf);
					_aid_calculator->cal_weight(all_edges[i][j],aid_weight,tmp_conf);
					aid_edge_conf.push_back(tmp_conf);
				}
			}
			std::sort(primary_edge_conf.begin(),primary_edge_conf.end());
			std::sort(aid_edge_conf.begin(),aid_edge_conf.end());
		}

		void cal_weight(const EDGE_P<REAL, DIM> edge,std::pair<REAL,REAL>& weight,REAL& conf){
			std::pair<REAL,REAL> primary_weight,aid_weight;
			REAL p_conf = 0,a_conf = 0;
			_primary_calculator->cal_weight(edge,primary_weight,p_conf);
			_aid_calculator->cal_weight(edge,aid_weight,a_conf);
			// 如果两种weight不一致，则将primary_weight调低
			if (primary_weight.first >= primary_weight.second != aid_weight.first >= aid_weight.second) {
				double e = 1e-5;
				if(_use_primary(p_conf,a_conf)){
					weight = std::pair<REAL,REAL>(primary_weight.first,primary_weight.second);
				}else{
					weight = std::pair<REAL,REAL>(aid_weight.first,aid_weight.second);
				}

				if(weight.first > 0.5){
					weight.first = 0.5 + e;
					weight.second = 0.5 - e;
				}else{
					assert(weight.second >= 0.5);
					weight.first = 0.5 - e;
					weight.second = 0.5 + e;
				}
				conf = 0.5+e;
			}
			else {
				weight = primary_weight;
				conf = p_conf;
			}
		}

		nlohmann::json get_config(){
			nlohmann::json j;
			j["name"] = "CoCalculator";
			j["primary_calculator"] = _primary_calculator->get_config();
			j["aid_calculator"] = _aid_calculator->get_config();
			j["gap"] = _gap;
			return j;
		}

		void visualization(const EDGE_P<REAL, DIM> edge,std::string floder){
			_primary_calculator->visualization(edge,floder + "/primary");
			_aid_calculator->visualization(edge,floder + "/aid");
			nlohmann::json j;
			j["config"] = get_config();
			j["cal_res"] = caculate_res(edge);

			std::pair<REAL,REAL> primary_weight,aid_weight;
			REAL p_conf = 0,a_conf = 0;

			_primary_calculator->cal_weight(edge,primary_weight,p_conf);
			_aid_calculator->cal_weight(edge,aid_weight,a_conf);
			nlohmann::json j1;
			j1["primary_weight"] = primary_weight;
			j1["aid_weight"] = aid_weight;
			j1["primary_conf"] = p_conf;
			j1["aid_conf"] = a_conf;
			j1["use_primary"] = _use_primary(p_conf,a_conf)?"true":"false";
			j["choose_log"] = j1;
			std::ofstream out(floder + "/log.json");
			out << j.dump(4);
			out.close();
		}

	};

	template<typename REAL, int DIM>
	edge_weight_calculator<REAL,DIM>* get_edge_caculator_from_json(nlohmann::json config_j,IPSR_Factory<REAL, DIM>* factory, const linked_graph<REAL, DIM>* lg){
		std::string name = config_j["name"];
		if(name == "overlap_diff_calculator"){
			return new overlap_diff_calculator<REAL,DIM>(config_j["use_op_confidence"] == "true");
		}else if(name == "consistency_calculator"){
			auto fun_config_j = config_j["consistency_calculator_config"];
			int width = fun_config_j["width"];
			int height = fun_config_j["height"];
			std::string projection_type = fun_config_j["projection_type"];
			int subdivision_times = fun_config_j["subdivision_times"];
			Rasterlization::PROJECTION_TYPE type = projection_type == "perspective" ? Rasterlization::PROJECTION_TYPE::PERSPECTIVE:Rasterlization::PROJECTION_TYPE::ORTHOGRAPHIC;
			VertexGetterP vg = get_vertex_getter_by_name<REAL,DIM>(fun_config_j["vertex_getter"]);
			std::string measurement_type = fun_config_j["measurement_type"];
			MeshConsistency<REAL,DIM>* measurer;

			if (measurement_type == "MaxAB") {
				measurer = new MaxAB<REAL, DIM>(vg, width, height, type);
			}
			else if (measurement_type == "OrientationConsistency") {
				measurer = new OrientationConsistency<REAL, DIM>(vg, width, height, type);
			}else{
				printf("measurement_type %s not supported\n",measurement_type.c_str());
				assert(false);
			}
			return new consistency_calculator<REAL, DIM>(factory, lg,measurer,subdivision_times);
		}else if(name == "dipole_calculator"){
			config_j = config_j["dipole_calculator_config"];
			return new dipole_calculator<REAL,DIM>(config_j["eps"],lg);
		}else if(name == "CoCalculator"){
			config_j = config_j["CoCalculator"];
			config_j = config_j["CoCalculator_config"];
			return new CoCalculator<REAL,DIM>(get_edge_caculator_from_json<REAL,DIM>(config_j["primary_calculator"],factory,lg),get_edge_caculator_from_json<REAL, DIM>(config_j["aid_calculator"],factory,lg),config_j["gap"],lg);
		}else if(name == "point_consistency_calculator"){
			return new point_consistency_calculator<REAL,DIM>(config_j["point_consistency_config"]["radius"],lg);
		}else if(name == "consistency_calculator_plus"){
			return new consistency_calculator_plus<REAL,DIM>(factory,lg);
		}else{
			assert(false);
		}
	}

	// 使用邻接链表存储的图
	template<typename REAL, int DIM>
	class linked_graph {
		std::map<VERTEX_P<REAL, DIM>, int> _vertex2idx; //_vertex2idx[_nodes[i]] = i。有时候有node指针，但不知道其在_nodes中的索引
		std::map<std::pair<int, int>, EDGE_P<REAL, DIM>> _idx2edge;

		edge_weight_calculator<REAL, DIM>* _calculator = NULL;

		nlohmann::json metric_j;
		nlohmann::json config_j;
		nlohmann::json log_j;

		// 重建顶点与索引的映射
		void _build_mapper(){
			_vertex2idx.clear();
			for(int i = 0;i<_nodes.size();i++){
				_vertex2idx[_nodes[i]] = i;
			}
			for(int i = 0;i<_nodes.size();i++){
				for(int j = 0;j<_edges[i].size();j++){
					int s = _vertex2idx[_edges[i][j]->start], e = _vertex2idx[_edges[i][j]->end];
					_idx2edge[std::pair<int,int>(s,e)] = _edges[i][j];
				}
			}
		}

		int find(std::vector<int>& parent, int x){
			if(parent[x] != x){
				parent[x] = find(parent,parent[x]);
			}
			return parent[x];
		}

		/**
		 * @brief 
		 * 初始状态下，每个顶点都是一个集合，每个集合中只有一个顶点。
		 * step1:对边随机排序，此后每次选择一条有向边，
		 * 	a.如果该边的两个顶点不在同一个集合中，则将这两个集合合并,
		 * 	  然后，将这条边及其逆边加入到最小生成树中
		 * 	b.如果该边的两个顶点在同一个集合中，则跳过这条边
		 * step2:此时的最小生成树为一个无向图，需要将其转为树
		 * 	    用visited数组记录已经访问的顶点，然后从任意一个顶点开始遍历
		 * 		对于每个顶点：
		 * 	   	   遍历顶点的边：
		 * 		      如果边的另一个顶点未被访问，则将这条边加入到最小生成树中
		 *            如果边的另一个顶点已经被访问，则跳过这条边
		 * 		树的构建结束条件为所有的顶点都被访问
		 * 			
		 * @return linked_graph 其中node部分为原始的node，edge部分为生成树的边
		 */
		linked_graph random_kruskal(){
			_build_mapper();
			std::vector<int> nodes;
			std::vector<std::pair<int,int>> sorted_edges;//排序后的边
			std::vector<std::vector<std::pair<int,int>>> temp_edges(_nodes.size());//生成树的中间结果，此时任意两个顶点之间有两条边
			
			std::vector<VERTEX_P<REAL,DIM>> res_nodes;
			std::vector<std::vector<EDGE_P<REAL,DIM>>> res_edges(_nodes.size());//生成树的边

			for(int i = 0;i<_nodes.size();i++){
				nodes.push_back(i);
				res_nodes.push_back(_nodes[i]);
				for(int j = 0;j<_edges[i].size();j++){
					int s = _vertex2idx[_edges[i][j]->start];
					int e = _vertex2idx[_edges[i][j]->end];

					double brief = _edges[i][j]->get_brief();
					brief = std::sqrt(brief);//防止brief过大,会很慢
					//global_var::config_j["graph_ipsr_config"]["linked_graph_config_random_kruskal_brief"] = "edge_brief";
					//double brief = 1
					//global_var::config_j["graph_ipsr_config"]["linked_graph_config"]["random_kruskal_brief"] = "always 1";

#ifdef _UNSERTAIN
					assert(s != e);
					assert(s == i);
#endif
					sorted_edges.push_back(std::pair<int,int>(s,e));
					for (int add = 1; add < brief; add++) {
						sorted_edges.push_back(std::pair<int, int>(s, e));
					}
				}
			}

			// 随机排序
			std::random_shuffle(sorted_edges.begin(),sorted_edges.end());

			// 并查集
			std::vector<int> parent(nodes.size());
			for(int i = 0;i<parent.size();i++){
				parent[i] = i;
			}	

			// step1：遍历边 构建生成树的中间结果
			for(int i = 0;i<sorted_edges.size();i++){
				int s = sorted_edges[i].first, e = sorted_edges[i].second;
				if(find(parent,s) != find(parent,e)){
					// 生成树的中间结果中，任意两个顶点之间有两条边
					temp_edges[s].push_back(std::pair<int,int>(s,e));
					temp_edges[e].push_back(std::pair<int,int>(e,s));
					parent[find(parent, e)] = find(parent, s);
				}
			}

			// 并查集的根节点的个数即为连通分量的个数
			int cluster_num = 0;
			std::vector<int> cluster_id(nodes.size(),0);
			for(int i=0;i<parent.size();i++){
				parent[i] = find(parent,i);
				if(cluster_id[parent[i]] == 0){
					cluster_id[parent[i]] = ++cluster_num;
				}
			}
			printf("cluster_num=%d\t", cluster_num);

			// step2：将生成树的中间结果转为生成树
			std::vector<int> visited(nodes.size(),0);
			std::queue<int> q;
			q.push(0);
			visited[0] = 1;
			while(!q.empty()){
				int cur = q.front();
				q.pop();
				for(int i = 0;i<temp_edges[cur].size();i++){
					int next = temp_edges[cur][i].second;
					// temp树中,可能有些边其实并不存在 这可能会导致原本联通的分量断连 FIXME
					if (_idx2edge.find(std::pair<int, int>(cur, next)) == _idx2edge.end()) {
						continue;
					}
					if(visited[next] == 0){
						visited[next] = 1;
						q.push(next);
						res_edges[cur].push_back(_idx2edge[std::pair<int, int>(cur, next)]);
					}
				}

				// 如果q为空,将未被访问的一个顶点(如果存在)加入到队列中.
				if (q.empty()) {
					for (int i = 0; i < nodes.size(); i++) {
						if (!visited[i]) {
							visited[i] = 1;
							q.push(i);
							break;
						}
					}
				}
			}
			return linked_graph(res_nodes,res_edges);
		}

		// 检查这个图是否是森林
		bool check_forest(){
			std::vector<int> visited(_nodes.size(),0);
			std::queue<int> q;
			q.push(0);
			visited[0] = 1;
			while(!q.empty()){
				int cur = q.front();
				q.pop();
				for(int i = 0;i<_edges[cur].size();i++){
					int next = _vertex2idx[_edges[cur][i]->end];
					if(visited[next] == 0){
						visited[next] = 1;
						q.push(next);
					}else{
						return false;
					}
				}
				if (q.empty()) {
					for (int i = 0; i < visited.size(); i++) {
						if (visited[i] == 0) {
							visited[i] = 1;
							q.push(i);
							break;
						}
					}
				}
			}
			return true;
		}

		// 撤销所有翻转 
		void _recover(){
			for(int i = 0;i<_nodes.size();i++){
				if(_nodes[i]->is_inved()){
					_nodes[i]->inv_op();
				}
			}
		}

		// 对齐两个翻转序列。比较二者的欧氏距离，如果大于一半，那么将后者翻转
		// 返回值为1表示翻转了，0表示没有翻转
		int align_flip(std::vector<int>& flip1, std::vector<int>& flip2){
			assert(flip1.size() == flip2.size());
			int dist = 0 , inv_dist = 0;
			for(int i = 0;i<flip1.size();i++){
				if(flip1[i] != flip2[i]){
					dist += _nodes[i]->cnt();
				}else{
					inv_dist += _nodes[i]->cnt();
				}
			}
			if(dist > inv_dist){
				for(int i = 0;i<flip2.size();i++){
					assert(flip2[i] == 1 || flip2[i] == 0);
					flip2[i] = 1 - flip2[i];
				}
				return 1;
			}
			return 0;
		}

		
		/**********************翻转算法 翻转node并返回翻转历史*************************************/
		/**
		 * @brief 
		 * 计算准确率。对于每个块，与_gt进行比较，
		 * 		若该翻转状态下与gt法向量差距相比另一种翻转状态下更小，则认为该翻转状态是正确的
		 * @param flag 
		 * @return float 
		 */
		float flip_acc(){
			int acc = 0;
			for(int i = 0;i<_nodes.size();i++){
				auto loss_pair = _nodes[i]->_cal_loss();
				if(loss_pair.first <= loss_pair.second){
					acc++;
				}
			}
			acc = std::max(acc, (int)_nodes.size() - acc);
			return (float)acc / _nodes.size();
		}

		/************************************翻转算法end*************************************************/

	public:

		enum FLIP_ARG{
			GREEDY,
			BRUTE_FORCE,
			SINGLE_TREE,
			MUTI_TREE,
			OPTIMAL_FILP,
			BEST_TREE
		};

		// 原始顶点
		//POINTS_NORMALS _ori_op;
		// 顶点列表
		std::vector<VERTEX_P<REAL,DIM>> _nodes;
		// 边列表的列表，每个顶点对应一个边列表
		std::vector<std::vector<EDGE_P<REAL, DIM>>> _edges;

		linked_graph(
			std::vector<VERTEX_P<REAL,DIM>> nodes, 
			std::vector<std::vector<EDGE_P<REAL,DIM>>> edges
			){
			assert(nodes.size() == edges.size());
			_nodes = nodes;
			_edges = edges;
			_build_mapper();
		}

		linked_graph(){
			_nodes = std::vector<VERTEX_P<REAL, DIM>>();
			_edges = std::vector<std::vector<EDGE_P<REAL, DIM>>>();
		}

		void update(IpsrController<REAL, DIM>* controller = NULL){
#pragma omp parallel for
			for (int i = 0; i < _nodes.size(); i++) {
				_nodes[i]->iter(controller);
			}
		}

		void parallel_update(int iter,IpsrController<REAL, DIM>* controller = NULL){
			std::mutex *mtx = new std::mutex[_nodes.size()];
			lzd_tools::thread_safe_int finished(0);
			int total_count = iter * _nodes.size();
			printf("start parallel update. %d nodes has total_iter%d\n",_(int)nodes.size(),total_count);
#pragma omp parallel for
			for(int i = 0;i<iter;i++){
				for (int j = 0; j < _nodes.size(); j++) {
					mtx[j].lock();
					_nodes[j]->iter(controller);
					mtx[j].unlock();
					++finished;
					printf("\r%d/%d", finished.get(), total_count);
				}
			}
			printf("\nparallel update done, saved %d iters by early stop\n", total_count - finished.get());
		}
		

		void set_edge_weight_calculator(edge_weight_calculator<REAL, DIM>* calculator){
			_calculator = calculator;
			config_j["edge_weight_calculator"] = calculator->get_config();
		}

		void update_all_op(){
			printf("update2op...\t");
#pragma omp parallel for
			for (int i = 0; i < _nodes.size();i++) {
				_nodes[i]->update2op();
			}
			printf("update2op end!\n");
		}

		/**
		 * @brief 
		 * 计算所有边的权重
		 */
		void cal_edges_weight() {
			assert(_calculator != NULL);
			printf("start calculateing weight...\n");
			lzd_tools::thread_safe_int cnt(0);
			int total_edge = 0;
			for (int i = 0; i < _nodes.size(); i++) {
				total_edge += _edges[i].size();
			}
#pragma omp parallel for
			for (int i = 0; i < _edges.size(); i++) {
				for (int j = 0; j < _edges[i].size(); j++) {
					std::pair<REAL, REAL> weight;
					REAL conf = 0;
					_calculator->cal_weight(_edges[i][j], weight,conf);
					if (weight.first == 0 || weight.second == 0) {
						printf("weight %d, %d", weight.first,weight.second);
					}
					if (conf < 0.5) {
						printf("conf=%f\n", conf);
						printf("weight=%f,%f\n", weight.first, weight.second);
						assert (false);
					}
					if(static_conf["linked_graph"]["use_strength"] == true) {
						weight.first = weight.first * _edges[i][j]->get_connection_strength();
						weight.second = weight.second * _edges[i][j]->get_connection_strength();
					}
					_edges[i][j]->set_weight(weight,conf);
					++cnt;
					printf("\rcal weight %d/%d ", cnt.get(), total_edge);
				}
			}
			printf("calculate weight done\n");
		}

		int cal_w_count(const std::vector<int>& flag) {
			assert(flag.size() == _nodes.size());
			int sum = 0;
			for (int i = 0; i < _nodes.size(); i++) {
				for (int j = 0; j < _edges[i].size(); j++) {
					int s = i, t = _vertex2idx[_edges[i][j]->end];
					auto weight = _edges[i][j]->get_weight(flag[s] == flag[t]);
					if (weight.first > weight.second)sum++;
				}
			}
			return sum;
		}

		REAL cal_all_weight_sum_according_to_flag(const std::vector<int>& flag){
			assert(flag.size() == _nodes.size());
			REAL sum = 0;
			for (int i = 0; i < _nodes.size(); i++) {
				for (int j = 0; j < _edges[i].size(); j++) {
					int s = i, t = _vertex2idx[_edges[i][j]->end];
					auto weight = _edges[i][j]->get_weight(flag[s] == flag[t]);
					sum += weight.first;
				}
			}
			return sum;
		}

		// 求所有边的权重和（要求所有边的权重已经计算）
		REAL cal_all_weight_sum() {
			REAL sum = 0;
			for (int i = 0; i < _nodes.size(); i++) {
				for (int j = 0; j < _edges[i].size(); j++) {
						sum += _edges[i][j]->get_weight().first;
				}
			}
			return sum;
		}

		std::vector<int> flip_graph(graph_arg::FlipGraph* fp){
			_recover();
			using namespace graph_arg;
			std::vector<int> res;
			graph_arg::FlipableGraph g(to_json());
			auto ttt = fp->flip(&g);
			res.resize(ttt.size());
			for (int i = 0; i < ttt.size(); i++)res[i] = ttt[i];
			nlohmann::json j;
			printf("weight_sum=%f \t edge_acc=%f \t flip_acc = %f \n", g.cal_current_weight_sum(), cal_edge_acc(g._edges), cal_flip_acc(g._nodes));
			j["config"] = fp->get_config();
			j["log"] = fp->get_log();
			log_j["flip_log"].push_back(j);
			for (int i = 0; i < _nodes.size(); i++) {
				if (res[i] != _nodes[i]->is_inved()) {
					_nodes[i]->inv_op();
				}
			}
			return res;
		}

		// comb: 将所有的顶点的点云合并成一个大的点云 
		// seg_id: 记录每个点云的来源
		POINTS_NORMALS comb(std::vector<int>& seg_id) {
			POINTS_NORMALS res;
			seg_id.clear();
			for(int i = 0;i<_nodes.size();i++){
				for(int j = 0;j<_nodes[i]->_op.size();j++){
					if (_nodes[i]->_op_type[j] == 0) {
						res.push_back(_nodes[i]->_op[j]);
						seg_id.push_back(i);
					}
				}
			}
			return res;
		}

		void comb(POINTS_NORMALS& op, POINTS_NORMALS& gt, std::vector<int>& seg_id) {
			op.clear();
			gt.clear();
			seg_id.clear();
			for (int i = 0; i < _nodes.size(); i++) {
				auto opi = _nodes[i]->get_self_op();
				auto gti = _nodes[i]->get_self_gt();
				for (int j = 0; j < opi.size(); j++) {
					op.push_back(opi[j]);
					gt.push_back(gti[j]);
					seg_id.push_back(i);
				}

			}
		}

		/***********************************visualization && log && metric*************************************/
		double edge_acc() {
			int edge_num = 0, correct_edge = 0;
			for (int i = 0; i < _nodes.size(); i++) {
				for (int j = 0; j < _edges[i].size(); j++) {
					edge_num++;
					if (!_edges[i][j]->is_wrong_edge()) {
						correct_edge++;
					}
				}
			}
			return double(correct_edge) / double(edge_num);	
		}

		/**
		 * @brief 
		 * // 将图转为json格式
		 * {
		 * 	"graph":[
		 * 		{
		 * 			"vertex":{...},
		 * 			"edges":[
		 * 				{...},...
		 * 			]
		 * 		},...	
		 * 	]
		 * }
		 * @return * nlohmann::json 
		 */
		nlohmann::json to_json(){
			nlohmann::json _j;
			_j["graph"] = nlohmann::json::array();
			_j["graph_metric"]["edge_acc"] = edge_acc();
			_j["graph_metric"]["flip_acc"] = flip_acc();
			for(int i = 0;i<_nodes.size();i++){
				nlohmann::json node;
				node["vertex"] = _nodes[i]->to_json();
				node["vertex"]["id"] = i;
				node["edges"] = nlohmann::json::array();
				for(int j = 0;j<_edges[i].size();j++){
					auto temp = _edges[i][j]->to_json();
					assert(i == _vertex2idx[_nodes[i]]);
					temp["start"] = i;
					temp["end"] = _vertex2idx[_edges[i][j]->end];
					node["edges"].push_back(temp);
				}
				_j["graph_topology"].push_back(node);
			}
			return _j;
		}

		std::string to_string(){
			return to_json().dump();
		}

		/**
		 * @brief 
		 * 打印图的拓扑结构
		 * @param adpath 例如"/graph_ipsr/"
		 * @param adname 用于给文件命名
		 */
		void print_topology(std::string adpath = "",std::string adname = "") {
			std::string path = _nodes[0]->_handle->get_out_put_base(adpath) + adname + "graph_topology.json";
			std::ofstream out(path);
			out << to_string();
			out.close();
		}

		nlohmann::json get_log() {
			nlohmann::json node_log;
			for (int i = 0; i < _nodes.size(); i++) {
				nlohmann::json nd;
				nd["id"] = i;
				nd["log"] = _nodes[i]->get_log();
				nd["metric"] = _nodes[i]->cal_metric().to_json();
				node_log.push_back(nd);
			}
			// 根据j["metric"]["avg_nd_loss"]对j["nodes"]进行排序
			std::sort(node_log.begin(), node_log.end(), [](nlohmann::json a, nlohmann::json b) {
				return a["metric"]["avg_nd_loss"].get<double>() < b["metric"]["avg_nd_loss"].get<double>();
				});
			log_j["node_log"] = node_log;
			return log_j;
		}
		
		nlohmann::json get_metric() {
			return metric_j;
		}

		nlohmann::json get_config() {
			return config_j;
		}

		/**
		 * @brief 
		 * 打印图中,各个节点的结果
		 */
		void print_log(std::string adpath = "",std::string adname = ""){
			std::string path = _nodes[0]->_handle->get_out_put_base(adpath) + adname + "graph_log.json";
			nlohmann::json j = get_log();
			std::ofstream out(path);
			out << j.dump(4);
			out.close();		
		}

		void draw_mesh(std::string adpath = "",std::string adname = ""){
			MESH mesh;
			using COLOR_NAME = lzd_tools::beautifulColorMapper<int>::COLOR_NAME;

			// 将所有顶点加入到mesh中
			std::vector<int> seg_id;
			POINTS_NORMALS op = comb(seg_id);
			// 打印op
			std::string path = _nodes[0]->_handle->get_out_put_base(adpath) + adname + "graph_op.ply";
			lzd_tools::InFiniteMap<int>* colormap = new lzd_tools::RandomColorMapper<int>(seg_id);
			lzd_tools::op2ply(op, path, XForm<REAL, DIM + 1>().Identity(), colormap);

			//for(int i = 0;i<op.size();i++){
			//	mesh.first.push_back(op[i].first);
			//}
			
			// 所有点的类别设为0
			seg_id = std::vector<int>();

			// 每个顶点上绘制一个球
			std::vector<Point<REAL,DIM>> centers;
			for(int i = 0;i<_nodes.size();i++){
				Point<REAL,DIM> c = _nodes[i]->center_op;
				auto sphere = lzd_tools::get_sphere(c,0.01);
				lzd_tools::add_topology(mesh, sphere);
				//int node_status = _nodes[i]->inved_status_according2gt()==0?COLOR_NAME::GRAY_BLUE:COLOR_NAME::SEAWEED_GREEN;
				int node_status = _nodes[i]->inved_status_according2gt()==0?0:1;
				for(int j = 0;j<sphere.first.size();j++){
					seg_id.push_back(node_status);
				}				
				centers.push_back(c);			
			}
			// 将所有的边加入到mesh中
			for(int i = 0;i<_nodes.size();i++){
				for(int j = 0;j<_edges[i].size();j++){
					int edge_status = _edges[i][j]->is_wrong_edge() ? 2:3; 
					int s = i, e = _vertex2idx[_edges[i][j]->end];
					auto arrow = lzd_tools::get_arrow(centers[s],centers[e],0.02);					
					for(int k = 0;k<arrow.first.size();k++){
						seg_id.push_back(edge_status);
					}
					lzd_tools::add_topology(mesh, arrow);
				}
			}
			path = _nodes[0]->_handle->get_out_put_base(adpath) + adname + "graph_topology.ply";
			lzd_tools::InFiniteMap<int>* graphcolormap = new lzd_tools::beautifulColorMapper<int>(seg_id);
			lzd_tools::mesh2ply(mesh,path,XForm<REAL,DIM+1>().Identity(),graphcolormap);
		}

		void draw_edge(std::string adpath = "", std::string adname = "",std::string measurement = "brief") {
			MESH all_edges;
			std::vector<REAL> edge_metric;
			// 所有的点加入到mesh中
			std::vector<Point<REAL, DIM>> centers;
			for (int i = 0; i < _nodes.size(); i++) {
				Point<REAL, DIM> c = _nodes[i]->center_op;
				auto sphere = lzd_tools::get_sphere(c, 0.01);
				lzd_tools::add_topology(all_edges, sphere);
				for (int j = 0; j < sphere.first.size(); j++) {
					edge_metric.push_back(_nodes[i]->get_vertex_conf());
				}
				centers.push_back(c);
			}
			// 将所有的边加入到mesh中
			for(int i = 0;i<_nodes.size();i++){
				for(int j = 0;j<_edges[i].size();j++){
					int s = i, e = _vertex2idx[_edges[i][j]->end];
					auto arrow = lzd_tools::get_arrow(_nodes[s]->center_op,_nodes[e]->center_op,0.02);
					lzd_tools::add_topology(all_edges, arrow);
					for(int k = 0;k<arrow.first.size();k++){
						if (measurement == "brief")edge_metric.push_back(_edges[i][j]->get_brief());
						else if (measurement == "connection_strength")edge_metric.push_back(_edges[i][j]->get_connection_strength());
					}
				}
			}
			std::string path = _nodes[0]->_handle->get_out_put_base(adpath) + adname + "graph_edge_" + measurement + ".ply";
			lzd_tools::InFiniteMap<REAL>* colormap = new lzd_tools::GradientMapper<REAL>(edge_metric);
			lzd_tools::mesh2ply(all_edges,path,XForm<REAL,DIM+1>().Identity(),colormap);
		}

		lzd_tools::PointNormalMetric<REAL,DIM> cal_metric(){
			POINTS_NORMALS gt;
			POINTS_NORMALS op;
			for(int i = 0;i<_nodes.size();i++){
				auto opi = _nodes[i]->get_self_op();
				auto gti = _nodes[i]->get_self_gt();
				gt.insert(gt.end(),gti.begin(),gti.end());
				op.insert(op.end(),opi.begin(),opi.end());
			}
			return lzd_tools::PointNormalMetric<REAL,DIM>(gt,op);
		}

		void query_res(const POINTS_NORMALS& gt, POINTS_NORMALS& res_op) {
			POINTS_NORMALS self_op;
			get_all_op(self_op);
			assert(gt.size() == res_op.size());
			// self_op's kdtree
			vector<kdt::KDTreePoint> vertices;
			for (int i = 0; i < self_op.size(); i++) {
				array<REAL, 3> p = { self_op[i].first[0],self_op[i].first[1],self_op[i].first[2] };
				vertices.push_back(kdt::KDTreePoint(p));
			}
			kdt::KDTree<kdt::KDTreePoint> kdtree(vertices);

			double total_dist = 0;
#pragma omp parallel for
			for (int i = 0; i < gt.size(); i++) {
				array<REAL, 3> p = { gt[i].first[0],gt[i].first[1],gt[i].first[2] };
				auto res = kdtree.knnSearch(p, 1);
				res_op[i].first = gt[i].first;
				res_op[i].second = self_op[res[0]].second;
				total_dist += Distance(gt[i].first, self_op[res[0]].first);
			}
			printf("query_test::dist between gt and self_op = %f\n", total_dist);
		}
		void save_res(std::string adpath, std::string adname, IPSR_Factory<REAL,DIM>* adhandle){
			nlohmann::json save_option = ConfigManager::get_common_config()["save_option"];
			std::string path = adhandle->get_out_put_base(adpath) + adname;
			POINTS_NORMALS op, gt;
			std::vector<int> seg_id;
			comb(op,gt,seg_id);
			std::string op_path = path + adhandle->get_resname() + "_orientedpoint.ply";
			std::string gt_path = path + adhandle->get_resname() + "_GT_samples.ply";
			lzd_tools::InFiniteMap<int>* colormap = new lzd_tools::RandomColorMapper<int>(seg_id);
			if(save_option["orientedpoint"] == true){
				lzd_tools::op2ply(op,op_path,XForm<REAL,DIM+1>().Identity(),colormap);
			}
			if(save_option["GT_samples"] == true){
				lzd_tools::op2ply(gt,gt_path,XForm<REAL,DIM+1>().Identity(),colormap);
			}

			if(save_option["cleaned_mesh"] == true){
				assert(adhandle != NULL);
				std::string surface_path = path + adhandle->get_resname() + "_surface.ply";
				auto mesh = adhandle->get_mesh_from_points_normals(op);
				std::vector<double> pmdist;
				auto cleanedmesh = shrink_boundary(mesh, op, pmdist,20);
				output_ply(surface_path, cleanedmesh, XForm<REAL, DIM + 1>().Identity());
			}
			print_topology(adpath+"/graph_ipsr/",adname);
			draw_edge(adpath+"/graph_ipsr/",adname,"brief");
			draw_edge(adpath+"/graph_ipsr/",adname,"connection_strength");
			draw_mesh(adpath+"/graph_ipsr/",adname);
		}

		void get_all_op(POINTS_NORMALS& op){
			op.clear();
			for(int i = 0;i<_nodes.size();i++){
				auto opi = _nodes[i]->get_self_op();
				op.insert(op.end(),opi.begin(),opi.end());
			}
		}

		void save_all_edge(std::string adpath = "/graph_ipsr/all_edge/", std::string adname = "") {
			assert(_calculator != NULL);
			std::string base_floder = _nodes[0]->_handle->get_out_put_base(adpath);
			rmkdir(adpath);
#pragma omp parallel for
			for (int i = 0; i < _edges.size(); i++) {
				for (int j = 0; j < _edges[i].size(); j++) {
					auto e = _edges[i][j];
					int s = _vertex2idx[e->start], t = _vertex2idx[e->end];
					printf("edge (%d,%d), weigth = (%f,%f)n", s, t, e->get_weight().first, e->get_weight().second);
					std::string floder = base_floder + adname + std::to_string(s) + "_" + std::to_string(t) + "/";
					rmkdir(floder);
					_calculator->visualization(e, floder);
					
				}
			}
			cout << "save all edge done!\n";
		}
		/**
		 * @brief 
		 * 绘制错误的边(即无法从权重中判断是否对齐)
		 * @param calculator  
		 */
		void save_wrong_edge(std::string adpath = "/graph_ipsr/wrong_edge/", std::string adname = ""){
			assert(_calculator != NULL);
			std::string base_floder = _nodes[0]->_handle->get_out_put_base(adpath);
			rmkdir(adpath);

#pragma omp parallel for
			for(int i = 0;i<_edges.size();i++){
				for(int j = 0;j<_edges[i].size();j++){
					auto e = _edges[i][j];
					if(e->is_wrong_edge()){
						int s = _vertex2idx[e->start], t = _vertex2idx[e->end];
						printf ("wrong edge (%d,%d), weigth = (%f,%f)n", s, t, e->get_weight().first, e->get_weight().second);
						std::string floder = base_floder + adname +  std::to_string(s) + "_" + std::to_string(t) + "/";
						rmkdir(floder);
						_calculator->visualization(e,floder);
					}
				}
			}
			cout << "save wrong edge done!\n";
		}
		
		void save_optim(std::string adpath = "/graph_ipsr/", std::string adname = "optimflip") {
			// save optimal flip
			std::string path = _nodes[0]->_handle->get_out_put_base(adpath) + adname + ".ply";
			POINTS_NORMALS op, gt;
			std::vector<int> seg_id;
			auto optim_filp = flip_graph(new graph_arg::OptimFlip());
			comb(op, gt, seg_id);
			lzd_tools::InFiniteMap<int>* colormap = new lzd_tools::RandomColorMapper<int>(seg_id);
			lzd_tools::op2ply(op, path, XForm<REAL, DIM + 1>().Identity(), colormap);
		}
		
		// 注意node的get_resname()和实际输出的是不一样的。
		void load_seg(std::string self_arg_flag,std::string refer_arg_flag, std::string adpath = "/graph_ipsr/lg_seg/", std::string adname = "res") {
			double total_dist = 0, max_dist=0, total_count = 0;
			std::string floder = _nodes[0]->_handle->get_out_put_base(adpath);
			floder = floder.replace(floder.find(self_arg_flag), self_arg_flag.size(), refer_arg_flag);
			printf("load seg from %s\n", floder.c_str());
			// 检查文件夹是否存在
			if(GetFileAttributesA(floder.c_str()) == INVALID_FILE_ATTRIBUTES){
				printf("floder %s not exist\n", floder.c_str());
				assert(false);
			}

			std::string log_path = floder + adname + "log.json";
			for (int i = 0; i < _nodes.size(); i++) {
				std::string path = floder + adname + std::to_string(i) + ".ply";
				// 读取ply文件
				POINTS_NORMALS op;
				ply_reader<REAL, DIM>(path, op);
				if (op.size() != _nodes[i]->_op.size()) {
					printf("op size not match, %d, %d\n", op.size(), _nodes[i]->_op.size());
					assert(false);
				}
				for (int j = 0; j < op.size(); j++) {
					double d = Distance(op[j].first, _nodes[i]->_op[j].first);
					_nodes[i]->_op[j].second = op[j].second;
					total_count += 1;
					max_dist = max(d, max_dist);
					total_dist += d;
				}
				_nodes[i]->iupdate2op();
			}
			printf("load seg done, total_count = %f,total_dist = %f, avg_dist = %f, max_dist = %f\n", total_count, total_dist, total_dist/total_count,max_dist);
		}


		/**
		 * @brief 
		 * 存储当前每个顶点的结果
		 * @param adpath 
		 * @param adname 例如结束阶段调用时,可以设置为"res"  
		 */
		void save_seg(std::string adpath = "/graph_ipsr/lg_seg/", std::string adname = "") {
			nlohmann::json nodes_j;
			for(int i = 0;i<_nodes.size();i++){
				nlohmann::json j = _nodes[i]->to_json();
				j["id"] = _vertex2idx[_nodes[i]];
				nodes_j.push_back(j);
			}
			// 保存seg
			std::string floder = _nodes[0]->_handle->get_out_put_base(adpath);
			rmkdir(floder);
			for(int i = 0;i<_nodes.size();i++){
				std::string path = floder + adname + std::to_string(_vertex2idx[_nodes[i]]) + ".ply";
				lzd_tools::op2ply(_nodes[i]->_op,path,XForm<REAL,DIM+1>().Identity(),std::make_pair(_nodes[i]->_op_type, lzd_tools::get_regular_colormap(0, 1)));
			}
			// 保存log
			std::string log_path = floder + adname + "log.json";
			std::ofstream out(log_path);
			out << get_log().dump(4);
			out.close();	
		}

		// TODO
		void save_gt(std::string adpath = "/graph_ipsr/lg_seg/", std::string adname = "gt") {
			// 保存seg
			std::string floder = _nodes[0]->_handle->get_out_put_base(adpath);
			rmkdir(floder);
			for (int i = 0; i < _nodes.size(); i++) {
				std::string path = floder + adname + std::to_string(_vertex2idx[_nodes[i]]) + ".ply";
				lzd_tools::op2ply(_nodes[i]->get_all_gt(), path, XForm<REAL, DIM + 1>().Identity(), std::make_pair(_nodes[i]->_op_type, lzd_tools::get_regular_colormap(0, 1)));
			}
		}

}; 

/******************************************顶点的组织********************************************************/
	// 可能考虑将顶点用简单的立方体，或者octree来组织。这部分完全面向顶点。但每个组织方式，需要提供一个转为图的接口。
	template<typename REAL, int DIM>
	class NodeAllocator {
	public:
		virtual linked_graph<REAL, DIM> op_orchestration(
			IPSR_Factory<REAL, DIM> ipsr_factory, // 用于从op生成ipsr
			int (*_update_plan)(Period& p), //ipsr更新方式
			NormalEstimation<REAL, DIM>* estimator //用来初始化法向量
		) = 0;

		virtual nlohmann::json get_config() = 0;

	};

	// 将立方体分割成一个立方体网格,每个小格为一个正方体
	template<typename REAL, int DIM>
	class simple_cube : public NodeAllocator<REAL,DIM>{
		typedef std::vector<int> CUBE_COORDINATE;
		// 立方体网格中的一个立方体的边界(矩形边界)
		typedef struct rectangle_bbox {
			static REAL INF;
			std::pair<REAL, REAL> boundaries[DIM];//每个维度的上/下界,first为下界, second为上界
			rectangle_bbox(){
				for(int d = 0;d<DIM;d++){
					boundaries[d].first = rectangle_bbox::INF;
					boundaries[d].second = -1*rectangle_bbox::INF;
				}
			}
			rectangle_bbox(POINTS_NORMALS op){
				for(int i = 0;i<op.size();i++){
					update(op[i].first);
				}
			}
			void update(Point<REAL,DIM> p){
				for(int d = 0;d<DIM;d++){
					if(p[d] < boundaries[d].first)boundaries[d].first = p[d];
					if(p[d] > boundaries[d].second)boundaries[d].second = p[d];
				}
			}

			bool check_isinside(Point<REAL, DIM> p) {
				// 允许一定误差
				for (int i = 0; i < DIM; i++) {
					REAL gap = 0.0001 * (boundaries[i].second - boundaries[i].first);
					if (p[i] < boundaries[i].first  - gap || p[i] > boundaries[i].second + gap) return false;
				}
				return true;
			}

		}rectangle_bbox;
		
		// 立方体网格中的一个单元
		typedef struct simple_cube_grid{
			static REAL overlap_rate;
			// 构造函数参数
			POINTS_NORMALS grid_op;//该小格内的点云
			rectangle_bbox op_bbox;//grid_op的边界
			rectangle_bbox grid_bbox;//该小格的边界
			
			int num;//该小格内的点云数量
			std::vector<CUBE_COORDINATE> neighbers;//邻居的3维坐标
			//std::vector<std::vector<int>> extand_overlap;//与邻居的重合部分在邻居的grid_op中的索引
			
			//first: 该维度上的边界小于overlap_rate的点; second: 该维度上的边界大于1-overlap_rate的点
			std::vector<std::pair<std::vector<int>, std::vector<int>>> border_idxs;

			simple_cube_grid(POINTS_NORMALS op, rectangle_bbox grid_bbox){
				grid_op = op;
				num = op.size();
				op_bbox = rectangle_bbox(op);
				this->grid_bbox = grid_bbox;
				neighbers.clear();

				// 计算border_idxs
				border_idxs.clear();
				for(int d = 0;d<DIM;d++){
					std::vector<int> pre;
					std::vector<int> latter;
					REAL pre_bound = grid_bbox.boundaries[d].first + overlap_rate * (grid_bbox.boundaries[d].second - grid_bbox.boundaries[d].first);
					REAL latter_bound = grid_bbox.boundaries[d].second - overlap_rate * (grid_bbox.boundaries[d].second - grid_bbox.boundaries[d].first);
					for(int i = 0;i<op.size();i++){
						if(op[i].first[d] < pre_bound){
							pre.push_back(i);
						}
						else if(op[i].first[d] > latter_bound){
							latter.push_back(i);
						}
					}
					border_idxs.push_back(std::pair<std::vector<int>, std::vector<int>>(pre, latter));
				}
			}
			// 判断一个点是否在该grid内部	
			bool check_isinside(Point<REAL,DIM> p){
				return grid_bbox.check_isinside(p);
			}
			// 判断该grid是否为有效的
			bool is_vaild(){
				return grid_op.size() > simple_cube<REAL, DIM>::_op_num_treshold;
			}
		}GRID;

		static int _op_num_treshold;//每个立方体中的点云数量阈值 如果小于这个阈值，则不考虑
		
		// 由构造函数提供
		int min_grid;//指定后，每个维度上至少会切分min_grid刀；_stride = (最短的维度+eps）/min_grid
		int max_grid;//TODO 没法保证同时满足min_grid和max_grid
		POINTS_NORMALS _op;//原始点云
		
		// 在_spilt_op阶段计算
		rectangle_bbox _op_bbox;//原始点云的边界
		REAL _stride;// 划分的步长。由于是立方体，所以只需要一个步长
		T_CUBE<GRID*> _grids;//立方体网格
		int x_num, y_num, z_num;//立方体网格各个维度的长度


		void _spilt_op(){
			_op_bbox = rectangle_bbox(_op);
			//计算步长
			_stride = 0.0;
			REAL min_stride = 1e9;
			for(int i = 0;i<DIM;i++){
				REAL stride = _op_bbox.boundaries[i].second - _op_bbox.boundaries[i].first;
				if(stride < min_stride) min_stride = stride;
			}
			_stride = min_stride / min_grid;
			// 将包围盒的下界减少_stride,使得所有含有点云的grid的索引都是正数
			for(int i = 0;i<DIM;i++){
				_op_bbox.boundaries[i].first -= (_stride*1.001);
			}
			//计算网格数量
			x_num = floor((_op_bbox.boundaries[0].second - _op_bbox.boundaries[0].first) / _stride) + 1;
			y_num = floor((_op_bbox.boundaries[1].second - _op_bbox.boundaries[1].first) / _stride) + 1;
			z_num = floor((_op_bbox.boundaries[2].second - _op_bbox.boundaries[2].first) / _stride) + 1;

			T_CUBE<POINTS_NORMALS> _grid_op;
			T_CUBE<rectangle_bbox> _grid_bboxs;// 各个grid的包围盒
			
			_grid_op.resize(x_num);
			_grid_bboxs.resize(x_num);
			for (int i = 0; i < x_num; i++) {
				_grid_op[i].resize(y_num);
				_grid_bboxs[i].resize(y_num);
				for (int j = 0; j < y_num; j++) {
					_grid_op[i][j].resize(z_num);
					_grid_bboxs[i][j].resize(z_num);
					for(int k = 0;k<z_num;k++){
						_grid_op[i][j][k] = POINTS_NORMALS();
						_grid_bboxs[i][j][k] = rectangle_bbox();
						unsigned int tdim[3] = { i,j,k };
						for(int l = 0;l<DIM;l++){
							_grid_bboxs[i][j][k].boundaries[l].first = _op_bbox.boundaries[l].first + tdim[l] * _stride;
							_grid_bboxs[i][j][k].boundaries[l].second = _op_bbox.boundaries[l].first + (tdim[l] + 1) * _stride;
						}
					}
				}
			}
			for(int i = 0; i < _op.size(); i++) {
				int x_index = floor((_op[i].first[0] - _op_bbox.boundaries[0].first) / _stride);
				int y_index = floor((_op[i].first[1] - _op_bbox.boundaries[1].first) / _stride);
				int z_index = floor((_op[i].first[2] - _op_bbox.boundaries[2].first) / _stride);
				assert(x_index > 0 && x_index < x_num);
				assert(y_index > 0 && y_index < y_num);
				assert(z_index > 0 && z_index < z_num);
				assert(_grid_bboxs[x_index][y_index][z_index].check_isinside(_op[i].first));
				_grid_op[x_index][y_index][z_index].push_back(_op[i]);
			}

			cube_resize(_grids, x_num, y_num, z_num);
			for (int i = 0; i < x_num; i++) {
				for (int j = 0; j < y_num; j++) {
					for (int k = 0; k < z_num; k++) {
						_grids[i][j][k] = new GRID(_grid_op[i][j][k], _grid_bboxs[i][j][k]);
					}  
				}
			}
			for (int i = 0; i < x_num; i++) {
				for (int j = 0; j < y_num; j++) {
					for (int k = 0; k < z_num; k++) {
						_grids[i][j][k]->neighbers = _get_vaild_neighbor_3idx({ i,j,k });
					}
				}
			}
		}

		// 返回i，j，k坐标的邻居的三维坐标
		std::vector<CUBE_COORDINATE> _get_neighbor_3idx(CUBE_COORDINATE _3idx){
			int i = _3idx[0], j = _3idx[1], k = _3idx[2];
			assert(i >= 0 && i < x_num);
			assert(j >= 0 && j < y_num);
			assert(k >= 0 && k < z_num);
			std::vector<CUBE_COORDINATE> neighbors;
			neighbors.clear();
			if(i - 1 >= 0) neighbors.push_back(CUBE_COORDINATE({i - 1, j, k}));
			if(i + 1 < x_num) neighbors.push_back(CUBE_COORDINATE({i + 1, j, k}));
			if(j - 1 >= 0) neighbors.push_back(CUBE_COORDINATE({i, j - 1, k}));
			if(j + 1 < y_num) neighbors.push_back(CUBE_COORDINATE({i, j + 1, k}));
			if(k - 1 >= 0) neighbors.push_back(CUBE_COORDINATE({i, j, k - 1}));
			if(k + 1 < z_num) neighbors.push_back(CUBE_COORDINATE({i, j, k + 1}));
			return neighbors;
		}

		// 返回i，j，k坐标的合法邻居的三维坐标
		std::vector<CUBE_COORDINATE> _get_vaild_neighbor_3idx(CUBE_COORDINATE _3idx){
			std::vector<CUBE_COORDINATE> neighbors = _get_neighbor_3idx(_3idx);
			std::vector<CUBE_COORDINATE> vaild_neighbors;
			vaild_neighbors.clear();
			for(int i = 0; i < neighbors.size(); i++){
				if(_grids[neighbors[i][0]][neighbors[i][1]][neighbors[i][2]]->is_vaild()) {
					vaild_neighbors.push_back(neighbors[i]);
				}
			}
			return vaild_neighbors;
		}
		
		/**
		 * @brief 
		 * 完成grid到vertex以及edge的一些准备工作。
		 * @param op vertex的点云
		 * @param flag vertex中的点的属性。0表示来自grid，1表示来自邻居grid
		 * @param correspon 1点一定会同时出现在两个vertex中，需要维护一一对应关系以便计算overlap。因此是一个pair数组的列表
		 * 					每个pair数组记录了1点的其中一个来源
		 * 					数组中的每个pair记录了1点在本身中的索引（first）以及在邻居中的索引（second）。
		 * 					
		 */
		void _grid2op(CUBE_COORDINATE _3idx, POINTS_NORMALS& op, std::vector<int>& flag, std::vector<std::vector<std::pair<int,int>>>& correspon){
			op.clear();
			flag.clear();
			GRID* grid = _grids[_3idx[0]][_3idx[1]][_3idx[2]];
			// 复制op
			for(int i = 0;i<grid->grid_op.size();i++){
				op.push_back(grid->grid_op[i]);
				flag.push_back(0);
			}
			correspon.clear();

			for(auto n:grid->neighbers){
				assert(_grids[n[0]][n[1]][n[2]]->is_vaild());
				POINTS_NORMALS& neibor_op = _grids[n[0]][n[1]][n[2]]->grid_op;
				std::vector<int> idx;//点在邻居中的索引
				for(int d = 0;d<DIM;d++){
					if(_3idx[d]!=n[d]){
						assert(_3idx[(d+1)%3] == n[(d+1)%3] && _3idx[(d+2)%3] == n[(d+2)%3]);//另外两个维度相同
						if(_3idx[d] == n[d] - 1){
							idx = _grids[n[0]][n[1]][n[2]]->border_idxs[d].first;
						}else if(_3idx[d] == n[d] + 1){
							idx = _grids[n[0]][n[1]][n[2]]->border_idxs[d].second;
						}else{
							assert(0);
						}
					}
				}
				std::vector<std::pair<int,int>> idx_pair;
				for(int i = 0;i<idx.size();i++){
					op.push_back(neibor_op[idx[i]]);
					flag.push_back(1);
					idx_pair.push_back(std::make_pair(op.size()-1,idx[i])); 
				}
				correspon.push_back(idx_pair);
			}
		}


public:
		simple_cube(const POINTS_NORMALS op, int min_grid = 4):_op(op),min_grid(min_grid){
			_spilt_op();
		}

		nlohmann::json get_config(){
			nlohmann::json j;
			j["name"] = "simple_cube";
			j["min_grid"] = min_grid;
			return j;
		}

		/**
		 * @brief 
		 * 
		 * @return linked_graph<REAL, DIM> 
		 */
		linked_graph<REAL, DIM> op_orchestration(
			IPSR_Factory<REAL, DIM> ipsr_factory,
			int (*_update_plan)(Period& p),
			NormalEstimation<REAL,DIM>* estimator	
		){
			linked_graph<REAL, DIM> ihp_graph;
			T_CUBE<VERTEX_P<REAL,DIM>> vertexp_cube;
			T_CUBE<std::vector<EDGE_P<REAL, DIM>>> edgelist_cube;

			//初始化顶点和边的指针数组
			vertexp_cube.resize(x_num);
			for(int i = 0;i<x_num;i++){
				vertexp_cube[i].resize(y_num);
				for(int j = 0;j<y_num;j++){
					vertexp_cube[i][j].resize(z_num);
					for(int k = 0;k<z_num;k++){
						vertexp_cube[i][j][k] = new VERTEX<REAL,DIM>();
					}
				}
			}
			edgelist_cube.resize(x_num);
			for(int i = 0;i<x_num;i++){
				edgelist_cube[i].resize(y_num);
				for(int j = 0;j<y_num;j++){
					edgelist_cube[i][j].resize(z_num);
					for(int k = 0;k<z_num;k++){
						if(!_grids[i][j][k]->is_vaild())continue;
						for(auto n:_grids[i][j][k]->neighbers){
							assert(_grids[n[0]][n[1]][n[2]]->is_vaild());	
							auto edgep = new EDGE<REAL, DIM>(vertexp_cube[i][j][k], vertexp_cube[n[0]][n[1]][n[2]]);
							edgelist_cube[i][j][k].push_back(edgep);
						}
					}
				}
			}			

			// 填充vertex和edge
			for (int i = 0; i < x_num; i++) {
				for (int j = 0; j < y_num; j++) {
					for (int k = 0; k < z_num; k++) {
						if(!_grids[i][j][k]->is_vaild())continue;
						// 填充vertex
						POINTS_NORMALS tempop;
						std::vector<std::vector<std::pair<int,int>>> tempcor;
						_grid2op({i,j,k},tempop,vertexp_cube[i][j][k]->_op_type,tempcor);
						vertexp_cube[i][j][k]->_handle = ipsr_factory.create_ipsr_from_op(tempop,_update_plan,estimator);
						vertexp_cube[i][j][k]->_op = tempop;
						vertexp_cube[i][j][k]->center_op = lzd_tools::get_center_on_op(tempop);
						
						// 填充edge
						for (int t = 0; t < _grids[i][j][k]->neighbers.size();t++){
							edgelist_cube[i][j][k][t]->overlap_idxs = tempcor[t];
						}
					}
				}
			}

			// 构建图
			for(int i = 0;i<x_num;i++){
				for(int j = 0;j<y_num;j++){
					for(int k = 0;k<z_num;k++){
						if(!_grids[i][j][k]->is_vaild())continue;
						ihp_graph._nodes.push_back(vertexp_cube[i][j][k]);
						ihp_graph._edges.push_back(edgelist_cube[i][j][k]);
					}
				}
			}
			return ihp_graph;
		}

		// 获得所有立方体的中心顶点(注意不是坐标)
		std::vector<Point<REAL, DIM>> get_centers() {
			std::vector<Point<REAL, DIM>> centers;
			for (int i = 0; i < x_num; i++) {
				for (int j = 0; j < y_num; j++) {
					for (int k = 0; k < z_num; k++) {
						if(!_grids[i][j][k]->is_vaild())continue;
						centers.push_back(lzd_tools::get_center_on_op(_grids[i][j][k]->grid_op));
					}
				}
			}
			return centers;
		}



}; 

	template<typename REAL, int DIM>
	class SeedGrowing:public NodeAllocator<REAL, DIM> {
		using pcl_op_cloud = pcl::PointCloud<pcl::PointNormal>;
		using pcl_kdtree = pcl::search::KdTree<pcl::PointNormal>;
		
		typedef struct SinglePoint {
			int id;
			//std::pair<Point<REAL, DIM>, Normal<REAL, DIM>> op;
			SinglePoint(int id) :id(id) {
				//op = _points_normals[id];
			}
		}SinglePoint;

		typedef struct PointEdge {
			int s;
			int e;
			REAL distance;

			PointEdge(int s, int e, REAL distance) :s(s), e(e),distance(distance) {
				assert(distance >= 0); // 有的模型存在重合的点
			}
		}PointEdge;

		typedef struct Cluster {
			int cluster_id;

			std::vector<int> op_type;//0表示local，1表示share
			std::vector<SinglePoint> points;

			std::map<int, std::vector<std::pair<int, int>>> share_points_cor;//share_points在cluster中的索引和在其他cluster中的索引。后续map的每个元素作为edge的overlap_idxs

			Cluster(int cluster_id) :cluster_id(cluster_id) {
				points.clear();
				op_type.clear();
				share_points_cor.clear();
			}

			int add_point(SinglePoint p, int type) {
				points.push_back(p);
				op_type.push_back(type);
				return points.size() - 1;
			}

			//此时扩张以及结束，根据点出现过的簇，计算本簇所有的点的原始id，以及在原始簇中的索引
			void build_cor(const std::vector<std::vector<std::pair<int, int>>>& idx_in_clusters,const std::vector<int>& _belonging) {
				for (int i = 0; i < points.size(); i++) {
					int pidx = points[i].id;//在全局中的索引
					for (int j = 0; j < idx_in_clusters[pidx].size(); j++) {
						int cidx = idx_in_clusters[pidx][j].first;//重合点所属的cluster_id
						if (cidx == cluster_id)continue;
						if (cidx != _belonging[pidx])continue;
						int pidx_in_c = idx_in_clusters[pidx][j].second;
						share_points_cor[cidx].push_back(std::make_pair(i, pidx_in_c));
					}
				}
			}
		}Cluster;

		typedef struct KnnGraph {
			// input
			const POINTS_NORMALS _points_normals;
			float _radius;//搜索半径上限
			int _k;//邻居数量上限

			// data
			pcl_kdtree::Ptr global_tree;
			pcl_op_cloud::Ptr cloud;
			std::vector<SinglePoint> points;
			std::vector<std::vector<PointEdge>> edges;// edges[i]表示第i个点的所有邻居

			// ************************************private**************************************************

			/**
			 * @brief 
			 * 维护一个点集。
			 * 开始时，点集内只有一个点；将点集中所有邻居中，距离最近的点加入点集。
			 * 重复上述过程，直到点集中的点数量达到k
			 */
			std::vector<int> __prim_neighbers_search(const std::vector<std::vector<float>>& dist_mat,int start_idx){
				int visited_num = 0, n = dist_mat.size();
				std::vector<float> distance(n,1e9);
				std::vector<bool> visited(n,false);
				
				distance[start_idx] = 0;
				
				for(int i = 0;i<n;i++){
					int t = -1;
					for(int j = 0;j<n;j++){
						if(visited[j])continue;
						if(t == -1||distance[j] < distance[t]){
							t = j;
						}
					}
					visited[t] = true;
					visited_num++;
					for(int j = 0;j<n;j++)distance[j] = std::min(distance[j],dist_mat[t][j]);
					if(visited_num >= _k)break;
				}
				std::vector<int> neighbers;
				for(int i = 0;i<n;i++){
					if(visited[i])neighbers.push_back(i);
				}
				return neighbers;			
			}

			// /**
			//  * @brief 
			//  * 
			//  * @param neighbers 
			//  * @param idx
			//  */
			// std::vector<int> _get_neighbers(int idx) {
			// 	assert(false);//这个方法用不了.
			// 	std::vector<int> neighbers;
			// 	std::vector<int> temp_idx;
			// 	std::vector<float> temp_distance;
			// 	global_tree->radiusSearch(cloud->points[idx], _radius, temp_idx, temp_distance);
				
			// 	//**********prim搜索**********
			// 	assert(temp_idx[0] == idx);
			// 	int idx_in_temp = 0;
				
			// 	std::vector<std::vector<float>> dist_mat;
			// 	dist_mat.resize(temp_idx.size(),std::vector<float>(temp_idx.size(),0));
			// 	for(int i = 0;i<temp_idx.size();i++){
			// 		for(int j = i + 1;j<temp_idx.size();j++){
			// 			REAL d = Distance(_points_normals[temp_idx[i]].first,_points_normals[temp_idx[j]].first);
			// 			dist_mat[i][j] = d;
			// 			dist_mat[j][i] = d;
			// 		}
			// 	}
			// 	neighbers = __prim_neighbers_search(dist_mat,idx_in_temp);	
			// 	for(int i = 0;i<neighbers.size();i++){
			// 		neighbers[i] = temp_idx[neighbers[i]];
			// 	}
			// 	return neighbers;
			// }

			/**
			 * @brief
			 * 从seed出发，计算所有点到seed的最短路径长
			 * @param distance
			 * @param seed
			 * @param early_stop 增长到一定点的数量后停止
			 * @return int 可达的点的数量
			 */
			int _seed_bfs(std::vector<REAL>& distance, int seed, int early_stop = -1) {
				lzd_tools::AcumutableTimer seed_bfs("seed_bfs");
				int num = 0;
				distance.clear();
				distance.resize(_points_normals.size(), 1e9);
				std::queue<int> q;
				q.push(seed);
				num++;
				distance[seed] = 0;
				while (!q.empty()) {
					int idx = q.front();
					q.pop();
					for (auto e : edges[idx]) {
						// 如果e.e没有被访问过 将e.e加入队列
						if (distance[e.e] == 1e9) {
							q.push(e.e);
							num++;
						}
						// 更新e.e的distance
						if (distance[e.e] > distance[e.s] + e.distance) {
							distance[e.e] = distance[e.s] + e.distance;
						}
					}
					if (num > early_stop && early_stop >0) {
						return num;
					}
				}
				return num;
			}

			// /**
			//  * @brief
			//  * 从seed出发，计算所有点到seed的最短路径长
			//  * @param distance
			//  * @param seed
			//  * @return int 可达的点的数量
			//  */
			// int _seed_bfs_discrete(std::vector<REAL>& distance, int seed) {
			// 	int num = 0;
			// 	distance.clear();
			// 	distance.resize(_points_normals.size(), 1e9);
			// 	std::queue<int> q;
			// 	q.push(seed);
			// 	num++;
			// 	distance[seed] = 0;
			// 	while (!q.empty()) {
			// 		int idx = q.front();
			// 		q.pop();
			// 		for (auto e : edges[idx]) {
			// 			// 如果e.e没有被访问过 将e.e加入队列
			// 			if (distance[e.e] == 1e9) {
			// 				q.push(e.e);
			// 				num++;
			// 			}
			// 			// 更新e.e的distance
			// 			if (distance[e.e] > distance[e.s] + 1) {
			// 				distance[e.e] = distance[e.s] + 1;
			// 			}
			// 		}
			// 	}
			// 	return num;
			// }

			// ************************************api**************************************************
			/**
			 * @brief Construct a new Knn Graph object
			 * 每个点的邻居数量上限为k
			 * @param ops 点云
			 * @param radius knn搜索半径上限
			 * @param k knn邻居数量上限
			 */
			KnnGraph(POINTS_NORMALS& ops, float radius, int k) :_radius(radius), _k(k),_points_normals(ops) {
				printf("init knngraph... \t");
				lzd_tools::AcumutableTimer initknn("init knngraph");
				cloud = op2pcl(ops);
				global_tree = pcl_kdtree::Ptr(new pcl_kdtree);
				global_tree->setInputCloud(cloud);

				//构建点
				printf("build points... \t");
				for (int i = 0; i < _points_normals.size(); i++) {
					points.push_back(SinglePoint(i));
				}
				//构建边
				printf("build edges... \t");
				edges.clear();
				edges.resize(_points_normals.size());
#pragma omp parallel for
				for (int i = 0; i < _points_normals.size(); i++) {
					//global_var::bar.set_progress((int)double(i) / _points_normals.size());
					std::vector<int> idx;
					std::vector<float> distance;
					global_tree->nearestKSearch(cloud->points[i], _k, idx, distance);
					for (int j = 0; j < min((int)idx.size(),_k); j++) {
						if (i == idx[j])continue;
						if (distance[j] >= _radius)break;
						if (edges[i].size() >= k)break;
						edges[i].push_back(PointEdge(i, idx[j], Distance(ops[i].first,ops[idx[j]].first)));
					}
				}
				printf("finish!\n");
			}

			/**
			 * @brief Get the neighbor object
			 * 得到idx的邻居列表
			 * @param idx  
			 */
			std::vector<int> get_neighbor(int idx) {
				std::vector<int> res;
				for (auto e : edges[idx]) {
					res.push_back(e.e);
				}
				return res;
			}

			/**
			 * @brief Get the nearest op idx object
			 * 得到距离p最近的点的索引
			 * @param p 
			 * @return int 最近点的索引 
			 */
			int get_nearest_op_idx(Point<REAL, DIM> p) {
				std::vector<int> idx;
				std::vector<float> distance;
				pcl::PointNormal pp;
				for (int i = 0; i < DIM; i++) {
					pp.data[i] = p[i];
					pp.normal[i] = 0;
				}
				global_tree->nearestKSearch(pp, 1, idx, distance);
				return idx[0];
			}

			// int get_nearest_op_idx(uint idx){
			// 	std::vector<PointEdge> temp = edges[idx];

			// }

		}KnnGraph;




// input
		POINTS_NORMALS const _points_normals;
		// REAL _radius_expand_rate;

// data
		KnnGraph _knn_graph; // 定义"联通"。

		/////////////种子//////////////////////
		std::vector<int> _seeds;
		/**
		 * @brief 
		 * 记录所有点到seed的最短路径长度
		 * _distance_matrix[i]的长度=点的数量,
		 * _distance_matrix[i][j]表示第j个点到第i个seed的最短路径长度
		 */
		std::vector<std::vector<REAL>> _distance_matrix;
		std::vector<int> _coverage_num;//记录每个seed的可达点的数量；如果某个种子的可达数量比较小，说明这个种子可能出现在噪点中

		/////////////簇///////////////////////
		std::vector<Cluster> _clusters;
		/**
		 * @brief 
		 * 记录所有点所出现的cluster以及在对应cluster中的索引; 
		 * len(_idx_in_cluster[i])表示第i个点所出现过的cluster的数量
		 * _idx_in_cluster[i][j].first表示第i个点出现的第j个cluster的id, _idx_in_cluster[i][j].second表示第i个点在第j个cluster中的索引;
		 */
		std::vector<std::vector<std::pair<int,int>>> _idx_in_cluster;
		/**
		 * @brief 
		 * 每个点只作为local_points出现一次。记录其作为local_points出现的cluster的id
		 */
		std::vector<int> _idx;
		
		//记录每个cluster之间的邻接关系。_cluster_graph[i][j] 中存储了cluster i和cluster j之间的所有边界点
		std::vector<std::vector<std::vector<int>>> _cluster_graph;

		std::vector<std::vector<double>> _edge_betweeness;
		std::vector<double> _node_betweeness;

		// 每个点的label
		aux_arg::LabelGetter<REAL, DIM>* _label_getter;
		std::vector<int> _label;

		nlohmann::json config_j;

// aid_func

		// 根据_distance_matrix，判断第[start:_points_normals.size())]的点是否都被访问过
		bool is_all_arrival(int &start){
			if(_distance_matrix.size() == 0)return false;
			for(int j = start;j<_points_normals.size();j++){
				bool visited = false;
				for(int i = 0;i<_seeds.size();i++){
					if(_distance_matrix[i][j] != 1e9){
						visited = true;
						break;
					}
				}
				if(!visited){
					start = j;
					return false;
				}
			}
			return true;
		}

/****************************************************种子的选择*******************************************************************************************/
		/**
		 * @brief
		 * 选择种子
		 * 任务1: 选择种子在点云中的索引 _seeds
		 * 任务2: 计算所有点到种子的最短路径长度 distance_matrix
		 * 任务3: 计算每个种子的可达点的数量 coverage_num
		*/
		

		// 从不可达点中随机选择种子 确保所有点都被访问到
		// 可以在其他种子选择策略使用前/后调用
		void __scan_unreachable(){
			int check_start = 0;
			while(!is_all_arrival(check_start)){
				int seed = check_start;
				std::vector<REAL> distance;
				int num = _knn_graph._seed_bfs(distance,seed);
				_seeds.push_back(seed);
				_distance_matrix.push_back(distance);
				_coverage_num.push_back(num);
				
			}
	
		}
		
		// 种子的策略可以改变，只需要保证种子可以到达所有点
		// 随机选择seed，并初始化_distance_matrix、_coverage_num。
		void _random_seed_selection(int seed_num = 100){
			assert(config_j.find("seed_selection_alg") == config_j.end());
			config_j["seed_selection_alg"] = "_random_seed_selection";
			_seeds.clear();
			_distance_matrix.clear();
			
			__scan_unreachable(); // TODO NOT TESTED

			for(int i = _seeds.size();i<seed_num;i++){
				int seed = rand() % _points_normals.size();
				while(std::find(_seeds.begin(),_seeds.end(),seed) != _seeds.end()){
					seed = rand() % _points_normals.size();
				}
				std::vector<REAL> distance;
				_seeds.push_back(seed);
				int num = _knn_graph._seed_bfs(distance,_seeds[i]);
				_distance_matrix.push_back(distance);
				_coverage_num.push_back(num);
			}
		}

		void _simple_cube_selection(int min_grid = 3){
			assert(config_j.find("seed_selection_alg") == config_j.end());
			nlohmann::json j;
			j["name"] = "simple_cube_selection";
			
			_seeds.clear();
			_distance_matrix.clear();
			__scan_unreachable(); // TODO NOT TESTED

			simple_cube<REAL, DIM> sc(_points_normals,min_grid);
			j["config"] = sc.get_config();
			config_j["seed_selection_alg"] = j;

			std::vector<Point<REAL,DIM>> centers = sc.get_centers();
			for(int i = 0;i<centers.size();i++){
				int seed = _knn_graph.get_nearest_op_idx(centers[i]);
				std::vector<REAL> distance;
				_seeds.push_back(seed);
				int num = _knn_graph._seed_bfs(distance,seed);
				_distance_matrix.push_back(distance);
				_coverage_num.push_back(num);
			}
		}

		/**
		 * @brief 
		 * 快速选择k个种子 保证每个种子所在的格子内大概为n/k个点
		 * @param k 种子大概数量
		 */
		void _k_grid_selection(int k=200){
			assert(config_j.find("seed_selection_alg") == config_j.end());
			config_j["seed_selection_alg"] = "_k_grid_selection";
			if(k > _points_normals.size()){
				k = _points_normals.size();
				printf("k is too large, set k = %d\n",k);
			}
			if (k <= 0) {
				printf("undified k, set k = 1\n");
				k = 1;
			}

			_seeds.clear();
			_distance_matrix.clear();
			__scan_unreachable(); // TODO NOT TESTED

			std::vector<int> idx = lzd_tools::k_select<REAL,DIM>(_points_normals,k);
			for(int i = 0;i<idx.size();i++){
				std::vector<REAL> distance;
				_seeds.push_back(idx[i]);
				int num = _knn_graph._seed_bfs(distance,idx[i]);
				_distance_matrix.push_back(distance);
				_coverage_num.push_back(num);
			}
		}

		void _adaptive_grid_selection(int point_per_grid = 2000){
			int k = _points_normals.size() / max(1,point_per_grid);
			_k_grid_selection(k);
		}

		void _select_by_label() {
			assert(config_j.find("seed_selection_alg") == config_j.end());
			config_j["seed_selection_alg"] = "_select_by_label";
			_seeds.clear();
			_label.resize(_points_normals.size());
			for (int i = 0; i < _label.size(); i++) {
				_label[i] = _label_getter->get_label(_points_normals[i].first);
			}
			std::map<int,std::vector<int>> label2idx;
			for (int i = 0; i < _label.size(); i++) {
				label2idx[_label[i]].push_back(i);
			}
			
//#pragma omp parallel for
			for (auto it = label2idx.begin(); it != label2idx.end(); it++) {
				int center = it->second[it->second.size()/2];
				std::vector<REAL> distance;
				_seeds.push_back(center);
				int num = _knn_graph._seed_bfs(distance, center);
				_distance_matrix.push_back(distance);
				_coverage_num.push_back(num);
			}
			// int seed_num = _seeds.size();
			// __scan_unreachable();
			// if(_seeds.size() != seed_num){
			// 	printf("add %d seeds because of unreachable points\n", _seeds.size() - seed_num);
			// }
		}


/****************************************************生成簇*************************************************************************************************/
		/**
		 * @brief 生成簇 (_clusters,每个簇的簇id从0~n)
		 * 任务1: 将每个点分配给一个种子 (_idx ∈ {簇id})
		 * 任务2: 分配的同时,计算点在簇中的索引(_idx_in_cluster)
		*/
		void _allocate_points_by_labels(){
			config_j["allocate_points_alg"] = "_allocate_points_by_labels";

			_idx.clear();
			_clusters.clear();
			_idx_in_cluster.clear();
			
			for(int i = 0;i<_points_normals.size();i++){
				_idx.push_back(_label[i]);
			}
			for(int i = 0;i<_seeds.size();i++){
				_clusters.push_back(Cluster(i));
			}

			for(int i = 0;i<_points_normals.size();i++){
				_idx_in_cluster.push_back(std::vector<std::pair<int, int>>());
				for(int j = 0;j<_seeds.size();j++){
					if(j==_idx[i]) {
						int t = _clusters[j].add_point(SinglePoint(i),0);
						_idx_in_cluster[i].push_back(std::make_pair(j,t));
					}
				}
			}
		}

		// 根据_distance_matrix，将每个点分配到距离最近的种子
		void _allocate_points_by_dist(){
			config_j["allocate_points_alg"] = "_allocate_points_by_dist";
			// 计算_idx 
			_idx.clear();
			//std::vector<int> cnt_idx(_seeds.size(),0);
			for(int i = 0;i<_points_normals.size();i++){
				REAL min_dist = 1e9;
				int min_idx = -1;
				for(int j = 0;j<_seeds.size();j++){
					if(_distance_matrix[j][i] < min_dist){
						min_dist = _distance_matrix[j][i];
						min_idx = j;
					}
				}
				_idx.push_back(min_idx);
				//cnt_idx[min_idx]++;
			}

			_clusters.clear();
			_idx_in_cluster.clear();
			for(int i = 0;i<_seeds.size();i++){
				_clusters.push_back(Cluster(i));
			}
			
			// 创建簇,并划分簇的local points
			for(int i = 0;i<_points_normals.size();i++){
				_idx_in_cluster.push_back(std::vector<std::pair<int, int>>());
				for(int j = 0;j<_seeds.size();j++){
					if(j==_idx[i]) {
						int t = _clusters[j].add_point(SinglePoint(i),0);
						_idx_in_cluster[i].push_back(std::make_pair(j,t));
					}
				}
			}
		}


/****************************************************生成图*************************************************************************************************/
		/**
		 * @brief 根据_idx,得到每个簇的邻居簇，根据簇的邻接关系得到_cluster_graph。
		 * 定义两个簇互为邻居簇的条件为：A簇中存在一个点，该点的邻域内存在点属于B簇。
		 * 定义边界点如下：某个点的邻域内的点，部分属于该点所在簇{A}，部分属于其他簇{B,C,D,...}，则该点为边界点，表示为A->B,C,D...
		 */
		void _cal_cluster_graph(){
			_cluster_graph.clear();
			_cluster_graph.resize(_seeds.size());
			for(int i = 0;i<_seeds.size();i++){
				_cluster_graph[i].resize(_seeds.size());
				for(int j = 0;j<_seeds.size();j++){
					_cluster_graph[i][j] = std::vector<int>();
				}
			}

			for(int i = 0;i<_points_normals.size();i++){
				std::vector<int> nei_ij = _knn_graph.get_neighbor(i);
				std::set<int> nei_ij_belonging;
				for(int j = 0;j<nei_ij.size();j++){
					nei_ij_belonging.insert(_idx[nei_ij[j]]);
				}
				for(auto it = nei_ij_belonging.begin();it != nei_ij_belonging.end();it++){
					if(*it == _idx[i])continue;
					_cluster_graph[_idx[i]][*it].push_back(i);
				}
			}
		}

		/**
		 * @brief 
	     * 将t簇并入s簇。需要修改：_idx、_idx_in_cluster、cluster、_cluster_graph、_distance_matrix
		 * _seed:不修改。将t中的并入s，无seed产生或者删除
		 * _idx[i]:如果_idx[i] == t,则 s -> _idx[i] 
		 * _idx_in_cluster[i]:if t in _idx_in_cluster[i] && s not in _idx_in_cluster[i],then replace t with s
		 * cluster:将t中的点加入s中。此时cluster中一定只有0点
		 * _cluster_graph:将t的边界点加入s的边界点中（s，t之间的边除外）
		 * _distance_matrix:_distance_matrix[s][:] = min(_distance_matrix[s][:],_distance_matrix[t][:])
		 */
		void __delete_edge(int s, int t) {
			int seedt = _seeds[t];
			int cnt_t = 0;
			for (int i = 0; i < _points_normals.size(); i++) {
				if (_idx[i] == t) {
					_idx[i] = s;
					assert(_idx_in_cluster[i].size() == 1);
					assert(_idx_in_cluster[i][0].first == t);
					int idx_in_cluster_t = _idx_in_cluster[i][0].second;
					assert(_clusters[t].op_type[idx_in_cluster_t] == 0);
					int idx_in_cluster_s = _clusters[s].add_point(_clusters[t].points[idx_in_cluster_t], 0);
					_idx_in_cluster[i][0] = std::make_pair(s, idx_in_cluster_s);
					cnt_t++;
				}
			}
			assert(cnt_t == _clusters[t].op_type.size());

			_cluster_graph[s][t].clear();
			_cluster_graph[t][s].clear();
			for (int i = 0; i < _seeds.size(); i++) {
				for (int j = 0; j < _cluster_graph[i][t].size(); j++) {
					_cluster_graph[i][s].push_back(_cluster_graph[i][t][j]);
				}
				_cluster_graph[i][t].clear();
				for (int j = 0; j < _cluster_graph[t][i].size(); j++) {
					_cluster_graph[s][i].push_back(_cluster_graph[t][i][j]);
				}
				_cluster_graph[t][i].clear();
			}
			for(int i = 0;i<_distance_matrix[s].size();i++){
				_distance_matrix[s][i] = std::min(_distance_matrix[s][i],_distance_matrix[t][i]);
			}

			_distance_matrix[t] = std::vector<REAL>(_distance_matrix[t].size(), 1e9);

/*			for (int i = 0; i < _cluster_graph.size(); i++) {
				_cluster_graph[i].erase(_cluster_graph[i].begin() + t);
			}*/					
			_clusters[t].points.clear();
			_clusters[t].op_type.clear();
			//_seeds.erase(_seeds.begin() + t);
		}

		// 计算edge_betweeness和node_betweeness，并正则化
		void _cal_betweeness() {
			aux_arg::VV<int> cluster_ad_graph;
			cluster_ad_graph.resize(_seeds.size());
			for (int i = 0; i < _seeds.size(); i++) {
				cluster_ad_graph[i].resize(_seeds.size(), 0);
				for (int j = 0; j < _seeds.size(); j++) {
					if (_cluster_graph[i][j].size() > 0)cluster_ad_graph[i][j] = 1;
					else cluster_ad_graph[i][j] = -1;
				}
			}
			_edge_betweeness = aux_arg::_boost_cal_edge_betweenness(cluster_ad_graph);
			// normalize _edge_betweeness
			REAL max_edge_betweeness = 0;
			for (int i = 0; i < _edge_betweeness.size(); i++) {
				for (int j = 0; j < _edge_betweeness[i].size(); j++) {
					if (_edge_betweeness[i][j] > max_edge_betweeness) {
						max_edge_betweeness = _edge_betweeness[i][j];
					}
				}
			}
			for (int i = 0; i < _edge_betweeness.size(); i++) {
				for (int j = 0; j < _edge_betweeness[i].size(); j++) {
					if (_edge_betweeness[i][j] < 0) {
						assert(_edge_betweeness[i][j] == -1);
					}
					_edge_betweeness[i][j] /= max_edge_betweeness;
				}
			}			
			
			_node_betweeness = aux_arg::_boost_cal_node_betweeness(cluster_ad_graph);
			REAL max_node_betweeness = 0;
			for (int i = 0; i < _node_betweeness.size(); i++) {
				assert(_node_betweeness[i] >= 0);
				if (_node_betweeness[i] > max_node_betweeness) {
					max_node_betweeness = _node_betweeness[i];
				}
			}
			for (int i = 0; i < _node_betweeness.size(); i++) {
				_node_betweeness[i] /= max_node_betweeness;
			}
		}

		// 检测cut edge,并消除第一个cutedge。返回cut edge的数量
		// int __minedege_comb(std::vector<int> sub_seeds) {
		// 	if (sub_seeds.size() <= 1)return 0;
		// 	bgl_api::UGRAPH g;
		// 	g.resize(sub_seeds.size());
		// 	for (int i = 0; i < sub_seeds.size(); i++) {
		// 		g[i].resize(sub_seeds.size());
		// 	}
		// 	// 从_cluster_graph中构建UGRAPH
		// 	// 只要_cluster_graph[i][j]或者_cluster_graph[j][i]中有边，就认为i和j之间有边
		// 	for (int i = 0; i < sub_seeds.size(); i++) {
		// 		for (int j = i + 1; j < sub_seeds.size(); j++) {
		// 			int s = sub_seeds[i], t = sub_seeds[j];
		// 			if (_cluster_graph[s][t].size() > 0 || _cluster_graph[t][s].size() > 0) {
		// 				double betweeness = max(_edge_betweeness[s][t], _edge_betweeness[t][s]);
		// 				assert(betweeness > 0 && betweeness <= 1);
		// 				double w = 1 / betweeness;
		// 				w *= 100;
		// 				w = min(w, 1e9);
		// 				g[i][j] = (int)w, g[j][i] = (int)w;
		// 				//_edge_betweeness[s][t] = w, _edge_betweeness[t][s] = w;
		// 			}
		// 			else {
		// 				g[i][j] = 0, g[j][i] = 0;
		// 			}
		// 		}
		// 	}
		// 	std::vector<std::pair<int, int>> res = bgl_api::bgl_stoer_wanger_mincut(g);
		// 	assert(res.size() != 0);
		// 	for (int i = 0; i < 1; i++) {
		// 		__delete_edge(sub_seeds[res[i].first], sub_seeds[res[i].second]);
		// 	}
		// 	return res.size();
		// }

	//***************一些块过滤器,用于合并块************************

		// int _low_minedgecut_filter_iter() {
		// 	bgl_api::UGRAPH g;
		// 	g.resize(_seeds.size());
		// 	for (int i = 0; i < _seeds.size(); i++) {
		// 		g[i].resize(_seeds.size());
		// 	}
		// 	// 从_cluster_graph中构建UGRAPH
		// 	// 只要_cluster_graph[i][j]或者_cluster_graph[j][i]中有边，就认为i和j之间有边
		// 	for (int i = 0; i < _seeds.size(); i++) {
		// 		for (int j = i+1; j < _seeds.size(); j++) {
		// 			if (_cluster_graph[i][j].size() > 0 || _cluster_graph[j][i].size() > 0) {
		// 				g[i][j] = 1, g[j][i] = 1;
		// 			}else{
		// 				g[i][j] = -1, g[j][i] = -1;
		// 			}
		// 		}
		// 	}

		// 	// 得到最大的联通子图
		// 	auto componets = bgl_api::bgl_connected_components(g);
		// 	int c = 0;
		// 	for (int i = 0; i < componets.size(); i++) {
		// 		c = max(c,__minedege_comb(componets[i]));
		// 	}
		// 	return c;
		// }

		// void _low_minedgecut_filter(int max_epoch = 10,bool if_save_mid = false){
		// 	nlohmann::json j;
		// 	j["name"] = ("_low_minedgecut_filter");
		// 	j["filter_times"] = max_epoch;
		// 	config_j["edge_filters"].push_back(j);
			
		// 	std::string path = global_var::data_output_base + "/graph_ipsr/temp/";
		// 	rmkdir(path);
		// 	_cal_betweeness();
		// 	_draw_topology(path + "filter_before");
		// 	for (int i = 0; i < max_epoch; i++) {
		// 		int c = _low_minedgecut_filter_iter();
		// 		if(if_save_mid)_draw_topology(path + "filter_after" + std::to_string(i));
		// 		_cal_betweeness();
		// 		if (c >= 4)break;
		// 	}

		// }

		// void _high_edge_betweenness_filter(){
		// 	config_j["edge_filter"].push_back("_high_edge_betweenness_filter");
		// 	std::string path = global_var::data_output_base + "/graph_ipsr/temp/";
		// 	rmkdir(path);
		// 	_cal_betweeness();
		// 	_draw_topology(path + "filter_before");

		// 	//global_var::config_j["graph_ipsr_config"]["SeedGrowingConfig"]["high_edge_betweenness_filter_delete_cnt"] = "const 3";

		// 	// 删除betweenness top10%的边
		// 	int cnt = 0;
		// 	while (cnt < 0) {
		// 		_cal_betweeness();
		// 		// 找到betweenness最大的边
		// 		std::pair<int, int> e;
		// 		REAL max_betweeness = 0;
		// 		for (int i = 0; i < _edge_betweeness.size(); i++) {
		// 			for (int j = 0; j < _edge_betweeness[i].size(); j++) {
		// 				if (_edge_betweeness[i][j] > max_betweeness) {
		// 					max_betweeness = _edge_betweeness[i][j];
		// 					e = std::make_pair(i, j);
		// 				}
		// 			}
		// 		}
		// 		__delete_edge(e.first, e.second);
		// 		cnt++;
		// 		_draw_topology(path + "filter_after" + std::to_string(cnt));
		// 	}
		// }

		void _small_cluster_filter() {
			nlohmann::json j;
			j["name"] = "_small_cluster_filter";
			int avg_num = 0;
			for(int i = 0;i<_clusters.size();i++){
				avg_num += _clusters[i].points.size();
			}
			avg_num /= _clusters.size();
			int num_to_delete = avg_num/2;
			j["num_to_delete"] = num_to_delete;

			for(int i = 0;i<_clusters.size();i++){
				if(_clusters[i].points.size() > num_to_delete)continue;
				// 查找所有邻居簇 找到 cluster.size / boundary.size 最小者
				int s = _clusters[i].cluster_id, t = -1, tsz_p_boundary_points = 1e9;
				assert(s == i);
				for(int j = 0;j<_clusters.size();j++){
					if(s == j)continue;
					if(_cluster_graph[s][j].size() == 0 || _cluster_graph[j][s].size() == 0)continue;
					int tpbp = _clusters[j].points.size() / _cluster_graph[s][j].size();
					if(tpbp < tsz_p_boundary_points){
						t = j, tsz_p_boundary_points = tpbp;
					}
				}
				if(t != -1){
					printf("merge cluster_%d's size = %d,  cluster_%d's size = %d, boundary points num between them = % d\n", s, _clusters[s].points.size(), t, _clusters[t].points.size(), _cluster_graph[s][t].size());
					__delete_edge(t,s);
					printf("after merge, cluster_%d's size = %d, cluster_%d's size = %d\n",s,_clusters[s].points.size(),t,_clusters[t].points.size());
				}else{
					printf("cluster %d,'s size = %d but can't find a cluster to merge\n",s,_clusters[s].points.size());
				}
			}
			config_j["node_filter"].push_back(j);
		}

		// // 从外而内扩张簇（不好用）
		// void _cluster_expand_by_other() {
		// 	for(int i = 0;i<_points_normals.size();i++){
		// 		for(int j = 0;j<_seeds.size();j++){
		// 			REAL dist_to_seed = _distance_matrix[_idx[i]][i];
		// 			// 条件:如果该点与簇k的距离<该点与本簇距离的1.25倍,且两个簇存在邻接关系,则将该点加入簇k
		// 			bool flag = true;
		// 			flag = flag && _distance_matrix[j][i] <= _radius_expand_rate * dist_to_seed;
		// 			flag = flag && _cluster_graph[_idx[i]][j].size() > 0;
		// 			if(!flag)continue;	
		// 			int t = _clusters[j].add_point(SinglePoint(i),1);
		// 			_idx_in_cluster[i].push_back(std::make_pair(j,t));
		// 		}
		// 	}			
		// }

		/**
		 * @brief 
		 * 从边界点遍历扩张簇。定义cluster <A,B>的boundary_point：属于A，但其k邻域内存在B簇的点
		 * for A in _clusters: // 
		 * 	for <A,B> in _cluster_graph[A]:	 // 
		 *	   初始化一个处理队列q,boundary point点加入队列,并标记为已访问;
		 *     初始化一个距离数组,表示与boundary point的距离，将B中的点距离设为1e9;其余点距离设为0，它们不会被"加入"到处理队列中
		 *     从q开始,bfs遍历,计算出队点的distance，直到q为空;将距离矩阵中大于零且小于阈值的点加入A簇,作为A关于B的外点
		 */
		void _cluster_expand_by_boundary() {
			int max_dist = seedgrowing_config_s["expand_max_iter"];
			nlohmann::json j;
			lzd_tools::AcumutableTimer expand("cluster_expand_by_boundary");
			j["name"] = "cluster_expand_by_boundary";
			j["max_dist"] = max_dist;
			config_j["cluster_expand_method"] = j;
			for(int i = 0;i<_clusters.size();i++){
				for(int j = 0;j<_clusters.size();j++){
					if(i == j)continue;
					if(_cluster_graph[i][j].size() == 0)continue;
					std::vector<int> boundary_points = _cluster_graph[i][j];
					std::vector<int> distance(_points_normals.size(),0);
					std::queue<int> q;
					for(int k = 0;k<boundary_points.size();k++){
						assert(_idx[boundary_points[k]] == i); // A->B的边界点一定是A的local_point
						q.push(boundary_points[k]);
					}
					for(int k = 0;k<_clusters[j].points.size();k++){
						if (_clusters[j].op_type[k] != 0)continue;
						distance[_clusters[j].points[k].id] = 1e9;
					}
					while(!q.empty()){
						int idx = q.front();
						q.pop();
						for(auto e:_knn_graph.edges[idx]){
							assert(idx == e.s);
							if(distance[e.e]  == 1e9){
								q.push(e.e);
								assert(distance[idx] != 1e9);
							}
							if(distance[e.e] > distance[e.s] + 1){
								distance[e.e] = distance[e.s] + 1;
							}
						}
					}
					for(int k = 0;k<_points_normals.size();k++){
						if(distance[k] > 0 && distance[k] < max_dist){
							int t = _clusters[i].add_point(SinglePoint(k),1);
							_idx_in_cluster[k].push_back(std::make_pair(i,t));
						}
					}
				}
			}
		}
				
		// 根据计算好的_distance_matrix，计算_idx，划分簇,扩大簇,并计算_idx_in_cluster。
		void _graphing(){
			_cal_cluster_graph();

			if(seedgrowing_config_s["small_cluster_filter"] == true)_small_cluster_filter();
			
			//_low_minedgecut_filter();
			// 计算当前每个边的betweenness，将betweenness很大的边两侧的点合并为一个点
			// _cluster_expand_by_other();
			_cluster_expand_by_boundary();
			for(int i = 0;i<_seeds.size();i++){
				_clusters[i].build_cor(_idx_in_cluster,_idx);
			}
		}

// help func

		// 打印点云和拓扑结构
		// 要求:_idx、_betweeness已经计算好
		// path即ply文件路径（无后缀名）
		void _draw_topology(std::string path){
			// 根据_idx绘制点云
			lzd_tools::InFiniteMap<int>* mapper = new lzd_tools::RandomColorMapper<int>(_idx);
			lzd_tools::op2ply(_points_normals,path+"op.ply", XForm<REAL, DIM + 1>().Identity(), mapper);

			std::vector<REAL> op_metric;
			MESH mesh;
			// 绘制种子点
			for(int i = 0;i<_clusters.size();i++){
				auto sphere = lzd_tools::get_sphere(_points_normals[_seeds[i]].first,0.01);
				lzd_tools::add_topology(mesh,sphere);
				for(int j = 0;j<sphere.first.size();j++){
					op_metric.push_back(_node_betweeness[i]);
				}
			}
			// 绘制边
			for(int i = 0;i<_cluster_graph.size();i++){
				for(int j = 0;j<_cluster_graph[i].size();j++){
					if(_cluster_graph[i][j].size() == 0)continue;
					auto arrow = lzd_tools::get_arrow(_points_normals[_seeds[i]].first,_points_normals[_seeds[j]].first,0.02);
					lzd_tools::add_topology(mesh,arrow);
					for(int k = 0;k<arrow.first.size();k++){
						op_metric.push_back(_edge_betweeness[i][j]);
					}
				}
			}
			lzd_tools::mesh2ply(mesh, path+"topology.ply", XForm<REAL, DIM + 1>().Identity(), new lzd_tools::GradientMapper<REAL>(op_metric));
		}
	
		// 保存每个簇的距离图（可视化每个点到簇中心的距离）
		void _print_distance(std::string path,int seed_idx){
			POINTS_NORMALS op;
			Cluster& cluster = _clusters[seed_idx];
			std::vector<REAL> op_dist;
			//std::vector<REAL> alldist;
			
			/*for (int i = 0; i < _clusters.size(); i++) {
				for (int j = 0; j < _distance_matrix[i].size(); j++) {
					if(_distance_matrix[i][j]!=1e9)alldist.push_back(_distance_matrix[i][j]);
				}
			}*/
			
			for (int i = 0; i < cluster.points.size(); i++) {
				if (cluster.op_type[i] == 1)continue;
				int t = cluster.points[i].id;
				op.push_back(_points_normals[t]);
				op_dist.push_back(_distance_matrix[seed_idx][t]);
			}
			
			path = path + std::to_string(seed_idx) + "dist.ply";
			lzd_tools::InFiniteMap<REAL>* mapper = new lzd_tools::GradientMapper<REAL>(op_dist);
			mapper->_flag = op_dist;
			lzd_tools::op2ply(op,path,XForm<REAL,DIM+1>().Identity(),mapper);
		}

		void _select_seed(){
			auto sp_config = seedgrowing_config_s["seed_selection_alg_spconfig"];
			lzd_tools::AcumutableTimer clock(sp_config["name"]);
			if(sp_config["name"] == "k_selection"){
				_k_grid_selection(sp_config["k"]);
			}else if(sp_config["name"] == "adaptive_grid_selection"){
				_adaptive_grid_selection(sp_config["point_per_grid"]);
			}else if(sp_config["name"] == "simple_cube_selection"){
				_simple_cube_selection(sp_config["cut_times"]);
			}else if(sp_config["name"] == "random_seed_selection"){
				_random_seed_selection(sp_config["seed_num"]);
			}else{
				printf("no such seed selection algorithm named %s\n",sp_config["name"].get<std::string>().c_str());
				exit(0);
			}
		}


public:
		// static config

		static nlohmann::json seedgrowing_config_s;

		SeedGrowing(POINTS_NORMALS& points_normals, aux_arg::LabelGetter<REAL,DIM>* label_getter):_points_normals(points_normals),_label_getter(label_getter),_knn_graph(points_normals,seedgrowing_config_s["knn_graph"]["radius"], seedgrowing_config_s["knn_graph"]["k"]) {
			_label_getter = label_getter;
			printf("SeedGrowing start init...\n");
			printf("start select seeds...\n");
			// 选择种子
			_select_by_label();
			int t = 0;
			// TODO 存在不可达点
			if(!(is_all_arrival(t))){
				printf("exist unreachable points\n");
			}

			printf("start clustering...\n");
			_allocate_points_by_labels();

			// 生成图
			_graphing();
			printf("SeedGrowing init done!\n");
			// 打印簇的数量
			printf("cluster num = %d\n",_clusters.size());
			config_j["name"] = "SeedGrowing";		
		}


		SeedGrowing(POINTS_NORMALS& points_normals):_points_normals(points_normals), _knn_graph(points_normals, seedgrowing_config_s["knn_graph"]["radius"], seedgrowing_config_s["knn_graph"]["k"]) {
			printf("SeedGrowing start init...\n");
			printf("start select seeds...\n");
			// 选择种子
			
			_select_seed();
			int t = 0;
			assert(is_all_arrival(t));

			printf("start clustering...\n");
			// 生成cluster
			_allocate_points_by_dist();

			// 生成图
			_graphing();

			printf("SeedGrowing init done!\n");
			config_j["name"] = "SeedGrowing";
		}

		nlohmann::json get_config() {
			return config_j;
		}

		linked_graph<REAL, DIM> op_orchestration(
			IPSR_Factory<REAL, DIM> ipsr_factory,
			int (*_update_plan)(Period& p),
			NormalEstimation<REAL,DIM>* estimator
		){
			std::vector<VERTEX_P<REAL,DIM>> vertexp;
			std::vector<std::vector<EDGE_P<REAL, DIM>>> edgelist;
			vertexp.clear();
			edgelist.clear();

			std::map<int,int> clusterId2Idx;
			std::map<int,int> vpIdx2clusterIdx;

			// 初始化顶点和边的指针数组
			for (int i = 0; i < _clusters.size(); i++)
			{
				assert(_clusters[i].cluster_id == i);
				if(_clusters[i].points.size() == 0){
					continue;
				}
				auto nvp = new VERTEX<REAL,DIM>();
				POINTS_NORMALS tempop(_clusters[i].points.size());
				std::vector<int> tempflag(_clusters[i].points.size());
				for(int j = 0;j<_clusters[i].points.size();j++){
					tempop[j] = _points_normals[_clusters[i].points[j].id];
					tempflag[j] = _clusters[i].op_type[j];
				}
				nvp->_handle = ipsr_factory.create_ipsr_from_op(tempop, _update_plan, new DoingNothing<REAL,DIM>());
				nvp->_op = tempop;
				nvp->_op_type = tempflag;
				nvp->center_op = _points_normals[_seeds[i]].first;
				clusterId2Idx[_clusters[i].cluster_id] = vertexp.size();
				vpIdx2clusterIdx[vertexp.size()] = _clusters[i].cluster_id;  
				vertexp.push_back(nvp);
				edgelist.push_back(std::vector<EDGE_P<REAL, DIM>>(0));
			}
			lzd_tools::thread_safe_int finished(0);

			printf("start init\v");
#pragma omp parallel for
			for(int i = 0;i<vertexp.size();i++){
				vertexp[i]->_handle->init_op_normal(estimator);
				printf("%d/%d\r",finished.get(), vertexp.size());
				++finished;
			}
			printf("\n");
			
			// 填充edge
			for(int i = 0;i<vertexp.size();i++){
				for(int j = 0;j<vertexp.size();j++){
					int clusteri = vpIdx2clusterIdx[i],clusterj = vpIdx2clusterIdx[j], s=i,e=j;
					if(clusteri == clusterj)continue;
					// if(_cluster_graph[clusteri][clusterj].size() == 0)continue;
					EDGE<REAL, DIM>* edgep;
					if(_clusters[clusteri].share_points_cor.find(clusterj) == _clusters[clusteri].share_points_cor.end()){
						edgep = new EDGE<REAL, DIM>(vertexp[s], vertexp[e]);
					}else{
						edgep = new EDGE<REAL, DIM>(vertexp[s], vertexp[e], _clusters[clusteri].share_points_cor[clusterj]);
					}
					edgep -> set_connection_strength((REAL)_cluster_graph[clusteri][clusterj].size());
					if(!edgep->is_vaild())continue;
					edgelist[s].push_back(edgep);
				}
			}

			std::string path = vertexp[0]->_handle->get_out_put_base("/graph_ipsr/seedcluster/");
			// for(int i = 0;i<vertexp.size();i++){
			// 	_print_distance(path,i);
			// }

			linked_graph<REAL, DIM> ihp_graph(vertexp, edgelist);	
			return ihp_graph;			
		}	
	};

	/**
	 * @brief 
	 * 根据label分配点到簇 再根据簇生成图
	 */
	template<typename REAL, int DIM>
	class AllocateByLabels:public NodeAllocator<REAL,DIM> {
		private:
			std::vector<int> _label;
			POINTS_NORMALS _points_normals;

		public:
			AllocateByLabels(std::vector<int> label, POINTS_NORMALS points_normals):_label(label),_points_normals(points_normals){
			}
			
			nlohmann::json get_config() {
				nlohmann::json j;
				j["name"] = "AllocateByLabels";
				return j;
			}

			//TODO 写一个自己的...
			linked_graph<REAL, DIM> op_orchestration(
				IPSR_Factory<REAL, DIM> ipsr_factory,
				int (*_update_plan)(Period& p),
				NormalEstimation<REAL,DIM>* estimator
			){
				aux_arg::LabelGetter<REAL,DIM>* label_getter = new aux_arg::LabelGetter<REAL,DIM>(_points_normals,_label);
				SeedGrowing<REAL,DIM> sg(_points_normals,label_getter);
				return sg.op_orchestration(ipsr_factory,_update_plan,estimator);
			}		
	};

	template<typename REAL, int DIM>
	NodeAllocator<REAL,DIM>* createAllocatorByJson(nlohmann::json j, POINTS_NORMALS points_normals){
		std::string name = j["name"];
		if(name == "AllocateByLabels"){
			std::vector<int> label;
			std::string ip = j["ip"];
			int port = j["port"];
			lzd_tools::Socket socket(ip,port);
			assert(socket.connectToServer());
			nlohmann::json req = j["method"];
			req["data_size"] = points_normals.size();
			socket.SendJson(req.dump());
			std::string resbuff;
			socket.ReceiveJson(resbuff);
			nlohmann::json res = nlohmann::json::parse(resbuff);
			assert (res["status"] == "OK");
			std::vector<REAL> data;
			for(int i =0;i<points_normals.size();i++){
				for (int j = 0; j < DIM; j++) {
					data.push_back(points_normals[i].first[j]);
				}
			}
			socket.SendDoubleArray(data);
			socket.ReceiveIntArray(label,points_normals.size());
			return new AllocateByLabels<REAL,DIM>(label,points_normals);
		}else if(name == "SeedGrowing"){
			return new SeedGrowing<REAL,DIM>(points_normals);
		}else if(name == "simple_cube"){
			return new simple_cube<REAL,DIM>(points_normals);
		}else{
			printf("error allocator name\n");
			return nullptr;
		}
		return nullptr;

	}

	typedef double REAL;
	const unsigned int DIM = 3U;
	int simple_cube<REAL, DIM>::_op_num_treshold = 10;
	REAL simple_cube<REAL, DIM>::rectangle_bbox::INF = 1e9;
	REAL simple_cube<REAL, DIM>::simple_cube_grid::overlap_rate = 0.25;
	int graph_edge<REAL, DIM>::_overlap_treshold = 0;
	nlohmann::json static_conf = nlohmann::json();
	nlohmann::json SeedGrowing<REAL, DIM>::seedgrowing_config_s = nlohmann::json();

	static void assign_static_conf(){
		static_conf = ConfigManager::get_config_in_config_floder("graph_ipsr.json")["static_conf"];
		SeedGrowing<REAL, DIM>::seedgrowing_config_s = static_conf["SeedGrowing"];
		graph_edge<REAL, DIM>::_overlap_treshold = static_conf["graph_edge"]["overlap_treshold"];
	}

}//namespace
