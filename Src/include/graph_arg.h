#pragma once
#include <nlohmann/json.hpp>
#include <vector>



namespace graph_arg {
	typedef double REAL;
	typedef unsigned int nid;
	typedef std::pair<REAL,REAL> weight_pair;

	struct Tree;
	class TreeVistor;

	class FlipableGraph;
	class FlipGraph;

	typedef struct FlipableNode {
		const nid id;
		int inv_time; // 被反转的次数
		
		bool same_to_gt;// 初始时(inv_time==0时)是否和gt一样;0表示不一样，1表示一样

		int _size;// 点云数量

		FlipableNode(nid id, bool if_gt = 1);

		FlipableNode(nlohmann::json j);

		bool is_inv() const;

		bool is_same() const; // 当前是否和gt一样

		int flip();
	}FlipableNode;

	typedef struct ConfidenceEdge {
		const FlipableNode* start;
		const FlipableNode* end;
		REAL weight;
		REAL inv_weight;
		REAL confidence;

		ConfidenceEdge(const FlipableNode* start, const FlipableNode* end, REAL weight, REAL inv_weight, REAL confidence = 1);

		//delete ConfidenceEdge(nlohmann::json j);// 无法从json构造
		
		weight_pair get_weight(int if_same = -1) const;

		bool is_good_edge();//1表示是好边，0表示是坏边

		// 布尔化权重
		void boolize_weight();

		// 乘以confidence
		void multiply_conf();

	}ConfidenceEdge;

	typedef struct Tree {
		FlipableNode* _root;
		std::vector<Tree*> _childrens;

		Tree(FlipableNode* root);
		void add_child(Tree* child);

		/**
		 * @brief 
		 * 后序遍历树
		 * @param vistor 
		 */
		void PostTraverse(TreeVistor* vistor);

		/**
		 * 获得当前树的所有节点 
		 */
		std::vector<FlipableNode*> get_all_nodes();

		// nlohmann::json to_json();
	}Tree;

	class TreeVistor {
	public:
		virtual void visit(Tree* t) = 0;
	};

	class align_to_root : public TreeVistor {
		const FlipableGraph* graph;
	public:
		align_to_root(FlipableGraph* graph);
		void visit(Tree* t);
	};

	class FlipGraph {
	public:
		virtual std::vector<bool> flip(FlipableGraph* g) = 0;
		virtual nlohmann::json get_config() = 0;	
		virtual nlohmann::json get_log() = 0;
		
		static FlipGraph* get_flip_arg_from_json(nlohmann::json config_j);
	
	};



	class OptimFlip : public FlipGraph {
	public:
		std::vector<bool> flip(FlipableGraph* g);
		nlohmann::json get_config();
		nlohmann::json get_log();
	};

	/**
	 * @brief 
	 * 暴力翻转，遍历所有可能的翻转序列
	 */
	class BruteForceFlip:public FlipGraph{
	public:
		std::vector<bool> flip(FlipableGraph* g);
		nlohmann::json get_config();
		nlohmann::json get_log();
	};

	/**
	 * @brief 
	 * 
	 */
	class ForestFlip:public FlipGraph{
	public:
		std::vector<bool> flip(FlipableGraph* g);
		nlohmann::json get_config();
		nlohmann::json get_log();
	};

	/**
	 * @brief 
	 * 随机进行times次翻转;将翻转序列对齐后,对每个点投票
	 */
	class VoteFlip :public FlipGraph {
		int _times;
		FlipGraph* _base_alg;
		nlohmann::json _log_j;
	public:
		VoteFlip(FlipGraph* base_alg,int times);
		std::vector<bool> flip(FlipableGraph* g);
		nlohmann::json get_config();
		nlohmann::json get_log();
	};
	
	/**
	 * @brief 
	 * 多个翻转序列中选择weight_sum最小的
	 */
	class ChoseBestFlip :public FlipGraph {
		int _times;
		FlipGraph* _base_alg;
		nlohmann::json _log_j;
	public:
		ChoseBestFlip(FlipGraph* base_alg,int times);
		std::vector<bool> flip(FlipableGraph* g);
		nlohmann::json get_config();
		nlohmann::json get_log();
	};

	/**
	 * @brief 
	 * 搜索weight_sum最小的翻转; 区别于ChoseBestFlip，ChoseBestFlip在每次尝试后会reset flip状态
	 */
	class SearchBestFlip :public FlipGraph {
		int _times;
		FlipGraph* _base_alg;
		nlohmann::json _log_j;
	public:
		SearchBestFlip(FlipGraph* base_alg,int times);
		std::vector<bool> flip(FlipableGraph* g);
		nlohmann::json get_config();
		nlohmann::json get_log();
	};

	class MIQPFlip :public FlipGraph {
		nlohmann::json _log_j;
		std::string _ip;
		int _port;
	public:
		MIQPFlip(nlohmann::json j);
		MIQPFlip(std::string ip, int port) :_ip(ip), _port(port) {};
		std::vector<bool> flip(FlipableGraph* g);
		nlohmann::json get_config();
		nlohmann::json get_log();
	};

	double cal_flip_acc(std::vector<FlipableNode*> _node);

	double cal_edge_acc(std::vector<std::vector<ConfidenceEdge*>> _edges);

    // 邻接矩阵存储的graph
	class FlipableGraph {
		nlohmann::json log_j;
		
		
	public:
		std::vector<FlipableNode*> _nodes;
		std::vector<std::vector<ConfidenceEdge*>> _edges; // 矩阵存储的边。其中edge[i][j]表示i->j的边。 NULL 表示不存在边

		// 对齐两个状态序列
		void _align_flip(const std::vector<bool>& flip1, std::vector<bool>& flip2);
		
		void reset_flip_status();

		/**
		 * @brief 
		 * @param node_status 0表示不翻转，1表示翻转 
		 */
		void apply_flip(const std::vector<bool>& node_status);

		FlipableGraph(std::vector<FlipableNode*> nodes, std::vector<std::vector<ConfidenceEdge*>> edges);

		FlipableGraph(nlohmann::json j);

		// nlohmann::json get_config();
		// nlohmann::json get_metric();
		// nlohmann::json get_log();
		// nlohmann::json to_json();
		/**
		 * @brief 计算_node1和_node2之间所有边的权重
		 * @return weight_pair 
		 */
		weight_pair cal_weight(const std::vector<FlipableNode*>& _node1, const std::vector<FlipableNode*>& _node2) const;
		
		/**
		 * @brief 得到当前所有边的权重和
		 * @return REAL 
		 */
		REAL cal_current_weight_sum() const;

		// 每个联通分量随机生成一颗树
		std::vector<Tree*> get_random_forest();

		nlohmann::json get_metric();
		
		// 将所有边的权重布尔化
		void boolize_weight();

		void multiply_conf();
	
	};

};