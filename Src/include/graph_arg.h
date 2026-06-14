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
		int inv_time;
		
		bool same_to_gt;

		int _size;

		FlipableNode(nid id, bool if_gt = 1);

		FlipableNode(nlohmann::json j);

		bool is_inv() const;

		bool is_same() const;

		int flip();
	}FlipableNode;

	typedef struct ConfidenceEdge {
		const FlipableNode* start;
		const FlipableNode* end;
		REAL weight;
		REAL inv_weight;
		REAL confidence;

		ConfidenceEdge(const FlipableNode* start, const FlipableNode* end, REAL weight, REAL inv_weight, REAL confidence = 1);

		
		weight_pair get_weight(int if_same = -1) const;

		bool is_good_edge();

		void boolize_weight();

		void multiply_conf();

	}ConfidenceEdge;

	typedef struct Tree {
		FlipableNode* _root;
		std::vector<Tree*> _childrens;

		Tree(FlipableNode* root);
		void add_child(Tree* child);

		/**
		 * @brief 
		 * @param vistor 
		 */
		void PostTraverse(TreeVistor* vistor);

		/**
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

	class FlipableGraph {
		nlohmann::json log_j;
		
		
	public:
		std::vector<FlipableNode*> _nodes;
		std::vector<std::vector<ConfidenceEdge*>> _edges;

		void _align_flip(const std::vector<bool>& flip1, std::vector<bool>& flip2);
		
		void reset_flip_status();

		/**
		 * @brief 
		 */
		void apply_flip(const std::vector<bool>& node_status);

		FlipableGraph(std::vector<FlipableNode*> nodes, std::vector<std::vector<ConfidenceEdge*>> edges);

		FlipableGraph(nlohmann::json j);

		// nlohmann::json get_config();
		// nlohmann::json get_metric();
		// nlohmann::json get_log();
		// nlohmann::json to_json();
		/**
		 * @return weight_pair 
		 */
		weight_pair cal_weight(const std::vector<FlipableNode*>& _node1, const std::vector<FlipableNode*>& _node2) const;
		
		/**
		 * @return REAL 
		 */
		REAL cal_current_weight_sum() const;

		std::vector<Tree*> get_random_forest();

		nlohmann::json get_metric();
		
		void boolize_weight();

		void multiply_conf();
	
	};

};