/*
Copyright (c) 2022, Fei Hou and Chiyu Wang, Institute of Software, Chinese Academy of Sciences.
All rights reserved.

The code can only be used for academic purpose, and cannot be used for commercial
purpose without written permission.

Redistribution and use in source for academic purpose, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <IPSR.h>
#include <graph_ipsr.h>
#include <ipsr_controller.h>
#include <o3d_api.h>
#include <nlohmann/json.hpp>
#include <configuration.h>
#include <graph_arg.h>



using namespace std;

unsigned int seed = 0;

double RESET_TRESHOLD = 0.1;
int POINT_RESET_ITER = 500;
int FACE_RESET_ITER = 500;

string input_data_root = "D:/WorkData/ipsr_explore/input/";
string out_data_root = "D:/WorkData/ipsr_explore/out/";
typedef double REAL;
const unsigned int DIM = 3U;
int psr_threads = 4;


void classic_ipsr(const std::string& input_name, const std::string& output_path, int iters, double pointweight, int depth, int k_neighbors);
void graph_ipsr(const std::string& input_name, const std::string& output_path, int iters, double pointweight, int depth, int k_neighbors);
void bitree_ipsr(const std::string& input_name, const std::string& output_path, int iters, double pointweight, int depth, int k_neighbors);
void test_graph(const std::string& input_name, const std::string& output_path, int iters, double pointweight, int depth, int k_neighbors);

std::map <IPSR_TYPE, IPSR_ENTRANCE> IPSR_ENTRANCE_MAP = {
	{IPSR_TYPE::CLASSIC_IPSR, classic_ipsr},
	{IPSR_TYPE::BITREE_IPSR, bitree_ipsr},
	{IPSR_TYPE::GRAPH_IPSR, graph_ipsr},
	{IPSR_TYPE::TEST_GRAPH, test_graph}
};

nlohmann::json time_cost;
lzd_tools::Timer timer(time_cost, "main");

static void output_logj(std::string path) {
	std::ofstream out(path);
	if (!out.is_open()) {
		std::cout << "ERROR: save log failed" << std::endl;
		return;
	}
	timer.Done();
	global_var::metric_j["time_cost"] = time_cost;
	global_var::metric_j["accumulat time cost"] = lzd_tools::AcumutableTimer::get_time_map();
	nlohmann::json log_j;
	log_j["config"] = global_var::config_j;
	log_j["metric"] = global_var::metric_j;
	out << log_j.dump(4);
	out.close();
	
	std::string floder = path.substr(0, path.find_last_of("/"));
	std::string name = path.substr(path.find_last_of("/") + 1);
	std::string log_path = floder + "/log/" + name;
	rmkdir(floder + "/log/");
	std::ofstream out_log(log_path);
	if (!out_log.is_open()) {
		std::cout << "ERROR: save log failed" << std::endl;
		return;
	}
	nlohmann::json log_j2;
	log_j2["config"] = global_var::config_j;
	log_j2["log"] = global_var::res_log_j;
	out_log << log_j2.dump(4);
	out_log.close();
}

static void save_config(std::string config_name){
	nlohmann::json config_j = ConfigManager::get_config_in_config_floder(config_name);
	rmkdir(global_var::data_output_base + "/config/");
	std::string path = global_var::data_output_base + "/config/" + config_name;
	std::ofstream out(path);
	if (!out.is_open()) {
		std::cout << "ERROR: save config failed" << std::endl;
		return;
	}
	out << config_j.dump(4);
}

void classic_ipsr(const string& input_name, const string& output_path, int iters, double pointweight, int depth, int k_neighbors) {
	global_var::config_j["mode"] = "classic_ipsr";
	nlohmann::json classic_ipsr_confg = ConfigManager::get_config_in_config_floder("classic_ipsr.json");
	save_config("classic_ipsr.json");
	classic_ipsr_confg["RandomInitConf"]["seed"] = seed;
	IPSR_Factory<REAL,DIM> ipsr_factory(input_name, output_path, iters, pointweight, depth, k_neighbors, seed, psr_threads);

	NormalEstimation<REAL, DIM>* estimator = get_estimator_from_json<REAL, DIM>(classic_ipsr_confg["estimator"]);
	
	IPSR_HANDLE_P ipsrhp = ipsr_factory.create_single_ipsr_handle(get_update_plan_by_name(classic_ipsr_confg["update_plan"]), estimator);

	IpsrController<REAL, DIM>* ipsr_controller  = new SingleIpsrController<REAL,DIM>(classic_ipsr_confg["save_op"], classic_ipsr_confg["save_mesh"]);
	
	int ep = -1;
	do {
		if (ConfigManager::get_common_config()["save_grid"]) {
			ipsr_factory.save_grid(ipsrhp->_points_normals,ipsrhp->get_out_put_base("/temp/") + to_string(ep) + ".grid");
		}
		ep = ipsr_controller->iter(ipsrhp);
		printf("\r %d", ep);
	} while (ep != -1); 
	ipsrhp->save_res();
	global_var::config_j["ipsr_controller_config"] = ipsr_controller->get_config();
	global_var::config_j["ipsr_config"] = ipsrhp->get_config();
	global_var::config_j["init_func"] = estimator->get_config();
	global_var::metric_j["res"] = ipsrhp->get_metric();
	global_var::res_log_j["ipsr_log"] = ipsrhp->get_log();
	output_logj(output_path+ipsrhp->get_resname()+"_log.json");
}

void bitree_ipsr(const string& input_name, const string& output_path, int iters, double pointweight, int depth, int k_neighbors) {
	printf("bitree_ipsr is not implemented in this build.\n");
}

void graph_ipsr(const string& input_name, const string& output_path, int iters, double pointweight, int depth, int k_neighbors) {
	GRAPH_IPSR::assign_static_conf();
	global_var::config_j["mode"] = "graph_ipsr";
	nlohmann::json sp_config;
	nlohmann::json com_conf = ConfigManager::get_common_config();
	nlohmann::json graph_ipsr_conf = ConfigManager::get_config_in_config_floder("graph_ipsr.json");
	save_config("graph_ipsr.json");
	graph_ipsr_conf["estimator"]["RandomInitConf"]["seed"] = seed;

	int seg_depth = graph_ipsr_conf["seg_sample_depth"];
	IPSR_Factory<REAL, DIM> seg_ipsr_factory(input_name, output_path, iters, pointweight, seg_depth, k_neighbors, seed, psr_threads);
	IPSR_Factory<REAL, DIM> global_ipsr_factory(input_name, output_path, iters, pointweight, depth, k_neighbors, seed, psr_threads);
	
	NormalEstimation<REAL, DIM>* estimator = get_estimator_from_json<REAL, DIM>(graph_ipsr_conf["estimator"]);
	
	auto update_plan = get_update_plan_by_name(graph_ipsr_conf["update_plan"]);
	
	vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>> raw_points_normals;
	string input_name_base = input_name.substr(0, input_name.find_last_of('.'));
	ply_reader<REAL, DIM>(input_name, raw_points_normals);

	POINTS_NORMALS unique_ori_op;
	int dpn = lzd_tools::duplicate_filter(raw_points_normals,unique_ori_op);
	printf("delect %d duplicate points\n", dpn);
	vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>> ori_points_normals = unique_ori_op;
	
	
	if (depth > 0) {
		printf("global downsampling...\n");
		ori_points_normals = global_ipsr_factory.down_sample_points(unique_ori_op, depth);
	}
	else {
		printf("skip global downsampling\n");
	}

	lzd_tools::Timer allocate_node(time_cost, "allocate_node");
	GRAPH_IPSR::NodeAllocator<REAL, DIM>* al = GRAPH_IPSR::createAllocatorByJson(graph_ipsr_conf["NodeAllocator"], ori_points_normals);
	allocate_node.Done();

	lzd_tools::Timer orchectration(time_cost, "orchectration and ipsr initialization");
	GRAPH_IPSR::linked_graph<REAL, DIM> ipsr_graph = al->op_orchestration(seg_ipsr_factory, update_plan, estimator);
	
	sp_config["ipsr_spilter"] = al->get_config();
	ipsr_graph.update_all_op();
	if (!ConfigManager::get_common_config()["no_save"]) {
		ipsr_graph.save_seg("/graph_ipsr/lg_seg/", "init");
		ipsr_graph.save_gt("/graph_ipsr/lg_seg/");
	}
	
	IpsrController<REAL, DIM>* ipsr_controller = new SingleIpsrController<REAL,DIM>(-1,-1);
	
	
	if (com_conf["save_option"]["init_normal"] == true){
		POINTS_NORMALS init_op;
		ipsr_graph.get_all_op(init_op);
		assert(init_op.size() == ori_points_normals.size());
		std::string init_op_path = global_ipsr_factory.get_out_put_base() + global_ipsr_factory.get_resname() + "_init_op.ply";
		lzd_tools::op2ply(init_op, init_op_path, XForm<REAL, DIM + 1>().Identity());		
	}
	orchectration.Done();

	
	lzd_tools::Timer graph_ipsr(time_cost, "graph ipsr ");
	if (iters < 0) {
		ipsr_graph.load_seg(global_ipsr_factory.get_resname(),graph_ipsr_conf["refer_arg_flag"]);
	}
	else{
		for (int i = 0; i < iters; i++) {
			ipsr_graph.update(ipsr_controller);
			printf("\r graph iters %d/%d", i,iters);
		}
		ipsr_graph.update_all_op();
	}

	if (!ConfigManager::get_common_config()["no_save"]) {
		ipsr_graph.save_seg("/graph_ipsr/lg_seg/", "res");
	}
	graph_ipsr.Done();


	lzd_tools::Timer cal_edge(time_cost, "cal edge");
	GRAPH_IPSR::edge_weight_calculator<REAL, DIM>* consist_calculator = GRAPH_IPSR::get_edge_caculator_from_json<REAL, DIM>(graph_ipsr_conf["edge_weight_calculator"], &seg_ipsr_factory,&ipsr_graph);
	ipsr_graph.set_edge_weight_calculator(consist_calculator);
	ipsr_graph.cal_edges_weight();
	if(graph_ipsr_conf["save_wrong_edge"] == true)
	{
		ipsr_graph.save_wrong_edge();
	}
	if (graph_ipsr_conf["save_all_edge"] == true)ipsr_graph.save_all_edge();
	cal_edge.Done();


	lzd_tools::Timer flip_nodes(time_cost, "flip nodes");
	std::string path;
	printf("start flip node...\n");
	std::string arg_name = graph_ipsr_conf["FlipArgName"];
	graph_arg::FlipGraph* flip_alg_opt = graph_arg::FlipGraph::get_flip_arg_from_json(graph_ipsr_conf["FlipArgConf"][arg_name]);
	ipsr_graph.flip_graph(flip_alg_opt);
	flip_nodes.Done();

	POINTS_NORMALS res_op(raw_points_normals.size());
	lzd_tools::Timer cal_res(time_cost, "output");

	if (!ConfigManager::get_common_config()["no_save"]) {
		ipsr_graph.save_res("", "", &global_ipsr_factory);
	}
	ipsr_graph.query_res(raw_points_normals,res_op);
	global_var::metric_j["res"] = lzd_tools::PointNormalMetric<REAL, DIM>(raw_points_normals, res_op).to_json();
	if (com_conf["save_option"]["ori_est"]) {
		lzd_tools::op2ply(res_op, global_ipsr_factory.get_out_put_base() + global_ipsr_factory.get_resname() + "_ori_est.ply");
	}
	if (!ConfigManager::get_common_config()["no_save"]) {
		ipsr_graph.save_optim();
		ipsr_graph.query_res(raw_points_normals, res_op);
		global_var::metric_j["optim_res"] = lzd_tools::PointNormalMetric<REAL, DIM>(raw_points_normals, res_op).to_json();
	}

	path = global_ipsr_factory.get_out_put_base() + global_ipsr_factory.get_resname() + "_log.json";
	sp_config["config"] = ipsr_graph.get_config();
	sp_config["init_func"] = estimator->get_config();
	global_var::config_j["graph_ipsr_config"] = sp_config;
	global_var::config_j["ipsr_controller_config"] = ipsr_controller->get_config();
	global_var::config_j["ipsr_config"] = ipsr_graph._nodes[0]->_handle->get_config();
	global_var::res_log_j["graph_ipsr_log"] = ipsr_graph.get_log();
	output_logj(path);
}

void test_graph(const string& input_name, const string& output_path, int iters, double pointweight, int depth, int k_neighbors) {
	printf("test_graph is not implemented in this build.\n");
}

void betweenness(const string& input_name, const string& output_path, int iters, double pointweight, int depth, int k_neighbors) {

}

int main(int argc, char *argv[])
{
	string input_name;
	int iters = 10;
	double pointweight = 10;
	int depth = 10;
	int k_neighbors = 10;
	IPSR_TYPE ipsr_type = IPSR_TYPE::CLASSIC_IPSR;
	for (int i = 1; i < argc; i += 2)
	{
		if (strcmp(argv[i], "--in") == 0)
		{
			input_name = argv[i + 1];
			}
		else if (strcmp(argv[i], "--iters") == 0)
		{
			long v = strtol(argv[i + 1], nullptr, 10);
			iters = static_cast<int>(v);
		}
		else if (strcmp(argv[i], "--pointWeight") == 0)
		{
			pointweight = strtod(argv[i + 1], nullptr);
			if (pointweight < 0.0 || pointweight == HUGE_VAL || pointweight == -HUGE_VAL)
			{
				printf("invalid value of --pointWeight");
				return 0;
			}
		}
		else if (strcmp(argv[i], "--depth") == 0)
		{
			long v = strtol(argv[i + 1], nullptr, 10);
			if (!valid_parameter(v))
			{
				printf("WARNING depth<0");
			}
			depth = static_cast<int>(v);
		}
		else if (strcmp(argv[i], "--neighbors") == 0)
		{
			long v = strtol(argv[i + 1], nullptr, 10);
			if (!valid_parameter(v))
			{
				printf("invalid value of --neighbors");
				return 0;
			}
			k_neighbors = static_cast<int>(v);
		}
		else if (strcmp(argv[i], "--input_data_root") == 0)
		{
			input_data_root = argv[i + 1];
		}
		else if (strcmp(argv[i], "--seed") == 0)
		{
			seed = strtol(argv[i + 1], nullptr, 10);
		}
		else if (strcmp(argv[i], "--out_data_root") == 0)
		{
			out_data_root = argv[i + 1];
		}
		else if (strcmp(argv[i], "--ipsr_type") == 0)
		{
			ipsr_type = static_cast<IPSR_TYPE>(strtol(argv[i + 1], nullptr, 10));
		}
		else if (strcmp(argv[i], "--psr_threads") == 0)
		{
			long v = strtol(argv[i + 1], nullptr, 10);
			if (v < 1)
			{
				printf("invalid value of --psr_threads");
				return 0;
			}
			psr_threads = static_cast<int>(v);
		}
		else if (strcmp(argv[i], "--config_folder") == 0)
		{
			std::string _config_folder = argv[i + 1];
			if (!_config_folder.empty() && _config_folder.back() != '/' && _config_folder.back() != '\\') {
				_config_folder += "/";
			}
			if(ConfigManager::check_config_floder(_config_folder)){
				ConfigManager::config_floder = _config_folder;
			}else{
				printf("invalid config path of %s\n", _config_folder.c_str());
				return 0;
			}
		}
		else
		{
			printf("unknown parameter of %s\n", argv[i]);
			return 0;
		}
	}

	if (argc <= 1 || input_name.empty())
	{
		// ask for parameters
		cout << "Please input the parameters:\n";
		cout << "--in                      input .ply model\n";
		cin >> input_name;
		cout << "--iters       maximum number of iterations\n";
		cin >> iters;
		cout << "--pointWeight screened weight of SPSR\n";
		cin >> pointweight;
		cout << "--depth       maximum depth of the octree\n";
		cin >> depth;
		cout << "--neighbors   number of the nearest neighbors to search\n";
		cin >> k_neighbors;
		cout << "--seed        random seed\n";
		cin >> seed;
		cout << "--ipsr_type   type of IPSR, \n0 for classic_ipsr\n1 for bitree_ipsr\n2 for graph_ipsr\n3 for test_graph\n";
		int temp;
		cin >> temp;
		ipsr_type = static_cast<IPSR_TYPE>(temp);
		cout << "--psr_threads  number of PoissonRecon threads\n";
		cin >> psr_threads;
		std::string config_floder;
		cout << "--config_floder config file floder\n";
		cin >> config_floder;
		if (!config_floder.empty() && config_floder.back() != '/' && config_floder.back() != '\\') {
			config_floder += "/";
		}
		if (ConfigManager::check_config_floder(config_floder))ConfigManager::config_floder = config_floder;
		else {
			cout << "invaild config_floder\n";
			assert(false);
		}
	}
	
	string input_path = input_data_root + input_name + ".ply";

	if (!std::filesystem::exists(input_path)) {
		string user_input_path = input_name;
		if (std::filesystem::exists(user_input_path)) {
			input_path = user_input_path;
			input_name = input_name.substr(input_name.find_last_of('/') + 1);
			input_name = input_name.substr(input_name.find_last_of('\\') + 1);
			input_name = input_name.substr(0, input_name.find_last_of('.')); 
			printf("use %s as input_name\n", input_name.c_str());
		}
		else {
			printf("\ninput file %s does not exist\n", input_path.c_str());
			printf("use %s as input file directly\n", user_input_path.c_str());
			printf("Current directory: %s\n", std::filesystem::current_path().string().c_str());
			return 0;
		}
	}

	IPSR_Factory<REAL, DIM> ipsr_factory(input_path, out_data_root, iters, pointweight, depth, k_neighbors, seed, psr_threads);

	string output_path = out_data_root + "/" + ipsr_factory.get_resname() + "/" + input_name + "/";
	rmkdir(output_path);

	printf("\nIterative Poisson Surface Reconstruction (iPSR)\n");
	printf("Parameters:\n");
	printf("--in          %s\n", input_path.c_str());
	printf("--out         %s\n", output_path.c_str());
	printf("--iters       %d\n", iters);
	printf("--pointWeight %f\n", pointweight);
	printf("--depth       %d\n", depth);
	printf("--neighbors   %d\n\n", k_neighbors);


	global_var::config_j["iters"] = iters;
	global_var::config_j["pointweight"] = pointweight;
	global_var::config_j["depth"] = depth;
	global_var::config_j["k_neighbors"] = k_neighbors;
	global_var::config_j["seed"] = seed;
	global_var::config_j["psr_threads"] = psr_threads;
	global_var::config_j["input_path"] = input_path;
	global_var::config_j["output_path"] = output_path;

	global_var::data_output_base = output_path;

	timer.Start();

	save_config("common.json");
	IPSR_ENTRANCE ipsr = IPSR_ENTRANCE_MAP[ipsr_type];
	ipsr(input_path, output_path, iters, pointweight, depth, k_neighbors);

	printf("%s", lzd_tools::AcumutableTimer::get_time_map().dump(4).c_str());
	printf("Main finish sussfully!!!");
	return 0;
}
