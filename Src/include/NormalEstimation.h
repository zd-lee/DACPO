#pragma once
#include <utility.h>
#include <tools.h>
#include <o3d_api.h>
#include <pcl_api.h>
#include <configuration.h>
#include <socket_api.h>
#include <graph_arg.h>
#include <dipole.h>

template<typename REAL,int DIM>
class NormalEstimation{
public:
    virtual int Estimate(ORIENTED_POINTS& op) = 0;
    virtual nlohmann::json get_config() = 0;
};
template<typename REAL, int DIM>
class DoingNothing : public NormalEstimation<REAL, DIM>{
public:
    nlohmann::json get_config(){
        nlohmann::json j;
        j["name"] = "DoingNothing";
        return j;
    }

    int Estimate(ORIENTED_POINTS& op){
        return 0;
    }
};

template<typename REAL, int DIM>
class RandomInit : public NormalEstimation<REAL, DIM>{
    int _seed;
public:
    RandomInit(int seed):_seed(seed){}
    nlohmann::json get_config(){
        nlohmann::json j;
        j["name"] = "RandomInit";
        j["seed"] = _seed;
        return j;
    }

    int Estimate(ORIENTED_POINTS& op){
        Normal<REAL, DIM> zero_normal(Point<REAL, DIM>(0, 0, 0));
        srand(_seed);
        for(int i = 0;i<op.size();i++){
            do{
                op[i].second = Point<REAL, DIM>(rand()%1001, rand()%1001 - 500.0, rand()%1001 - 500.0);

            }while(op[i].second == zero_normal);
            normalize(op[i].second);
        }
        return 0;
    }
};


template<typename REAL, int DIM>
class Hoppe1994 : public NormalEstimation<REAL, DIM> {
    int _K;
    REAL _lambda;
    REAL _cos_theta;
    NormalEstimation<REAL, DIM>* _pre_estimator;

public:
    Hoppe1994(int K, REAL lambda = 0, REAL cos_theta = 0 ) :
        //_pre_estimator(pre_estimator),
        _K(K), _lambda(lambda), _cos_theta(cos_theta) {}


    nlohmann::json get_config() {
        nlohmann::json j;
        j["name"] = "Hoppe1994";
        //j["pre_estimator"] = _pre_estimator->get_config();
        j["K"] = _K;
        j["lambda"] = _lambda;
        j["cos_theta"] = _cos_theta;
        return j;
    }

    int Estimate(ORIENTED_POINTS& points_normals) {
        //NormalEstimation<REAL, DIM>* pre_e = new RandomInit<REAL,DIM>(0);
        //pre_e->Estimate(points_normals);
        //_pre_estimator->Estimate(points_normals);// ΪOrientNormalsConsistentTangentPlane�ṩ������
        Normal<REAL, DIM> zero_normal(Point<REAL, DIM>(0, 0, 0));
        open3d::geometry::PointCloud pcd;
        pcd.points_.resize(points_normals.size());
        pcd.normals_.resize(points_normals.size());
        for (int i = 0; i < points_normals.size(); i++) {
            pcd.points_[i] = Eigen::Vector3d(points_normals[i].first[0], points_normals[i].first[1], points_normals[i].first[2]);
            //pcd.normals_[i] = Eigen::Vector3d(points_normals[i].second.normal[0], points_normals[i].second.normal[1], points_normals[i].second.normal[2]);
        }
        pcd.EstimateNormals(open3d::geometry::KDTreeSearchParamKNN());
        //o3d_norm_estimate(points_normals);;
        pcd.OrientNormalsConsistentTangentPlane(_K, _lambda, _cos_theta);
        for (int i = 0; i < points_normals.size(); i++) {
            Normal<REAL, DIM> n;
            for (int j = 0; j < DIM; j++)n.normal[j] = pcd.normals_[i][j];
            points_normals[i].second = n;
            assert(!(points_normals[i].second == zero_normal));
        }
        return 0;
    }

};


template<typename REAL, int DIM>
class PclNormalEstimation : public NormalEstimation<REAL, DIM> {
    int _k;
public:
    PclNormalEstimation(int k = 5) :_k(k) {}
    int Estimate(ORIENTED_POINTS& op) {
        NormalEstimation<REAL, DIM>* pre_e = new RandomInit<REAL, DIM>(0);
        pre_e->Estimate(op);
        Normal<REAL, DIM> zero_normal(Point<REAL, DIM>(0, 0, 0));
        pcl::PointCloud<pcl::PointNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointNormal>);
        cloud->points.resize(op.size());
        for (int i = 0; i < op.size(); i++) {
            cloud->points[i].x = op[i].first[0];
            cloud->points[i].y = op[i].first[1];
            cloud->points[i].z = op[i].first[2];
        }
        pcl::NormalEstimation<pcl::PointNormal, pcl::PointNormal> ne;
        ne.setInputCloud(cloud);
        pcl::search::KdTree<pcl::PointNormal>::Ptr tree(new pcl::search::KdTree<pcl::PointNormal>());
        ne.setSearchMethod(tree);
        pcl::PointCloud<pcl::PointNormal>::Ptr cloud_normals(new pcl::PointCloud<pcl::PointNormal>);
        ne.setKSearch(_k);
        ne.compute(*cloud_normals);
        for (int i = 0; i < op.size(); i++) {
            Normal<REAL, DIM> n;
            n.normal[0] = cloud_normals->points[i].normal_x;
            n.normal[1] = cloud_normals->points[i].normal_y;
            n.normal[2] = cloud_normals->points[i].normal_z;
            op[i].second = n;
            assert(!(op[i].second == zero_normal));
        }
        return 0;
    }

    nlohmann::json get_config() {
        nlohmann::json j;
        j["name"] = "PclNormalEstimation";
        j["k"] = _k;
        return j;
    }
};

template<typename REAL, int DIM>
class PCAEstimation : public NormalEstimation<REAL, DIM> {
    int _k;
public:
    PCAEstimation(int k = 10) :_k(k) {}
    int Estimate(ORIENTED_POINTS& points_normals) {
        Normal<REAL, DIM> zero_normal(Point<REAL, DIM>(0, 0, 0));
        open3d::geometry::PointCloud pcd;
        pcd.points_.resize(points_normals.size());
        pcd.normals_.resize(points_normals.size());
        for (int i = 0; i < points_normals.size(); i++) {
            pcd.points_[i] = Eigen::Vector3d(points_normals[i].first[0], points_normals[i].first[1], points_normals[i].first[2]);
        }
        pcd.EstimateNormals(open3d::geometry::KDTreeSearchParamKNN(_k));
        for (int i = 0; i < points_normals.size(); i++) {
            Normal<REAL, DIM> n;
            for (int j = 0; j < DIM; j++)n.normal[j] = pcd.normals_[i][j];
            points_normals[i].second = n;
            assert(!(points_normals[i].second == zero_normal));
        }
        return 0;
    }

    nlohmann::json get_config() {
        nlohmann::json j;
        j["name"] = "PCANormalEstimation";
        j["k"] = _k;
        return j;
    }
};


// TODO
template<typename REAL, int DIM>
class FieldPropagationEstimate: public NormalEstimation<REAL, DIM>{
    NormalEstimation<REAL, DIM>* _pre_estimator;
    FIELD_PROPAGATION::MaskFieldPropagation _propagator;
    
public:

    FieldPropagationEstimate(NormalEstimation<REAL, DIM>* pre_estimator, int C, double eps, int times, int mask_size, int knngraph_k, double knngraph_radius) :{
		_pre_estimator = pre_estimator;
		_propagator = FIELD_PROPAGATION::MaskFieldPropagation(FIELD_PROPAGATION::MaskField(C,eps), times, mask_size, knngraph_k, knngraph_radius);
    }

	FieldPropagationEstimate(nlohmann::json config) :_propagator(FIELD_PROPAGATION::MaskField(config["C"], config["eps"]), config["times"], config["knn_mask"], config["k_neighbors"], config["r"], config["diffuse"])  
    {
		_pre_estimator = get_estimator_from_json<REAL, DIM>(config["pre_estimator"]);
	}

    nlohmann::json get_config(){
        nlohmann::json j;
		j["name"] = "FieldPropagationEstimate";
		j["pre_estimator"] = _pre_estimator->get_config();
		j["propagator"] = _propagator.get_config();
        return j;
    }

    int Estimate(ORIENTED_POINTS& op){
		_pre_estimator->Estimate(op);
		Eigen::MatrixXd pts(op.size(), 6);
        // convert to Eigen::MatrixXd
#pragma omp parallel for
        for (int i = 0; i < op.size(); i++) {
            Eigen::RowVectorXd row(6);
            row << op[i].first[0], op[i].first[1], op[i].first[2], op[i].second.normal[0], op[i].second.normal[1], op[i].second.normal[2];
			pts.row(i) = row;
        }
		std::vector<bool> flip_status = _propagator.xie_propagation_points_on_bfs_tree(pts);
		for (int i = 0; i < op.size(); i++) {
			if (flip_status[i])op[i].second *= -1;
		}
		return 0;
    }
};

template<typename REAL, int DIM>
class SocketNormalEstimation : public NormalEstimation<REAL, DIM> {
    std::string _ip;
    int _port;

    std::string FUNCTION_NAME;
    nlohmann::json _func_config; 

    int _estimate(ORIENTED_POINTS& op) {
        lzd_tools::Socket socket(_ip, _port);//每次都要重新连接
		if (!socket.connectToServer()) {
			std::cerr << "Failed to connect to server" << std::endl;
			assert(false);
			return -1;
		}

        // 发送请求头
        nlohmann::json request;
        request["function_name"] = FUNCTION_NAME;
        request["function_config"] = _func_config;
        request["data_size"] = op.size();
        if (!socket.SendJson(request.dump())) {
			assert(false);
			std::cerr << "Failed to send request" << std::endl;
			return -1;
		}
        // 接受确认
        std::string response_buff;
        nlohmann::json response;
        if(!socket.ReceiveJson(response_buff)) {
			assert(false);
			std::cerr << "Failed to receive response" << std::endl;
			return -1;
		}
        response = nlohmann::json::parse(response_buff);
        if (response["status"] != "OK") {
			assert(false);
			std::cerr << "Error: " << response["error"] << std::endl;
			return -1;
		}
        // 发送数据
		std::vector<double> data;
		for (int i = 0; i < op.size(); i++) {
			for (int j = 0; j < DIM; j++) {
				data.push_back(op[i].first[j]);
			}
		}
		if (!socket.SendDoubleArray(data)) {
			assert(false);
			std::cerr << "Failed to send data" << std::endl;
			return -1;
		}
		// 接受数据
        data.clear();
        data.resize(op.size() * DIM * 2);
        if (!socket.ReceiveDoubleArray(data, op.size() * DIM * 2)) {
            assert(false);
            std::cerr << "Failed to receive data" << std::endl;
            return -1;
        }
        double offset = 0;
        for (int i = 0; i < op.size(); i++) {
            for (int j = 0; j < DIM; j++) {
                offset += abs(op[i].first[j] - data[i * DIM * 2 + j]);
                op[i].second.normal[j] = data[i * DIM * 2 + DIM + j];
            }
        }
        if (offset > 1e-3) {
            printf("offset = %f\n", offset);
            std::cerr << "ERROR!! offset too large" << std::endl;
        }
        // 将零向量替换为估计的法向量
        Normal<REAL, DIM> zero_normal(Point<REAL, DIM>(0, 0, 0));
        int zero_count = 0;
        for (int i = 0; i < op.size(); i++) {
            if (op[i].second == zero_normal) {
                // 随机初始化
                do {
                    op[i].second = Point<REAL, DIM>(rand() % 1001, rand() % 1001 - 500.0, rand() % 1001 - 500.0);
                } while (op[i].second == zero_normal);
                normalize(op[i].second);
                zero_count++;
            }
        }
        if(zero_count > 0) {
            printf("Warning: %d / %d zero normals are replaced by random normals\n", zero_count, op.size());
        }
        return 0;

        
    }
public:
    SocketNormalEstimation(std::string function_name, nlohmann::json func_config) :FUNCTION_NAME(function_name), _func_config(func_config) {
        _ip = ConfigManager::get_common_config()["estimator_server"]["ip"];
        _port = ConfigManager::get_common_config()["estimator_server"]["port"];
    }

    // 允许并发调用
    int Estimate(ORIENTED_POINTS& op) {
        return _estimate(op);
    }

    nlohmann::json get_config() {
        nlohmann::json j;
        j["name"] = "SocketNormalEstimation";
        j["function_name"] = FUNCTION_NAME;
        j["function_config"] = _func_config;
        return j;
    }

};
    



template<typename REAL, int DIM>
NormalEstimation<REAL, DIM>* get_estimator_from_json(nlohmann::json estimator_config) {
    NormalEstimation<REAL, DIM>* estimator = new RandomInit<REAL, DIM>(0);
    if (estimator_config["name"] == "Hoppe1994") {
        nlohmann::json cj = estimator_config["Hoppe1994Conf"];
        int K = cj["K"];
        int lambda = cj["lambda"];
        REAL alpha = cj["alpha"];
        estimator = new Hoppe1994<REAL, DIM>(K, alpha, lambda);
    }
    else if (estimator_config["name"] == "DoingNothing") {
        estimator = new DoingNothing<REAL, DIM>();
    }
    else if (estimator_config["name"] == "PclNormalEstimation") {
        nlohmann::json cj = estimator_config["PclNormalEstimationConf"];
        int k = cj["k"];
        estimator = new PclNormalEstimation<REAL, DIM>(k);
    }
    else if (estimator_config["name"] == "SocketNormalEstimation") {
        nlohmann::json cj = estimator_config["SocketNormalEstimationConf"];
        std::string function_name = cj["function_name"];
        nlohmann::json func_config = cj["function_config"];
        estimator = new SocketNormalEstimation<REAL, DIM>(function_name, func_config);
    }
    else if(estimator_config["name"] == "RandomInit"){
        seed = estimator_config["RandomInitConf"]["seed"];
        estimator = new RandomInit<REAL, DIM>(seed);
	}
	else if (estimator_config["name"] == "FieldPropagationEstimate") {
        estimator = new FieldPropagationEstimate<REAL,DIM>(estimator_config["FieldPropagationEstimateConf"]);
	}else if (estimator_config["name"] == "PCA") {
        estimator = new PCAEstimation<REAL, DIM>(estimator_config["PCAConf"]["k"]);
    }
    else {
        cout << "Error: estimator name not found" << endl;
        cout << "Using default estimator: RandomInit" << endl;
    }
    return estimator;
};

