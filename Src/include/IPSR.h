#pragma once 
#include <array>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include "pcl_api.h"
#include "o3d_api.h"
#include "kdtree.h"
#include "utility.h"
#include "tools.h"
#include "orientation_consistency.h"
#include "NormalEstimation.h"
#include "ipsr_op_filter.h"

#define INDIRECT_MAPPING 0
#define DIRECT_MAPPING 1

#define IPSR_HANDLE_P ipsr_handle<REAL,DIM>* 

using std::array;
using std::endl;
using std::to_string;

struct Period{
    int* p_ep;
    double * p_old_var;
    double * p_new_var;
    std::vector<std::string>* p_cmd;
    int* update_strategy;
    int print_op_freq = 5;
    int max_iters = 30;
    int print_mesh = 1000;

    Period();
};

extern std::map<int (*)(Period& p),std::string> update_plane_namelist;

int interact_update(Period& p);
int diri_ind_ipsr(Period& p);
int diri_dir_ipsr(Period& p);
int neumman_ind_ipsr(Period& p);
int neumman_dir_ipsr(Period& p);

int comb_plan5(Period& p);

int comb_update6(Period& p);

int comb_update7(Period& p);

int comb_update8(Period& p);

int specify_update(Period& p);

typedef int (*UpdatePlan)(Period&);

UpdatePlan get_update_plan_by_name(std::string name);

//        REAL theta = calculate_angle(points_normals[i].second, pcl_normals[i].second);
//        REAL bias = 20;
//            points_normals[i].second = get_rand_norm<REAL, DIM>(0);
//            count_fix++;
//            flag[i] = 1;
lzd_tools::FiniteMap get_mincut_colormap();

template<typename REAL, unsigned int DIM>
MESH shrink_boundary(const MESH& mesh, const POINTS_NORMALS& op, std::vector<REAL>& pmdist, int k = 20);

/****************************************
 * ******************************************/
template<typename REAL, unsigned int DIM>
class ipsr_handle {
    nlohmann::json _config_j;
    // nlohmann::json _log_j;

    std::vector<Update_filter<REAL, DIM>*> _filters;

public:
    std::vector<std::string> _cmd;
    int _max_iters;
    int _k_neighbors;
    int _seed;
    int (*update_plan)(Period& p);
    POINTS_NORMALS _ori_points_normals;
public:
    POINTS_NORMALS _points_normals;
    POINTS_NORMALS _gt_points_normals;
    POINTS_NORMALS _res_points_normals;
    POINTS_NORMALS _init_normals;
    int _epoch;
    Period _p_;
    int _update_strategy;
    REAL avg_max_diff;
    REAL old_avg_max_diff;
    XForm<REAL, DIM + 1> ixform;
    MESH mesh;
    std::vector<Normal<REAL, DIM>> avg_normal;
    std::vector<float> i_nearestSample_var;
    std::vector<float> _points_normals_diff;
    std::vector<double> _weight_samples;
    kdt::KDTree<kdt::KDTreePoint> input_point_tree;
    std::vector<std::vector<int>> nearestSamples;

    std::ostringstream _out;
    std::vector <REAL> avg_max_var_log;
    std::vector <REAL> sample_var_log;
    std::vector <REAL> loss_log;
    std::vector <REAL> avg_pmdist_log;
    std::vector <int> mesh_size_log;
    std::vector <std::string> behavior_log;

    nlohmann::json log_j;
    nlohmann::json mj;

    std::vector<IPSR_HANDLE_P> spilt_hps;
    std::vector<std::vector<int>> _index_mapper;
    bool _if_spilt = false;

    std::vector<int> dlabel;

public:
    ipsr_handle<REAL, DIM>(
        const std::vector<std::string> cmd,
        int max_iters,
        int k_neighbors,
        int (*update_plan)(Period& p),
        const POINTS_NORMALS& points_normals,
        int seed = 0
    ) {
        this->_cmd = std::vector<std::string>(cmd.size());
        for (int i = 0; i < cmd.size(); i++) {
            this->_cmd[i] = cmd[i];
            
        }
        this->_max_iters = max_iters;
        this->update_plan = update_plan;
        _k_neighbors = k_neighbors;
        _seed = seed;


        this->_p_ = Period();
        this->_p_.p_ep = &this->_epoch;
        this->_p_.p_old_var = &this->old_avg_max_diff;
        this->_p_.p_new_var = &this->avg_max_diff;
        this->_p_.p_cmd = &(this->_cmd);
        this->_p_.update_strategy = &this->_update_strategy;
        this->_p_.max_iters = this->_max_iters;

        assert (points_normals.size() > 0);
        _ori_points_normals = points_normals;

        this->_epoch = -1;//epoch == -1 before init
        this->avg_max_diff = 0;
        this->old_avg_max_diff = 0;
    }

    ipsr_handle<REAL, DIM>(const ipsr_handle<REAL, DIM>& ipsr) {
        this->_cmd = std::vector<std::string>(ipsr._cmd.size());
        for (int i = 0; i < ipsr._cmd.size(); i++) {
            this->_cmd[i] = ipsr._cmd[i];
        }
        this->_max_iters = ipsr._max_iters;
        this->update_plan = ipsr.update_plan;
        this->_ori_points_normals = ipsr._ori_points_normals;
        this->_points_normals = ipsr._points_normals;
        this->_res_points_normals = ipsr._res_points_normals;
        this->_k_neighbors = ipsr._k_neighbors;
        this->_seed = ipsr._seed;

        this->_p_ = Period();
        this->_p_.p_ep = &this->_epoch;
        this->_p_.p_old_var = &this->old_avg_max_diff;
        this->_p_.p_new_var = &this->avg_max_diff;
        this->_p_.p_cmd = &(this->_cmd);
        this->_p_.update_strategy = &this->_update_strategy;
        this->_p_.max_iters = this->_max_iters;
        this->_epoch = ipsr._epoch;// 
        this->_update_strategy = ipsr._update_strategy;
        this->avg_max_diff = ipsr.avg_max_diff;
        this->old_avg_max_diff = ipsr.old_avg_max_diff;
        this->ixform = ipsr.ixform;
        this->mesh = ipsr.mesh;
        this->avg_normal = ipsr.avg_normal;
        this->i_nearestSample_var = ipsr.i_nearestSample_var;
        this->_points_normals_diff = ipsr._points_normals_diff;
        this->_weight_samples = ipsr._weight_samples;
        this->_gt_points_normals = ipsr._gt_points_normals;


        std::vector<kdt::KDTreePoint> vertices;
        vertices.reserve(this->_points_normals.size());
        for (size_t i = 0; i < this->_points_normals.size(); ++i)
        {
            array<REAL, 3> _p_{this->_points_normals[i].first[0], this->_points_normals[i].first[1], this->_points_normals[i].first[2]};
            vertices.push_back(kdt::KDTreePoint(_p_));
        }
        input_point_tree.build(vertices);
    }

    // init_type: 0:nothing 1:random 2:pca 3:pcl


    /**
     * @brief 
     * @return int 
     */
    int ipsr_initialization(){
        _epoch = 0;
        std::vector<char* > argv_str(this->_cmd.size());
        for (int i = 0; i < this->_cmd.size(); i++) {
            argv_str[i] = &this->_cmd[i][0];
        }

        this->_points_normals = sample_points_entrance<REAL, DIM>((int)argv_str.size(), argv_str.data(), _ori_points_normals, this->ixform, &this->_weight_samples);
        if (_points_normals.size() == 0) {
            printf("WARNING::ipsr_handle::init::sample_points.size() == 0\n");
            this->ixform = XForm<REAL, DIM + 1>().Identity();
            this->_points_normals = _ori_points_normals;
            _weight_samples = std::vector<double>(_points_normals.size(), 1.0);
        }

        assert(ixform.determinant() != 0);
        assert(ixform.determinant() != NAN);
        _gt_points_normals.reserve(_points_normals.size());
        for (size_t i = 0; i < _points_normals.size(); i++) _gt_points_normals.push_back(_points_normals[i]);
        std::vector<kdt::KDTreePoint> vertices;
        vertices.reserve(this->_points_normals.size());
        for (size_t i = 0; i < this->_points_normals.size(); ++i)
        {
            array<REAL, 3> _p_{ this->_points_normals[i].first[0], this->_points_normals[i].first[1], this->_points_normals[i].first[2] };
            vertices.push_back(kdt::KDTreePoint(_p_));
        }
        input_point_tree.build(vertices);
        _out = std::ostringstream();

        nlohmann::json filter_conf = ConfigManager::get_common_config()["update_filter"];
        nlohmann::json namelist = filter_conf["namelist"];
        for (auto name : namelist) {
            if(name == "NonZero"){
                _filters.push_back(new NonZeroFilter<REAL, DIM>());
            }
            else if(name == "InitFilter"){
                _filters.push_back(new InitFilter<REAL, DIM>(filter_conf["config_table"]["InitFilter"],_points_normals));
            }
            else {
                printf("ipsr_handle::init::unknow filter name\n");
            }
        }

        for (auto f:_filters)
        {
            _out<<f->get_config()<<endl;
        }
        return 0;
    }
    
    /**
     * */
    int init(NormalEstimation<REAL, DIM>* estimator)
    {
        ipsr_initialization();
        init_op_normal(estimator);
        return 0;
    }

    void init_op_normal(NormalEstimation<REAL,DIM>* estimator){
        estimator->Estimate(_points_normals);
        log_j["init_log"].push_back(estimator->get_config());
    }

    nlohmann::json get_config(){
        nlohmann::json j;
        j["max_iters"] = _max_iters;
        j["k_neighbors"] = _k_neighbors;
        j["seed"] = _seed;
        j["update_plan"] = update_plane_namelist[update_plan];
        return j;
    }

    nlohmann::json get_log(){
        nlohmann::json tlog;
        for(int i = 0;i<avg_max_var_log.size();i++){
            tlog[to_string(i)] = to_string(avg_max_var_log[i]);
        }
        log_j["avg_max_var_log"] = tlog;
        log_j["behavior"] = _out.str();
        return log_j;
    }

    nlohmann::json get_metric(){
        nlohmann::json j;
        j = lzd_tools::PointNormalMetric<REAL,DIM>(_points_normals, _gt_points_normals).to_json();
        return j;
    }
    REAL cal_avg_loss(bool permid_nd = true) {
        REAL avg_loss = 0;
        for (size_t i = 0; i < _points_normals.size(); i++) {
            avg_loss += calculate_angle(_gt_points_normals[i].second, _points_normals[i].second) * 180 / M_PI;
        }
        avg_loss /= _points_normals.size();
        if (permid_nd)avg_loss = std::min(avg_loss, 180 - avg_loss);
        return avg_loss;
    }

    int cal_and_log_matric(int top_k = -1) {
        lzd_tools::AcumutableTimer cal_matric("ipsr_handle::cal_and_log_matric");

        if (top_k <= 0)top_k = _points_normals.size();
        REAL avg_loss = cal_avg_loss();
        loss_log.push_back(avg_loss);

        mesh_size_log.push_back(mesh.second.size());

        REAL avg_pmdist;
        std::vector<REAL> pmdists;
        shrink_boundary(this->mesh,this->_points_normals,pmdists,1);
        if (pmdists.size() == 0) {
            avg_pmdist = 0;
        }
        else {
            sort(pmdists.begin(), pmdists.end(), std::greater<REAL>());
            avg_pmdist = std::accumulate(pmdists.begin(), pmdists.begin() + top_k, 0.0) / top_k;
        }
        avg_pmdist_log.push_back(avg_pmdist);
        return 0;
    }

    IPSR_HANDLE_P inherit_handle_with_different_op(POINTS_NORMALS op,NormalEstimation<REAL,DIM>* estimator = new DoingNothing<REAL,DIM>()) {
        IPSR_HANDLE_P p = new ipsr_handle<REAL, DIM>(
            _cmd,
            _max_iters,
            _k_neighbors,
            update_plan,
            op
        );
        p->init(estimator);
        // 
        *(p->_p_.p_ep) = *(_p_.p_ep);
        p->_update_strategy = _update_strategy;
        return p;
    }

    void get_mesh(MESH& target_mesh) {
        lzd_tools::AcumutableTimer get_mesh_clock("ipsr_handle::get_mesh");
        if (avg_max_var_log.size() == 0) {
            /**************************reconstruction*******************************************/
            std::vector<char*> argv_str(_cmd.size());

            for (size_t i = 0; i < _cmd.size(); i++) {
                char* p = new char[_cmd[i].size() + 1];
                strcpy(p, _cmd[i].c_str());
                argv_str[i] = p;
                // argv_str[i] = &_cmd[i][0];
            }
            target_mesh = poisson_reconstruction_entrance<REAL, DIM>((int)argv_str.size(), argv_str.data(), this->_points_normals, &this->_weight_samples);
        }
        else {
            target_mesh = mesh;
        }
    }

    void update_mesh() {
        lzd_tools::AcumutableTimer update_mesh_clock("ipsr_handle::update_mesh");
        /**************************reconstruction*******************************************/
        std::vector<char*> argv_str(_cmd.size());
        
        for (size_t i = 0; i < _cmd.size(); i++){
            char* p = new char[_cmd[i].size() + 1];
            strcpy(p, _cmd[i].c_str());
            argv_str[i] = p;
            // argv_str[i] = &_cmd[i][0];
        }
        this->mesh = poisson_reconstruction_entrance<REAL, DIM>((int)argv_str.size(), argv_str.data(), this->_points_normals, &this->_weight_samples);    
        if (_cmd[find_arg(_cmd, "--bType") + 1] == "3") {
            std::vector<REAL> pmdists;
            int K = ConfigManager::get_common_config()["shrink_boundary_config"]["K"];
            int start_iter = ConfigManager::get_common_config()["shrink_boundary_config"]["start_iter"];
            if (_epoch >= start_iter && start_iter >= 0) {
                this->mesh = shrink_boundary<REAL, DIM>(this->mesh, this->_points_normals, pmdists, K);
                _out << "clean mesh\n";
            }
        }
    }

    //void update_mesh2() {
    //    
    //    this->spilt_hps.clear();
    //    this->_index_mapper.clear();

    //    _index_mapper.resize(max_label + 2);

    //        segs[idx].push_back(_points_normals[i]);
    //        _index_mapper[idx].push_back(i);
    //        this->spilt_hps.push_back(inherit_handle_with_different_op(segs[i]));
    //        _index_mapper.pop_back();
    //        this->spilt_hps.push_back(inherit_handle_with_different_op(segs.back()));
    //    

    //        this->spilt_hps[i]->update_mesh();
    //    this->mesh.first.clear();
    //    this->mesh.second.clear();
    //    vector<int> tri_label;
    //            this->mesh.first.push_back(spilt_hps[i]->ixform * seg_mesh.first[j]);
    //            this->mesh.second.push_back({ j_tri[0] + start_idx,j_tri[1] + start_idx,j_tri[2] + start_idx });
    //            tri_label.push_back(i);
    //        start_idx += seg_mesh.first.size();
    int iter() {
        if(avg_max_var_log.size() > 0&&_epoch<_max_iters){
            _out << "\n\n***************************************************";
            _out << "iter::epoch: " << _epoch;
            _out << " end***************************************************\n";
            std::string _out_str = _out.str();
            behavior_log.push_back(_out_str);
            _epoch++;
        }

        _out.flush();
        _out << "\n\n--------------------------------------------------";
        _out << "iter::epoch: " << _epoch;
        _out << " start--------------------------------------------------\n";

        /**************************update_plan************************************/
        if (_epoch < 0) {
            _out << "iter:: you should init the hander before iter\n";
        }

        if (_epoch >= _max_iters) {
            return -1;
        }
        if (update_plan(this->_p_) == -1) {
            return -1;
        }
        if (_cmd[find_arg(_cmd, "--bType") + 1] == "3") {
            _out << "iter::use Neumann  ";
        }
        else if (_cmd[find_arg(_cmd, "--bType") + 1] == "2") {
            _out << "iter::use Dirichlet\n";
        }
        else {
            _out << "iter::error boundary case\n";
        }
        // clear the mesh
        std::vector<Point<REAL, DIM>>().swap(this->mesh.first);
        std::vector<std::vector<int>>().swap(this->mesh.second);

        update_mesh();
        

        lzd_tools::AcumutableTimer update_normal_t("ipsr_handle::iter_update_normal");
        /**************************recompute the face normals********************************/
        // store the tri_face_normals of each point corresponding to the nearestSamples
        std::vector<Normal<REAL, DIM>> tri_face_normals(mesh.second.size());
        // nearestSamples[i][j] implies the index of the j-th nearest sample point of the i-th face center in mesh.second 
        nearestSamples.clear();
        nearestSamples.resize(mesh.second.size());
#pragma omp parallel for
        for (int i = 0; i < (int)mesh.second.size(); i++) {
            if (mesh.second[i].size() != 3) {
                _out << "Error: mesh.second[i].size() != 3" << endl;
            }
            Point<REAL, DIM> c = mesh.first[mesh.second[i][0]] + mesh.first[mesh.second[i][1]] + mesh.first[mesh.second[i][2]];
            c /= 3;
            array<REAL, 3> a{c[0], c[1], c[2]};
            tri_face_normals[i] = Point<REAL, DIM>::CrossProduct(mesh.first[mesh.second[i][1]] - mesh.first[mesh.second[i][0]], mesh.first[mesh.second[i][2]] - mesh.first[mesh.second[i][0]]);
            nearestSamples[i] = input_point_tree.knnSearch(kdt::KDTreePoint(a), _k_neighbors);
        }

        Normal<REAL, DIM> zero_normal(Point<REAL, DIM>(0, 0, 0));
        std::vector<Normal<REAL, DIM>> projective_normals(_points_normals.size(), zero_normal);


        std::vector<std::vector<int>> sample_tri_idx(_points_normals.size());

        /****************************************
        * indirect mapping sampling
        ******************************************/
        if (_update_strategy == INDIRECT_MAPPING) {
            sample_tri_idx = inverse_map<REAL, DIM>(nearestSamples, _points_normals.size());
            _out << "iter::indirectly mapping..." << endl;
        }
        /****************************************
         * direct mapping sampling
         * ******************************************/
        else if (_update_strategy == DIRECT_MAPPING) {
            _out << "iter::directly mapping..." << endl;
            // rebuild the kdtree of the tri_face_normals in every iter. that's the main reason why direct mapping is slower than indirect mapping
            kdt::KDTree<kdt::KDTreePoint> tri_center_tree;
            std::vector<kdt::KDTreePoint> tri_centers;
            tri_centers.reserve(mesh.second.size());
#pragma omp parallel for
            for (int i = 0; i < mesh.second.size(); i++) {
                Point<REAL, DIM> c = mesh.first[mesh.second[i][0]] + mesh.first[mesh.second[i][1]] + mesh.first[mesh.second[i][2]];
                c /= 3;
                array<REAL, 3> a{c[0], c[1], c[2]};
                tri_centers.push_back(kdt::KDTreePoint(a));
            }
            tri_center_tree.build(tri_centers);

#pragma omp parallel for
            for (int i = 0; i < (int)_points_normals.size(); i++) {
                array<REAL, 3> a{_points_normals[i].first[0], _points_normals[i].first[1], _points_normals[i].first[2]};
                auto idxs = tri_center_tree.knnSearch(a, _k_neighbors);
                sample_tri_idx[i].resize(idxs.size());
                for (int j = 0; j < idxs.size(); j++) {
                    sample_tri_idx[i][j] = idxs[j];
                }
            }
        }
        else {
            _out << "iter::Error: unknown update strategy" << endl;
            return -1;
        }


#pragma omp parallel for
        for (int i = 0; i < _points_normals.size(); i++) {
            for (int j = 0; j < sample_tri_idx[i].size(); j++) {
                projective_normals[i] += tri_face_normals[sample_tri_idx[i][j]];
            }
            normalize<REAL, DIM>(projective_normals[i]);
        }
        update_normal_t.Done();

        lzd_tools::AcumutableTimer log("ipsr_handle::iter::cal_metric and log");
        // update the normals of points and compute the average difference between the normals of points before and after the update
        old_avg_max_diff = avg_max_diff;
        avg_max_diff = 0;
        size_t heap_size = static_cast<size_t>(ceil(_points_normals.size() / 1000.0));
        size_t tmp = std::min(_points_normals.size(), static_cast <size_t>(10));
        heap_size = std::max(heap_size, tmp);
        assert(heap_size != 0);
        avg_max_diff = cal_points_normals_diff(projective_normals, _points_normals, _points_normals_diff, heap_size);
        
        std::vector<int> if_update(_points_normals.size(), 1);
        int not_update_count = 0;        
        for(auto filter:_filters){
            filter->And(if_update,projective_normals);
        }
        for (int i = 0; i < _points_normals.size(); i++) {
            if (if_update[i] == 0) {
                not_update_count++;
            }
            else {
                _points_normals[i].second = projective_normals[i];
            }
        }
        _out<<("not_update_count: %d\n", not_update_count);

        avg_max_var_log.push_back(avg_max_diff);
        _out << "iter::avg_max_diff: " << avg_max_diff << endl;

        if (!ConfigManager::get_common_config()["no_save"]) {
            cal_i_avg_var(sample_tri_idx, tri_face_normals);
            std::vector<std::pair<int, REAL>> var_idx;
            for (int i = 0; i < i_nearestSample_var.size(); i++) {
                var_idx.push_back(std::make_pair(i, i_nearestSample_var[i]));
            }
            std::sort(var_idx.begin(), var_idx.end(), [](std::pair<int, REAL> a, std::pair<int, REAL> b) {return a.second > b.second; });
            REAL avg_var = 0;
            for (size_t i = 0; i < heap_size; i++) {
                avg_var += var_idx[i].second;
            }
            avg_var /= heap_size;
            sample_var_log.push_back(avg_var);
            _out << "iter::avg_max_samplevar: " << avg_var << endl;
            cal_and_log_matric(heap_size);
            _out << "iter::avg_max_pmdist: " << avg_pmdist_log[avg_pmdist_log.size() - 1] << endl;
            _out << "iter::mesh_size: " << mesh_size_log[mesh_size_log.size() - 1] << endl;
            _out << "iter::avg_loss: " << loss_log[loss_log.size() - 1] << endl;
        }
        

        /**************************************************************
         *  last iteration
         **************************************************************/
        if (_epoch == _max_iters-1 || avg_max_diff < 0.1) {
            if (ConfigManager::get_common_config()["PCA_refine"]) {
				_out << "iter::use PCA to refine the normals\n";
                NormalEstimation<REAL, DIM>* estimator = new PCAEstimation<REAL, DIM>();
                auto pca_normal = POINTS_NORMALS(_points_normals);
                estimator->Estimate(pca_normal);
#pragma omp parallel for
                for (int i = 0; i < _points_normals.size(); i++) {
                    if (dot<REAL,DIM>(_points_normals[i].second, pca_normal[i].second) < 0) {
                        pca_normal[i].second *= -1;
                    }
                    _points_normals[i].second = pca_normal[i].second;
                }
                update_mesh();
            }
            int zero_c = cal_i_avg_var(sample_tri_idx, tri_face_normals);
            _epoch = _max_iters-1;
            _out << "iter::" << zero_c << " points haven't been selected by any face(" << _points_normals.size() << " totally)" << endl;
        }  
        return _epoch;
    }
    void end_iter() {
        _out << "iter " << _epoch << " is end_iter...\n";
        // map the _points_normals to the _res_points_normals
        _res_points_normals.resize(_ori_points_normals.size());
        for (int i = 0; i < _ori_points_normals.size(); ++i)
        {
            auto p = _ori_points_normals[i].first;
            _res_points_normals[i].first = p;
            p = ixform * p;
            auto nearest_op = input_point_tree.knnSearch(kdt::KDTreePoint({ p[0], p[1], p[2] }), 1);
            _res_points_normals[i].second = _points_normals[nearest_op[0]].second;
        }
        _out << "end_iter done\n";
    }

    int cal_i_avg_var(const MULTI_MAP mapper, const std::vector<Normal<REAL, DIM>>& sample_normals) {
        lzd_tools::AcumutableTimer clock("ipsr_handle::cal_i_avg_var");
        Normal<REAL, DIM> zero_normal(Point<REAL, DIM>(0, 0, 0));
        int zero_c = 0;
        i_nearestSample_var.resize(_points_normals.size());
        avg_normal.resize(_points_normals.size());
        // calculate the average normal of each point
        for (size_t i = 0; i < _points_normals.size(); i++) {
            avg_normal[i] = zero_normal;
            if (mapper[i].size() == 0) {
                zero_c++;
                continue;
            }

            for (size_t j = 0; j < mapper[i].size(); j++) {
                avg_normal[i] += sample_normals[mapper[i][j]];
            }
            avg_normal[i] /= mapper[i].size();
        }

        // calculate the variance of the tri_face_normals
        for (size_t i = 0; i < _points_normals.size(); ++i)
        {
            REAL variance = 0;
            for (size_t j = 0; j < mapper[i].size(); ++j)
            {
                // calculate the angle between the normal and the average normal as the variance
                REAL diff = abs(calculate_angle(Normal<REAL, DIM>(sample_normals[mapper[i][j]]), avg_normal[i]));
                variance += diff;
            }
            if (mapper[i].size() != 0) {
                variance /= mapper[i].size();
                i_nearestSample_var[i] = variance;
            }
            else {
                i_nearestSample_var[i] = 0;
            }
        }
        return zero_c;
    }

    // @return
    float cal_points_normals_diff(const std::vector<Normal<REAL, DIM>>& new_normals, const POINTS_NORMALS& old_op, std::vector<float>& points_normals_diff, int K) {
        assert(new_normals.size() == old_op.size());
        points_normals_diff.resize(new_normals.size());
        Normal<REAL, DIM> zero_normal(Point<REAL, DIM>(0, 0, 0));
#pragma omp parallel for
        for (int i = 0; i < new_normals.size(); i++) {
            if (new_normals[i] == zero_normal) {
                points_normals_diff[i] = 0;
                continue;
            }
            else {
                points_normals_diff[i] = Point<REAL, DIM>::SquareNorm((new_normals[i] - old_op[i].second).normal);
            }
        }
        std::vector<float> temp(points_normals_diff);
        std::sort(temp.begin(), temp.end(), std::greater<REAL>());
        float avg = 0;
        assert(K <= temp.size());
        for (int i = 0; i < K; i++) {
            avg += sqrt(temp[i]);
        }
        avg /= float(K);
        return avg;
    }

    int save_current_ops(std::string addition_folder = "/temp/", std::string addition_name = "",
        XForm<REAL, DIM + 1> transform = XForm<REAL, DIM + 1>().Identity(),
        std::pair<std::vector<int>, lzd_tools::FiniteMap> flagandColorMap = std::make_pair(std::vector<int>(), lzd_tools::FiniteMap())) 
    
    {
            std::string path = get_out_put_base(addition_folder);
            std::string name = addition_name + "_op_it" + std::to_string(*(this->_p_.p_ep)) + ".ply";
            // output_sample_points_and_normals<REAL, DIM>(path + name, _points_normals, transform * ixform, color_flag);
            lzd_tools::op2ply(_points_normals, path + name, transform * ixform, flagandColorMap);
            if (avg_max_var_log.size() == 0 || sample_var_log.size() == 0)return 0;
            _save_log(path+name);
            return 0; 
    }

    int save_current_mesh(std::string addition_folder = "/temp/", std::string addition_name = "",
        XForm<REAL, DIM + 1> transform = XForm<REAL, DIM + 1>().Identity()
     )
    {   
        std::string path = get_out_put_base(addition_folder);
        std::string name = addition_name + "_mesh_it" + std::to_string(*(this->_p_.p_ep)) + ".ply";
        if (this->mesh.first.size() == 0) {
            printf("WRONG! NO MESH TO SAVE!");
            return -1;
        }
        output_ply(path + name, this->mesh, transform * ixform);
        return 0;

    }

    double cal_save_mesh_consistency(std::string adname = "", VertexGetterP func = get_tetrahedron_vertex) {
        std::string path = _cmd[find_arg(this->_cmd, "--out") + 1] + "/consistency/";
        rmkdir(path);
        std::string name = "mesh_op_it" + std::to_string(*(this->_p_.p_ep));
        std::vector<REAL> pmdist;
        std::vector<char* > argv_str(_cmd.size());
        for (int i = 0; i < argv_str.size(); i++) {
            char* p = new char[_cmd[i].size() + 1];
            strcpy(p, _cmd[i].c_str());
            argv_str[i] = p;
            // argv_str[i] = &_cmd[i][0];
        }
        // this->mesh = poisson_reconstruction_entrance<REAL, DIM>((int)argv_str.size(), argv_str.data(), this->_points_normals, &this->_weight_samples);
        auto _mesh = poisson_reconstruction_entrance<REAL, DIM>((int)argv_str.size(), argv_str.data(), _points_normals, &_weight_samples);



        auto cleaned_mesh = shrink_boundary(_mesh,this->_points_normals,pmdist);
        
        std::string pname = path + name;
        if (adname == "")pname = adname;
        else pname += adname;
        MaxAB<REAL,DIM> oc(func);
        std::vector<double> temp;
        double consistency = oc.cal_consistency(cleaned_mesh,temp,pname);
        
        printf("consistency = %f\n", consistency);

        return consistency;
    }

    std::string get_out_put_base(std::string addition_folder = "") {
        std::string res = _cmd[find_arg(this->_cmd, "--out") + 1] +  addition_folder;
        rmkdir(res);
        return res;
    }

    std::string get_resname() {
        std::string depth = _cmd[find_arg(this->_cmd, "--depth") + 1];
        std::string point_weight = _cmd[find_arg(this->_cmd, "--pointWeight") + 1];
        std::string filename = "it_" + to_string(_max_iters);
        filename += "_dp_" + depth;
        filename += "_nb_" + to_string(_k_neighbors);
        filename += "_sd_" + to_string(_seed);
        filename += "_pt_" + point_weight;
        return filename;
    }

    int save_res(std::string adname="") {
        std::string output_path = _cmd[find_arg(this->_cmd, "--out") + 1];
        std::string depth = _cmd[find_arg(this->_cmd, "--depth") + 1];

        rmkdir(output_path);
        std::string res_name = output_path + get_resname() + adname;
        
        std::vector<REAL> pmdist;
        if(this->mesh.first.size() == 0) {
            get_mesh(this->mesh);
        }
        nlohmann::json j = ConfigManager::get_common_config();
        if(j["save_option"]["uncleaned_mesh"] == true) {         
            output_ply(res_name + "_surface.ply",this->mesh,ixform);
        }
        if(j["save_option"]["cleaned_mesh"] == true) {
            MESH cmesh = clean_mesh(this->mesh, _points_normals, pmdist, 20);
            output_ply(res_name + "_cleaned_surface.ply",cmesh,ixform);
        }
        if(j["save_option"]["orientedpoint"] == true) {
            lzd_tools::op2ply(_points_normals, res_name + "_orientedpoint.ply", ixform);
        }
        if(j["save_option"]["GT_samples"] == true) {
            // gt
            lzd_tools::op2ply(_gt_points_normals, res_name + "_GT_samples.ply", ixform);
        }
        if(j["save_option"]["ori_est"] == true) {
            POINTS_NORMALS ori_est = POINTS_NORMALS(_ori_points_normals);
            auto unused_mapper = std::vector<int>();
            printf("update to origin points\n"); 
            XForm <REAL, DIM + 1> xform = ixform.inverse();
            auto dist = update_target(ori_est, xform, unused_mapper);
            printf("update to origin points done, point pair total dist = %f\n", dist);
            lzd_tools::op2ply(ori_est, res_name + "_ori_est.ply", XForm<REAL, DIM + 1>().Identity());
        }

        //string var_path = res_name + "_variance.txt";
        //ofstream var_file(var_path);
        //    var_file << i_nearestSample_var[i] << endl;
        //var_file.close();
        //
        //string component_path = res_name + "_component.txt";
        //ofstream component_file(component_path);
        //        component_file << component_index[i][j] << " ";
        //    component_file << endl;
        //component_file.close();
        //
        //string diff_path = res_name + "_points_normals_diff.txt";
        //ofstream diff_file(diff_path);
        //    diff_file << _points_normals_diff[i] << endl;
        //diff_file.close();
        //
        //string pmdist_path = res_name + "_pmdist.txt";
        //ofstream pmdist_file(pmdist_path);
        //    pmdist_file << pmdist[i] << endl;
        //pmdist_file.close();

        std::string log_path = output_path + "iter_logs/";
        rmkdir(log_path);
        _save_log(log_path + get_resname() + "_");
        return 0;
    }

    bool if_end() {
        return _epoch == _max_iters;
    }    

    POINTS_NORMALS get_res(){
        POINTS_NORMALS res = POINTS_NORMALS(_ori_points_normals);
        auto trasform = ixform.inverse();
        for(int i = 0; i < res.size(); i++) {
            auto c = trasform * res[i].first;
            std::array<REAL, 3> a{c[0], c[1], c[2]};
            int n = input_point_tree.nnSearch(kdt::KDTreePoint(a));
            res[i].second = _points_normals[n].second;
        }   
        return res;
    }

    /**
     * @brief 
     */
    REAL update_target(POINTS_NORMALS& target, XForm<REAL, DIM + 1> transform, std::vector<int>& mapper) {
        REAL dist = 0;
        mapper.resize(target.size());
        if(_points_normals.size()<1) {
			printf("update_target::WARNING _points_normals.size()<1\n");
			return 1e9;
		}

#pragma omp parallel for
        for(int i = 0; i < target.size(); i++) {
            auto c = transform * target[i].first;
            std::array<REAL, 3> a{c[0], c[1], c[2]};
            int n = input_point_tree.nnSearch(kdt::KDTreePoint(a));
            target[i].second = _points_normals[n].second;
            mapper[i] = n;
        }   
        return dist;
    }

private:
    int _save_log(std::string path) {
        std::ofstream out(path + "avg_max_var_log.txt");
        for (int i = 0; i < avg_max_var_log.size(); i++) {
            out << avg_max_var_log[i] << std::endl;
        }
        out.close();
        std::ofstream out2(path + "sample_var_log.txt");
        for (int i = 0; i < sample_var_log.size(); i++) {
            out2 << sample_var_log[i] << std::endl;
        }
        out2.close();

        std::ofstream out3(path + "loss_log.txt");
        for (int i = 0; i < loss_log.size(); i++) {
            out3 << loss_log[i] << std::endl;
        }
        out3.close();

        std::ofstream out4(path + "avg_pmdist_log.txt");
        for (int i = 0; i < avg_pmdist_log.size(); i++) {
            out4 << avg_pmdist_log[i] << std::endl;
        }
        out4.close();

        std::ofstream out5(path + "mesh_size_log.txt");
        for (int i = 0; i < mesh_size_log.size(); i++) {
            out5 << mesh_size_log[i] << std::endl;
        }
        out5.close();

        std::ofstream out6(path + "behavior_log.txt");
        for (int i = 0; i < behavior_log.size(); i++) {
            out6 << behavior_log[i] << std::endl;
        }
        out6.close();
        return 0;
    }
};

template<typename REAL, unsigned int DIM>
void _select_tri(const MESH& mesh, const POINTS_NORMALS& op,
                std::vector<REAL>& pmdist,std::vector<int>& if_select,
                int k = 20)
{
    if (mesh.first.size() == 0) {
        printf("select_mesh::WARNING mesh.first.size()==%zu, op.size()==%zu", mesh.first.size(), op.size());
        return;
    }
    if_select = std::vector<int>(mesh.second.size(), 0);
    std::vector<Point<REAL, DIM>> center_points(mesh.second.size(), Point<REAL, DIM>({0,0,0}));
    std::vector<kdt::KDTreePoint> vertices(mesh.second.size(), kdt::KDTreePoint({0,0,0}));
    kdt::KDTree<kdt::KDTreePoint> kdtree;
    for (int i = 0; i < mesh.second.size(); ++i)
    {
        auto v = (mesh.first[mesh.second[i][0]] + mesh.first[mesh.second[i][1]] + mesh.first[mesh.second[i][2]]) / 3;
        center_points[i] = v;
        array<REAL, 3> _p_{v[0], v[1], v[2]};
        vertices[i] = (kdt::KDTreePoint(_p_));
    }
    kdtree.build(vertices);
    pmdist.resize(op.size());
    for(int i=0;i<op.size();i++){
        Point<REAL, DIM> p = op[i].first;
        auto idxs = kdtree.knnSearch(kdt::KDTreePoint({p[0], p[1], p[2]}), k);
        pmdist[i] = distance(p, center_points[idxs[0]]);
        for(int j=0;j<idxs.size();j++){
            if_select[idxs[j]] = 1;
        }
    }
}


/**
 * @brief 
 */
template<typename REAL, unsigned int DIM>
MESH clean_mesh(const MESH& mesh, const POINTS_NORMALS& op,std::vector<REAL>& pmdist,int k = 20){
    if (mesh.first.size() == 0) {
        printf("clean_mesh::WARNING mesh.first.size()==%zu, op.size()==%zu", mesh.first.size(), op.size());
        return mesh;
    }
    std::vector<int> if_select(mesh.second.size(), 0);
    _select_tri(mesh,op,pmdist,if_select,k);    

    MESH cleaned_mesh;
    cleaned_mesh.first = std::vector<Point<REAL, DIM>>(mesh.first);
    for(int i=0;i<mesh.second.size();i++){
        if(if_select[i] == 1){
            cleaned_mesh.second.push_back(mesh.second[i]);
        }
    }
    return cleaned_mesh;
}

/**
 * @brief 
 */
template<typename REAL, unsigned int DIM>
MESH shrink_boundary(const MESH& mesh, const POINTS_NORMALS& op,std::vector<REAL>& pmdist,int k){
    unsigned int n = mesh.first.size();
    unsigned int m = mesh.second.size();
    std::vector<int> t;
    std::vector<bool> is_near_tri(m, false);
    std::vector<bool> is_near_point(n, false);
    std::unordered_map<unsigned int, unsigned int> edge_cnt;
    std::vector<bool> is_boundary(n, true);
    std::vector<std::vector<unsigned int>> edge_list(n, std::vector<unsigned int>());
    
    _select_tri(mesh,op,pmdist, t,k);
    for (int i = 0; i < m; i++)is_near_tri[i] = t[i] > 0;

    unsigned int tn = n;
    unsigned int cnt =0;
    while(tn>>=1)cnt++;
    auto edge_hash = [&cnt](const unsigned int& a, const unsigned int& b)->unsigned int {
        return a < b ? (a << cnt) + b : (b << cnt) + a;
    };
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < 3; j++) {
            int a = mesh.second[i][j], b = mesh.second[i][(j + 1) % 3];
            is_near_point[a] = is_near_point[a] || is_near_tri[i];
            is_near_point[b] = is_near_point[b] || is_near_tri[i];
        }
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < 3; j++) {
            int a = mesh.second[i][j], b = mesh.second[i][(j + 1) % 3];
            auto ehs = edge_hash(a,b);
            if (edge_cnt.find(ehs) != edge_cnt.end()) {
                is_boundary[a] = is_boundary[b] = false;
            }
            else {
                edge_cnt.insert(std::make_pair(ehs, 1));
            }
        }
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < 3; j++) {
            int a = mesh.second[i][j], b = mesh.second[i][(j + 1) % 3];
            assert(a < n && b < n);
            if (is_near_point[a] || is_near_point[b])continue;
            edge_list[a].push_back(b), edge_list[b].push_back(a);
        }
    }

    std::queue<unsigned int> q;
    for (int i = 0; i < n; i++) {
        if (!is_near_point[i] && is_boundary[i]) {
            q.push(i);
        }
    }
    while (!q.empty()) {
        auto font = q.front();
        q.pop();
        for (auto nei : edge_list[font]) {
            if (!is_near_point[nei] && !is_boundary[nei]) {
                is_boundary[nei] = true;
                q.push(nei);
            }
        }
    }

    std::vector<bool> is_select_point(is_near_point);
    std::vector<bool> is_select_tri(is_near_tri);
    for (int i = 0; i < n; i++)is_select_point[i] = is_select_point[i] || !is_boundary[i];

#pragma omp parallel for
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < 3; j++)is_select_tri[i] = is_select_tri[i] || is_select_point[mesh.second[i][j]];
    }
    MESH cleaned_mesh;
    cleaned_mesh.first = std::vector<Point<REAL, DIM>>(mesh.first);
    for (int i = 0; i < mesh.second.size(); i++) {
        if (is_select_tri[i])cleaned_mesh.second.push_back(mesh.second[i]);
    }
    return cleaned_mesh;
}


template<typename REAL, unsigned int DIM>
int fix_update(IPSR_HANDLE_P p_ipsr, int (*fix_plan)(Period& p) = diri_dir_ipsr, int fix_times = 3) {
    
    printf("fix_update...\n");
    p_ipsr->_max_iters += fix_times;
    auto tp = p_ipsr->update_plan;
    p_ipsr->update_plan = fix_plan;
    for (int i = 0; i < fix_times; i++) {
        p_ipsr->iter();
    }
    p_ipsr->update_plan = tp;
    p_ipsr->_max_iters -= 3;
    printf("fix_update done\n");
    return 0;
}
