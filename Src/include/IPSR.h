#pragma once 
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include "pcl_api.h"
#include "o3d_api.h"
#include "kdtree.h"
#include "utility.h"
#include "tools.h"
//#include <mincut_api.h>
#include "orientation_consistency.h"
#include "NormalEstimation.h"
#include "ipsr_op_filter.h"

#define INDIRECT_MAPPING 0
#define DIRECT_MAPPING 1

#define IPSR_HANDLE_P ipsr_handle<REAL,DIM>* 

// 用于管理当前epoch的行为(参数、更新策略、是否打印)
struct Period{
    // 当前epoch
    int* p_ep;
    // 上一个ep的var
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

// 交互式更新
int interact_update(Period& p);
// dirichlet 边界条件 + 间接映射
int diri_ind_ipsr(Period& p);
// dirichlet 边界条件 + 直接映射
int diri_dir_ipsr(Period& p);
// neumman 边界条件 + 间接映射
int neumman_ind_ipsr(Period& p);
// neumman 边界条件 + 直接映射
int neumman_dir_ipsr(Period& p);

int comb_plan5(Period& p);

int comb_update6(Period& p);

int comb_update7(Period& p);

int comb_update8(Period& p);

int specify_update(Period& p);

typedef int (*UpdatePlan)(Period&);

UpdatePlan get_update_plan_by_name(std::string name);

//template<typename REAL, int DIM>
//int pcl_fix(POINTS_NORMALS& points_normals, POINTS_NORMALS& pcl_normals){
//    std::vector<int> flag(pcl_normals.size(),0);
//    if(points_normals.size() != pcl_normals.size()){
//        printf("pcl_fix::points_normals.size() != pcl_normals.size()\n");
//        return -1;
//    }
//    int count_fix = 0;
//    for(int i=0;i<points_normals.size();i++){
//        // 计算两个法向量的夹角，并转为角度制
//        REAL theta = calculate_angle(points_normals[i].second, pcl_normals[i].second);
//        theta = 180 * (theta/M_PI); 
//        REAL bias = 20;
//        // 如果夹角在80-100度之间,则points_normals[i]的法向量替换为pcl_normals[i]的法向量
//        if(theta>90-bias&&theta<90+bias){
//            points_normals[i].second = get_rand_norm<REAL, DIM>(0);
//            count_fix++;
//            flag[i] = 1;
//        }
//    }
//    return count_fix;
//}

lzd_tools::FiniteMap get_mincut_colormap();

/****************************************
 * 面向points_normals的ipsr
 * ******************************************/
template<typename REAL, int DIM>
class ipsr_handle {
    nlohmann::json _config_j;
    // nlohmann::json _log_j;

    std::vector<Update_filter<REAL, DIM>*> _filters;

public:
    //需要初始化时指定
    std::vector<std::string> _cmd;
    int _max_iters;
    int _k_neighbors;
    int _seed;
    int (*update_plan)(Period& p);
    POINTS_NORMALS _ori_points_normals;// 原始点云,不会归一化或者下采样
    // 运行期间更新
public:
    POINTS_NORMALS _points_normals;//使用对_ori_points_normals进行采样并移动后得到的点。psr会在采样阶段建立octree树，进行采样并生成权重，这个步骤没法跳过
    POINTS_NORMALS _gt_points_normals;//采样后的groud_truth，即采样后的_ori_points_normals，同样经过的采样与偏移
    POINTS_NORMALS _res_points_normals;//最终结果 (弃用)
    POINTS_NORMALS _init_normals; // 初始化的法向量
    int _epoch; // 当前迭代次数
    Period _p_; // 保存当前的一些状态
    int _update_strategy;// 更新策略 
    REAL avg_max_diff;// 衡量本次迭代的角度变化量  
    REAL old_avg_max_diff;
    XForm<REAL, DIM + 1> ixform; //_points_normals * ixform = _ori_points_normals在原空间中的采样
    MESH mesh;
    std::vector<Normal<REAL, DIM>> avg_normal;// 每个点与所有选择了该点的三角面片的法向量的均值（需要用来计算方差）
    std::vector<float> i_nearestSample_var;// 每个点与所有选择了该点的三角面片的法向量的方差
    std::vector<float> _points_normals_diff; // 每个点的法向量的更新幅度。
    std::vector<double> _weight_samples;// psr采样阶段生成的
    kdt::KDTree<kdt::KDTreePoint> input_point_tree;//由迭代期间的点建立的kdtree
    std::vector<std::vector<int>> nearestSamples;// 三角面片的k个最近点在_point_normals索引

    // *******************************以下部分在拷贝时不需要复制***************************
    std::ostringstream _out;// 输出流
    // 日志
    std::vector <REAL> avg_max_var_log;// 点变换量的topk均值
    std::vector <REAL> sample_var_log;// 采样方差的topk均值
    std::vector <REAL> loss_log;// 每次迭代后的loss
    std::vector <REAL> avg_pmdist_log;// 点到mesh距离的topk均值
    std::vector <int> mesh_size_log;// 每次迭代后mesh的大小 (三角形个数)
    std::vector <std::string> behavior_log;// 迭代期间的行为(TODO)
    // std::vector <std::pair<int,int>> init_log;// 记录在哪个epoch使用了哪种初始化方式

    nlohmann::json log_j;
    nlohmann::json mj;

    // 供spile_iter使用
    std::vector<IPSR_HANDLE_P> spilt_hps;
    std::vector<std::vector<int>> _index_mapper;//spilte_res的点在原点云中的索引
    bool _if_spilt = false;

    // 其他的一些数据
    std::vector<int> dlabel;// 每个点所属的簇

public:
    // 构造函数
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

    // 析构函数
    // ~ipsr_handle<REAL,DIM>(){
    //     //printf("ipsr_handle::~ipsr_handle()\n");
    //     // TODO
    // }

    // 拷贝构造函数 将ipsr_handle<REAL,DIM> ipsr的值拷贝到当前对象
    // 更新阶段也维持原样。
    ipsr_handle<REAL, DIM>(const ipsr_handle<REAL, DIM>& ipsr) {
        this->_cmd = std::vector<std::string>(ipsr._cmd.size());
        for (int i = 0; i < ipsr._cmd.size(); i++) {
            this->_cmd[i] = ipsr._cmd[i];
        }
        // 参数以及_point_normals不变
        this->_max_iters = ipsr._max_iters;
        this->update_plan = ipsr.update_plan;
        this->_ori_points_normals = ipsr._ori_points_normals;
        this->_points_normals = ipsr._points_normals;
        this->_res_points_normals = ipsr._res_points_normals;
        this->_k_neighbors = ipsr._k_neighbors;
        this->_seed = ipsr._seed;

        this->_p_ = Period(); // 重新赋值 _p_
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


        vector<kdt::KDTreePoint> vertices;
        vertices.reserve(this->_points_normals.size());
        for (size_t i = 0; i < this->_points_normals.size(); ++i)
        {
            array<REAL, 3> _p_{this->_points_normals[i].first[0], this->_points_normals[i].first[1], this->_points_normals[i].first[2]};
            vertices.push_back(kdt::KDTreePoint(_p_));
        }
        input_point_tree.build(vertices);
    }

    // sample_points、init_normal、build_kdtree
    // init_type: 0:nothing 1:random 2:pca 3:pcl
    // gt_points_normals: 从_point_normals中采样得到的点云


    /**
     * @brief 
     * step1: octree下采样,并拷贝一份到_gt_points_normals
     * step2: 建立_points_normals的kdtree
     * step3: 初始化_out;初始化纠正器
     * @return int 
     */
    int ipsr_initialization(){
        _epoch = 0;
        // 采样点
        std::vector<char* > argv_str(this->_cmd.size());
        for (int i = 0; i < this->_cmd.size(); i++) {
            argv_str[i] = &this->_cmd[i][0];
        }

        this->_points_normals = sample_points_entrance<REAL, DIM>((int)argv_str.size(), argv_str.data(), _ori_points_normals, this->ixform, &this->_weight_samples);
        if (_points_normals.size() == 0) {
            printf("WARNING::ipsr_handle::init::sample_points.size() == 0\n");
            // 使用原始点云作为采样点
            // TODO 可能会有问题
            this->ixform = XForm<REAL, DIM + 1>().Identity();
            this->_points_normals = _ori_points_normals;
            _weight_samples = std::vector<double>(_points_normals.size(), 1.0);
        }

        assert(ixform.determinant() != 0);
        assert(ixform.determinant() != NAN);
        // 将_points_normals值拷贝到_gt_points_normals
        _gt_points_normals.reserve(_points_normals.size());
        for (size_t i = 0; i < _points_normals.size(); i++) _gt_points_normals.push_back(_points_normals[i]);
        // 建立points_normals的kdtree
        std::vector<kdt::KDTreePoint> vertices;
        vertices.reserve(this->_points_normals.size());
        for (size_t i = 0; i < this->_points_normals.size(); ++i)
        {
            array<REAL, 3> _p_{ this->_points_normals[i].first[0], this->_points_normals[i].first[1], this->_points_normals[i].first[2] };
            vertices.push_back(kdt::KDTreePoint(_p_));
        }
        input_point_tree.build(vertices);
        //初始化_out为一个字符串流
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
     * 建议使用ipsr_initialization() + init_op_normal()
     * */
    int init(NormalEstimation<REAL, DIM>* estimator)
    {
        ipsr_initialization();
        // 初始化点云法向量
        init_op_normal(estimator);
        return 0;
    }

    void init_op_normal(NormalEstimation<REAL,DIM>* estimator){
        estimator->Estimate(_points_normals);
        // std::pair<int, std::string> p = std::make_pair(_epoch, init_type_namelist[init_type]);
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

    // 得到当前_points_normals的metric
    nlohmann::json get_metric(){
        nlohmann::json j;
        // if(loss_log.size() == 0){
        //     return j;
        // }
        j = lzd_tools::PointNormalMetric<REAL,DIM>(_points_normals, _gt_points_normals).to_json();
        return j;
    }
    // 计算当前_points_normals的loss(假定真值存在)
    REAL cal_avg_loss(bool permid_nd = true) {
        REAL avg_loss = 0;
        for (size_t i = 0; i < _points_normals.size(); i++) {
            avg_loss += calculate_angle(_gt_points_normals[i].second, _points_normals[i].second) * 180 / M_PI;
        }
        avg_loss /= _points_normals.size();
        if (permid_nd)avg_loss = std::min(avg_loss, 180 - avg_loss);
        return avg_loss;
    }

    // 计算一些迭代期间用不到的matric
    // 可指定pmdist取topmaxk的均值
    int cal_and_log_matric(int top_k = -1) {
        lzd_tools::AcumutableTimer cal_matric("ipsr_handle::cal_and_log_matric");

        if (top_k <= 0)top_k = _points_normals.size();
        // 计算当前epoch的avg_loss = avg(gt_points_normals - points_normals))
        REAL avg_loss = cal_avg_loss();
        loss_log.push_back(avg_loss);
        // printf("avg_loss: %f\n", avg_loss); //不应该在这里打印

        // 记录当前的mesh_size
        mesh_size_log.push_back(mesh.second.size());

        // 计算当前epoch的avg_distance
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
        // printf("avg_max_pmdist: %f\n", avg_pmdist);
        return 0;
    }

    // 得到一个除了点云数据以外，其他都和当前对象相同的ipsr_handle对象
    // init_type seed为init阶段参数
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

    // 得到一个曲面
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
        // 如果是Neumann边界条件 需要清理mesh
        if (_cmd[find_arg(_cmd, "--bType") + 1] == "3") {
            std::vector<REAL> pmdists;
            int K = ConfigManager::get_common_config()["shrink_boundary_config"]["K"];
            int start_iter = ConfigManager::get_common_config()["shrink_boundary_config"]["start_iter"];// start_iter<0表示永不清理
            if (_epoch >= start_iter && start_iter >= 0) {
                this->mesh = shrink_boundary<REAL, DIM>(this->mesh, this->_points_normals, pmdists, K);
                _out << "clean mesh\n";
            }
        }
    }

    // 先对目标进行聚类 然后再分别更新
    //void update_mesh2() {
    //    double _eps = 0.02;
    //    int _min_points = 10;
    //    
    //    this->spilt_hps.clear();
    //    this->_index_mapper.clear();
    //    std::vector<POINTS_NORMALS> segs;
    //    std::vector<int> label;
    //    o3d_dbscan_extraction(_points_normals, label, _eps, _min_points);
    //    int max_label = *max_element(label.begin(), label.end());

    //    segs.resize(max_label + 2);// segs[max_label+1]中存放other
    //    _index_mapper.resize(max_label + 2);

    //    for (int i = 0; i < _points_normals.size(); i++) {
    //        int idx = label[i] >= 0 ? label[i] : (max_label + 1);
    //        segs[idx].push_back(_points_normals[i]);
    //        _index_mapper[idx].push_back(i);
    //    }
    //    for (int i = 0; i < segs.size() - 1; i++) {
    //        assert(segs[i].size() > 0);
    //        this->spilt_hps.push_back(inherit_handle_with_different_op(segs[i]));
    //    }
    //    if (segs.back().size() == 0) {
    //        _index_mapper.pop_back();
    //    }
    //    else {
    //        this->spilt_hps.push_back(inherit_handle_with_different_op(segs.back()));
    //    }
    //    _if_spilt = true;
    //    

    //    for (int i = 0; i < this->spilt_hps.size(); i++) {
    //        this->spilt_hps[i]->update_mesh();
    //    }
    //    int start_idx = 0;
    //    this->mesh.first.clear();
    //    this->mesh.second.clear();
    //    vector<int> tri_label;
    //    for (int i = 0; i < this->spilt_hps.size(); i++) {
    //        MESH seg_mesh = this->spilt_hps[i]->mesh;
    //        for (int j = 0; j < seg_mesh.first.size(); j++) {
    //            this->mesh.first.push_back(spilt_hps[i]->ixform * seg_mesh.first[j]);
    //        }
    //        for (int j = 0; j < seg_mesh.second.size(); j++) {
    //            std::vector<int> j_tri = seg_mesh.second[j];
    //            this->mesh.second.push_back({ j_tri[0] + start_idx,j_tri[1] + start_idx,j_tri[2] + start_idx });
    //            tri_label.push_back(i);
    //        }
    //        start_idx += seg_mesh.first.size();
    //    }
    //    _if_spilt = false;
    //    std::string name = std::to_string((int)(this));
    //    lzd_tools::mesh2ply(mesh, get_out_put_base("/update2/") + name + "mesh.ply", this->ixform, std::make_pair(tri_label, lzd_tools::get_regular_colormap(0, max_label+2)));
    //    lzd_tools::op2ply(_points_normals,get_out_put_base("/update2/")+ name + "op.ply", this->ixform, std::make_pair(label, lzd_tools::get_regular_colormap(-1, max_label + 2)));
    //}

    // 迭代一次;当avg_max_diff小于阈值时停止迭代,返回-1,否则返回epoch
    // _out: 输出流。默认新建一个字符串流
    int iter() {
        // 上一个迭代（如果存在）已经结束，将_out中的内容输出到字符串流中,并更新_epoch值。注意_out是ofstream类型的，所以需要先转换为ostringstream类型
        if(avg_max_var_log.size() > 0&&_epoch<_max_iters){
            // 绘制分割线，长度为50个*，正中间为iter::epoch
            _out << "\n\n***************************************************";
            _out << "iter::epoch: " << _epoch;
            _out << " end***************************************************\n";
            std::string _out_str = _out.str();
            behavior_log.push_back(_out_str);
            _epoch++;
        }

        // 清空_out
        _out.flush();
        // 绘制分割线，长度为50，正中间为iter::epoch
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
        //update_mesh2();
        

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
        // **temporary store** the tri_face_normals of each point in point_normals for the following for-loop
        // projective_normals[i] store the normal of the i-th point in point_normals after indirect mapping or direct mapping
        std::vector<Normal<REAL, DIM>> projective_normals(_points_normals.size(), zero_normal);


        // 不论是直接映射还是间接映射,都需要在surface上采样若干个面片,并根据面片的法向量计算点的法向量。记录每个点采样的面索引
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
            // 计算三角面片的中心点。
#pragma omp parallel for
            for (int i = 0; i < mesh.second.size(); i++) {
                Point<REAL, DIM> c = mesh.first[mesh.second[i][0]] + mesh.first[mesh.second[i][1]] + mesh.first[mesh.second[i][2]];
                c /= 3;
                array<REAL, 3> a{c[0], c[1], c[2]};
                tri_centers.push_back(kdt::KDTreePoint(a));
            }
            tri_center_tree.build(tri_centers);

#pragma omp parallel for
            // 用knn找出每个点的采样面片
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


        // 根据sample_tri_idx对每个点的法向量更新
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

        // 计算每个点的sample_var
        if (!ConfigManager::get_common_config()["no_save"]) {
            cal_i_avg_var(sample_tri_idx, tri_face_normals);
            std::vector<std::pair<int, REAL>> var_idx;
            for (int i = 0; i < i_nearestSample_var.size(); i++) {
                var_idx.push_back(std::make_pair(i, i_nearestSample_var[i]));
            }
            std::sort(var_idx.begin(), var_idx.end(), [](std::pair<int, REAL> a, std::pair<int, REAL> b) {return a.second > b.second; });
            // 统计采样方差较大的点的均值
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
                // 使用pca法向量来修正
                NormalEstimation<REAL, DIM>* estimator = new PCAEstimation<REAL, DIM>();
                auto pca_normal = POINTS_NORMALS(_points_normals);
                estimator->Estimate(pca_normal);
#pragma omp parallel for
                for (int i = 0; i < _points_normals.size(); i++) {
                    // 如果法向量与pca法向量的夹角大于90度，则取反
                    if (dot<REAL,DIM>(_points_normals[i].second, pca_normal[i].second) < 0) {
                        pca_normal[i].second *= -1;
                    }
                    _points_normals[i].second = pca_normal[i].second;
                }
                update_mesh();
            }
            int zero_c = cal_i_avg_var(sample_tri_idx, tri_face_normals);
            _epoch = _max_iters-1;
            // printf("%d points haven't been selected by any face(%d totally)\n", zero_c, _points_normals.size());
            _out << "iter::" << zero_c << " points haven't been selected by any face(" << _points_normals.size() << " totally)" << endl;
        }  
        return _epoch;
    }
    // 分别取方差较大的点与方差较小的点为source和sink
    //int chose_source_sink3(int& source, int& sink) {
    //    source = -1, sink = -1;
    //    // 求方差的topk
    //    int k = i_nearestSample_var.size() / 2;
    //    if (k == 0) {
    //        k = 1;
    //        _out << "chose_source_sink3::WARNNING points less than 10!!\n";
    //    }
    //    std::vector<std::pair<int, REAL>> var_idx;
    //    for (int i = 0; i < i_nearestSample_var.size(); i++) {
    //        var_idx.push_back(std::make_pair(i, i_nearestSample_var[i]));
    //    }
    //    std::sort(var_idx.begin(), var_idx.end(), [](std::pair<int, REAL> a, std::pair<int, REAL> b) {return a.second > b.second; });
    //    int tidx = rand() % k;
    //    // 在topk中随机选择一个点作为圆心
    //    if (var_idx.size() == 0) {
    //        source = rand() % _points_normals.size();
    //        sink = rand() % _points_normals.size();
    //    }
    //    else {
    //        source = var_idx[tidx].first;
    //        sink = var_idx[var_idx.size() - tidx - 1].first;
    //    }

    //    return 1;
    //}

    //// 多次mincut,如果flow小于0则返回 (之前试过取最小者,会导致多次在相同的地方cut)
    //int multi_cut(int& source, int& sink, REAL& flow, std::vector<int>& labels, REAL radius = 0.03, int chances = 3) {
    //    source = -1, sink = -1;
    //    flow = 0;
    //    for (int i = 0; i < chances; i++) {
    //        int tsource = -1, tsink = -1;
    //        chose_source_sink3(tsource, tsink);
    //        std::vector<int> tlabels(this->_points_normals.size());// -1表示other，
    //        // 求最小流
    //        auto tflow = min_cut(_points_normals, tlabels, tsource, tsink, radius);
    //        if (tflow <0) {
    //            flow = tflow;
    //            source = tsource;
    //            sink = tsink;
    //            labels = tlabels;
    //            return 0;
    //        }
    //        else {
    //            flow = min(flow, tflow);
    //        }
    //    }
    //    return 0;
    //}

    //int chose_source_sink2(std::pair<int, int>& source_sink, REAL rmax = 0.05, REAL rmin = 0.048) {
    //    int source = -1, sink = -1;
    //    // 求方差的topk
    //    int k = i_nearestSample_var.size() / 2;
    //    std::vector<std::pair<int, REAL>> var_idx;
    //    for (int i = 0; i < i_nearestSample_var.size(); i++) {
    //        var_idx.push_back(std::make_pair(i, i_nearestSample_var[i]));
    //    }
    //    std::sort(var_idx.begin(), var_idx.end(), [](std::pair<int, REAL> a, std::pair<int, REAL> b) {return a.second > b.second; });
    //    int tidx = rand() % k;
    //    // 在topk中随机选择一个点作为圆心
    //    int center = var_idx[tidx].first;

    //    // int max_var_id = 0;
    //    // REAL max_var = 0;
    //    // for (int i = 0; i < i_nearestSample_var.size(); i++) {
    //    //     if (i_nearestSample_var[i] > max_var) {
    //    //         max_var = i_nearestSample_var[i];
    //    //         max_var_id = i;
    //    //     }
    //    // }
    //    // int center = max_var_id;

    //    // 以center为圆心,以rmax为半径得到一个大圆,以rmin为半径得到一个小圆
    //    auto rmax_idx = input_point_tree.radiusSearch(kdt::KDTreePoint({ _points_normals[center].first[0], _points_normals[center].first[1], _points_normals[center].first[2] }), rmax);
    //    auto rmin_idx = input_point_tree.radiusSearch(kdt::KDTreePoint({ _points_normals[center].first[0], _points_normals[center].first[1], _points_normals[center].first[2] }), rmin);

    //    // 排除一些边界case
    //    if (rmax_idx.size() < rmin_idx.size() + 2) {
    //        // printf("chose_source_sink2::WARNING: rmax_idx.size() < 2\n");
    //        _out << "chose_source_sink2::WARNING: rmax_idx.size() < 2\n";
    //        return -1;
    //    }

    //    // 大圆内部的点与小圆内部的点的差集
    //    std::vector<int> ring_idx;
    //    // 对两个集合排序,然后求差集
    //    std::sort(rmax_idx.begin(), rmax_idx.end());
    //    std::sort(rmin_idx.begin(), rmin_idx.end());
    //    std::set_difference(rmax_idx.begin(), rmax_idx.end(), rmin_idx.begin(), rmin_idx.end(), std::inserter(ring_idx, ring_idx.begin()));
    //    // 取第一个点为基准点pivot,连接pivot于center的向量,然后依次计算ring_idx中的点与pivot的向量的夹角,根据夹角对ring_idx中的点进行排序
    //    std::vector<std::pair<int, REAL>> ring_idx_angle;
    //    auto center_pivot = _points_normals[center].first - _points_normals[ring_idx[0]].first;
    //    for (int i = 0; i < ring_idx.size(); i++) {
    //        auto center_pi = _points_normals[center].first - _points_normals[ring_idx[i]].first;
    //        REAL angle = calculate_angle(Normal<REAL, DIM>(center_pi), Normal<REAL, DIM>(center_pivot));
    //        ring_idx_angle.push_back(std::make_pair(ring_idx[i], angle));
    //    }
    //    std::sort(ring_idx_angle.begin(), ring_idx_angle.end(), [](std::pair<int, REAL> a, std::pair<int, REAL> b) {return a.second < b.second; });
    //    // 为了简化问题,将ring的个数限制为偶数
    //    if (ring_idx_angle.size() % 2 == 1)ring_idx_angle.pop_back();

    //    int ring_count = ring_idx_angle.size(), half = ring_count / 2;
    //    int start_pos = 0;
    //    // 遍历start_pos = 0,1,2,...,half-1, 将ring_idx_angle分为等大的两部分:[start_pos, (start_pos+half-1)%count]和[(start_pos+half)%count, (start_pos+ring_count-1)%count]
    //    // 选择使得两部分的角度差最大的start_pos
    //    // 这是一个动态规划问题,可以优化为O(n)的复杂度
    //    // 第一步,先求出当start_pos=0时,两个部分的角度差
    //    REAL max_diff = 0;
    //    for (int i = 0; i < half; i++) {
    //        max_diff += abs(ring_idx_angle[(i + half) % ring_count].second - ring_idx_angle[i].second);
    //    }
    //    REAL temp_diff = max_diff;
    //    for (int i = 1; i < half; i++) {
    //        temp_diff += abs(ring_idx_angle[(i + half) % ring_count].second - ring_idx_angle[i].second);
    //        temp_diff -= ring_idx_angle[i - 1].second - ring_idx_angle[(i + half - 1) % ring_count].second;
    //        if (temp_diff > max_diff) {
    //            max_diff = temp_diff;
    //            start_pos = i;
    //        }
    //    }
    //    source = ring_idx_angle[start_pos].first;
    //    sink = ring_idx_angle[(start_pos + half) % ring_count].first;
    //    source_sink = std::make_pair(source, sink);
    //    return 0;
    //}

    //// 向外提供对ipsr实例进行切割的函数
    //int min_cut_self()
    //{
    //    // printf("*************min_cut_start*************\n");
    //    _out << "*************min_cut_start*************\n";
    //    // 使用diri+dir修复
    //    fix_update(this,diri_dir_ipsr,3);

    //    mj["iter"] = _epoch;
    //    std::vector<int> labels(this->_points_normals.size());// -1表示other，
    //    std::string path = _cmd[find_arg(this->_cmd, "--out") + 1] + "/min_cut/";
    //    if (GetFileAttributesA(path.c_str()) == INVALID_FILE_ATTRIBUTES) {
    //        CreateDirectoryA(path.c_str(), NULL);
    //    }

    //    string label_path = path + "min_cut" + std::to_string(*(this->_p_.p_ep)) + ".txt";

    // 
    //    // 方案一 随机选择两个点
    //    // int source = -1, sink = -1;

    //    // 方案二 多次mincut取flow较小者
    //    int source = -1, sink = -1;
    //    REAL flow = 0;
    //    multi_cut(source, sink, flow, labels);
    //    _out << "maxflow_before_inv" << flow << endl;
    //    mj["maxflow_before_inv"] = flow;
    //    REAL nei_radius = 0.01;
    //    if (dlabel.size() != _points_normals.size()) {
    //        std::cout << "cluster before mincut\n";
    //        dlabel.resize(this->_points_normals.size(), 0);
    //        o3d_dbscan_extraction<REAL, DIM>(this->_points_normals, dlabel, nei_radius, 1);
    //    }
    //    // 打印labels信息
    //    int source_count = 0, sink_count = 0, other_count = 0;
    //    for (int i = 0; i < labels.size(); i++) {
    //        if (labels[i] == SOURCE_PART_LABEL || labels[i] == SOURCE_LABEL) {
    //            if (dlabel[i] == dlabel[source]) {
    //                source_count++;
    //            }
    //            else {
    //                labels[i] = -1;
    //            }
    //        }
    //        else if (labels[i] == SINK_PART_LABEL || labels[i] == SINK_LABEL) {
    //            if (dlabel[i] == dlabel[sink]) {
    //                sink_count++;
    //            }
    //            else {
    //                labels[i] = -1;
    //            }
    //        }
    //        else {
    //            other_count++;
    //        }
    //    }
    //    // printf("source_count = %d, sink_count = %d, other_count = %d\n", source_count, sink_count, other_count);
    //    _out << "source_count = " << source_count << ", sink_count = " << sink_count << ", other_count = " << other_count << endl;
    //    

    //    save_current_ops("/min_cut/","ori", XForm<REAL, DIM + 1>().Identity(), std::make_pair(labels, get_mincut_colormap()));

    //    //// 按照目前的权值定义，认为flow<0是有效的分割
    //    // 用一致性判别
    //    if (true) {
    //        // printf("max_cut_self::inverse %d points\n", -1 * std::max(-1 * source_count, -1 * ((int)labels.size() - source_count)));
    //        double before_consistency = cal_save_mesh_consistency("before",get_dodecahedron_vertex);
    //        
    //        _out << "max_cut_self::inverse " << -1 * std::max(-1 * source_count, -1 * ((int)labels.size() - source_count)) << " points\n";
    //        // 选择较小的part 翻转法向量
    //        int inv_type = source_count < sink_count ? SOURCE_PART_LABEL : SINK_PART_LABEL;
    //        int inv_s_idx = source_count < sink_count ? source : sink;
    //        for (int i = 0; i < labels.size(); i++) {
    //            // 
    //            if (labels[i] == inv_type)_points_normals[i].second.normal = _points_normals[i].second.normal * -1;
    //        }

    //        // 翻转后检查
    //        double after_consistency = cal_save_mesh_consistency("after",get_dodecahedron_vertex);
    //        std::vector<int> tlabes;
    //        REAL tflow = min_cut(_points_normals, tlabes, source, sink, nei_radius);
    //        mj["maxflow_after_inv"] = tflow;
    //        mj["consistence_after_inv"] = after_consistency;
    //        mj["consistency_before_inv"] = before_consistency;

    //        if (after_consistency<before_consistency) {
    //            //printf("max_cut_self::flow = %f  -->  tflow = %f , perfer to cancel invertion\n",flow,tflow);
    //            //_out << "max_cut_self::flow = " << flow << "  -->  tflow = " << tflow << " , perfer to cancel invertion\n";
    //            for (int i = 0; i < labels.size(); i++) {
    //                if (labels[i] == inv_type)_points_normals[i].second.normal = _points_normals[i].second.normal * -1; // 使用原来的label
    //            }
    //            mj["if_flip_after_mincut"] = "false";
    //        }
    //        else {
    //            save_current_ops("/min_cut/", "after_inv", XForm<REAL, DIM + 1>().Identity(), std::make_pair(tlabes, get_mincut_colormap()));
    //            mj["if_flip_after_mincut"] = "true";
    //        }
    //    }
    //    log_j["min_cut_self_log"].push_back(mj);
    //    _out << "*************min_cut_end*************\n";
    //    return 0;
    //}

    //int min_edge_cut(){
    //    
    //}

    // 完成最后一个iter
    void end_iter() {
        // printf("iter %d is end_iter...\n", _epoch);
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
        // printf("end_iter done\n");
        _out << "end_iter done\n";
    }

    // int cal_avg_diff(float sample_rate=0.001){
    //     //TODO
    // }

    /// @brief 计算每个点所采样的法向量的均值以及方差
    /// @param mapper 采样的索引
    /// @param sample_normals 采样点的法向量 
    /// @return 返回无法映射的点的个数 
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

    // 计算每个点法向量的变化量
    // @return
    float cal_points_normals_diff(const std::vector<Normal<REAL, DIM>>& new_normals, const POINTS_NORMALS& old_op, std::vector<float>& points_normals_diff, int K) {
        assert(new_normals.size() == old_op.size());
        points_normals_diff.resize(new_normals.size());
        Normal<REAL, DIM> zero_normal(Point<REAL, DIM>(0, 0, 0));
#pragma omp parallel for
        for (int i = 0; i < new_normals.size(); i++) {
            // 由于ipsr在遇到无法映射的点时,会将其法向量保留，因此该点处的diff为0
            if (new_normals[i] == zero_normal) {
                points_normals_diff[i] = 0;
                continue;
            }
            else {
                points_normals_diff[i] = Point<REAL, DIM>::SquareNorm((new_normals[i] - old_op[i].second).normal);
            }
        }
        // 返回topk的平均值
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

    // 默认打印为{--out}/temp/it_{current_iter}.ply
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
        //std::vector<double>t;
        //output_ply(path + name, shrink_boundary(this->mesh,this->_points_normals,t,20), transform * ixform);
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
      //  auto visiable_flag = Rasterlization::random_view(cleaned_mesh);
      //  int min_val = 1e9, max_val = -1e9;
      //  for(int i=0;i<visiable_flag.size();i++){
      //      if(visiable_flag[i]>10)visiable_flag[i] = 10;
      //      min_val = std::min(min_val, visiable_flag[i]);
      //      max_val = std::max(max_val, visiable_flag[i]);
      //  }
      //  auto flagandcolormap = std::make_pair(visiable_flag,lzd_tools::get_regular_colormap(min_val,max_val));
        
        std::string pname = path + name;
        if (adname == "")pname = adname;
        else pname += adname;
        MaxAB<REAL,DIM> oc(func);
        std::vector<double> temp;
        double consistency = oc.cal_consistency(cleaned_mesh,temp,pname);
        
        printf("consistency = %f\n", consistency);

        //lzd_tools::mesh2ply(cleaned_mesh,path + name, this->ixform,flagandcolormap);
        return consistency;
    }

    // 返回存储路径根目录；会递归地创建目录
    std::string get_out_put_base(std::string addition_folder = "") {
        std::string res = _cmd[find_arg(this->_cmd, "--out") + 1] +  addition_folder;
        rmkdir(res);
        return res;
    }

    // 生成后缀名(it_xx_..)
    std::string get_resname() {
        std::string depth = _cmd[find_arg(this->_cmd, "--depth") + 1];
        std::string point_weight = _cmd[find_arg(this->_cmd, "--pointWeight") + 1];
        std::string filename = "it_" + to_string(_max_iters);
        filename += "_dp_" + depth;
        filename += "_nb_" + to_string(_k_neighbors);
        filename += "_sd_" + to_string(_seed);
        filename += "_pt_" + point_weight;
        //std::string res_name = filename;
        return filename;
    }

    int save_res(std::string adname="") {
        std::string output_path = _cmd[find_arg(this->_cmd, "--out") + 1];
        std::string depth = _cmd[find_arg(this->_cmd, "--depth") + 1];
        //std::string k_neighbors = to_string(k_neighbors);
        //std::string seed = _cmd[find_arg(this->_cmd, "--seed") + 1];

        rmkdir(output_path);
        std::string res_name = output_path + get_resname() + adname;
        
        std::vector<REAL> pmdist;
        // 最终面片
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
            // 最终点
            lzd_tools::op2ply(_points_normals, res_name + "_orientedpoint.ply", ixform);
        }
        if(j["save_option"]["GT_samples"] == true) {
            // gt
            lzd_tools::op2ply(_gt_points_normals, res_name + "_GT_samples.ply", ixform);
        }
        // 更新到原始点
        if(j["save_option"]["ori_est"] == true) {
            POINTS_NORMALS ori_est = POINTS_NORMALS(_ori_points_normals);
            auto unused_mapper = std::vector<int>();
            printf("update to origin points\n"); 
            XForm <REAL, DIM + 1> xform = ixform.inverse();
            auto dist = update_target(ori_est, xform, unused_mapper);
            printf("update to origin points done, point pair total dist = %f\n", dist);
            lzd_tools::op2ply(ori_est, res_name + "_ori_est.ply", XForm<REAL, DIM + 1>().Identity());
        }

         //// 打印variance
        //string var_path = res_name + "_variance.txt";
        //ofstream var_file(var_path);
        //for (size_t i = 0; i < i_nearestSample_var.size(); i++)
        //{
        //    var_file << i_nearestSample_var[i] << endl;
        //}
        //var_file.close();
        //
        //// 打印component_index
        //auto component_index = spilt_mesh<REAL, DIM>(cmesh).com_points;
        //string component_path = res_name + "_component.txt";
        //ofstream component_file(component_path);
        //for (size_t i = 0; i < component_index.size(); i++)
        //{
        //    for (size_t j = 0; j < component_index[i].size(); j++)
        //    {
        //        component_file << component_index[i][j] << " ";
        //    }
        //    component_file << endl;
        //}
        //component_file.close();
        //
        //// 打印points_normals_diff
        //string diff_path = res_name + "_points_normals_diff.txt";
        //ofstream diff_file(diff_path);
        //for (size_t i = 0; i < _points_normals_diff.size(); i++)
        //{
        //    diff_file << _points_normals_diff[i] << endl;
        //}
        //diff_file.close();
        //
        //// 打印pmdist
        //string pmdist_path = res_name + "_pmdist.txt";
        //ofstream pmdist_file(pmdist_path);
        //for (size_t i = 0; i < pmdist.size(); i++)
        //{
        //    pmdist_file << pmdist[i] << endl;
        //}
        //pmdist_file.close();

        // 打印log
        // 创建log文件夹
        string log_path = output_path + "iter_logs/";
        rmkdir(log_path);
        _save_log(log_path + get_resname() + "_");
        return 0;
    }

    bool if_end() {
        return _epoch == _max_iters;
    }    

    // 用目前的_points_normals来更新原始点云的法向量,返回带法向量的原始点云
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
     * @param target 待更新的点云 
     * @param transform 待更新的点云到_points_normals的变换 
     * @param mapper 用于记录target中的点在_points_normals中的索引
     * @return * REAL 用于记录每个点对的距离之和。如果很大，说明变换不合适
     */
    REAL update_target(POINTS_NORMALS& target, XForm<REAL, DIM + 1> transform, std::vector<int>& mapper) {
        REAL dist = 0;// TODO
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

template<typename REAL, int DIM>
void _select_tri(const MESH& mesh, const POINTS_NORMALS& op,
                std::vector<REAL>& pmdist,std::vector<int>& if_select,
                int k = 20)
{
    if (mesh.first.size() == 0) {
        printf("select_mesh::WARNING mesh.first.size()==%d, op.size()==%d", mesh.first.size(), op.size());
        return;
    }
    if_select = std::vector<int>(mesh.second.size(), 0);
    // 由三角面片中心构建kdtree
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
    // 遍历op中的点，找到距离op中的点最近的k个三角面片，设置if_select
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
 * 从mesh中提取距离op近的三角面片
 * @param mesh 原始点云 不作改变 
 * @param op 参考点云
 * @param pm_dist 返回op中的点到mesh中的点的平均距离 
 * @param k knn的k值
 * @return 提取后的MESH 
 */
template<typename REAL, int DIM>
MESH clean_mesh(const MESH& mesh, const POINTS_NORMALS& op,std::vector<REAL>& pmdist,int k = 20){
    if (mesh.first.size() == 0) {
        printf("clean_mesh::WARNING mesh.first.size()==%d, op.size()==%d", mesh.first.size(), op.size());
        return mesh;
    }
    std::vector<int> if_select(mesh.second.size(), 0);
    _select_tri(mesh,op,pmdist,if_select,k);    

    // 直接删除距离远的面片
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
 * clean_mesh进阶版,只从边界收缩
 * clean_mesh对于迭代来说有时候是有害的，因为会破坏原本的双层模型。
 * 我们需要clean，是因为Neumann边界条件下产生的飞边会导致法向量出错。
 * 至少在迭代期间，这个clean过程不应该产生更多的边界点
 */
template<typename REAL, int DIM>
MESH shrink_boundary(const MESH& mesh, const POINTS_NORMALS& op,std::vector<REAL>& pmdist,int k = 20){
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
            is_near_point[a] = is_near_point[a] || is_near_tri[i];// 只要所在的三角形被选中，该点即被选中
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
        for (int j = 0; j < 3; j++)is_select_tri[i] = is_select_tri[i] || is_select_point[mesh.second[i][j]]; // 这里扩大了一圈
    }
    MESH cleaned_mesh;
    cleaned_mesh.first = std::vector<Point<REAL, DIM>>(mesh.first);
    for (int i = 0; i < mesh.second.size(); i++) {
        if (is_select_tri[i])cleaned_mesh.second.push_back(mesh.second[i]);
    }
    return cleaned_mesh;
}


// 使用指定次计划外的更新方式。如果此时ipsr已经接近end（_epoch>_max_iters-fix_times），那么会导致最后实际的迭代次数大于_max_iter;
template<typename REAL, int DIM>
int fix_update(IPSR_HANDLE_P p_ipsr, int (*fix_plan)(Period& p) = diri_dir_ipsr, int fix_times = 3) {
    //printf("fix_update::记得该回去\n");
    //return 0;
    
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
