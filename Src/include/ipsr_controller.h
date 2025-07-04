#pragma once 
#include <IPSR.h>


/***********************************ipsr status filter***********************************************/
/**
 * @brief 判断ipsr是否处于某种状态
 */
template<typename REAL, int DIM>
class BasicIPSRFilter {
public:
    virtual bool check(const IPSR_HANDLE_P handle) = 0;
    virtual nlohmann::json get_config() = 0;
};

template<typename REAL, int DIM>
class ModFilter : public BasicIPSRFilter<REAL,DIM> {
    int mod;
public:
    bool check(const IPSR_HANDLE_P handle) {
        // 有时候_epoch与迭代的次数会不一致
        return handle->avg_max_var_log.size() % mod == 0;
    }

    ModFilter(int mod = 5) {
        this->mod = mod;
    }

    nlohmann::json get_config() {
        nlohmann::json j;
        j["name"] = "ModFilter";
        j["mod"] = mod;
        return j;
    }
};

template<typename REAL, int DIM>
class ShockingIPSRFilter : public BasicIPSRFilter<REAL,DIM> {
    int window_size;//窗口大小 查看最近window_size个epoch的信息    
public:
    ShockingIPSRFilter(int window_size = 5) {
        this->window_size = window_size;
    }

    bool check(const IPSR_HANDLE_P handle) {
        bool flag = true;
        BasicIPSRFilter<REAL, DIM>* f = new ModFilter<REAL, DIM>(5);// 防止不断reinit。
        flag = flag && f->check(handle);
        int count_inc = 0, count_dec = 0;
        auto var_log = handle->avg_max_var_log;
        if (var_log.size() < window_size) {
            return false;
        }
        for (int i = var_log.size() - window_size; i < var_log.size(); i++) {
            if (var_log[i] > var_log[i - 1]) {
                count_inc++;
            }
            else if (var_log[i] < var_log[i - 1]) {
                count_dec++;
            }
        }
        flag = flag && (count_dec < count_inc);
        return flag;
    }

    nlohmann::json get_config() {
        nlohmann::json j;
        j["name"] = "ShockingIPSRFilter";
        j["window_size"] = window_size;
        return j;
    }
};

template<typename REAL, int DIM>
class SufficientEpochFilter : public BasicIPSRFilter<REAL, DIM> {
    int _at_least_epoch_count;

public:
    SufficientEpochFilter(int at_least_epoch_count = 10) {
        this->_at_least_epoch_count = at_least_epoch_count;
    }

    bool check(const IPSR_HANDLE_P handle) {
        return handle->_epoch <= handle->_max_iters - _at_least_epoch_count;
    }

    nlohmann::json get_config() {
        nlohmann::json j;
        j["name"] = "SufficientEpochFilter";
        j["at_least_epoch_count"] = _at_least_epoch_count;
        return j;
    }
};

template<typename REAL, int DIM>
class ForbidFilter : public BasicIPSRFilter<REAL, DIM> {
public:
    bool check(const IPSR_HANDLE_P handle) {
        return false;
    }

    nlohmann::json get_config() {
        nlohmann::json j;
        j["name"] = "ForbidFilter";
        return j;
    }
};


template<typename REAL, int DIM>
bool filter_ipsr(const IPSR_HANDLE_P handle,const std::vector<BasicIPSRFilter<REAL,DIM>*> filters) {
    for (auto filter : filters) {
        if (!filter->check(handle)) {
            return false;
        }
    }
    return true;
};


/***********************************ipsr vistor***********************************************/
template<typename REAL, int DIM>
class IpsrController{
public:
    virtual int iter(IPSR_HANDLE_P _handle) = 0;
    virtual nlohmann::json get_log() = 0;
    virtual nlohmann::json get_config() = 0;
};

/**
 * @brief 
 * 初步规划: 负责ipsr_hp的reinit、mincut等行为的调用
 */
template<typename REAL, int DIM>
class SingleIpsrController : public IpsrController<REAL,DIM> {
    //std::vector<BasicIPSRFilter<REAL,DIM>*> _reinit_filters;
    //std::vector<BasicIPSRFilter<REAL,DIM>*> _mincut_filters;
    int _minimum_to_iter; // op数量小于这个值时不执行。
    // 增加配置项时，记得在get_config中添加
    int _save_op;// 每隔_save_op轮存一个op
    int _save_mesh;
    NormalEstimation<REAL, DIM>* _estimator; //如不为空,则第一次调用iter时使用

public:
    //
    SingleIpsrController(
        //std::vector<BasicIPSRFilter<REAL, DIM>*> reinit_filters = std::vector<BasicIPSRFilter<REAL, DIM>*>(),
        //std::vector<BasicIPSRFilter<REAL, DIM>*> mincut_filters = std::vector<BasicIPSRFilter<REAL, DIM>*>(),
        int save_op = -1,
        int save_mesh = -1,
        NormalEstimation<REAL, DIM>* estimator = nullptr
        //,int minimum_to_iter = 100
        ):_estimator(estimator) {
        //this->_reinit_filters = reinit_filters;
        //this->_mincut_filters = mincut_filters;
        //_minimum_to_iter = minimum_to_iter;
        _save_op = save_op;
        _save_mesh = save_mesh;
    }

    // static SingleIpsrController* classic_ipsr_controller()

    nlohmann::json get_config() {
        nlohmann::json j;
        //j["minimum_to_iter"] = _minimum_to_iter;
        //for (auto filter : _reinit_filters) {
        //    j["reinit_filters"].push_back(filter->get_config());
        //}
        //for (auto filter : _mincut_filters) {
        //    j["mincut_filters"].push_back(filter->get_config());
        //}
        return j;
    }

    int iter(IPSR_HANDLE_P _handle) {
        //if (_handle->_points_normals.size() < _minimum_to_iter) {
        //    printf("_points_normals size = %d,too few points to iter, skip!!\n",_handle->_points_normals.size());
        //    return 0;
        //}
        if ((_handle->_epoch == -1 || _handle->_epoch == 0) && _estimator != nullptr)_handle->init_op_normal(_estimator);
        int ep = _handle->iter();
         if (_save_op >0 && (ep % _save_op) == 0) {
         	_handle->save_current_ops("/temp/");
         }
         if (_save_mesh > 0 && (ep % _save_mesh) == 0) {
             _handle->save_current_mesh("/temp/");
         }
        //if (filter_ipsr(_handle, _reinit_filters)) {
        //    int seed = rand();
        //     _handle->init_op_normal(new RandomInit<REAL, DIM>(seed));
        //}

        //if (filter_ipsr(_handle, _mincut_filters)) {
        //    _handle->min_cut_self();
        //}
        return ep;
    }

    nlohmann::json get_log() {
        assert(false);
        // controller的log,定义不明晰。因为ipsr_handle也有log。
        nlohmann::json j;
        j["reinit_times"] = -1;
        j["mincut_times"] = -1;
        return j;
    }

};


/**
 * @brief 
 * 存储原void ipsr(...)的参数 统一管理
 * @tparam REAL 
 * @tparam DIM 
 */
template<typename REAL, int DIM>
class IPSR_Factory{
    std::vector<std::string> _cmd;
    std::string _input_name;
    std::string _output_path;
    int _iters;
    int _k_neighbors;
    int _seed; // 历史遗留问题;
    double _pointweight;
    
    // Factory只是统一存储参数, 便于在不同场景创建ipsrhp;
    nlohmann::json get_config() {
        assert(false);
    }

public:

    // 使用原入口参数
    IPSR_Factory(const std::string& input_name, const std::string& output_path,
    int iters, double pointweight, int depth, int k_neighbors,int seed = 0
    ){
        std::string command = "PoissonRecon --in " + input_name + " --out " + output_path + "  --bType 2 --depth " + to_string(depth) + " --pointWeight " + to_string(pointweight);
        _cmd = split(command);
        _input_name = input_name;
        _output_path = output_path;
        _iters = iters;
        _k_neighbors = k_neighbors;
        _seed = seed;
        _pointweight = pointweight;
    }

    IPSR_Factory(const IPSR_Factory<REAL,DIM>& factory){
        _cmd = factory._cmd;
        _input_name = factory._input_name;
        _output_path = factory._output_path;
        _iters = factory._iters;
        _k_neighbors = factory._k_neighbors;
        _seed = factory._seed;
    }

    // 获得IPSR对象; op对象根据input_name读取
    IPSR_HANDLE_P create_single_ipsr_handle(
        int (*_update_plan)(Period& p),
        NormalEstimation<REAL, DIM>* estimator
    ){
        POINTS_NORMALS op;
        ply_reader<REAL,DIM>(_input_name,op);
        return create_ipsr_from_op(op,_update_plan,estimator);
    }

    // 获得以任意op为输入的IPSR对象
    IPSR_HANDLE_P create_ipsr_from_op(
        POINTS_NORMALS op,
        int (*_update_plan)(Period& p),
        NormalEstimation<REAL, DIM>* estimator
    ){
        IPSR_HANDLE_P handle = new ipsr_handle<REAL,DIM>(_cmd,_iters,_k_neighbors,_update_plan,op,_seed);
        handle->init(estimator);
        return handle;
    }
    
    void change_self_args(std::string arg, std::string value) {
        change_args(_cmd, arg, value);
    }


/******************************************************************一些常用的ipsr handle****************************************************************************/

    IPSR_HANDLE_P get_classic_neumman_ipsr_handle(int seed = _seed){
        return create_single_ipsr_handle(neumman_ind_ipsr,new RandomInit<REAL, DIM>(seed));
    }
    
    IPSR_HANDLE_P get_classic_dirichlet_ipsr_handle(int seed = _seed) {
        return create_single_ipsr_handle(diri_ind_ipsr,new RandomInit<REAL, DIM>(seed));
    }
    

/******************************************************************一些依赖这组参数的help func****************************************************************************/
    
    // 单纯的降采样,同时变换回原空间
    // 注意，psr的_sample_points实质是完成了三件事：0.归一化点 1.使用八叉树降采样 2.给降采样后的点赋权。这个方法相当于只用到了1的结果 
    // @depth 降采样的层数。如果不指定，则不进行下采样
    POINTS_NORMALS down_sample_points(const POINTS_NORMALS& ori, int depth = -1){
        if (depth == -1) {
            return POINTS_NORMALS(ori);
        }
        XForm<REAL,DIM+1> iXForm;
        std::vector<double> weight_sample;
        std::vector<char*> argv_str(_cmd.size());
        std::vector<std::string> tcmd(_cmd);
        change_args(tcmd, "--depth", to_string(depth));
        for (int i = 0; i < tcmd.size(); i++) {
            char* p = new char[tcmd[i].size() + 1];
            strcpy(p, tcmd[i].c_str());
            argv_str[i] = p;
        }
        POINTS_NORMALS op = sample_points_entrance<REAL,DIM>((int)argv_str.size(),argv_str.data(),ori,iXForm,&weight_sample);
        // 将op变换回原空间
        for (int i = 0; i < op.size(); i++) {
            op[i].first = iXForm * op[i].first;
        }
        printf("down sample points from %d to %d\n",(int)ori.size(), (int)op.size());
        return op;
    }

    // 从点云中重建mesh
    // 如果depth==-1，则使用默认的depth=10来重建
    MESH get_mesh_from_points_normals(const POINTS_NORMALS& op){
        XForm<REAL, DIM + 1> iXForm;
        std::vector<double> weight_sample;
        std::vector<char*> argv_str(_cmd.size());
        std::vector<std::string> temp_cmd(_cmd);

        if (temp_cmd[find_arg(temp_cmd, "--depth") + 1] == "-1") {
            printf("WARNING::--depth == -1, use 10 for alternate\n");
            change_args(temp_cmd, "--depth", to_string(10));
        }
        for (int i = 0; i < temp_cmd.size(); i++) {
            argv_str[i] = new char[temp_cmd[i].size() + 1];
            strcpy(argv_str[i], temp_cmd[i].c_str());
        }
        POINTS_NORMALS sampled_op = sample_points_entrance<REAL,DIM>((int)argv_str.size(),argv_str.data(),op,iXForm,&weight_sample);
        auto mesh = poisson_reconstruction_entrance<REAL,DIM>((int)argv_str.size(),argv_str.data(),sampled_op,&weight_sample);
        // 变换回原空间
        for (int i = 0; i < mesh.first.size(); i++) {
            mesh.first[i] = iXForm * mesh.first[i];
        }
        return mesh;
    } 

	void save_grid(const POINTS_NORMALS& op, std::string filename) {
        XForm<REAL, DIM + 1> iXForm;
        std::vector<double> weight_sample;
        std::vector<std::string> temp_cmd(_cmd);
        temp_cmd.push_back("--voxel");
        temp_cmd.push_back(filename);
        
        //temp_cmd.push_back(filename);

        if (temp_cmd[find_arg(temp_cmd, "--depth") + 1] == "-1") {
            printf("WARNING::--depth == -1, use 10 for alternate\n");
            change_args(temp_cmd, "--depth", to_string(10));
        }
        std::vector<char*> argv_str(temp_cmd.size());
        for (int i = 0; i < temp_cmd.size(); i++) {
            argv_str[i] = new char[temp_cmd[i].size() + 1];
            strcpy(argv_str[i], temp_cmd[i].c_str());
        }
        POINTS_NORMALS sampled_op = sample_points_entrance<REAL, DIM>((int)argv_str.size(), argv_str.data(), op, iXForm, &weight_sample);
        poisson_reconstruction_entrance<REAL, DIM>((int)argv_str.size(), argv_str.data(), sampled_op, &weight_sample);
        return;
    }

    std::string get_out_put_base(std::string addition_folder = "") {
        std::string res = _cmd[find_arg(this->_cmd, "--out") + 1] + addition_folder;
        rmkdir(res);
        return res;
    }

    std::string get_resname() {
        std::string depth = _cmd[find_arg(this->_cmd, "--depth") + 1];
        std::string filename = "it_" + to_string(_iters);
        filename += "_dp_" + depth;
        filename += "_nb_" + to_string(_k_neighbors);
        filename += "_sd_" + to_string(_seed);
        filename += "_pt_" + to_string(_pointweight);
        //std::string res_name = filename;
        return filename;
    }

};