#pragma once 
#include <IPSR.h>


/***********************************ipsr status filter***********************************************/
/**
 * @brief Checks whether an IPSR handle is in a specific state.
 */
template<typename REAL, unsigned int DIM>
class BasicIPSRFilter {
public:
    virtual bool check(const IPSR_HANDLE_P handle) = 0;
    virtual nlohmann::json get_config() = 0;
};

template<typename REAL, unsigned int DIM>
class ModFilter : public BasicIPSRFilter<REAL,DIM> {
    int mod;
public:
    bool check(const IPSR_HANDLE_P handle) {
        // _epoch can differ from the iteration count.
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

template<typename REAL, unsigned int DIM>
class ShockingIPSRFilter : public BasicIPSRFilter<REAL,DIM> {
    int window_size;
public:
    ShockingIPSRFilter(int window_size = 5) {
        this->window_size = window_size;
    }

    bool check(const IPSR_HANDLE_P handle) {
        bool flag = true;
        BasicIPSRFilter<REAL, DIM>* f = new ModFilter<REAL, DIM>(5);
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

template<typename REAL, unsigned int DIM>
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

template<typename REAL, unsigned int DIM>
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


template<typename REAL, unsigned int DIM>
bool filter_ipsr(const IPSR_HANDLE_P handle,const std::vector<BasicIPSRFilter<REAL,DIM>*> filters) {
    for (auto filter : filters) {
        if (!filter->check(handle)) {
            return false;
        }
    }
    return true;
};


/***********************************ipsr vistor***********************************************/
template<typename REAL, unsigned int DIM>
class IpsrController{
public:
    virtual int iter(IPSR_HANDLE_P _handle) = 0;
    virtual nlohmann::json get_log() = 0;
    virtual nlohmann::json get_config() = 0;
};

/**
 * @brief 
 * Coordinates actions such as reinitialization and min-cut updates.
 */
template<typename REAL, unsigned int DIM>
class SingleIpsrController : public IpsrController<REAL,DIM> {
    int _minimum_to_iter;
    int _save_op;
    int _save_mesh;
    NormalEstimation<REAL, DIM>* _estimator;

public:
    //
    SingleIpsrController(
        int save_op = -1,
        int save_mesh = -1,
        NormalEstimation<REAL, DIM>* estimator = nullptr
        //,int minimum_to_iter = 100
        ):_estimator(estimator) {
        //this->_reinit_filters = reinit_filters;
        //this->_mincut_filters = mincut_filters;
        _save_op = save_op;
        _save_mesh = save_mesh;
    }

    // static SingleIpsrController* classic_ipsr_controller()

    nlohmann::json get_config() {
        nlohmann::json j;
        //j["minimum_to_iter"] = _minimum_to_iter;
        //    j["reinit_filters"].push_back(filter->get_config());
        //    j["mincut_filters"].push_back(filter->get_config());
        return j;
    }

    int iter(IPSR_HANDLE_P _handle) {
        if ((_handle->_epoch == -1 || _handle->_epoch == 0) && _estimator != nullptr)_handle->init_op_normal(_estimator);
        int ep = _handle->iter();
         if (_save_op >0 && (ep % _save_op) == 0) {
         	_handle->save_current_ops("/temp/");
         }
         if (_save_mesh > 0 && (ep % _save_mesh) == 0) {
             _handle->save_current_mesh("/temp/");
         }
        //     _handle->init_op_normal(new RandomInit<REAL, DIM>(seed));
        //    _handle->min_cut_self();
        return ep;
    }

    nlohmann::json get_log() {
        assert(false);
        // Controller-level log; IPSR handles also keep their own logs.
        nlohmann::json j;
        j["reinit_times"] = -1;
        j["mincut_times"] = -1;
        return j;
    }

};


/**
 * @brief 
 * Stores the former ipsr(...) parameters in one place.
 * @tparam REAL 
 * @tparam DIM 
 */
template<typename REAL, unsigned int DIM>
class IPSR_Factory{
    std::vector<std::string> _cmd;
    std::string _input_name;
    std::string _output_path;
    int _iters;
    int _k_neighbors;
    int _seed;
    int _psr_threads;
    double _pointweight;
    
    nlohmann::json get_config() {
        assert(false);
    }

public:

    IPSR_Factory(const std::string& input_name, const std::string& output_path,
    int iters, double pointweight, int depth, int k_neighbors, int seed = 0, int psr_threads = 4
    ){
        if (psr_threads < 1) {
            psr_threads = 1;
        }
        std::string command = "PoissonRecon --in " + input_name + " --out " + output_path + "  --bType 2 --depth " + to_string(depth) + " --pointWeight " + to_string(pointweight) + " --threads " + to_string(psr_threads);
        _cmd = split(command);
        _input_name = input_name;
        _output_path = output_path;
        _iters = iters;
        _k_neighbors = k_neighbors;
        _seed = seed;
        _psr_threads = psr_threads;
        _pointweight = pointweight;
    }

    IPSR_Factory(const IPSR_Factory<REAL,DIM>& factory){
        _cmd = factory._cmd;
        _input_name = factory._input_name;
        _output_path = factory._output_path;
        _iters = factory._iters;
        _k_neighbors = factory._k_neighbors;
        _seed = factory._seed;
        _psr_threads = factory._psr_threads;
    }

    IPSR_HANDLE_P create_single_ipsr_handle(
        int (*_update_plan)(Period& p),
        NormalEstimation<REAL, DIM>* estimator
    ){
        POINTS_NORMALS op;
        ply_reader<REAL,DIM>(_input_name,op);
        return create_ipsr_from_op(op,_update_plan,estimator);
    }

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


/****************************************************************** Common IPSR handles ****************************************************************************/

    IPSR_HANDLE_P get_classic_neumman_ipsr_handle(int seed = -1){
        if (seed < 0) seed = _seed;
        return create_single_ipsr_handle(neumman_ind_ipsr,new RandomInit<REAL, DIM>(seed));
    }
    
    IPSR_HANDLE_P get_classic_dirichlet_ipsr_handle(int seed = -1) {
        if (seed < 0) seed = _seed;
        return create_single_ipsr_handle(diri_ind_ipsr,new RandomInit<REAL, DIM>(seed));
    }
    

/****************************************************************** Helpers using the factory parameters ***********************************************************/
    
    // Down-sample points and transform them back to the original space.
    // PSR sampling also normalizes and weights points; this helper only keeps the octree sampling result.
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
        // Transform sampled points back to the original space.
        for (int i = 0; i < op.size(); i++) {
            op[i].first = iXForm * op[i].first;
        }
        printf("down sample points from %d to %d\n",(int)ori.size(), (int)op.size());
        return op;
    }

    // Reconstruct a mesh from oriented points. If depth == -1, use depth 10.
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
        // Transform the mesh back to the original space.
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
        return filename;
    }

};
