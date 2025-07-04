#include "tools.h"
#include "NormalEstimation.h"

/**
 * @brief 
 * 判断是否更新某个点
 */
template<typename REAL, int DIM>
class Update_filter {
public:
    // 请使用And
    virtual std::vector<int> if_update(const std::vector<Normal<REAL, DIM>>& project_normal) = 0;
    virtual nlohmann::json get_config() = 0;

    /**
     * @brief 
     * @param update_flag 原本的更新flag 
     * @param project_normal 
     * @return int 本次被拒绝的点的数量
     */
    int And(std::vector<int>& update_flag,const std::vector<Normal<REAL, DIM>>& project_normal){
        std::vector<int> new_update_flag = if_update(project_normal);
        int refusenum = 0;
        for(int i = 0;i<update_flag.size();i++){
            if(update_flag[i] == 1 && new_update_flag[i] == 0){
                refusenum++;
            }
            update_flag[i] = update_flag[i] & new_update_flag[i];
        }
        return refusenum;
    }
};

template<typename REAL, int DIM>
class NonZeroFilter : public Update_filter<REAL, DIM> {
    nlohmann::json _config_j;
    nlohmann::json _log_j;
public:
    NonZeroFilter() {
        _config_j["name"] = "NonZeroFilter";
    }

    nlohmann::json get_config() {
        return _config_j;
    }

    nlohmann::json get_log() {
        return _log_j;
    }
    std::vector<int> if_update(const std::vector<Normal<REAL, DIM>>& project_normal) {
        std::vector<int> res(project_normal.size(), 1);
        int refusenum = 0;
        for (int i = 0; i < project_normal.size(); i++) {
            if (project_normal[i] == Normal<REAL, DIM>(Point<REAL, DIM>(0, 0, 0))) {
                res[i] = 0;
                refusenum++;
            }
        }
        _log_j["refusenum"] = refusenum;
        return res;
    }
};


/**
 * @brief 
 * 对于与init结果的差距大于minthreshold且小于maxthreshold的点,拒绝更新
 */
template<typename REAL, int DIM>
class InitFilter : public Update_filter<REAL, DIM> {
    ORIENTED_POINTS init_res;
    nlohmann::json _config_j;
    nlohmann::json _log_j;
    REAL minthreshold;
    REAL maxthreshold;

public:

    InitFilter(NormalEstimation<REAL, DIM>* estimator,
        REAL minthreshold, REAL maxthreshold,
        const ORIENTED_POINTS& op) :
        minthreshold(minthreshold), maxthreshold(maxthreshold) {
        init_res = op;
        estimator->Estimate(init_res);
        _config_j["estimator"] = estimator->get_config();
    }

    InitFilter(nlohmann::json config, const ORIENTED_POINTS& op) {
		minthreshold = config["minthreshold"];
		maxthreshold = config["maxthreshold"];
		init_res = op;
        NormalEstimation<REAL,DIM>* estimator = get_estimator_from_json<REAL,DIM>(config["estimator"]);
        estimator->Estimate(init_res);
        _config_j["estimator"] = estimator->get_config();
	}

    // 直接使用现成的initop
    InitFilter(ORIENTED_POINTS& inited_op, REAL minthreshold, REAL maxthreshold) {
        init_res = inited_op;
        minthreshold = minthreshold;
        maxthreshold = maxthreshold;
        _config_j["estimator"] = "directly use inited_op";
    }


    nlohmann::json get_config() {
        _config_j["name"] = "InitFilter";
        _config_j["minthreshold"] = minthreshold;
        _config_j["maxthreshold"] = maxthreshold;
        return _config_j;
    }

    nlohmann::json get_log() {
        return _log_j;
    }

    std::vector<int> if_update(const std::vector<Normal<REAL, DIM>>& project_normal) {
        int refusenum = 0;
        assert(project_normal.size() == init_res.size());
        std::vector<int> res(project_normal.size(), 1);
        for (int i = 0; i < project_normal.size(); i++) {
            REAL rangle = calculate_angle(init_res[i].second, project_normal[i]);
            REAL angle = rangle * 180 / M_PI;
            if (angle > minthreshold && angle < maxthreshold) {
                res[i] = 0;
                refusenum++;
            }
        }
        _log_j["refusenum"] = refusenum;
        return res;
    }
};

