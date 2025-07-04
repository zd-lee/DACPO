#pragma once
#include <Eigen/Dense>
#include <tools.h>

template<typename REAL>
using E_FIELD_GRAD = std::vector<std::vector<REAL>>;
template<typename REAL>
using INTERSACTION_MATRIX = std::vector<std::vector<REAL>>;


template<typename REAL, int DIM>
E_FIELD_GRAD<REAL> field_grad(const ORIENTED_POINTS& sources, const ORIENTED_POINTS& means, REAL _eps) {
    E_FIELD_GRAD<REAL> E_total(means.size(), std::vector<REAL>(DIM, 0));
#pragma omp parallel for
    for (int i = 0; i < sources.size(); i++) {
        for (int j = 0; j < means.size(); j++) {
            auto p = sources[i].second.normal;
            REAL R[DIM];
            for (int k = 0; k < DIM; k++) {
                R[k] = sources[i].first[k] - means[j].first[k];
            }
            REAL R_norm = 0;
            for (int k = 0; k < DIM; k++) {
                R_norm += R[k] * R[k];
            }
            R_norm = sqrt(R_norm);
            if (R_norm <= 0) {
                continue; // 不考虑重合的点
            }
            REAL R_unit[DIM];
            for (int k = 0; k < DIM; k++) {
                R_unit[k] = R[k] / R_norm;
            }
            REAL R_norm_cube = R_norm * R_norm * R_norm;
            for (int k = 0; k < DIM; k++) E_total[j][k] += (3 * p[k] * R_unit[k] * R_unit[k] - p[k]) / (R_norm_cube + _eps);
        }
    }
    return E_total;
};

template<typename REAL, int DIM>
REAL cal_field_intersaction(const ORIENTED_POINTS& sources, const ORIENTED_POINTS& means, REAL eps) {
    E_FIELD_GRAD<REAL> st_E = field_grad(sources, means, eps);
    REAL st_interaction = 0;
    for (int i = 0; i < st_E.size(); i++) {
        for (int j = 0; j < DIM; j++) {
            st_interaction += st_E[i][j] * means[i].second.normal[j];
        }
    }
    E_FIELD_GRAD<REAL> ts_E = field_grad(means, sources, eps);
    REAL ts_interaction = 0;
    for (int i = 0; i < ts_E.size(); i++) {
        for (int j = 0; j < DIM; j++) {
            ts_interaction += ts_E[i][j] * sources[i].second.normal[j];
        }
    }
    return st_interaction + ts_interaction;
};

namespace FIELD_PROPAGATION {
    class MaskField {
    public:
        int C;
        double eps;

        MaskField(int C_,double eps_)
            : C(C_), eps(eps_) {}

        /**
         * @brief
         * @param source N*6
         * @param target N*6
         * @param mask N*K
         * @return Eigen::MatrixXd N*K. @return[i][j] = source[mask[i][j]]对target[i]产生的影响
         */
        Eigen::MatrixXd mask_xie_intersection(const Eigen::MatrixXd& source,
            const Eigen::MatrixXd& target,
            const Eigen::MatrixXi& mask);

		nlohmann::json get_config() {
			nlohmann::json j;
			j["C"] = C;
			j["eps"] = eps;
			return j;
		}

    private:

        /**
         * @brief
         * @param source N*6
         * @param target 1*6
         * @param mask 1*K
		 * @return Eigen::MatrixXd K*3
         */
        Eigen::MatrixXd mask_xie_field(const Eigen::MatrixXd& source,
            const Eigen::RowVectorXd& target,
            const Eigen::VectorXi& knn_mask);
    };


    struct KnnGraph {
        int _k;
        double _radius;
        std::vector<std::vector<int>> _graph_idx;
        std::vector<std::vector<double>> _graph_weight;
        const Eigen::MatrixXd& _pts;

        KnnGraph(const Eigen::MatrixXd& pts, int k, double radius);

        /**
         * @brief Construct a new Knn Graph without radius and weight
         * @param pts
         * @param k
         */
        KnnGraph(const Eigen::MatrixXd& pts, int k);
        /**
         * @brief
         * @param startpoints
         * @return std::vector<int> 遍历的顺序
         */
        std::vector<int> get_BFS_route(int startpoints);

        /**
         * @brief 得到邻居最少的点的邻居数
         * @return int 
         */
        int get_shortest_nei();

	};


    class MaskFieldPropagation {
        MaskField _field;
        int _times;
        bool _diffuse;
        int mask_size;
        int knngraph_k;
        double knngraph_radius;
        std::vector<bool> _xie_propagation_points_in_order(const Eigen::MatrixXd& pts, const std::vector<int>& order, 
            const Eigen::MatrixXd& pre_caled_field_intersaction_mat, const Eigen::MatrixXi& KnnMask);



    public:
        MaskFieldPropagation(MaskField field_, int times_, int mask_size_, int knngraph_k_, double knngraph_radius_, bool diffuse)
            : _field(field_), _times(times_), mask_size(mask_size_), knngraph_k(knngraph_k_), knngraph_radius(knngraph_radius_),_diffuse(diffuse) {
        }

        std::vector<bool> xie_propagation_points_on_bfs_tree(const Eigen::MatrixXd& pts);
        
		nlohmann::json get_config() {
			nlohmann::json j;
			j["field"] = _field.get_config();
			j["times"] = _times;
			j["mask_size"] = mask_size;
			j["knngraph_k"] = knngraph_k;
			j["knngraph_radius"] = knngraph_radius;
            j["diffuse"] = _diffuse;
			return j;
		}
    };



}

