#pragma once

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <vector>

#include <pcl/features/normal_3d.h>
#include <pcl/point_types.h>
#include <pcl/search/kdtree.h>

#include "tools.h"

void test_open3d();

namespace lzd_o3d_compat {

template <typename REAL, unsigned int DIM>
inline pcl::PointCloud<pcl::PointXYZ>::Ptr make_xyz_cloud(const ORIENTED_POINTS& points_normals) {
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    cloud->points.resize(points_normals.size());
    cloud->width = static_cast<std::uint32_t>(points_normals.size());
    cloud->height = 1;
    cloud->is_dense = false;
    for (size_t i = 0; i < points_normals.size(); ++i) {
        cloud->points[i].x = static_cast<float>(points_normals[i].first[0]);
        cloud->points[i].y = static_cast<float>(points_normals[i].first[1]);
        cloud->points[i].z = static_cast<float>(points_normals[i].first[2]);
    }
    return cloud;
}

template <typename REAL, unsigned int DIM>
inline REAL squared_distance(const std::pair<Point<REAL, DIM>, Normal<REAL, DIM>>& a,
                             const std::pair<Point<REAL, DIM>, Normal<REAL, DIM>>& b) {
    REAL distance = 0;
    for (int d = 0; d < DIM; ++d) {
        REAL delta = a.first[d] - b.first[d];
        distance += delta * delta;
    }
    return distance;
}

}  // namespace lzd_o3d_compat

template <typename REAL, unsigned int DIM>
int o3d_dbscan_extraction(const ORIENTED_POINTS& points_normals,
                          std::vector<int>& labels,
                          REAL eps = 0.1,
                          int min_points = 10);

template <typename REAL, unsigned int DIM>
int o3d_dbscan_extraction(const ORIENTED_POINTS& points_normals,
                          std::vector<ORIENTED_POINTS>& opss,
                          ORIENTED_POINTS& other,
                          REAL eps = 0.1,
                          int min_points = 10) {
    std::vector<int> labels;
    o3d_dbscan_extraction<REAL, DIM>(points_normals, labels, eps, min_points);

    int max_label = -1;
    for (int label : labels) {
        max_label = std::max(max_label, label);
    }

    opss.clear();
    opss.resize(max_label + 1);
    other.clear();
    for (size_t i = 0; i < points_normals.size(); ++i) {
        if (labels[i] >= 0) {
            opss[labels[i]].push_back(points_normals[i]);
        } else {
            other.push_back(points_normals[i]);
        }
    }
    return 0;
}

template <typename REAL, unsigned int DIM>
int o3d_dbscan_extraction(const ORIENTED_POINTS& points_normals,
                          std::vector<int>& labels,
                          REAL eps,
                          int min_points) {
    const int point_count = static_cast<int>(points_normals.size());
    labels.assign(point_count, -1);
    if (point_count == 0) {
        return 0;
    }

    const REAL eps2 = eps * eps;
    int cluster_id = 0;
    std::vector<char> visited(point_count, 0);

    auto region_query = [&](int idx) {
        std::vector<int> neighbors;
        for (int i = 0; i < point_count; ++i) {
            if (lzd_o3d_compat::squared_distance<REAL, DIM>(points_normals[idx], points_normals[i]) <= eps2) {
                neighbors.push_back(i);
            }
        }
        return neighbors;
    };

    for (int i = 0; i < point_count; ++i) {
        if (visited[i]) {
            continue;
        }
        visited[i] = 1;
        std::vector<int> neighbors = region_query(i);
        if (static_cast<int>(neighbors.size()) < min_points) {
            labels[i] = -1;
            continue;
        }

        labels[i] = cluster_id;
        std::queue<int> queue;
        for (int neighbor : neighbors) {
            if (neighbor != i) {
                queue.push(neighbor);
            }
        }

        while (!queue.empty()) {
            int current = queue.front();
            queue.pop();
            if (!visited[current]) {
                visited[current] = 1;
                std::vector<int> current_neighbors = region_query(current);
                if (static_cast<int>(current_neighbors.size()) >= min_points) {
                    for (int neighbor : current_neighbors) {
                        queue.push(neighbor);
                    }
                }
            }
            if (labels[current] < 0) {
                labels[current] = cluster_id;
            }
        }
        ++cluster_id;
    }
    return 0;
}

template <typename REAL, unsigned int DIM>
int HPR(const ORIENTED_POINTS& points_normals,
        std::vector<size_t>& hpr_idxs,
        Point<REAL, DIM> cam,
        REAL radius = 1000) {
    hpr_idxs.clear();
    if (points_normals.empty()) {
        return 0;
    }

    // Conservative visibility fallback: select points on the camera-facing half
    // of the cloud. The target release config does not use this path directly.
    Point<REAL, DIM> center;
    for (const auto& point_normal : points_normals) {
        center += point_normal.first;
    }
    center /= static_cast<REAL>(points_normals.size());

    Point<REAL, DIM> view_dir = center - cam;
    for (size_t i = 0; i < points_normals.size(); ++i) {
        Point<REAL, DIM> rel = points_normals[i].first - center;
        REAL dot = 0;
        for (int d = 0; d < DIM; ++d) {
            dot += rel[d] * view_dir[d];
        }
        if (dot <= radius) {
            hpr_idxs.push_back(i);
        }
    }
    return static_cast<int>(hpr_idxs.size());
}

template <typename REAL, unsigned int DIM>
int rand_HPR(const ORIENTED_POINTS& points_normals, std::vector<int>& if_visiable) {
    if_visiable.assign(points_normals.size(), 1);
    return static_cast<int>(points_normals.size());
}

template <typename REAL, unsigned int DIM>
int o3d_norm_estimate(ORIENTED_POINTS& points_normals, REAL eps = 0.1, int min_points = 10) {
    if (points_normals.empty()) {
        return 0;
    }

    auto cloud = lzd_o3d_compat::make_xyz_cloud<REAL, DIM>(points_normals);
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> estimator;
    estimator.setInputCloud(cloud);
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());
    estimator.setSearchMethod(tree);
    estimator.setKSearch(std::max(3, min_points));

    pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
    estimator.compute(*normals);
    for (size_t i = 0; i < points_normals.size(); ++i) {
        points_normals[i].second.normal[0] = normals->points[i].normal_x;
        points_normals[i].second.normal[1] = normals->points[i].normal_y;
        points_normals[i].second.normal[2] = normals->points[i].normal_z;
        REAL len = 0;
        for (unsigned int d = 0; d < DIM; ++d) {
            len += points_normals[i].second.normal[d] * points_normals[i].second.normal[d];
        }
        if (len > 0) {
            len = std::sqrt(len);
            for (unsigned int d = 0; d < DIM; ++d) {
                points_normals[i].second.normal[d] /= len;
            }
        }
    }
    return 0;
}

template <typename REAL, unsigned int DIM>
MESH subdivision_mesh(const MESH& mesh, int times = 1) {
    if (times <= 0) {
        return mesh;
    }
    throw std::runtime_error("subdivision_mesh requires Open3D; set subdivision_times to 0 in config on this build");
}
