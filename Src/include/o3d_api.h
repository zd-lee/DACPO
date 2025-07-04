#pragma once
#include <open3d/geometry/PointCloud.h>
#include "Open3D/Open3D.h"
#include "tools.h"

void test_open3d();

/// @brief 对points_normals进行DBSCAN聚类，将聚类结果存储在opss中，将噪声点存储在other中
/// @param eps 查找半径 
/// @param min_points 用于判定核心点的最小邻域点数。当min_points=1时，算法退化为传统的基于密度的聚类算法
/// @return 
template <typename REAL, int DIM>
int o3d_dbscan_extraction(const ORIENTED_POINTS& points_normals,std::vector<ORIENTED_POINTS>& opss,ORIENTED_POINTS& other,REAL eps = 0.1, int min_points = 10){
    // Create a point cloud
    open3d::geometry::PointCloud pcd;
    pcd.points_.resize(points_normals.size());
    pcd.normals_.resize(points_normals.size());
    for(int i=0;i<points_normals.size();i++){
        pcd.points_[i] = Eigen::Vector3d(points_normals[i].first[0], points_normals[i].first[1], points_normals[i].first[2]);
        pcd.normals_[i] = Eigen::Vector3d(points_normals[i].second.normal[0], points_normals[i].second.normal[1], points_normals[i].second.normal[2]);
    }

    // DBSCAN
    std::vector<int> labels = pcd.ClusterDBSCAN(eps, min_points, true);
    int max_label = *max_element(labels.begin(), labels.end());
    
    // Extract clusters into opss
    opss.clear();
    opss.resize(max_label+1);
    other.clear();

    for(int i=0;i<points_normals.size();i++){
        if (labels[i] >= 0) {
            opss[labels[i]].push_back(points_normals[i]);
        }
        else {
            other.push_back(points_normals[i]);
        }
        
    }
    return 0;
}

/// @brief 对points_normals进行DBSCAN聚类，将聚类结果存储在labels中。其中，labels[i]表示第i个点所属的类别，-1表示噪声点
/// @param labels 存储分类后的标签
/// @param eps Density parameter that is used to find neighbouring points（来自pcd库）
/// @param min_points 用于判定核心点的最小邻域点数。当min_points=1时，算法退化为传统的基于密度的聚类算法
/// @return 
template <typename REAL, int DIM>
int o3d_dbscan_extraction(const ORIENTED_POINTS& points_normals,std::vector<int>& labels,REAL eps = 0.1, int min_points = 10){
    // Create a point cloud
    open3d::geometry::PointCloud pcd;
    pcd.points_.resize(points_normals.size());
    for(int i=0;i<points_normals.size();i++){
        pcd.points_[i] = Eigen::Vector3d(points_normals[i].first[0], points_normals[i].first[1], points_normals[i].first[2]);
    }
    // DBSCAN
    labels = pcd.ClusterDBSCAN(eps, min_points, true);
    return 0;
}


/**
 * @brief 
 * 对点云进行HPR，将可见的点的索引放在hop_idxs中
 * @param points_normals 点云
 * @param hpr_idxs 可见点的索引
 * @param cam  相机位置
 * @param radius HPR算法的参数，越大则可见点越多 
 * @return int 可见点的数量
 */
template <typename REAL, int DIM>
int HPR(const ORIENTED_POINTS& points_normals, std::vector<size_t>& hop_idxs, Point<REAL,DIM> cam, REAL radius = 1000){
    open3d::geometry::PointCloud pcd;
    pcd.points_.resize(points_normals.size());
    for (int i = 0; i < points_normals.size(); i++)
    {
        pcd.points_[i] = Eigen::Vector3d(points_normals[i].first[0], points_normals[i].first[1], points_normals[i].first[2]);
    }
    //REAL diameter = pcd.GetMaxBound() - pcd.GetMinBound();
    Eigen::Vector3d camera_pos(cam[0], cam[1], cam[2]);
    auto res = pcd.HiddenPointRemoval(camera_pos, radius);
    hop_idxs = std::get<1>(res);
    return hop_idxs.size();
}

/**
 * @brief 
 * 对点云用一个随机的视点进行HPR. 
 * @param points_normals 
 * @param hop_idxs 
 * @return int 
 */
template <typename REAL, int DIM>
int rand_HPR(const ORIENTED_POINTS& points_normals, std::vector<int>& if_visiable) {
    printf("random hpr...\n");
    open3d::geometry::PointCloud pcd;
    pcd.points_.resize(points_normals.size());
    REAL maxlen = 0;
    for (int i = 0; i < points_normals.size(); i++)
    {
        pcd.points_[i] = Eigen::Vector3d(points_normals[i].first[0], points_normals[i].first[1], points_normals[i].first[2]);
        //maxlen = max(maxlen, pcd.points_[i].norm());
    }
    REAL diameter = (pcd.GetMaxBound() - pcd.GetMinBound()).norm();
    Eigen::Vector3d camera_pos = Eigen::Vector3d::Random();
    camera_pos.normalize();
    camera_pos *= diameter;
    camera_pos = camera_pos + pcd.GetCenter();
    auto res = pcd.HiddenPointRemoval(camera_pos, diameter*1000);
    auto residx = std::get<1>(res);
    if_visiable.resize(points_normals.size());
    for(auto idx:residx){
        if_visiable[idx] = 1;
    }
    printf("random hpr done\n");
    return residx.size();
}

template <typename REAL, int DIM>
int o3d_norm_estimate(ORIENTED_POINTS& points_normals, REAL eps = 0.1, int min_points = 10) {
    // Create a point cloud
    open3d::geometry::PointCloud pcd;
    pcd.points_.resize(points_normals.size());
    pcd.normals_.resize(points_normals.size());
    for (int i = 0; i < points_normals.size(); i++) {
        pcd.points_[i] = Eigen::Vector3d(points_normals[i].first[0], points_normals[i].first[1], points_normals[i].first[2]);
    }
    pcd.EstimateNormals();
    for (int i = 0; i < points_normals.size(); i++) {
        Normal<REAL, DIM> n;
        for(int j = 0;j<DIM;j++)n.normal[j] = pcd.normals_[i][j];
        points_normals[i].second = n;
    }
    return 0;
}

template <typename REAL, int DIM>
open3d::t::geometry::TriangleMesh convertMesh2o3dMesh(const MESH& mesh) {
    open3d::core::Tensor vertices({ (long long)mesh.first.size(), DIM }, open3d::core::Dtype::Float32);
    open3d::core::Tensor triangles({ (long long)mesh.second.size(), 3 }, open3d::core::Dtype::Int32);
    
    for (int i = 0; i < mesh.first.size(); i++) {
        for (int j = 0; j < DIM; j++) {
            vertices[i][j] = mesh.first[i][j];
        }
    }
    for (int i = 0; i < mesh.second.size(); i++) {
        for (int j = 0; j < 3; j++) {
            triangles[i][j] = mesh.second[i][j];
        }
    }
    return open3d::t::geometry::TriangleMesh(vertices, triangles);
}

template <typename REAL, int DIM>
MESH subdivision_mesh(const MESH& mesh,int times = 1) {
    open3d::geometry::TriangleMesh o3d_mesh;
    for (int i = 0; i < mesh.first.size(); i++) {
        o3d_mesh.vertices_.push_back(Eigen::Vector3d(mesh.first[i][0], mesh.first[i][1], mesh.first[i][2]));
    }
    for (int i = 0; i < mesh.second.size(); i++) {
        o3d_mesh.triangles_.push_back(Eigen::Vector3i(mesh.second[i][0], mesh.second[i][1], mesh.second[i][2]));
    }
    auto nmesh = o3d_mesh.SubdivideMidpoint(times);
    MESH res;
    for (int i = 0; i < nmesh->vertices_.size(); i++) {
        Point<REAL, DIM> p;
        for (int j = 0; j < DIM; j++) {
            p[j] = nmesh->vertices_[i][j];
        }
        res.first.push_back(p);
    }
    for (int i = 0; i < nmesh->triangles_.size(); i++) {
        std::vector<int> t;
        for (int j = 0; j < 3; j++) {
            t.push_back(nmesh->triangles_[i][j]);
        }
        res.second.push_back(t);
    }
    return res;
}