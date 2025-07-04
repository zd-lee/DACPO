#pragma once
#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/normal_3d.h>
#include <pcl/search/kdtree.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>
#include <iomanip> // for setw, setfill
#include "tools.h"

int pcd_seg(char * path);

//#define ORIENTED_POINTS std::vector<std::pair<Point<REAL, DIM>, Normal<REAL, DIM>>>
template<typename REAL, int DIM>
int est_norm(const std::vector <Point<REAL, DIM>>& points, std::vector <Normal<REAL, DIM>>& normals, int k_neighbor = 5) {
    auto point_cloud = xyz2pcl(points);
    // create the normal estimation class, and pass the input dataset to it
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
    ne.setInputCloud(point_cloud);

    // pass the original data (before downsampling) as the search surface
    ne.setSearchSurface(point_cloud);

    // create an empty kdtree representation, and pass it to the normal estimation object.
    // its content will be filled inside the object, based on the given surface dataset.
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());
    ne.setSearchMethod(tree);

    // output datasets
    pcl::PointCloud<pcl::Normal>::Ptr cloud_normals(new pcl::PointCloud<pcl::Normal>);


    // TODO 注意： radius过小会导致出现nan
    // use all neighbors in a sphere of radius 3cm
    //ne.setRadiusSearch(0.02);
    
    ne.setKSearch(10);

    // compute the features
    ne.compute(*cloud_normals);

    // copy the normals to normals
    normals.resize(cloud_normals->size());
    for (int i = 0; i < cloud_normals->size(); i++) {
        for (int j = 0; j < DIM; j++) {
            normals[i].normal[j] = cloud_normals->points[i].normal[j];
        }
    }
    return 0;
}


template <typename REAL, int DIM>
pcl::PointCloud<pcl::PointXYZ>::Ptr xyz2pcl(const std::vector<Point<REAL, DIM>>& points){
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
    for(int i=0;i<points.size();i++){
        pcl::PointXYZ p;
        for(int j=0;j<DIM;j++){
            p.data[j] = points[i][j];
        }
        cloud->push_back(p);
    }
    return cloud;
}

template <typename REAL, int DIM>
std::vector<Point<REAL, DIM>> pcl2xyz(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud){
    std::vector<Point<REAL, DIM>> points;
    for(int i=0;i<cloud->size();i++){
        Point<REAL, DIM> p;
        for(int j=0;j<DIM;j++){
            p[j] = cloud->points[i].data[j];
        }
        points.push_back(p);
    }
    return points;
}

template <typename REAL, int DIM>
pcl::PointCloud<pcl::PointNormal>::Ptr op2pcl(const std::vector<std::pair<Point<REAL, DIM>, Normal<REAL, DIM>>>& points_normals){
    pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>);
    for(int i=0;i<points_normals.size();i++){
        pcl::PointNormal pn;
        for(int j=0;j<DIM;j++){
            pn.data[j] = points_normals[i].first[j];
            pn.normal[j] = points_normals[i].second.normal[j];
        }
        cloud->push_back(pn);
    }
    return cloud;
}

template <typename REAL, int DIM>
std::vector<std::pair<Point<REAL, DIM>, Normal<REAL, DIM>>> pcl2op(pcl::PointCloud<pcl::PointNormal>::Ptr cloud){
    std::vector<std::pair<Point<REAL, DIM>, Normal<REAL, DIM>>> points_normals;
    for(int i=0;i<cloud->size();i++){
        std::pair<Point<REAL, DIM>, Normal<REAL, DIM>> pn;
        for(int j=0;j<DIM;j++){
            pn.first[j] = cloud->points[i].data[j];
            pn.second.normal[j] = cloud->points[i].normal[j];
        }
        points_normals.push_back(pn);
    }
    return points_normals;
}

/**
 * @brief 
 * 使用plane segmentation算法将点云分割为多个平面和剩余点
 * @param other_points 
 * @param threshold 
 * @return int 
 */
template <typename REAL, int DIM>
int pcd_plane_seg(const std::vector<Point<REAL, DIM>>& points, std::vector<std::vector<Point<REAL, DIM>>>& planes, std::vector<Point<REAL, DIM>>& other_points, REAL threshold = 0.02){
    // 将点云转换为pcl格式
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = xyz2pcl(points);
    std::cout << "PointCloud before filtering has: " << cloud->size () << " data points." << std::endl; //*

    // Create the filtering object: downsample the dataset using a leaf size of 1cm
    pcl::VoxelGrid<pcl::PointXYZ> vg;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered (new pcl::PointCloud<pcl::PointXYZ>);
    vg.setInputCloud (cloud);
    vg.setLeafSize (threshold, threshold, threshold);
    vg.filter (*cloud_filtered);
    std::cout << "PointCloud after filtering has: " << cloud_filtered->size ()  << " data points." << std::endl; //*

    // Create the segmentation object for the planar model and set all the parameters
    pcl::SACSegmentation<pcl::PointXYZ> seg;
    pcl::PointIndices::Ptr inliers (new pcl::PointIndices);
    pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_plane (new pcl::PointCloud<pcl::PointXYZ> ());

    seg.setOptimizeCoefficients (true);
    seg.setModelType (pcl::SACMODEL_PLANE); // 设置模型类型，检测平面
    seg.setMethodType (pcl::SAC_RANSAC); // 设置方法，随机采样一致性方法
    seg.setMaxIterations (100); // 设置最大迭代次数
    seg.setDistanceThreshold (threshold); // 设置距离阈值，点到平面的距离小于该阈值，认为该点属于该平面

    int nr_points = (int) cloud_filtered->size ();

    while(cloud_filtered->size () > 0.3 * nr_points){
        // Segment the largest planar component from the remaining cloud
        seg.setInputCloud (cloud_filtered); // 设置输入点云
        seg.segment (*inliers, *coefficients); // 分割点云，获得平面内点的索引和平面参数

        if (inliers->indices.size () == 0){
            std::cout << "Could not estimate a planar model for the given dataset." << std::endl;
            break;
        }

        // Extract the planar inliers from the input cloud
        pcl::ExtractIndices<pcl::PointXYZ> extract; // 创建点云提取对象
        extract.setInputCloud (cloud_filtered); // 设置输入点云
        extract.setIndices (inliers); // 设置分割出来的平面内点的索引
        extract.setNegative (false); // 设置提取内点还是外点，true表示提取外点，false表示提取内点

        // Get the points associated with the planar surface
        extract.filter (*cloud_plane); // 提取平面内点
        std::cout << "PointCloud representing the planar component: " << cloud_plane->size () << " data points." << std::endl;

        // 将平面内点转换为自定义的点格式
        std::vector<Point<REAL, DIM>> plane = pcl2xyz<REAL,DIM>(cloud_plane);
        planes.push_back(plane);

        // Remove the planar inliers, extract the rest
        extract.setNegative (true); // 设置提取外点
        extract.filter (*cloud_filtered); // 提取外点
    }

    // 将剩余的点转换为自定义的点格式
    for(int i=0;i<cloud_filtered->size();i++){
        Point<REAL, DIM> p;
        for(int j=0;j<DIM;j++){
            p[j] = cloud_filtered->points[i].data[j];
        }
        other_points.push_back(p);
    }

    return 0;
}

template <typename REAL, int DIM>
int pcl_euclidean_cluster_extraction(const std::vector<Point<REAL, DIM>>& points, std::vector<std::vector<Point<REAL, DIM>>>& clusters, std::vector<Point<REAL, DIM>>& other_points,REAL threshold = -1.0){
    // 将点云转换为pcl格式
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = xyz2pcl(points);
    std::cout << "PointCloud before filtering has: " << cloud->size () << " data points." << std::endl; //*

    bool if_downsample = false;
    auto cloud_filtered(cloud);
    // Create the filtering object: downsample the dataset using a leaf size of 1cm
    if (threshold > 0) {
        pcl::VoxelGrid<pcl::PointXYZ> vg;
        vg.setInputCloud(cloud);
        vg.setLeafSize(threshold, threshold, threshold);
        vg.filter(*cloud);
        std::cout << "PointCloud after filtering has: " << cloud_filtered->size() << " data points." << std::endl; //*
    }

    // Creating the KdTree object for the search method of the extraction
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);

    // 设置欧式聚类对象的参数
    pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
    ec.setClusterTolerance (threshold); // 设置近邻搜索的搜索半径
    //ec.setMinClusterSize (100); // 设置一个聚类需要的最少的点数目
    //ec.setMaxClusterSize (25000); // 设置一个聚类需要的最大点数目
    ec.setSearchMethod (tree); // 设置点云的搜索机制
    ec.setInputCloud (cloud_filtered); // 输入点云
    std::vector<pcl::PointIndices> cluster_indices; // 保存聚类的索引结果
    ec.extract (cluster_indices); // 执行欧式聚类，保存聚类结果到cluster_indices中

    // 遍历每一个聚类
    for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices.begin (); it != cluster_indices.end (); ++it){
        std::vector<Point<REAL, DIM>> cluster;
        // 遍历每一个聚类中的点
        for (std::vector<int>::const_iterator pit = it->indices.begin (); pit != it->indices.end (); ++pit){
            Point<REAL, DIM> p;
            for(int j=0;j<DIM;j++){
                p[j] = cloud_filtered->points[*pit].data[j];
            }
            cluster.push_back(p);
        }
        clusters.push_back(cluster);
    }
    return 0;
}

template <typename REAL, int DIM>
int pcl_euclidean_cluster_extraction(ORIENTED_POINTS points_normals,std::vector<ORIENTED_POINTS>& opss,ORIENTED_POINTS& other,REAL threshold = -1.0){
    std::vector<Point<REAL, DIM>> raw_points;
	std::vector<std::vector<Point<REAL,DIM>>> raw_planes;
	std::vector<Point<REAL,DIM>> rest_point;
	for (auto& p : points_normals) {
		raw_points.push_back(p.first);
	}
	int code = pcl_euclidean_cluster_extraction<REAL, DIM>(raw_points, raw_planes, rest_point);

    // 由raw_planes建立kd树
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = xyz2pcl(raw_points);
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
    tree->setInputCloud (cloud);

    // 对聚类中的每个点,在raw_points得到的kd树中找到最近的点,并将在points_normals中对应的点的法向量赋值给该点
    for(int i=0;i<raw_planes.size();i++){
        std::vector<std::pair<Point<REAL, DIM>, Normal<REAL, DIM>>> ops;
        for(auto& p:raw_planes[i]){
            std::vector<int> pointIdxNKNSearch(1);
            std::vector<float> pointNKNSquaredDistance(1);
            //convert p to pcl::PointXYZ
            pcl::PointXYZ p_pcl;
            for(int j=0;j<DIM;j++){
                p_pcl.data[j] = p[j];
            }
            tree->nearestKSearch (p_pcl, 1, pointIdxNKNSearch, pointNKNSquaredDistance);
            ops.push_back(std::make_pair(p,points_normals[pointIdxNKNSearch[0]].second));
        }
        opss.push_back(ops);
    }
    // 对rest_point中的每个点,在raw_points得到的kd树中找到最近的点,并将在points_normals中对应的点的法向量赋值给该点
    if(rest_point.size() <= 0)return 0;
    for(auto& p:rest_point){
        std::vector<int> pointIdxNKNSearch(1);
        std::vector<float> pointNKNSquaredDistance(1);
        //convert p to pcl::PointXYZ
        pcl::PointXYZ p_pcl;
        for(int j=0;j<DIM;j++){
            p_pcl.data[j] = p[j];
        }
        tree->nearestKSearch (p_pcl, 1, pointIdxNKNSearch, pointNKNSquaredDistance);
        other.push_back(points_normals[pointIdxNKNSearch[0]]);
    }
    return 0;
}


// 参考:https://pcl.readthedocs.io/projects/tutorials/en/master/min_cut_segmentation.html#min-cut-based-segmentation
// 