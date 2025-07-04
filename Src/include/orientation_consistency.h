#pragma once
#include "o3d_api.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>

#define NOT_SEEN -1 // 面片不可见 
#define SEE_INNER 1 // 看到面片内侧
#define SEE_OUTER 2 // 看到面片外侧
#define VIEW_POINT 3 // 观察点
#define LOOK_AT_POINT 4 // 观察目标点
#define PARALLEL 5 // 切面片

namespace Rasterlization{
   /**
     给定points, triangles, projection_matrix, view_matrix, image_width, image_height
     0. 初始化z_buffer,pixel_source.其中pixel_source[i,j]表示像素(i,j)对应的三角形的索引;
        初始化visibility_list,其中visibility_list[i]表示三角形i的可见性
        初始化triangle_depths,其中triangle_depths[i]表示三角形i的深度
     1. 将points投影到观察空间 view * xyzd,得到vxyzd
     2. 遍历三角形,设三角形的深度为三角形的三个顶点的最小z值
     3. 将vxyz通过正交/透视投影到裁剪空间 cxyzd = projection_mat * vxyzd,其中cx,cy,cz的范围是[-1, 1]
     4. 将cxyd通过视口变换,映射到像素空间,得到pxy
        像素空间的原点在图像的左上角,x轴向右,y轴向下,范围是[0, image_width-1]和[0, image_height-1]
        通过下面的变换将cxy映射到像素空间
         px = (cx + 1) * (image_width - 1) / 2
         py = (cy + 1) * (image_height - 1) / 2
     5. 遍历三角形,根据三角形的三个顶点的像素坐标,计算三角形的边界,得到min_x, max_x, min_y, max_y
        遍历边界内的像素,判断像素是否在三角形内部
        如果在三角形内部,计算像素的深度d,如果d小于z_buffer[y,x],则更新z_buffer[y,x] = d,并更新pixel_source[y,x] = i
     6. 遍历pixel_source,根据pixel_source更新visibility_list
    *
    */
   // 定义三角形顶点
   struct Triangle {
      int v1, v2, v3;
   };

   // 定义顶点
   struct Vertex {
      float x, y, z;
   };

   enum PROJECTION_TYPE {
       ORTHOGRAPHIC = 0,
       PERSPECTIVE   = 1
   };

   struct Image {
      enum BoundaryType {
            NOT_BOUNDARY = 0,
            DEPTH_BOUNDARY = 1,
            MESH_BOUNDARY = 2
      };

      int width, height;
      Eigen::MatrixXd zBuffer;
      Eigen::MatrixXi is_inside; //视点是否看向三角形内侧 1表示看向内侧 0表示看向外侧 -1表示无内容
      Eigen::MatrixXi triangle_id; // 三角形的索引 -1表示没有内容
      Eigen::MatrixXi cluster_id; // 簇的索引 -1表示无内容
      Eigen::MatrixXi is_boundary; // 是否是边界像素 0表示不是边界像素 1表示是depth边界，2表示是朝向边界且不是mesh边界
      
      double max_depth_gap;

      Image(int width, int height) : width(width), height(height) {
            zBuffer = Eigen::MatrixXd::Constant(height, width, 1e9);
            triangle_id = Eigen::MatrixXi::Constant(height, width, -1);
            cluster_id = Eigen::MatrixXi::Constant(height, width, -1);
            is_boundary = Eigen::MatrixXi::Constant(height, width, -1);
            is_inside = Eigen::MatrixXi::Constant(height, width, -1);
            max_depth_gap = 0.0f;
      }
      int getWidth() const { return width; }
      int getHeight() const { return height; }

      std::vector<std::pair<int,int>> get_neibor_pixel_idx(int i, int j); 
      std::vector<std::pair<int,int>> get_not_empty_neibor_pixel_idx(int i, int j);
      bool is_empty(int i, int j);

      void cal_boundary(double eps);

      /**
       * @brief 
       * 将由不同tri_type的三角形得到的像素分开。
       * boundary = len(unique(neighbor(tri_type[triangle_id[y,x]]))) == 2
       * @param tri_type 
       * @param flag boundary == true时，将is_boundary设为flag
       * @return * void 
       */
      void cal_boundary(std::vector<int>& tri_type,int flag = MESH_BOUNDARY);

      void reset() {
            zBuffer.setConstant(1e9);
            triangle_id.setConstant(-1);
            cluster_id.setConstant(-1);
            is_boundary.setConstant(-1);
            is_inside.setConstant(-1);
      }
      void draw(std::string path);

      void draw_cluster(std::string path);

   };

   class HasEdge {
   public:
       //    virtual bool operator()(const PixelData& a, const PixelData& b) = 0;
       virtual bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b, const Image* img) = 0;
   };


   /**
    * 定义联通: 两个像素点的深度值差值小于eps,且两个像素点的is_inside相同
    */
   class ContinousPixel : public HasEdge {
   private:
       double _eps;
   public:
       ContinousPixel(double eps) :_eps(eps) {}
       //   bool operator()(const PixelData& a, const PixelData& b){
       //         return abs(a.depth() - b.depth()) < _eps && a.is_inside() == b.is_inside();
       //   }
       bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b, const Image* img) {
            return abs(img->zBuffer(a.first, a.second) - img->zBuffer(b.first, b.second)) < _eps && img->is_inside(a.first, a.second) == img->is_inside(b.first, b.second);
       }
   };

   std::vector<std::vector<std::pair<int,int>>> ImageBfsOnNotEmptyPixel(Image* img, HasEdge* has_edge);

   // 视口变换 将裁剪空间的平面坐标映射到像素空间   
   Eigen::Vector2d viewportTransform(const Eigen::Vector2d& cxy, int imageWidth, int imageHeight);

   // 判断点是否在三角形内部
   bool isInsideTriangle(const Eigen::Vector2d& pxyP, const Eigen::Vector2d& pxyA, const Eigen::Vector2d& pxyB, const Eigen::Vector2d& pxyC);

   // 根据三角形的三个顶点的像素坐标更新像素源和深度缓冲区
   void updatePixelDepthBuffer(Eigen::MatrixXd& zBuffer, Eigen::MatrixXi& pixelSource, const Eigen::Vector2d& pxyA, const Eigen::Vector2d& pxyB, const Eigen::Vector2d& pxyC, double d, int i);

   // 计算视角矩阵
   // 要求up与viewPoints-lookat不平行 如果平行,则会随机生成一个up
   Eigen::Matrix4d calculateViewMatrix(const Eigen::Vector3d& viewPoints, const Eigen::Vector3d& lookat, const Eigen::Vector3d& up);

   // 计算正交投影矩阵
   Eigen::Matrix4d getOrthoMatrix(double lx, double rx, double by, double ty, double n, double f);
   
   /**
    * @brief 
    * @param triangles 索引数据
    * @param vertices 顶点数据
    * @param viewPoints 视点
    * @param lookat 注视点
    * @param up 上方
    * @param imageWidth 分辨率-宽度
    * @param imageHeight 分辨率-高度
    * @return * std::vector<int> 
    */
   std::vector<int> _computeVisibility(const std::vector<Triangle>& triangles, const std::vector<Vertex>& vertices, const Eigen::Vector3d& viewPoints,
                           const Eigen::Vector3d& lookat, const Eigen::Vector3d& up, int imageWidth, int imageHeight, PROJECTION_TYPE projection);
   
   void _rasterlization(Image* img, 
                        const std::vector<Triangle>& triangles, const std::vector<Vertex>& vertices,
                        const Eigen::Vector3d& viewPoints,const Eigen::Vector3d& lookat, const Eigen::Vector3d& up,
                        std::vector<int> &visibilityList, std::vector<double> &triangles_depths, PROJECTION_TYPE projection);

   /**
    * @brief 
    * 计算给定mesh的可视性
    * @param mesh surface 
    * @param viewPoints 相机位置
    * @param lookat 相机看向的点
    * @param up 用于辅助求view坐标系的一个正交轴
    * @param imageWidth 分辨率参数 越大越精确 
    * @param imageHeight 分辨率参数 越大越精确
    */
   template<typename REAL, int DIM>
   std::vector<int> computeVisibility(const MESH& mesh, const Point<float, 3>& viewPoints, const Point<float, 3>& lookat, const Point<float, 3>& up, int imageWidth = 1920, int imageHeight = 1080) {
       printf("caculating visibility for mesh...\n");
       std::vector<Rasterlization::Triangle> triangles;
       std::vector<Rasterlization::Vertex> vertices;
       for (auto& face : mesh.second)
       {
           Rasterlization::Triangle triangle;
           triangle.v1 = face[0];
           triangle.v2 = face[1];
           triangle.v3 = face[2];
           triangles.push_back(triangle);
       }
       for (auto& vertex : mesh.first)
       {
           Rasterlization::Vertex v;
           v.x = vertex[0];
           v.y = vertex[1];
           v.z = vertex[2];
           vertices.push_back(v);
       }
       Eigen::Vector3d viewPoints_eigen(viewPoints[0], viewPoints[1], viewPoints[2]);
       Eigen::Vector3d lookat_eigen(lookat[0], lookat[1], lookat[2]);
       Eigen::Vector3d up_eigen(up[0], up[1], up[2]);
       printf("start...\n");
       auto res =  _computeVisibility(triangles, vertices, viewPoints_eigen, lookat_eigen, up_eigen, imageWidth, imageHeight);
       printf("calculate visibility done!\n");
       return res;
   }


   /**
    * @brief 
    * 随机选取视角 返回该视角下三角形的可视性
    * @param mesh 
    * @return std::vector<int> 
    */
   template<typename REAL, int DIM>
   std::vector<int> random_view(const MESH& mesh,int image_width = 1920, int image_height = 1080) {
       // 随机生成观察点
       int randidx = rand() % mesh.first.size();
       auto viewPoints = mesh.first[randidx];
       // 随机生成观察目标点
       randidx = rand() % mesh.first.size();
       auto lookat = mesh.first[randidx];
       // 向上方向向量为(0, 1, 0)
       Point<float, 3> up(0.0, 1.0, 0.0);
       return computeVisibility(mesh, viewPoints, lookat, up, image_width, image_height);   
   }

   // TODO 可以改成矩阵运算
   // abandon 这只适用于透视投影
   // bool view_orientation(Eigen::Vector3d view_point, Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c); 

   double  _orientation_consistence(Image *img,const std::vector<Triangle>& triangles, const std::vector<Vertex>& vertices, const Eigen::Vector3d& viewPoints,
                           const Eigen::Vector3d& lookat, const Eigen::Vector3d& up, std::vector<int>& orient_flag,PROJECTION_TYPE projection);


   /**
    * @brief 
    * 将mesh转换为triangle和vertex 
    */
   template<typename REAL, int DIM>
   void mesh2tv(const MESH& mesh, std::vector<Triangle>& triangles, std::vector<Vertex>& vertices){
      triangles.clear();
      vertices.clear();
      for (auto& face : mesh.second)
      {
          Rasterlization::Triangle triangle;
          triangle.v1 = face[0];
          triangle.v2 = face[1];
          triangle.v3 = face[2];
          triangles.push_back(triangle);
      }
      for (auto& vertex : mesh.first)
      {
          Rasterlization::Vertex v;
          v.x = vertex[0];
          v.y = vertex[1];
          v.z = vertex[2];
          vertices.push_back(v);
      }
   }

}

/**
 * 随机选取若干个点,计算这些点的几何均值,作为view_point(防止view_point和lookat重合)
 */
template<typename REAL,int DIM>
Point<REAL, DIM> innner_random_view_point(const MESH& mesh,int K = 5) {
    Point<REAL, DIM> view_point;
    for (int i = 0; i < DIM; i++) {
        view_point[i] = 0;
    }
    for (int i = 0; i < K; i++) {
        int randidx = rand() % mesh.first.size();
        auto viewPoints = mesh.first[randidx];
        for (int j = 0; j < DIM; j++) {
            view_point[j] += viewPoints[j];
        }
    }
    for (int i = 0; i < DIM; i++) {
        view_point[i] /= K;
    }
    return view_point;
}

typedef std::vector<Eigen::Vector3d> (*VertexGetterP) ();

// 生成正四面体的顶点
std::vector<Eigen::Vector3d> get_tetrahedron_vertex();

// 生成正十二面体的面中心点
std::vector<Eigen::Vector3d> get_dodecahedron_vertex();

/**
 * @brief Get the viewpoints lookat object
 * @param mesh 待观察的mesh
 * @param view_points 返回的观察点 
 * @param lookat 返回的观察方向
 * @param up 返回的辅助向量
 * @param vertex_getter 生成观察点的基函数 
 * @return int 
 */
template<typename REAL,int DIM>
int get_viewpoints_lookat(const MESH& mesh,
    std::vector<Eigen::Vector3d>& view_points, std::vector<Eigen::Vector3d>& lookat,std::vector<Eigen::Vector3d>& up,
    VertexGetterP vertex_getter = get_tetrahedron_vertex) 
{
    // 计算mesh的中心点
    Point<REAL, DIM> mesh_center;
    for (int i = 0; i < DIM; i++) {
        mesh_center[i] = 0;
    }
    for (auto& vertex : mesh.first) {
        mesh_center += vertex;
    }
    mesh_center /= mesh.first.size();
    
    // 计算mesh的半径
    REAL mesh_radius = 0;
    for (auto& vertex : mesh.first) {
        mesh_radius = max(mesh_radius, Point<REAL, DIM>::SquareNorm((vertex - mesh_center)));
    }
    mesh_radius = std::sqrt(mesh_radius);
    Eigen::Vector3d emesh_center({ mesh_center[0],mesh_center[1],mesh_center[2] });

    // 生成顶点
    auto vertices = vertex_getter();
    // 将正四面体的四个顶点移动到mesh_center
    for (auto& vertex : vertices) {
        vertex = vertex * mesh_radius + emesh_center;
    }
    // 从centerpoint到vertices的向量,以及vertices到centerpoint的向量
    view_points.clear();
    lookat.clear();
    up.clear();
    for (auto& vertex : vertices) {
        view_points.push_back(vertex);
        lookat.push_back(emesh_center);
        up.push_back(Eigen::Vector3d(0, 0, 1));

        // FIXME 由于是直接用z值来判断远近,这部分有问题
        //view_points.push_back(emesh_center);
        //lookat.push_back(vertex);
        //up.push_back(Eigen::Vector3d(0, 0, 1));
    }
    return view_points.size();
}

// 已经被取代。请使用类OrientationConsistency
// /**
//  * @brief
//  * 随机选取若干个视角
//  * 对于每个视角,计算可见部分面片的内外侧一致性
//  *    内外侧一致性 consitency = max(内侧面片数, 外侧面片数) / (内侧面片数 + 外侧面片数)
//  *    判断view_point 所见的面片为内侧/外侧的方式:
//  *       连接view_point和面片的重心,如果该向量与面片法向量的夹角小于90度,则该面片为内侧面片,否则为外侧面片
//  * 返回所有视角的内外侧一致性均值
//  * @param mesh
//  * @param path 如果path不为空,则将每个视角下的可见面片写入ply文件
//  * @return double 
//  */
// template<typename REAL,int DIM>
// double orientation_consistence(const MESH& mesh, std::string path = "", VertexGetterP func = get_tetrahedron_vertex) {
//     printf("start caculate orientation_consistency...\n");
//     double sum_consistency = 0;
//     std::vector<Rasterlization::Triangle> triangles;
//     std::vector<Rasterlization::Vertex> vertices;
//     Rasterlization::mesh2tv(mesh, triangles, vertices);
//     std::vector<int> orient_flag;
//     std::vector<Eigen::Vector3d> view_points;
//     std::vector<Eigen::Vector3d> lookat;
//     std::vector<Eigen::Vector3d> up;
//     int times = get_viewpoints_lookat(mesh, view_points, lookat, up, func);
//     for(int i=0;i<view_points.size();i++){
//         auto consistency = Rasterlization::_orientation_consistence(triangles, vertices, view_points[i], lookat[i], up[i],orient_flag,1920, 1080);
//         sum_consistency += consistency;
//         if (path != "") {
//             auto vp = Point<REAL, DIM>({ view_points[i][0],view_points[i][1],view_points[i][2] });
//             auto lkp = Point<REAL, DIM>({ lookat[i][0] ,lookat[i][1] ,lookat[i][2] });
//             auto tmesh1 = lzd_tools::add_sphere(mesh, vp,0.01);
//             auto tmesh2 = lzd_tools::add_sphere(tmesh1, lkp,0.01);
//             // 同时改变orient_flag
//             for(int j = mesh.second.size();j<tmesh1.second.size();j++){
//                 orient_flag.push_back(VIEW_POINT);
//             }
//             for(int j = tmesh1.second.size();j<tmesh2.second.size();j++){
//                 orient_flag.push_back(LOOK_AT_POINT);
//             }            
//             //lzd_tools::mesh2ply(tmesh2, path +"_time_" + std::to_string(i) + "_" + std::to_string(consistency) + ".ply", XForm<REAL, DIM + 1>().Identity(), std::make_pair(orient_flag, get_orient_color()));
//             lzd_tools::mesh2ply(tmesh2, path + "_time_" + std::to_string(i) + ".ply", XForm<REAL, DIM + 1>().Identity(), std::make_pair(orient_flag, get_orient_color()));
//         }
//     }
//     printf("caculate orientation_consistency done!\n");
//     return sum_consistency / times;
// }
// template<typename REAL,int DIM>


lzd_tools::FiniteMap get_orient_color();
lzd_tools::FiniteMap get_100_color();


template<typename REAL, int DIM>
static std::string get_vertex_getter_name(VertexGetterP _func){
    if(_func == get_tetrahedron_vertex){
        return "get_tetrahedron_vertex";
    }else if(_func == get_dodecahedron_vertex){
        return "get_dodecahedron_vertex";
    }else{
        assert(false);
        return "unknown";
    }
}

template<typename REAL,int DIM>
static VertexGetterP get_vertex_getter_by_name(std::string name){
    if(name == "get_tetrahedron_vertex"){
        return get_tetrahedron_vertex;
    }else if(name == "get_dodecahedron_vertex"){
        return get_dodecahedron_vertex;
    }else{
        assert(false);
        return nullptr;
    }
}

/**
 * @brief 
 * 判断视线与点法向量的夹角
 * @param op 
 * @param view_point 
 * @return std::vector<bool> 
 */
template<typename REAL, int DIM>
std::vector<bool> is_inside(const ORIENTED_POINTS& op, const Point<REAL, DIM>& view_point){
    std::vector<bool> res(op.size());
#pragma omp parallel for
    for(int i = 0;i<op.size();i++){
        Point<REAL, DIM> normal = op[i].second.normal;
        Point<REAL, DIM> view_dir = op[i].first - view_point;
        res[i] = Point<REAL, DIM>::Dot(normal, view_dir) > 0;
    }
    return res;
}


template<typename REAL, int DIM>
class PointOrientationConsistency {
private:
    VertexGetterP _func;
    double _radius;//hpr算法中的参数，越大越精确

public:
    PointOrientationConsistency(double radius = 0.01,VertexGetterP func = get_dodecahedron_vertex) : _func(func), _radius(radius) {
        
    }
    nlohmann::json get_config(){
        nlohmann::json j;
        j["name"] = "PointOrientationConsistency";
        j["radius"] = _radius;
        j["vertex_getter"] = get_vertex_getter_name<REAL,DIM>(_func);
        return j;
    }

    // 使用open3d的HiddenPointRemoval算法计算点云的可见性
    double cal_consistency(ORIENTED_POINTS& op, std::vector<double>& consistencies, std::string path = "") {
        auto view_points = _func();
        int max_vp = 0;
        consistencies.push_back(0);
#pragma omp parallel
        for (int i = 0; i < view_points.size(); i++){
            auto vp = view_points[i];
            std::vector<size_t> hop_idxs;
            Point<REAL, DIM> cam = { vp[0],vp[1],vp[2] };
            int t = HPR<REAL,DIM>(op,hop_idxs, cam, _radius);
            assert(t > 0);
            if (t < max_vp) {
                continue;
            }
            else {
                max_vp = t;
            }
            std::vector<bool> _is_inside = is_inside(op, cam);
            int inside = 0, outside = 0;
            for (int i = 0; i < hop_idxs.size(); i++) {
                int idx = hop_idxs[i];
                if (_is_inside[idx]) {
                    inside++;
                }
                else {
                    outside++;
                }
            }
            double consistency = (double)max(inside, outside) / (double)(inside + outside);
            consistencies[0] = consistency;
            if (path != "") {
                std::vector<int> flag(op.size(), NOT_SEEN);
                for (int i = 0; i < hop_idxs.size(); i++) {
                    int idx = hop_idxs[i];
                    if (_is_inside[idx]) {
                        flag[idx] = SEE_INNER;
                    }
                    else {
                        flag[idx] = SEE_OUTER;
                    }
                }
                lzd_tools::op2ply(op, path + "_time_" + std::to_string(i) + ".ply", XForm<REAL, DIM + 1>().Identity(), std::make_pair(flag, get_orient_color()));
            }
        }
        double sum_consistency = 0;
        for (auto& c : consistencies) {
            sum_consistency += c;
        }
        // 打印结果
        printf("\n");
        return sum_consistency / view_points.size();
    }
};

template<typename REAL, int DIM>
void DrawTriangleMesh(const MESH& mesh, const std::vector<int>& orient_flag,
    std::string& path, const Eigen::Vector3d viewPoints,const Eigen::Vector3d lookat
) {
    auto vp = Point<REAL, DIM>({ viewPoints[0],viewPoints[1],viewPoints[2] });
    auto lkp = Point<REAL, DIM>({ lookat[0] ,lookat[1] ,lookat[2] });
    auto tmesh1 = lzd_tools::add_sphere(mesh, vp, 0.01);
    auto tmesh2 = lzd_tools::add_sphere(tmesh1, lkp, 0.01);
    std::vector<int> t_flag(orient_flag);
    // 同时改变orient_flag
    for (int j = mesh.second.size(); j < tmesh1.second.size(); j++) {
        t_flag.push_back(VIEW_POINT);
    }
    for (int j = tmesh1.second.size(); j < tmesh2.second.size(); j++) {
        t_flag.push_back(LOOK_AT_POINT);
    }
    lzd_tools::mesh2ply(tmesh2, path, XForm<REAL, DIM + 1>().Identity(), std::make_pair(t_flag, get_orient_color()));
}

template<typename REAL, int DIM>
void DrawTriangleCluster(const MESH& mesh, const std::vector<int>& cluster_id,
    std::string& path, const Eigen::Vector3d viewPoints, const Eigen::Vector3d lookat
) {
    auto vp = Point<REAL, DIM>({ viewPoints[0],viewPoints[1],viewPoints[2] });
    auto lkp = Point<REAL, DIM>({ lookat[0] ,lookat[1] ,lookat[2] });
    auto tmesh1 = lzd_tools::add_sphere(mesh, vp, 0.01);
    auto tmesh2 = lzd_tools::add_sphere(tmesh1, lkp, 0.01);

    std::vector<int> color_id(cluster_id);
    for(int i = 0;i<color_id.size();i++){
        assert(color_id[i] >= 0);
        color_id[i] = color_id[i] % 100;
    }

    std::vector<int> t_flag(color_id);
    // 同时改变orient_flag
    for (int j = mesh.second.size(); j < tmesh1.second.size(); j++) {
        t_flag.push_back(100);
    }
    for (int j = tmesh1.second.size(); j < tmesh2.second.size(); j++) {
        t_flag.push_back(100);
    }
    lzd_tools::mesh2ply(tmesh2, path, XForm<REAL, DIM + 1>().Identity(), std::make_pair(t_flag,get_100_color()));
}


template<typename REAL, int DIM>
class MeshConsistency {
public:
    virtual double cal_consistency(const MESH& mesh, std::vector<double>& consistencies, std::string path = "", std::vector<int> tri_type = std::vector<int>()) = 0;
    virtual nlohmann::json get_config() = 0;
};

/**
 * @brief 
 * 计算多视角下可见面片的一致性
 * MaxAB 的正确性比OrientationConsistency更高，这个方法后续将只作为对比使用
 * @tparam REAL 
 * @tparam DIM 
 */
template<typename REAL, int DIM>
class OrientationConsistency:public MeshConsistency<REAL,DIM>{
private:
    VertexGetterP _func;
    unsigned int _width;//分辨率宽度
    unsigned int _height;//分辨率高度
    Rasterlization::PROJECTION_TYPE _projection;
public:
    // static std::string get_vertex_getter_name(VertexGetterP _func){
    //     if(_func == get_tetrahedron_vertex){
    //         return "get_tetrahedron_vertex";
    //     }else if(_func == get_dodecahedron_vertex){
    //         return "get_dodecahedron_vertex";
    //     }else{
    //         assert(false);
    //         return "unknown";
    //     }
    // }
    OrientationConsistency(
        VertexGetterP func = get_dodecahedron_vertex,
        unsigned int width = 1920,unsigned int height = 1080, Rasterlization::PROJECTION_TYPE projection = Rasterlization::PROJECTION_TYPE::ORTHOGRAPHIC
    ):_func(func),_width(width),_height(height),_projection(projection){
    }
    /**
     * @brief 
     * @param mesh 
     * @param path 例如./../a, 会打印出./../a_time_0.ply,./../a_time_1.ply...
     * @param tri_type 如果指定 在绘制时会将不同tri_type的三角形得到的面片分开
     * @return double 表示看到的所有的朝向一致的面片数量 / 视点个数
     */
    double cal_consistency(const MESH& mesh, std::vector<double>& consistencies,
        std::string path = "", std::vector<int> tri_type = std::vector<int>()
        ){
        //printf("start caculate orientation_consistency...\n");
        double sum_consistency = 0;
        std::vector<Rasterlization::Triangle> triangles;
        std::vector<Rasterlization::Vertex> vertices;
        Rasterlization::mesh2tv(mesh, triangles, vertices);
        std::vector<int> orient_flag(mesh.second.size());
        std::vector<Eigen::Vector3d> view_points;
        std::vector<Eigen::Vector3d> lookat;
        std::vector<Eigen::Vector3d> up;
        int times = get_viewpoints_lookat(mesh, view_points, lookat, up, _func);
        Rasterlization::Image* img = new Rasterlization::Image(_width, _height);
        // 此处不可并行
        for(int i=0;i<view_points.size();i++){
            auto consistency = Rasterlization::_orientation_consistence(img,triangles, vertices, view_points[i], lookat[i], up[i],orient_flag,_projection);
            consistencies.push_back(consistency);
            sum_consistency += consistency;
            if (path != "") {
                if (tri_type.size() == mesh.second.size())img->cal_boundary(tri_type);
                DrawTriangleMesh(mesh, orient_flag, path + "_time_" + std::to_string(i) + ".ply", view_points[i], lookat[i]);
                img->draw(path + "_time_" + std::to_string(i) + ".jpg");
            }
            img->reset();
        }
        delete img;
        assert(sum_consistency != 0);
        //printf("caculate orientation_consistency done!\n");
        return sum_consistency / times;
    }
    nlohmann::json get_config(){
        nlohmann::json j;
        j["name"] = "OrientationConsistency";
        j["width"] = _width;
        j["height"] = _height;
        j["vertex_getter"] = get_vertex_getter_name<REAL,DIM>(_func);
        return j;
    }
};


template<typename REAL, int DIM>
class MaxAB:public MeshConsistency<REAL,DIM>{
    VertexGetterP _func;
    unsigned int _width;//分辨率宽度
    unsigned int _height;//分辨率高度
    Rasterlization::PROJECTION_TYPE _projection;

public:

    MaxAB(
        VertexGetterP func = get_dodecahedron_vertex,
        unsigned int width=400,unsigned int height=400, Rasterlization::PROJECTION_TYPE projection=Rasterlization::PROJECTION_TYPE::ORTHOGRAPHIC
    ):_func(func),_width(width),_height(height),_projection(projection){
    }

    double cal_consistency(const MESH& mesh, std::vector<double>& consistencies,
        std::string path = "", std::vector<int> tri_type = std::vector<int>()
        ){
        assert(tri_type.size() == mesh.second.size());
        //printf("start caculate orientation_consistency...\n");
        double sum_consistency = 0;
        std::vector<Rasterlization::Triangle> triangles;
        std::vector<Rasterlization::Vertex> vertices;
        Rasterlization::mesh2tv(mesh, triangles, vertices);
        std::vector<int> orient_flag(mesh.second.size());

        std::vector<Eigen::Vector3d> view_points;
        std::vector<Eigen::Vector3d> lookat;
        std::vector<Eigen::Vector3d> up;
        int times = get_viewpoints_lookat(mesh, view_points, lookat, up, _func);

        Rasterlization::Image* img = new Rasterlization::Image(_width, _height);
        std::vector<long long> ABSizes(times,0);
        long long sum_AB_size = 0;

        // 此处不可并行
        for(int i=0;i<view_points.size();i++){
            // TODO 这里暂时只能用_orientaion_consistence而不是_rasterlization 因为_orientaion_consistence会计算像素朝向
            Rasterlization::_orientation_consistence(img,triangles, vertices, view_points[i], lookat[i], up[i],orient_flag,_projection);
            Rasterlization::HasEdge* has_edge = new Rasterlization::ContinousPixel(img->max_depth_gap);
            std::vector<std::vector<std::pair<int, int>>> clusters = Rasterlization::ImageBfsOnNotEmptyPixel(img, has_edge);
            for(int cid=0; cid <clusters.size();cid++){
                auto& cluster = clusters[cid];
                long long A = 0, B = 0;
                for (auto& pixel : cluster) {
                    img->cluster_id(pixel.first, pixel.second) = cid;               
                    int tri_id = img->triangle_id(pixel.first, pixel.second);
                    assert(tri_id >= 0 && tri_id < tri_type.size());
                    if (tri_type[tri_id] == 1)A++;
                    else if(tri_type[tri_id] == 2)B++;
                    else assert(false);   
                }
                ABSizes[i] += A*B;
            }

            if (path != "") {
                img->cal_boundary(tri_type);
                DrawTriangleMesh(mesh, orient_flag, path + "_time_" + std::to_string(i) + ".ply", view_points[i], lookat[i]);
                img->draw(path + "_time_" + std::to_string(i) + ".jpg");
                img->draw_cluster(path + "_time_" + std::to_string(i) + "_cluster.jpg");
                
            }
            img->reset();
            consistencies.push_back((double)ABSizes[i]);
            sum_AB_size += ABSizes[i];
        }
        delete img;
        assert(sum_AB_size >= 0);
        assert(static_cast<double>(sum_AB_size)>=0);
        return static_cast<double>(sum_AB_size) / times;
    }

    nlohmann::json get_config(){
        nlohmann::json j;
        j["name"] = "MaxAB";
        j["width"] = _width;
        j["height"] = _height;
        j["vertex_getter"] = get_vertex_getter_name<REAL,DIM>(_func);
        j["projection"] = _projection==Rasterlization::PROJECTION_TYPE::PERSPECTIVE?"PERSPECTIVE":"ORTHOGRAPHIC";
        return j;
    }
};


template<typename REAL, int DIM>
class MaxAB_Mesh{
    VertexGetterP _func;

    typedef struct TriangleList{
        const std::vector<Point<REAL, DIM>>& _vertices; 
           
        std::vector<std::vector<int>> _triangles; // 所有的三角形网格共用一个顶点列表,所以只需要记录三角形的索引
        TriangleList(const MESH& mesh):_vertices(mesh.first){
            for(auto& face:mesh.second){
                _triangles.push_back(face);
            }
        }
        TriangleList(const std::vector<Point<REAL, DIM>>& vertices):_vertices(vertices){
            _triangles = std::vector<std::vector<int>>(0);
        }

        TriangleList() = delete;


        //重载[]运算符
        std::vector<int>& operator[](int idx){
            return _triangles[idx];
        }
        const std::vector<int>& operator[](int idx) const{
            return _triangles[idx];
        }

        int num_vertex(){
            return _vertices.size();
        }
        int num_triangles(){
            return _triangles.size();
        }

        void add_triangle(const std::vector<int>& face){
            _triangles.push_back(face);
        }

        const std::vector<Point<REAL, DIM>>& get_vertex_list(){
            return _vertices;
        }

    }TriangleList;

    // typedef std::pair<triangles,std::vector<int>> ABMesh; // 一个mesh和其对应的tri_type
    typedef struct ABMesh{
        TriangleList _mesh;
        std::vector<int> _tri_type;
        ABMesh(const MESH& ori_mesh,std::vector<int>& tri_type):_mesh(ori_mesh),_tri_type(tri_type){
            assert(tri_type.size() == _mesh.num_triangles());
        }

        ABMesh() = delete; // 
        ABMesh(const std::vector<Point<REAL, DIM>>& vertices):_mesh(vertices){
            assert(_mesh.num_triangles() == 0);
            _tri_type = std::vector<int>(0);
        }
        
        int num_vertex(){
            return _mesh.num_vertex();
        }
        int num_triangles(){
            return _mesh.num_triangles();
        }
        void add_triangle(const std::vector<int>& face,int tri_type){
            _mesh.add_triangle(face);
            _tri_type.push_back(tri_type);
        }

        std::vector<int>& get_triangle(int idx){
            return _mesh[idx];
        }

        std::pair<std::vector<int>,int> operator[](int idx){
            return std::make_pair(_mesh[idx],_tri_type[idx]);
        }

        void add_triangle(std::pair<std::vector<int>,int> face){
            _mesh.add_triangle(face.first);
            _tri_type.push_back(face.second);
        }

        const std::vector<Point<REAL, DIM>>& get_vertex_list() {
            return _mesh.get_vertex_list();
        }

        long long cal_size(){
            long long A = 0,B = 0, boundary_point = 0;
            // 如果有顶点属于两个不同类型的三角形,则该顶点是边界点
            std::vector<std::set<int>> vertex_2_tri_type(_mesh.num_vertex());
            for(int i = 0;i<_tri_type.size();i++){
                auto& face = _mesh[i];
                for(auto& v:face){
                    if(vertex_2_tri_type[v].find(_tri_type[i]) != vertex_2_tri_type[v].end()){
                        continue;
                    }
                    if(vertex_2_tri_type[v].size() == 1){
                        boundary_point++;
                    }
                    vertex_2_tri_type[v].insert(_tri_type[i]);
                }
                if(_tri_type[i] == 1){
                    A++;
                }else if(_tri_type[i] == 2){
                    B++;
                }
            }
            return A*B*boundary_point;
        }

    }ABMesh;




    // 计算每个面片的orient_flag 1表示朝向内侧 0表示朝向外侧;
    void cal_orientation(const MESH& mesh,Eigen::Vector3d view_point,std::vector<int>& orient_flag){
        lzd_tools::AcumutableTimer timer("MaxAB_Mesh::cal_orientation");
        orient_flag.resize(mesh.second.size());
        // TODO 有空换成矩阵相乘,可能会快一点
        // Eigen::MatrixXd view_ray(mesh.second.size(),3);
        // Eigen::MatrixXd tri_normals(mesh.second.size(),3);
        // Eigen::MatrixXd ray_dot_normals(mesh.second.size(),1);
      
        auto p2eip = [](const Point<REAL, DIM>& p) {
            Eigen::Vector3d res;
            for (int i = 0; i < DIM; i++) {
                res[i] = p[i];
            }
            return res;
        };

        for(int i = 0;i<mesh.second.size();i++){
            auto& face = mesh.second[i];
            auto& v1 = p2eip(mesh.first[face[0]]);
            auto& v2 = p2eip(mesh.first[face[1]]);
            auto& v3 = p2eip(mesh.first[face[2]]);
            Eigen::Vector3d normal = (v2 - v1).cross(v3 - v1);
            Eigen::Vector3d ray = view_point - v1; // WARNING 只用第一个顶点的法向量 简化计算了 
            orient_flag[i] = ray.dot(normal) > 0 ? 1 : 0;
        }
    }

    

    // 从一个三角列表中得到联通分量.两个三角联通 <=> 两个三角存在公共边
    std::vector<ABMesh> _get_components_mesh(ABMesh& abmesh){
        lzd_tools::AcumutableTimer timer("MaxAB_Mesh::_get_components_mesh");   
        std::vector<ABMesh> res;
        std::vector<int> vertex_root(abmesh.num_vertex());
        std::vector<int> cluster_size(abmesh.num_vertex(), 1);
        for(int i = 0;i<vertex_root.size();i++){
            vertex_root[i] = i;
        }

        // 计算每个三角形所属的簇
        std::function<int(int)> find_root = [&](int idx){
            if(vertex_root[idx] == idx){
                return idx;
            }else{
                vertex_root[idx] = find_root(vertex_root[idx]);
                return find_root(vertex_root[idx]);
            }
        };
        auto union_root = [&](int idx1,int idx2){
            int root1 = find_root(idx1);
            int root2 = find_root(idx2);
            if(root1 != root2){
                vertex_root[root2] = root1;
                cluster_size[root1] += cluster_size[root2];
                cluster_size[root2] = 0;
            }
        };        
        for(int i = 0;i<abmesh.num_triangles();i++){
            auto& face = abmesh.get_triangle(i);
            union_root(face[0],face[1]);
            union_root(face[1],face[2]);
        }

        // 初始化res
        std::map<int,int> root_2_clusterid;
        for(int i = 0;i<vertex_root.size();i++){
            if (cluster_size[i] < 3) {
                if (!(cluster_size[i] == 1 || cluster_size[i] == 0))printf("%d\t", cluster_size[i]); // 这里报错说明mesh内有单个边
                continue;
            }
            if(root_2_clusterid.find(vertex_root[i]) != root_2_clusterid.end())continue;
            root_2_clusterid[vertex_root[i]] = res.size();
            res.push_back(ABMesh(abmesh.get_vertex_list()));          
        }

        // 分配三角形
        for(int i = 0;i<abmesh.num_triangles();i++){
            auto& face = abmesh.get_triangle(i);
            assert(find_root(face[0]) == find_root(face[1]));
            assert (find_root(face[1]) == find_root(face[2]));
            int cluster_id = root_2_clusterid[find_root(face[0])];
                res[cluster_id].add_triangle(abmesh[i]);
        }
        return res;          
    }

    // 根据flag将mesh分成两部分
    std::pair<ABMesh, ABMesh> _spilt_mesh(ABMesh& abmesh,const std::vector<int>& flag){
        lzd_tools::AcumutableTimer timer("MaxAB_Mesh::_spilt_mesh");
        ABMesh res0(abmesh.get_vertex_list()),res1(abmesh.get_vertex_list());
        for(int i = 0;i<abmesh.num_triangles();i++){
            if(flag[i]){
                res0.add_triangle(abmesh[i]);
            }else{
                res1.add_triangle(abmesh[i]);
            }
        }
        return std::make_pair(res0,res1);
    }


public:
    MaxAB_Mesh(VertexGetterP func = get_dodecahedron_vertex):_func(func){
    }
    double cal_consistency(const MESH& mesh, std::vector<double>& consistencies,
        std::string path = "", std::vector<int> tri_type = std::vector<int>()
    ){
        assert(tri_type.size() == mesh.second.size());
        //printf("start caculate orientation_consistency...\n");
        double sum_consistency = 0;
        std::vector<Eigen::Vector3d> view_points;
        std::vector<Eigen::Vector3d> lookat;
        std::vector<Eigen::Vector3d> up;
        int times = get_viewpoints_lookat(mesh, view_points, lookat, up, _func);
        std::vector<long long> ABSizes(times,0);
        long long sum_AB_size = 0;

        for(int i = 0;i<times;i++){
            std::vector<int> orient_flag;
            long long absizei = 0;
            cal_orientation(mesh,view_points[i],orient_flag);
            ABMesh abmesh(mesh, tri_type);
            std::pair<ABMesh, ABMesh> inside_and_outside = _spilt_mesh(abmesh,orient_flag);
            std::vector<ABMesh> components = _get_components_mesh(inside_and_outside.first);
            auto components2 = _get_components_mesh(inside_and_outside.second);
            for(auto &component:components)absizei += component.cal_size();
            for (auto& component : components2)absizei += component.cal_size();


            // 如果不存在大于零的ABMesh,可以考虑早停
            // if(absizei == 0){
            //     break;
            // }
            sum_AB_size += absizei;
            consistencies.push_back(absizei);
            if (path != "") {
                MESH res_mesh(mesh.first, std::vector<std::vector<int>>(0));
                std::vector<int> cluster_id;
                int cid = 0;
                for(auto& component:components){
                    cid++;
                    if (component.cal_size() == 0)continue;
                    for(int j = 0;j<component.num_triangles();j++){
                        res_mesh.second.push_back(component[j].first);
                        cluster_id.push_back(cid);
                    }
                }
                for (auto& component : components2) {
                    cid++;
                    if (component.cal_size() == 0)continue;
                    for (int j = 0; j < component.num_triangles(); j++) {
                        res_mesh.second.push_back(component[j].first);
                        cluster_id.push_back(cid);
                    }
                }
                DrawTriangleCluster(res_mesh, cluster_id, path + "_time_cluster" + std::to_string(i) + ".ply", view_points[i], lookat[i]);
                DrawTriangleCluster(mesh, tri_type, path + "_time_" + std::to_string(i) + ".ply", view_points[i], lookat[i]);
            }
        
        }
        
        assert(sum_AB_size >= 0);
        assert(static_cast<double>(sum_AB_size)>=0);
        return static_cast<double>(sum_AB_size) / times;

    }
    
    nlohmann::json get_config(){
        nlohmann::json j;
        j["name"] = "MaxAB_Mesh";
        j["vertex_getter"] = get_vertex_getter_name<REAL,DIM>(_func);
        return j;
    }

};