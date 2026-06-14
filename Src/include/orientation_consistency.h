#pragma once
#include "o3d_api.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <set>

enum OrientationVisualFlag {
    NOT_SEEN = -1,
    SEE_INNER = 1,
    SEE_OUTER = 2,
    VIEW_POINT = 3,
    LOOK_AT_POINT = 4,
    PARALLEL = 5
};


namespace Rasterlization{
   /**
    * Rasterizes mesh triangles into an image-space z-buffer and marks visible
    * triangles from the current viewpoint.
    */
   struct Triangle {
      int v1, v2, v3;
   };

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
      Eigen::MatrixXi is_inside;
      Eigen::MatrixXi triangle_id;
      Eigen::MatrixXi cluster_id;
      Eigen::MatrixXi is_boundary;
      
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
       * boundary = len(unique(neighbor(tri_type[triangle_id[y,x]]))) == 2
       * @param tri_type 
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
    */
   class ContinousPixel : public HasEdge {
   private:
       double _eps;
   public:
       ContinousPixel(double eps) :_eps(eps) {}
       bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b, const Image* img) {
            return abs(img->zBuffer(a.first, a.second) - img->zBuffer(b.first, b.second)) < _eps && img->is_inside(a.first, a.second) == img->is_inside(b.first, b.second);
       }
   };

   std::vector<std::vector<std::pair<int,int>>> ImageBfsOnNotEmptyPixel(Image* img, HasEdge* has_edge);

   Eigen::Vector2d viewportTransform(const Eigen::Vector2d& cxy, int imageWidth, int imageHeight);

   bool isInsideTriangle(const Eigen::Vector2d& pxyP, const Eigen::Vector2d& pxyA, const Eigen::Vector2d& pxyB, const Eigen::Vector2d& pxyC);

   void updatePixelDepthBuffer(Eigen::MatrixXd& zBuffer, Eigen::MatrixXi& pixelSource, const Eigen::Vector2d& pxyA, const Eigen::Vector2d& pxyB, const Eigen::Vector2d& pxyC, double d, int i);

   Eigen::Matrix4d calculateViewMatrix(const Eigen::Vector3d& viewPoints, const Eigen::Vector3d& lookat, const Eigen::Vector3d& up);

   Eigen::Matrix4d getOrthoMatrix(double lx, double rx, double by, double ty, double n, double f);
   
   /**
    * @brief 
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
    * @param mesh surface 
    */
   template<typename REAL, unsigned int DIM>
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
       auto res =  _computeVisibility(triangles, vertices, viewPoints_eigen, lookat_eigen, up_eigen, imageWidth, imageHeight, Rasterlization::ORTHOGRAPHIC);
       printf("calculate visibility done!\n");
       return res;
   }


   /**
    * @brief 
    * @param mesh 
    * @return std::vector<int> 
    */
   template<typename REAL, unsigned int DIM>
   std::vector<int> random_view(const MESH& mesh,int image_width = 1920, int image_height = 1080) {
       int randidx = rand() % mesh.first.size();
       auto viewPoints = mesh.first[randidx];
       randidx = rand() % mesh.first.size();
       auto lookat = mesh.first[randidx];
       Point<float, 3> up(0.0, 1.0, 0.0);
       return computeVisibility(mesh, viewPoints, lookat, up, image_width, image_height);   
   }


   double  _orientation_consistence(Image *img,const std::vector<Triangle>& triangles, const std::vector<Vertex>& vertices, const Eigen::Vector3d& viewPoints,
                           const Eigen::Vector3d& lookat, const Eigen::Vector3d& up, std::vector<int>& orient_flag,PROJECTION_TYPE projection);


   /**
    * @brief 
    */
   template<typename REAL, unsigned int DIM>
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
 */
template<typename REAL,unsigned int DIM>
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

std::vector<Eigen::Vector3d> get_tetrahedron_vertex();

std::vector<Eigen::Vector3d> get_dodecahedron_vertex();

/**
 * @brief Get the viewpoints lookat object
 * @return int 
 */
template<typename REAL,unsigned int DIM>
int get_viewpoints_lookat(const MESH& mesh,
    std::vector<Eigen::Vector3d>& view_points, std::vector<Eigen::Vector3d>& lookat,std::vector<Eigen::Vector3d>& up,
    VertexGetterP vertex_getter = get_tetrahedron_vertex) 
{
    Point<REAL, DIM> mesh_center;
    for (int i = 0; i < DIM; i++) {
        mesh_center[i] = 0;
    }
    for (auto& vertex : mesh.first) {
        mesh_center += vertex;
    }
    mesh_center /= mesh.first.size();
    
    REAL mesh_radius = 0;
    for (auto& vertex : mesh.first) {
        mesh_radius = std::max(mesh_radius, Point<REAL, DIM>::SquareNorm((vertex - mesh_center)));
    }
    mesh_radius = std::sqrt(mesh_radius);
    Eigen::Vector3d emesh_center({ mesh_center[0],mesh_center[1],mesh_center[2] });

    auto vertices = vertex_getter();
    for (auto& vertex : vertices) {
        vertex = vertex * mesh_radius + emesh_center;
    }
    view_points.clear();
    lookat.clear();
    up.clear();
    for (auto& vertex : vertices) {
        view_points.push_back(vertex);
        lookat.push_back(emesh_center);
        up.push_back(Eigen::Vector3d(0, 0, 1));

    }
    return view_points.size();
}


lzd_tools::FiniteMap get_orient_color();
lzd_tools::FiniteMap get_100_color();


template<typename REAL, unsigned int DIM>
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

template<typename REAL,unsigned int DIM>
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
 * @param op 
 * @param view_point 
 * @return std::vector<bool> 
 */
template<typename REAL, unsigned int DIM>
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


template<typename REAL, unsigned int DIM>
class PointOrientationConsistency {
private:
    VertexGetterP _func;
    double _radius;

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
            double consistency = (double)std::max(inside, outside) / (double)(inside + outside);
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
        printf("\n");
        return sum_consistency / view_points.size();
    }
};

template<typename REAL, unsigned int DIM>
void DrawTriangleMesh(const MESH& mesh, const std::vector<int>& orient_flag,
    const std::string& path, const Eigen::Vector3d viewPoints,const Eigen::Vector3d lookat
    ) {
    auto vp = Point<REAL, DIM>({ viewPoints[0],viewPoints[1],viewPoints[2] });
    auto lkp = Point<REAL, DIM>({ lookat[0] ,lookat[1] ,lookat[2] });
    auto tmesh1 = lzd_tools::add_sphere(mesh, vp, 0.01);
    auto tmesh2 = lzd_tools::add_sphere(tmesh1, lkp, 0.01);
    std::vector<int> t_flag(orient_flag);
    for (int j = mesh.second.size(); j < tmesh1.second.size(); j++) {
        t_flag.push_back(VIEW_POINT);
    }
    for (int j = tmesh1.second.size(); j < tmesh2.second.size(); j++) {
        t_flag.push_back(LOOK_AT_POINT);
    }
    lzd_tools::mesh2ply(tmesh2, path, XForm<REAL, DIM + 1>().Identity(), std::make_pair(t_flag, get_orient_color()));
}

template<typename REAL, unsigned int DIM>
void DrawTriangleCluster(const MESH& mesh, const std::vector<int>& cluster_id,
    const std::string& path, const Eigen::Vector3d viewPoints, const Eigen::Vector3d lookat
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
    for (int j = mesh.second.size(); j < tmesh1.second.size(); j++) {
        t_flag.push_back(100);
    }
    for (int j = tmesh1.second.size(); j < tmesh2.second.size(); j++) {
        t_flag.push_back(100);
    }
    lzd_tools::mesh2ply(tmesh2, path, XForm<REAL, DIM + 1>().Identity(), std::make_pair(t_flag,get_100_color()));
}


template<typename REAL, unsigned int DIM>
class MeshConsistency {
public:
    virtual double cal_consistency(const MESH& mesh, std::vector<double>& consistencies, std::string path = "", std::vector<int> tri_type = std::vector<int>()) = 0;
    virtual nlohmann::json get_config() = 0;
};

/**
 * @brief 
 * Computes visible-triangle consistency from multiple viewpoints.
 * MaxAB is more accurate; keep this implementation for comparison.
 * @tparam REAL 
 * @tparam DIM 
 */
template<typename REAL, unsigned int DIM>
class OrientationConsistency:public MeshConsistency<REAL,DIM>{
private:
    VertexGetterP _func;
    unsigned int _width;
    unsigned int _height;
    Rasterlization::PROJECTION_TYPE _projection;
public:
    // static std::string get_vertex_getter_name(VertexGetterP _func){
    //     }else if(_func == get_dodecahedron_vertex){
    //     }else{
    OrientationConsistency(
        VertexGetterP func = get_dodecahedron_vertex,
        unsigned int width = 1920,unsigned int height = 1080, Rasterlization::PROJECTION_TYPE projection = Rasterlization::PROJECTION_TYPE::ORTHOGRAPHIC
    ):_func(func),_width(width),_height(height),_projection(projection){
    }
    /**
     * @brief 
     * @param mesh 
     * @param path Output prefix, e.g. `../a` writes `../a_time_0.ply`.
     * @param tri_type Optional triangle classes used for visualization.
     * @return Average visible consistent triangle count per viewpoint.
     */
    double cal_consistency(const MESH& mesh, std::vector<double>& consistencies,
        std::string path = "", std::vector<int> tri_type = std::vector<int>()
        ){
        double sum_consistency = 0;
        std::vector<Rasterlization::Triangle> triangles;
        std::vector<Rasterlization::Vertex> vertices;
        Rasterlization::mesh2tv<REAL, DIM>(mesh, triangles, vertices);
        std::vector<int> orient_flag(mesh.second.size());
        std::vector<Eigen::Vector3d> view_points;
        std::vector<Eigen::Vector3d> lookat;
        std::vector<Eigen::Vector3d> up;
        int times = get_viewpoints_lookat<REAL, DIM>(mesh, view_points, lookat, up, _func);
        Rasterlization::Image* img = new Rasterlization::Image(_width, _height);
        // Image state is reused in this loop.
        for(int i=0;i<view_points.size();i++){
            auto consistency = Rasterlization::_orientation_consistence(img,triangles, vertices, view_points[i], lookat[i], up[i],orient_flag,_projection);
            consistencies.push_back(consistency);
            sum_consistency += consistency;
            if (path != "") {
                if (tri_type.size() == mesh.second.size())img->cal_boundary(tri_type);
                DrawTriangleMesh<REAL, DIM>(mesh, orient_flag, path + "_time_" + std::to_string(i) + ".ply", view_points[i], lookat[i]);
                img->draw(path + "_time_" + std::to_string(i) + ".jpg");
            }
            img->reset();
        }
        delete img;
        assert(sum_consistency != 0);
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


template<typename REAL, unsigned int DIM>
class MaxAB:public MeshConsistency<REAL,DIM>{
    VertexGetterP _func;
    unsigned int _width;
    unsigned int _height;
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
        double sum_consistency = 0;
        std::vector<Rasterlization::Triangle> triangles;
        std::vector<Rasterlization::Vertex> vertices;
        Rasterlization::mesh2tv<REAL, DIM>(mesh, triangles, vertices);
        std::vector<int> orient_flag(mesh.second.size());

        std::vector<Eigen::Vector3d> view_points;
        std::vector<Eigen::Vector3d> lookat;
        std::vector<Eigen::Vector3d> up;
        int times = get_viewpoints_lookat<REAL, DIM>(mesh, view_points, lookat, up, _func);

        Rasterlization::Image* img = new Rasterlization::Image(_width, _height);
        std::vector<long long> ABSizes(times,0);
        long long sum_AB_size = 0;

        // Image state is reused in this loop.
        for(int i=0;i<view_points.size();i++){
            // _orientation_consistence computes per-pixel orientation.
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
                DrawTriangleMesh<REAL, DIM>(mesh, orient_flag, path + "_time_" + std::to_string(i) + ".ply", view_points[i], lookat[i]);
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


template<typename REAL, unsigned int DIM>
class MaxAB_Mesh{
    VertexGetterP _func;

    typedef struct TriangleList{
        const std::vector<Point<REAL, DIM>>& _vertices; 
           
        std::vector<std::vector<int>> _triangles;
        TriangleList(const MESH& mesh):_vertices(mesh.first){
            for(auto& face:mesh.second){
                _triangles.push_back(face);
            }
        }
        TriangleList(const std::vector<Point<REAL, DIM>>& vertices):_vertices(vertices){
            _triangles = std::vector<std::vector<int>>(0);
        }

        TriangleList() = delete;


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




    void cal_orientation(const MESH& mesh,Eigen::Vector3d view_point,std::vector<int>& orient_flag){
        lzd_tools::AcumutableTimer timer("MaxAB_Mesh::cal_orientation");
        orient_flag.resize(mesh.second.size());
      
        auto p2eip = [](const Point<REAL, DIM>& p) {
            Eigen::Vector3d res;
            for (int i = 0; i < DIM; i++) {
                res[i] = p[i];
            }
            return res;
        };

        for(int i = 0;i<mesh.second.size();i++){
            auto& face = mesh.second[i];
            Eigen::Vector3d v1 = p2eip(mesh.first[face[0]]);
            Eigen::Vector3d v2 = p2eip(mesh.first[face[1]]);
            Eigen::Vector3d v3 = p2eip(mesh.first[face[2]]);
            Eigen::Vector3d normal = (v2 - v1).cross(v3 - v1);
            Eigen::Vector3d ray = view_point - v1;
            orient_flag[i] = ray.dot(normal) > 0 ? 1 : 0;
        }
    }

    

    std::vector<ABMesh> _get_components_mesh(ABMesh& abmesh){
        lzd_tools::AcumutableTimer timer("MaxAB_Mesh::_get_components_mesh");   
        std::vector<ABMesh> res;
        std::vector<int> vertex_root(abmesh.num_vertex());
        std::vector<int> cluster_size(abmesh.num_vertex(), 1);
        for(int i = 0;i<vertex_root.size();i++){
            vertex_root[i] = i;
        }

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

        std::map<int,int> root_2_clusterid;
        for(int i = 0;i<vertex_root.size();i++){
            if (cluster_size[i] < 3) {
                if (!(cluster_size[i] == 1 || cluster_size[i] == 0))printf("%d\t", cluster_size[i]);
                continue;
            }
            if(root_2_clusterid.find(vertex_root[i]) != root_2_clusterid.end())continue;
            root_2_clusterid[vertex_root[i]] = res.size();
            res.push_back(ABMesh(abmesh.get_vertex_list()));          
        }

        for(int i = 0;i<abmesh.num_triangles();i++){
            auto& face = abmesh.get_triangle(i);
            assert(find_root(face[0]) == find_root(face[1]));
            assert (find_root(face[1]) == find_root(face[2]));
            int cluster_id = root_2_clusterid[find_root(face[0])];
                res[cluster_id].add_triangle(abmesh[i]);
        }
        return res;          
    }

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
        double sum_consistency = 0;
        std::vector<Eigen::Vector3d> view_points;
        std::vector<Eigen::Vector3d> lookat;
        std::vector<Eigen::Vector3d> up;
        int times = get_viewpoints_lookat<REAL, DIM>(mesh, view_points, lookat, up, _func);
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


            // A zero AB mesh can be used as an early-stop signal.
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
                DrawTriangleCluster<REAL, DIM>(res_mesh, cluster_id, path + "_time_cluster" + std::to_string(i) + ".ply", view_points[i], lookat[i]);
                DrawTriangleCluster<REAL, DIM>(mesh, tri_type, path + "_time_" + std::to_string(i) + ".ply", view_points[i], lookat[i]);
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
