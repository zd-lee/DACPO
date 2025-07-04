#include "orientation_consistency.h"
//#include <png_pp/png.hpp>
#include <queue>

std::vector<std::vector<std::pair<int, int>>> Rasterlization::ImageBfsOnNotEmptyPixel(Rasterlization::Image* img, Rasterlization::HasEdge* has_edge)
{
    std::vector<std::vector<std::pair<int, int>>> res;
    std::vector<std::vector<bool>> visited(img->getHeight(), std::vector<bool>(img->getWidth(), false));

    for (int i = 0; i < img->getHeight(); i++)
    {
        for (int j = 0; j < img->getWidth(); j++)
        {
            if (img->is_empty(i, j) || visited[i][j])continue;
            std::vector<std::pair<int, int>> cluster;
            std::queue<std::pair<int, int>> q;
            q.push(std::make_pair(i, j));
            visited[i][j] = true;
            while (!q.empty())
            {
                auto p = q.front();
                q.pop();
                cluster.push_back(p);
                auto neibors = img->get_not_empty_neibor_pixel_idx(p.first, p.second);
                for (auto &neibor : neibors)
                {
                    if (!visited[neibor.first][neibor.second] && has_edge->operator()(std::make_pair(p.first, p.second), neibor, img))
                    {
                        q.push(neibor);
                        visited[neibor.first][neibor.second] = true;
                    }
                }
            }
            res.push_back(cluster);
        }
    }
    return res;
}

Eigen::Vector2d Rasterlization::viewportTransform(const Eigen::Vector2d &cxy, int imageWidth, int imageHeight)
{
    Eigen::Vector2d pxy;
    pxy.x() = (cxy.x() + 1) * (imageWidth - 1) / 2;
    pxy.y() = (cxy.y() + 1) * (imageHeight - 1) / 2;
    return pxy;
}

bool Rasterlization::isInsideTriangle(const Eigen::Vector2d &pxyP, const Eigen::Vector2d &pxyA, const Eigen::Vector2d &pxyB, const Eigen::Vector2d &pxyC)
{
    //// 以A为原点,AB、AC为坐标轴，计算P的坐标(u,v)
    //double u = (pxyP - pxyA).dot(pxyB - pxyA) / (pxyB - pxyA).dot(pxyB - pxyA);
    //double v = (pxyP - pxyA).dot(pxyC - pxyA) / (pxyC - pxyA).dot(pxyC - pxyA);
    //// 判断P是否在三角形内部
    //return u >= 0 && v >= 0 && u + v <= 1;
    return true;

}

inline double minin3(double a, double b, double c)
{
    return min(min(a, b), c);
}

void Rasterlization::updatePixelDepthBuffer(Eigen::MatrixXd &zBuffer, Eigen::MatrixXi &pixelSource, const Eigen::Vector2d &pxyA, const Eigen::Vector2d &pxyB, const Eigen::Vector2d &pxyC, double d, int i)
{
    // 计算三角形的边界
    int minX = std::floor(minin3(pxyA.x(), pxyB.x(), pxyC.x()));
    int minY = std::floor(minin3(pxyA.y(), pxyB.y(), pxyC.y()));
    int maxX = std::ceil(max(max(pxyA.x(), pxyB.x()), pxyC.x()));
    int maxY = std::ceil(max(max(pxyA.y(), pxyB.y()), pxyC.y()));

    // 剔除边界外的像素
    if (minX < 0 || minY < 0 || maxX >= zBuffer.cols() || maxY >= zBuffer.rows())
    {
        return;
    }
    // 遍历边界内的像素
    for (int y = minY; y <= maxY; ++y)
    {
        for (int x = minX; x <= maxX; ++x)
        {
            // 判断像素是否在三角形内部
            Eigen::Vector2d pixel(x + 0.5, y + 0.5);
            if (isInsideTriangle(pixel, pxyA, pxyB, pxyC))
            {
                // 计算像素的深度d
                // 如果d小于zBuffer[y][x],则更新zBuffer[y][x] = d，并更新pixelSource[y][x] = i
                if (d < zBuffer(y, x))
                {
                    zBuffer(y, x) = d;
                    pixelSource(y, x) = i;
                }
            }
        }
    }
}

Eigen::Matrix4d Rasterlization::calculateViewMatrix(const Eigen::Vector3d &eye, const Eigen::Vector3d &center, const Eigen::Vector3d &up)
{
    Eigen::Vector3d z = (eye - center).normalized();
    Eigen::Vector3d x = up.cross(z).normalized();
    
    if(x.norm() == 0){
        auto nup = Eigen::Vector3d(up.x()+0.1300123,up.y()+0.0004211,up.z()+0.321098);
        x = nup.cross(z).normalized();
    }
    if (x.norm() == 0) {
        auto nup = Eigen::Vector3d(up.x() + 0.0001323, up.y() + 0.04211, up.z() + 0.0321098);
        x = nup.cross(z).normalized();
    }


    Eigen::Vector3d y = z.cross(x).normalized();
    Eigen::Matrix4d viewMatrix;
    viewMatrix << x.x(), x.y(), x.z(), -x.dot(eye),
        y.x(), y.y(), y.z(), -y.dot(eye),
        z.x(), z.y(), z.z(), -z.dot(eye),
        0, 0, 0, 1;
    return viewMatrix;
}

Eigen::Matrix4d Rasterlization::getOrthoMatrix(double left, double right, double bottom, double top, double _near, double _far)
{
    Eigen::Matrix4d orthoMatrix;
    orthoMatrix << 2 / (right - left), 0, 0, -(right + left) / (right - left),
        0, 2 / (top - bottom), 0, -(top + bottom) / (top - bottom),
        0, 0, -2 / (_far - _near), -(_far + _near) / (_far - _near),
        0, 0, 0, 1;
    return orthoMatrix;
}

/**
 * @brief Get the Perspective Matrix object
 * @param left 近裁剪面的左边界
 * @param right 近裁剪面的右边界
 * @param bottom 近裁剪面的下边界
 * @param top 近裁剪面的上边界
 * @param _near 近裁剪面的距离
 * @param _far 远裁剪面的距离
 * @return Eigen::Matrix4d 
 */
static Eigen::Matrix4d getPerspectiveMatrix(double left, double right, double bottom, double top, double _near, double _far)
{
    assert(_near * _far >= 0);
    Eigen::Matrix4d perspectiveMatrix;
    perspectiveMatrix << 2 * _near / (right - left), 0, (right + left) / (right - left), 0,
        0, 2 * _near / (top - bottom), (top + bottom) / (top - bottom), 0,
        0, 0, -(_far + _near) / (_far - _near), -2 * _far * _near / (_far - _near),
        0, 0, -1, 0;
    return perspectiveMatrix;
}


// 不生成图片,只计算三角网格的可见性 (TODO 其实还是计算了图片,只是没有保存  可以考虑优化)
std::vector<int> Rasterlization::_computeVisibility(
    const std::vector<Rasterlization::Triangle> &triangles, const std::vector<Rasterlization::Vertex> &vertices, 
    const Eigen::Vector3d &viewPoints, const Eigen::Vector3d &lookat, const Eigen::Vector3d &up, 
    int imageWidth, int imageHeight, Rasterlization::PROJECTION_TYPE projection)
{
    assert(triangles.size() > 0);
    std::vector<int> visibilityList(triangles.size(), 0);
    std::vector<double> triangles_depths(triangles.size(), 0.0f);
    Rasterlization::Image img(imageWidth, imageHeight);
    Rasterlization::_rasterlization(&img,triangles, vertices, viewPoints, lookat, up, visibilityList, triangles_depths,projection);
    return visibilityList;
}

void Rasterlization::_rasterlization(
    Rasterlization::Image* img,
    // 数据
    const std::vector<Rasterlization::Triangle> &triangles, const std::vector<Rasterlization::Vertex> &vertices, 
    // 视角相关
    const Eigen::Vector3d &viewPoints, const Eigen::Vector3d &lookat, const Eigen::Vector3d &up,
    // output 供三角网格使用
    std::vector<int> &visibilityList, std::vector<double> &triangles_depths,
    Rasterlization::PROJECTION_TYPE projection
)
{ 
    assert(triangles.size() > 0);    
    // 计算视角矩阵，将vertex投影到view space
    Eigen::Matrix4d viewMatrix = calculateViewMatrix(viewPoints, lookat, up);
    Eigen::MatrixX4d xyzd(vertices.size(), 4);
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        xyzd.row(i) << vertices[i].x, vertices[i].y, vertices[i].z, 1.0f;
    }
    Eigen::MatrixX4d vxyzd = (viewMatrix * xyzd.transpose()).transpose();   
 
    // 计算正交投影矩阵 (TODO 由于bbox操作麻烦 暂时先不优化)
    double rangeScale = 1.1f;
    Eigen::Matrix4d (*getProjectionMatrix)(double, double, double, double, double, double) = nullptr;
    if(projection == 0){
        getProjectionMatrix = getOrthoMatrix;
    }else{
        getProjectionMatrix = getPerspectiveMatrix;
    }

    Eigen::MatrixX4d project = getProjectionMatrix(
        vxyzd.col(0).minCoeff() * rangeScale, vxyzd.col(0).maxCoeff() * rangeScale,
        vxyzd.col(1).minCoeff() * rangeScale, vxyzd.col(1).maxCoeff() * rangeScale,
        vxyzd.col(2).minCoeff(), vxyzd.col(2).maxCoeff()
    );
    // 将vertex投影到裁剪空间
    Eigen::MatrixX4d cxyzd = (project * vxyzd.transpose()).transpose();
    Eigen::MatrixX2d cxy = cxyzd.leftCols(2);
    
    // 遍历三角形，更新zBuffer, triangle_id
    double max_depth_gap = 0; // 单个三角形内最大深度差
    for (size_t i = 0; i < triangles.size(); ++i)
    {
        Eigen::Vector2d pxyA = viewportTransform(cxy.row(triangles[i].v1).transpose(), img->width, img->height);
        Eigen::Vector2d pxyB = viewportTransform(cxy.row(triangles[i].v2).transpose(), img->width, img->height);
        Eigen::Vector2d pxyC = viewportTransform(cxy.row(triangles[i].v3).transpose(), img->width, img->height);
        double zA = cxyzd.row(triangles[i].v1).z(), zB = cxyzd.row(triangles[i].v2).z(), zC = cxyzd.row(triangles[i].v3).z();
        triangles_depths[i] = (zA + zB + zC) / 3.0f;
        max_depth_gap = max(max_depth_gap, abs(zA - zB));
        max_depth_gap = max(max_depth_gap, abs(zA - zC));
        max_depth_gap = max(max_depth_gap, abs(zB - zC));
        double d = triangles_depths[i];
        updatePixelDepthBuffer(img->zBuffer, img->triangle_id, pxyA, pxyB, pxyC, d, i);
    }

    img->cal_boundary(max_depth_gap);
    img->max_depth_gap = max_depth_gap;

    // 遍历pixelSource，根据pixelSource更新visibilityList
    for (int y = 0; y < img->height; ++y)
    {
        for (int x = 0; x < img->width; ++x)
        {
            int i = img->triangle_id(y, x);
            if (i != -1)
            {
                visibilityList[i] = 1;
            }
        }
    }
}

Eigen::Vector3d get_face_normal(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c) {
    Eigen::Vector3d ab = b - a;
    Eigen::Vector3d ac = c - a;
    Eigen::Vector3d n = ab.cross(ac);
    return n;
}

/**
*@brief 获得可见面片的一致性 consistency = max(内侧,外侧)
*/
double Rasterlization::_orientation_consistence(
    Rasterlization::Image *img,
    const std::vector<Triangle> &triangles, const std::vector<Vertex> &vertices, 
    const Eigen::Vector3d &viewPoints, const Eigen::Vector3d &lookat, const Eigen::Vector3d &up, 
    std::vector<int>& orient_flag,Rasterlization::PROJECTION_TYPE projection
    )
{
    // TODO 如果重复调用,可以将visibilityList缓存起来
    auto triangles_depths = std::vector<double>(triangles.size(), 0.0f);
    auto visibilityList = std::vector<int>(triangles.size(), 0);
    Rasterlization::_rasterlization(img, triangles, vertices, viewPoints, lookat, up,visibilityList, triangles_depths, projection);

    int inner_count = 0;
    int outer_count = 0;
    assert(triangles.size() == orient_flag.size());
    auto viewPoints_eigen = Eigen::Vector3d(viewPoints[0], viewPoints[1], viewPoints[2]);
    for (size_t i = 0; i < triangles.size(); ++i)
    {
        if (visibilityList[i] > 0)
        {
            auto a = Eigen::Vector3d(vertices[triangles[i].v1].x, vertices[triangles[i].v1].y, vertices[triangles[i].v1].z);
            auto b = Eigen::Vector3d(vertices[triangles[i].v2].x, vertices[triangles[i].v2].y, vertices[triangles[i].v2].z);
            auto c = Eigen::Vector3d(vertices[triangles[i].v3].x, vertices[triangles[i].v3].y, vertices[triangles[i].v3].z);            
            auto n = get_face_normal(a, b, c);
            Eigen::Vector3d view2face;
            if (projection == Rasterlization::PROJECTION_TYPE::PERSPECTIVE)
            {
                view2face = (a + b + c) / 3 - viewPoints_eigen;
            }
            else if (projection == Rasterlization::PROJECTION_TYPE::ORTHOGRAPHIC)
            {
                view2face = lookat - viewPoints_eigen;
            }else{
                assert(false);
            }
            auto val = n.dot(view2face);
            if(val < 0){
                inner_count++;
                orient_flag[i] = SEE_OUTER;
            }else if(val > 0){
                outer_count++;
                orient_flag[i] = SEE_INNER;
            }else{
                orient_flag[i] = PARALLEL;
            }
        }
        else {
            orient_flag[i] = NOT_SEEN;
        }
    }

    int NO_orientaion = 0;
    for (int i = 0;i<img->width;i++){
        for (int j = 0;j<img->height;j++){
            if(img->triangle_id(j,i) == -1){
                continue;
            }
            if(orient_flag[img->triangle_id(j,i)] == SEE_INNER){
                img->is_inside(j,i) = 1;
            }else if (orient_flag[img->triangle_id(j,i)] == SEE_OUTER){
                img->is_inside(j,i) = 0;
            }else{
                NO_orientaion++;
            }
        }
    }
  //  if (NO_orientaion) {
   //     printf("%d pixel has undefined orientation\n", NO_orientaion);
 //   }
    // return (double)max(inner_count,outer_count)/ (double)(inner_count+outer_count);// 比例会导致一些看得面片比较少的视角得分过高
    return (double)max(inner_count,outer_count);
}

/**
 * @brief Get the tetrahedron vertex object
 * @return std::vector<Eigen::Vector3d> 
 */
std::vector<Eigen::Vector3d> get_tetrahedron_vertex()
{
    std::vector<Eigen::Vector3d> vertices(4);
    vertices[0] = Eigen::Vector3d(1, 1, 1);
    vertices[1] = Eigen::Vector3d(1, -1, -1);
    vertices[2] = Eigen::Vector3d(-1, 1, -1);
    vertices[3] = Eigen::Vector3d(-1, -1, 1);
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        vertices[i] = vertices[i] / vertices[i].norm();
    }
    return vertices;
}

std::vector<Eigen::Vector3d> get_dodecahedron_vertex()
{
    std::vector<Eigen::Vector3d> vertices(12);
    double phi = (1 + sqrt(5)) / 2;
    vertices[0] = Eigen::Vector3d(0, 0, phi);
    vertices[1] = Eigen::Vector3d(0, 0, -phi);
    vertices[2] = Eigen::Vector3d(0, phi, 0);
    vertices[3] = Eigen::Vector3d(0, -phi, 0);
    vertices[4] = Eigen::Vector3d(phi, 0, 0);
    vertices[5] = Eigen::Vector3d(-phi, 0, 0);
    vertices[6] = Eigen::Vector3d(1, 1, 1);
    vertices[7] = Eigen::Vector3d(1, 1, -1);
    vertices[8] = Eigen::Vector3d(1, -1, 1);
    vertices[9] = Eigen::Vector3d(1, -1, -1);
    vertices[10] = Eigen::Vector3d(-1, 1, 1);
    vertices[11] = Eigen::Vector3d(-1, 1, -1);
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        vertices[i] = vertices[i] / vertices[i].norm();
    }
    return vertices;
}

lzd_tools::FiniteMap get_orient_color()
{
    lzd_tools::FiniteMap orient_color;
    orient_color.insert({ SEE_INNER, { 0, 44, 90} });//蓝色
    orient_color.insert({ SEE_OUTER, { 175, 102, 23 } });//橙色
    orient_color.insert({ NOT_SEEN, { 176, 168, 166 } });//灰色
    orient_color.insert({VIEW_POINT, { 219, 103, 76 } });//红色
    orient_color.insert({LOOK_AT_POINT, { 146, 50, 53 } });//红色
    orient_color.insert({ PARALLEL,{0,0,0} });// 黑色
    return orient_color;
}

lzd_tools::FiniteMap get_100_color()
{
    lzd_tools::FiniteMap color;
    for(int i = 0;i<global_var::color_map_100.size();i++){
        auto tmp = global_var::color_map_100[i];
        std::vector<int> tmp_color = {(int)tmp[0],(int)tmp[1],(int)tmp[2]};
        color.insert({i,tmp_color});
    }
    return color;
}

std::vector<std::pair<int, int>> Rasterlization::Image::get_neibor_pixel_idx(int i, int j)
{
    std::vector<std::pair<int,int>> res;
    if(i - 1 >= 0)res.push_back(std::make_pair(i-1,j));
    if(i + 1 < height)res.push_back(std::make_pair(i+1,j));
    if(j - 1 >= 0)res.push_back(std::make_pair(i,j-1));
    if(j + 1 < width)res.push_back(std::make_pair(i,j+1));
    return res;
}

std::vector<std::pair<int, int>> Rasterlization::Image::get_not_empty_neibor_pixel_idx(int i, int j)
{
    std::vector<std::pair<int,int>> res;
    if(i - 1 >= 0 && !is_empty(i-1,j))res.push_back(std::make_pair(i-1,j));
    if(i + 1 < height && !is_empty(i+1,j))res.push_back(std::make_pair(i+1,j));
    if(j - 1 >= 0 && !is_empty(i,j-1))res.push_back(std::make_pair(i,j-1));
    if(j + 1 < width && !is_empty(i,j+1))res.push_back(std::make_pair(i,j+1));
    return res;
}

// std::vector<Rasterlization::PixelData> Rasterlization::Image::get_neibor_pixel_data(int i, int j)
// {
//     std::vector<PixelData> res;
//     auto neibor_idx = get_neibor_pixel_idx(i,j);
//     for(auto idx:neibor_idx){
//         res.push_back(get_pixel_data(idx.first,idx.second));
//     }
//     return res;
// }

// std::vector<PixelData> Rasterlization::Image::get_not_empty_neibor_pixel_data(int i, int j)
// {
//     std::vector<PixelData> res;
//     auto neibor_idx = get_not_empty_neibor_pixel_idx(i,j);
//     for(auto idx:neibor_idx){
//         res.push_back(get_pixel_data(idx.first,idx.second));
//     }
//     return res;
// }

bool Rasterlization::Image::is_empty(int i, int j)
{
    return zBuffer(i,j) >= 1e9;
}

void Rasterlization::Image::cal_boundary(double eps)
{
    // 如果像素点的上下左右四个点中，存在深度值差值大于eps的，则认为是边界像素
    for (int y = 1; y < height-1; ++y)
    {
        for (int x = 1; x < width-1; ++x)
        {
            if (zBuffer(y, x) == 1e9)continue;            
            double d = zBuffer(y, x);
            if (abs(zBuffer(y - 1, x) - d) > eps || abs(zBuffer(y + 1, x) - d) > eps || abs(zBuffer(y, x - 1) - d) > eps || abs(zBuffer(y, x + 1) - d) > eps)
            {
                is_boundary(y, x) = DEPTH_BOUNDARY;
            }
        }
    }
}

void Rasterlization::Image::cal_boundary(std::vector<int> &tri_type,int flag)
{
    for (int y = 1; y < height - 1; ++y)
    {
        for (int x = 1; x < width - 1; ++x)
        {
            if (is_boundary(y,x) != -1 ||  zBuffer(y, x) == 1e9)continue;
            auto _tri_type = tri_type[triangle_id(y, x)];   
            if (zBuffer(y+1, x) != 1e9 && tri_type[triangle_id(y+1,x)] != _tri_type)
            {
                is_boundary(y, x) = flag;
                continue;
            }
            if (zBuffer(y, x+1) != 1e9 && tri_type[triangle_id(y,x+1)] != _tri_type )
            {
                is_boundary(y, x) = flag;
                continue;
            }
            if (zBuffer(y-1, x) != 1e9 && tri_type[triangle_id(y-1,x)] != _tri_type )
            {
                is_boundary(y, x) = flag;
                continue;
            }
            if (zBuffer(y-1, x) != 1e9 &&tri_type[triangle_id(y,x-1)] != _tri_type)
            {
                is_boundary(y, x) = flag;
                continue;
            }
        }
    }

}

void Rasterlization::Image::draw(std::string path)
{
//    // 生成图片
   printf("install pngpp to draw\n");
//    png::image<png::rgb_pixel> image(width, height);
//    // inside绘制成蓝色，outside绘制成橙色，边界绘制成红色
//    for (int y = 0; y < height; ++y)
//    {
//        for (int ix = width-1; ix >0; --ix)
//        {
//            int x = width - ix -1;
//            if (is_inside(y, x) == 1)
//            {
//                image[y][ix] = png::rgb_pixel(0, 44, 90);
//            }
//            else if (is_inside(y, x) == 0)
//            {
//                image[y][ix] = png::rgb_pixel(175, 102, 23);
//            }
//            else
//            {
//                image[y][ix] = png::rgb_pixel(176, 168, 166);
//            }

//            if (is_boundary(y, x) == DEPTH_BOUNDARY)
//            {
//                image[y][ix] = png::rgb_pixel(255, 0, 0);
//            }
//            if (is_boundary(y, x) == MESH_BOUNDARY)
//            {
//                image[y][ix] = png::rgb_pixel(0, 255, 255);
//            }
//        }
//    }
//    image.write(path);

}

void Rasterlization::Image::draw_cluster(std::string path)
{
   // 生成图片
   printf("install pngpp to draw cluster\n");
//    png::image<png::rgb_pixel> image(width, height);
//    for (int y = 0; y < height; ++y)
//    {
//        for (int ix = width - 1; ix > 0; --ix)
//        {
//            int x = width - ix - 1;
           
//            if (is_boundary(y, x) == DEPTH_BOUNDARY){
//                image[y][ix] = png::rgb_pixel(255, 0, 0);
//                continue;
//            }
//            if (is_boundary(y, x) == MESH_BOUNDARY){
//                image[y][ix] = png::rgb_pixel(0, 255, 255);
//                continue;
//            }

//            if (cluster_id(y, x) == -1)
//            {
//                image[y][ix] = png::rgb_pixel(255, 255, 255);
//            }
//            else
//            {
//                int id = cluster_id(y, x) % 100;
//                std::vector<unsigned int> color = global_var::color_map_100[id];
//                image[y][ix] = png::rgb_pixel(color[0], color[1], color[2]);
//            }
//        }
//    }
//    image.write(path);
}
