#include <o3d_api.h>
#include <string>



void test_open3d() {
    auto sphere = open3d::geometry::TriangleMesh::CreateSphere(1.0);
    sphere->ComputeVertexNormals();
    sphere->PaintUniformColor({ 0.0, 2.0, 0.0 });
    open3d::visualization::DrawGeometries({ sphere });
    return ;
}

