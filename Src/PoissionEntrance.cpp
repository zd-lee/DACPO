#include <PoissionEntrance.h>
#include <PoissonRecon.h>


POINTS_NORMALS_D sample_points_d(int argc, char* argv[], const POINTS_NORMALS_D& points_normals, XForm<double, 3 + 1>& iXForm, std::vector<double>* weight_samples)
{
    return sample_points(argc,argv,points_normals,iXForm,weight_samples);
}

MESH_D poisson_reconstruction_d(int argc, char* argv[], const POINTS_NORMALS_D& points_normals, const std::vector<double>* weight_samples)
{
    //static lzd_tools::resource_controller controller(10);
    //controller.P();
    MESH_D res =  poisson_reconstruction(argc, argv, points_normals, weight_samples);
    //controller.V();
    return res;
}
