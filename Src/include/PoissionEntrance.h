#pragma once 
#include "PointStream.h"
#include "PointStreamData.h"

#define MESH std::pair<std::vector<Point<REAL, DIM>>,std::vector<std::vector<int>>>   
// same as MESH, but use double
#define MESH_D std::pair<std::vector<Point<double, 3>>,std::vector<std::vector<int>>>

// first: points; second: normals
#define POINTS_NORMALS std::vector<std::pair<Point<REAL, DIM>, Normal<REAL, DIM>>>
// same as POINTS_NORMALS, but use double
#define POINTS_NORMALS_D std::vector<std::pair<Point<double, 3>, Normal<double, 3>>>

#define ORIENTED_POINTS std::vector<std::pair<Point<REAL, DIM>, Normal<REAL, DIM>>>


POINTS_NORMALS_D sample_points_d(int argc, char* argv[], const POINTS_NORMALS_D& points_normals, XForm<double, 3 + 1>& iXForm, std::vector<double>* weight_samples);
MESH_D poisson_reconstruction_d(int argc, char* argv[], const POINTS_NORMALS_D& points_normals, const std::vector<double>* weight_samples);



template <class REAL, unsigned int DIM>
POINTS_NORMALS sample_points_entrance(int argc, char *argv[], const POINTS_NORMALS &points_normals, XForm<REAL, DIM + 1> &iXForm, std::vector<double> *weight_samples) {
    return sample_points_d(argc,argv,points_normals,iXForm,weight_samples);
}

template <class REAL, unsigned int DIM>
MESH poisson_reconstruction_entrance(int argc, char *argv[], const POINTS_NORMALS &points_normals, const std::vector<double> *weight_samples) {
    return poisson_reconstruction_d(argc, argv, points_normals, weight_samples);
}