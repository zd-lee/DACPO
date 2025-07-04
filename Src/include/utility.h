/*
Copyright (c) 2022, Fei Hou and Chiyu Wang, Institute of Software, Chinese Academy of Sciences.
All rights reserved.

The code can only be used for academic purpose, and cannot be used for commercial
purpose without written permission.

Redistribution and use in source for academic purpose, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include <string>
#include <climits>
#include <fstream>
#include "happly.h"
#include "kdtree.h"
#include "PoissionEntrance.h"

template <class Real, unsigned int Dim>
void transform(std::vector<std::pair<Point<Real, Dim>, Normal<Real, Dim>>> &points_normals, const XForm<Real, Dim + 1> &iXForm)
{
	for (size_t i = 0; i < points_normals.size(); ++i)
	{
		points_normals[i].first = iXForm * points_normals[i].first;
	}
}

template <class Real, unsigned int Dim>
void ply_reader(const std::string& file, std::vector<std::pair<Point<Real, Dim>, Normal<Real, Dim>>>& points_normals)
{
	happly::PLYData plyIn(file);
	// Get data from the object

	std::vector<Real> x = plyIn.getElement("vertex").getProperty<Real>("x");
	std::vector<Real> y = plyIn.getElement("vertex").getProperty<Real>("y");
	std::vector<Real> z = plyIn.getElement("vertex").getProperty<Real>("z");
	
	bool exist_normal = plyIn.getElement("vertex").hasProperty("nx") && plyIn.getElement("vertex").hasProperty("ny") && plyIn.getElement("vertex").hasProperty("nz");
	std::vector<Real>nx, ny, nz;

	// if exist normal
	if(exist_normal){
		nx = plyIn.getElement("vertex").getProperty<Real>("nx");
		ny = plyIn.getElement("vertex").getProperty<Real>("ny");
		nz = plyIn.getElement("vertex").getProperty<Real>("nz");

	}
	Point<Real, Dim> p, n;

	for (size_t i = 0; i < x.size(); ++i) {
		p[0] = x[i];
		p[1] = y[i];
		p[2] = z[i];
		if (exist_normal) {
			n[0] = nx[i];
			n[1] = ny[i];
			n[2] = nz[i];
		}
		else
		{
			n[0] = n[1] = n[2] = 0.01;
		}
		Normal<Real, Dim> normal(n);
		points_normals.push_back(std::make_pair(p, normal));
	}

	// PLYInputPointStream<Real, Dim> ply(file.c_str());
	// Point<Real, Dim> data;
	// while (ply.nextPoint(data)) {
	// 	for (unsigned int i = 0; i < Dim; i++) {
	// 		p[i] = data[i];
	// 		n[i] = data[i];
	// 	}
	// 	Normal<Real, Dim> normal(n);
	// 	points_normals.push_back(std::make_pair(p,normal));
	// }


}

template <class Real, unsigned int Dim>
bool output_ply(const std::string &outFile, const std::pair<std::vector<Point<Real, Dim>>, std::vector<std::vector<int>>> &mesh, const XForm<Real, Dim + 1> &iXForm)
{
	const std::vector<Point<Real, Dim>> &points = mesh.first;
	const std::vector<std::vector<int>> &faces = mesh.second;

	std::ofstream plyfile;
	plyfile.open(outFile, std::ofstream::out);
	if (!plyfile)
	{
		printf("Cannot save result file %s\n", outFile.c_str());
		return false;
	}
	printf("writing to %s\n", outFile.c_str());

	plyfile << "ply\nformat ascii 1.0\n";
	plyfile << "element vertex " << points.size() << std::endl;
	plyfile << "property float x" << std::endl
		<< "property float y" << std::endl
		<< "property float z" << std::endl;
	plyfile << "element face " << faces.size() << std::endl;
	plyfile << "property list uchar int vertex_index" << std::endl;
	plyfile << "end_header" << std::endl;

	for (size_t i = 0; i < points.size(); ++i)
	{
		Point<Real, Dim> p = iXForm * points[i];
		plyfile << p[0] << " " << p[1] << " " << p[2] << std::endl;
	}

	for (size_t i = 0; i < faces.size(); ++i)
	{
		plyfile << faces[i].size();
		for (size_t j = 0; j < faces[i].size(); ++j)
			plyfile << " " << faces[i][j];
		plyfile << std::endl;
	}
	plyfile.close();
	return true;
}

/// <summary>
/// 根据颜色映射表 绘制ply。注意颜色表的个数等于flag中不同的值的个数
/// </summary>
/// <typeparam name="Real"></typeparam>
/// <typeparam name="Dim"></typeparam>
/// <param name="outFile"></param>
/// <param name="points_normals"></param>
/// <param name="iXForm"></param>
/// <param name="flag"></param>
/// <returns></returns>
template <class Real, unsigned int Dim>
bool output_sample_points_and_normals(const std::string& outFile, const std::vector<std::pair<Point<Real, Dim>,
	Normal<Real, Dim>>>& points_normals, const XForm<Real, Dim + 1>& iXForm,
	std::pair<std::vector<int>, std::vector<std::vector<int>>> flag = std::make_pair(std::vector<int>(), std::vector<std::vector<int>>()))
{
	std::ofstream plyfile;
	plyfile.open(outFile, std::ofstream::out);
	if (!plyfile)
	{
		printf("Cannot save result file %s\n", outFile.c_str());
		return false;
	}
	printf("writing to %s\n", outFile.c_str());

	bool has_label = flag.first.size() == points_normals.size();

	plyfile << "ply\nformat ascii 1.0\n";
	plyfile << "element vertex " << points_normals.size() << std::endl;
	plyfile << "property float x" << std::endl
		<< "property float y" << std::endl
		<< "property float z" << std::endl;
	plyfile << "property float nx" << std::endl
		<< "property float ny" << std::endl
		<< "property float nz" << std::endl;
	if (has_label){
		plyfile << "property uchar red" << std::endl;
		plyfile << "property uchar green" << std::endl;
		plyfile << "property uchar blue" << std::endl;
	}

	plyfile << "element face 0" << std::endl;
	plyfile << "property list uchar int vertex_index" << std::endl;
	plyfile << "end_header" << std::endl;

	for (size_t i = 0; i < points_normals.size(); ++i)
	{
		auto p = iXForm * points_normals[i].first;
		plyfile << p[0] << " " << p[1] << " " << p[2] << " ";
		plyfile << points_normals[i].second.normal[0] << " " << points_normals[i].second.normal[1] << " " << points_normals[i].second.normal[2];
		if(has_label){
			plyfile << " ";
			auto type = flag.first[i];
			if(type <=0){
				plyfile << "0 0 0" << std::endl;
			}
			else{
				int r = flag.second[type][0];
				int g = flag.second[type][1];
				int b = flag.second[type][2];
				plyfile << r << " " << g << " " << b << std::endl;
			}
		}else{
			plyfile << std::endl;
		}
	}
	plyfile.close();
	return true;
}

template <class Real, unsigned int Dim>
bool output_all_points_and_normals(const std::string& outFile, const std::vector<std::pair<Point<Real, Dim>,
	Normal<Real, Dim>>>& points_normals, const XForm<Real, Dim + 1>& iXForm,
	std::vector<int> labels) {
	


}



template <class Real, unsigned int Dim>
bool output_all_points_and_normals(const std::string &outFile, const std::string &input_name, 
	const std::vector<std::pair<Point<Real, Dim>, Normal<Real, Dim>>> &points_normals, 
	const kdt::KDTree<kdt::KDTreePoint> &tree, 
	const XForm<Real, Dim + 1> &iXForm)
{
	std::vector<std::pair<Point<Real, Dim>, Normal<Real, Dim>>> points_normals_all;
	ply_reader<Real, Dim>(input_name, points_normals_all);
	auto inv_iXForm = iXForm.inverse();
	for (size_t i = 0; i < points_normals_all.size(); ++i)
	{
		auto c = inv_iXForm * points_normals_all[i].first;
		std::array<Real, Dim> a{ c[0], c[1], c[2] };
		int n = tree.nnSearch(kdt::KDTreePoint(a));
		points_normals_all[i].second = points_normals[n].second;
	}

	std::ofstream plyfile;
	plyfile.open(outFile, std::ofstream::out);
	if (!plyfile)
	{
		printf("Cannot save result file %s\n", outFile.c_str());
		return false;
	}
	printf("writing to %s\n", outFile.c_str());

	plyfile << "ply\nformat ascii 1.0\n";
	plyfile << "element vertex " << points_normals_all.size() << std::endl;
	plyfile << "property float x" << std::endl
		<< "property float y" << std::endl
		<< "property float z" << std::endl;
	plyfile << "property float nx" << std::endl
		<< "property float ny" << std::endl
		<< "property float nz" << std::endl;
	plyfile << "element face 0" << std::endl;
	plyfile << "property list uchar int vertex_index" << std::endl;
	plyfile << "end_header" << std::endl;

	for (size_t i = 0; i < points_normals_all.size(); ++i)
	{
		const auto& p =  points_normals_all[i].first;
		plyfile << p[0] << " " << p[1] << " " << p[2] << " ";
		plyfile << points_normals_all[i].second.normal[0] << " " << points_normals_all[i].second.normal[1] << " " << points_normals_all[i].second.normal[2] << std::endl;
	}
	plyfile.close();
	return true;
}

template <class Real, int Dim>
bool operator==(const Normal<Real, Dim> &n1, const Normal<Real, Dim> &n2)
{
	for (int i = 0; i < Dim; ++i)
		if (n1.normal[i] != n2.normal[i])
			return false;
	return true;
}

template <class Real, unsigned int Dim>
void normalize(Normal<Real, Dim> &n)
{
	Real len = 0;
	for (unsigned int i = 0; i < Dim; ++i)
		len += n.normal[i] * n.normal[i];
	if (len != 0)
	{
		len = sqrt(len);
		for (unsigned int i = 0; i < Dim; ++i)
			n.normal[i] /= len;
	}
}

template <class Real, unsigned int Dim>
Real dot(const Normal<Real, Dim> &n1, const Normal<Real, Dim> &n2)
{
	Real d = 0;
	for (unsigned int i = 0; i < Dim; ++i)
		d += n1.normal[i] * n2.normal[i];
	return d;
}

// 返回弧度制夹角
template <class Real, unsigned int Dim>
Real calculate_angle(const Normal<Real, Dim> &n1, const Normal<Real, Dim> &n2)
{
	// normalize
	Normal<Real, Dim> n1n = n1;
	normalize(n1n);
	Normal<Real, Dim> n2n = n2;
	normalize(n2n);
	Real d = dot(n1n, n2n);
	if (d > 1)
		d = 1;//防止arcos报错
	else if (d < -1)
		d = -1;
	return acos(d);
}


inline std::vector<std::string> split(const std::string &s, char c = ' ')
{
	std::vector<std::string> str;
	unsigned pos = 0;
	while (pos < s.size())
	{
		while (pos < s.size() && s[pos] == c)
			++pos;

		unsigned int end = pos;

		do
		{
			++end;
		} while (end < s.size() && s[end] != c);

		if (pos < s.size())
			str.push_back(s.substr(pos, end - pos));
		pos = end;
	}

	return str;
}

inline bool valid_parameter(long v)
{
	return v >= 0 && v < INT_MAX;
}

#endif
