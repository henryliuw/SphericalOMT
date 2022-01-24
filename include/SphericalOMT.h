
#pragma once

#include <CGAL/number_utils.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_adaptation_traits_2.h>
#include <CGAL/Regular_triangulation_adaptation_policies_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>	
#include <vector>
#include <array>
#include <math.h>
#include <unordered_map>
#include <string>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>
#include <CGAL/Delaunay_triangulation_on_sphere_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>
#include <string>

#define _EPS_OMT 1e-6
#define _PI 3.14159265358979323846

typedef CGAL::Simple_cartesian<double>             Kernel;
typedef Kernel::Point_3                            Point;
typedef CGAL::Surface_mesh<Point>                  Mesh;
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Regular_triangulation_2<K> RT;
typedef CGAL::Regular_triangulation_adaptation_traits_2<RT>         AT;
typedef CGAL::Regular_triangulation_degeneracy_removal_policy_2<RT> DRP;
typedef CGAL::Voronoi_diagram_2<RT, AT, DRP>  VD;
typedef CGAL::Exact_predicates_inexact_constructions_kernel   KS2;
typedef CGAL::Delaunay_triangulation_on_sphere_traits_2<KS2>  TraitsS2;
typedef CGAL::Delaunay_triangulation_on_sphere_2<TraitsS2>    DToS2;
typedef TraitsS2::Point_3                                     Point_3;
//definition
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> readFile(std::string file);
std::pair<Eigen::MatrixXd, Eigen::VectorXd> SphereTriangulation(Eigen::MatrixXd vertex_positions); // return face connectivity matrix and face area


struct HASH{
    size_t operator()(const K::Point_2& p) const {
        double x = (double)p.x().exact();
        double y = (double)p.y().exact();
        std::size_t hi = std::hash<double>{}(x);
        std::size_t hj = std::hash<double>{}(y);
        return hi ^ (hj << 32); //maybe too big ? i dont know
    }
};

struct HASH_P3{
    size_t operator()(const Point_3& p) const {
        double x = (double)p.x();
        double y = (double)p.y();
        double z = (double)p.z();
        std::size_t hi = std::hash<double>{}(x);
        std::size_t hj = std::hash<double>{}(y);
        std::size_t hk = std::hash<double>{}(z);
        return hi ^ (hj << 16) ^ (hk << 32); //maybe too big ? i dont know
    }
};

// construct a OMT class with N vertices
class SphericalOMT
{
	// data members
    size_t size;
    std::unordered_map<K::Point_2, size_t, HASH> point2idx; // i dont know how to get the index of VD cell so maintain a hashtable here for matching
	std::vector<K::Point_2> points; // vertice coordinate in 2D complex plane
	std::vector<K::Weighted_point_2> wpoints; // store h in wpoints[i].weight()
    Eigen::VectorXd h; // wpoints is immutable so we maintain an h
    RT rt;
	VD vd; 
    std::vector<VD::Face_handle> faces;
    std::vector<CGAL::Polygon_2<K>> polys; // polygon representation of cells, should be updated each time
    Eigen::VectorXd nu;
    //Eigen::Matrix<double, N, 1> grad; // gradient of E(h)
	//Eigen::SparseMatrix<double> hess; // hessian of E(h)

    bool constructDiagram(bool verbose); // construct Power diagram and everything from current wpoints, rewrite rt, vd, faces, polys, return false if PD degenerate
	inline K::Point_2 sgProj(std::array<double, 3> vertex_position); // stereographic projection
    inline std::array<double, 3> isgProj(K::Point_2 point2); // inverse stereographic projection
    Eigen::SparseMatrix<double> hessian(); // computes the hessian matrix
    Eigen::SparseMatrix<double> hessByEdge();
    double hij(int i, int j);
    double lineIntegral(K::Segment_2 seg); // perform closed form integral calculation for hessian denominator
    double lineIntegral(VD::Edge_iterator ei);
    double curvatureIntegral(K::Segment_2 seg); // perform geodesic curvature integral for a segment 
    double polygonArea(CGAL::Polygon_2<K> pgon);
    Eigen::VectorXd Wh();
    void setHeight(Eigen::VectorXd h);
    void Nu(); // set nu
public:
	SphericalOMT(std::vector<std::array<double, 3>> vertex_positions); // construct the VD for the first time using inputs of 3D points
    SphericalOMT(Eigen::MatrixXd mat); // usinmg matrix 
	void printPoly(CGAL::Polygon_2<K> pgon);
    void printPolys();
    void debug();
    double step(double lbda); // perform a single step optimization and return the norm of gradient after optimization
    Eigen::MatrixXi faceMat(); // return the connectivity matrix for testing
    Eigen::MatrixXd mapBack();
    void setNu(std::string filepath);
    void printDistortion();
	//void gradient(); 
};
