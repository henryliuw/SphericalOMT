
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

#define _EPS_OMT 1e-6
#define _PI 3.14159265358979323846

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Regular_triangulation_2<K> RT;
typedef CGAL::Regular_triangulation_adaptation_traits_2<RT>         AT;
typedef CGAL::Regular_triangulation_degeneracy_removal_policy_2<RT> DRP;
typedef CGAL::Voronoi_diagram_2<RT, AT, DRP> VD;

//definition

struct HASH{
    size_t operator()(const K::Point_2& p) const {
        double x = (double)p.x().exact();
        double y = (double)p.y().exact();
        std::size_t hi = std::hash<double>{}(x);
        std::size_t hj = std::hash<double>{}(y);
        return hi ^ (hj << 32); //maybe too big ? i dont know
    }
};

template<int N> // construct a OMT class with N vertices
class SphericalOMT
{
	// data members
	int size;
    std::unordered_map<K::Point_2, double, HASH> point2idx; // i dont know how to get the index of VD cell so maintain a hashtable here for matching
	std::array<K::Point_2, N> points; // vertice coordinate in 2D complex plane
	std::array<K::Weighted_point_2, N> wpoints; // store h in wpoints[i].weight()
    Eigen::VectorXd h; // wpoints is immutable so we maintain an h
    RT rt;
	VD vd; 
    std::array<VD::Face_handle,N> faces;
    std::array<CGAL::Polygon_2<K>,N> polys; // polygon representation of cells, should be updated each time
    Eigen::VectorXd nu;
    //Eigen::Matrix<double, N, 1> grad; // gradient of E(h)
	//Eigen::SparseMatrix<double> hess; // hessian of E(h)

    bool constructDiagram(bool verbose); // construct Power diagram and everything from current wpoints, rewrite rt, vd, faces, polys, return false if PD degenerate
	inline K::Point_2 sgProj(std::array<double, 3> vertex_position); // stereographic projection
    Eigen::SparseMatrix<double> hessian(); // computes the hessian matrix
    double hij(int i, int j);
    double lineIntegral(K::Segment_2 seg); // perform closed form integral calculation for hessian denominator
    double curvatureIntegral(K::Segment_2 seg); // perform geodesic curvature integral for a segment 
    double polygonArea(CGAL::Polygon_2<K> pgon);
    Eigen::VectorXd Wh();
    void setHeight(Eigen::VectorXd h);
    void Nu(); // set nu
public:
	SphericalOMT(std::vector<std::array<double, 3>> vertex_positions); // construct the VD for the first time using inputs of 3D points
	void printPoly(CGAL::Polygon_2<K> pgon);
    void printPolys();
    void debug();
    double step(double lbda); // perform a single step optimization and return the norm of gradient after optimization
	//void gradient(); 
};

// implementation start
template<int N>
void SphericalOMT<N>::debug() {
    //std::cout << polygonArea(polys[4]);
    //for (auto v : vd.dual().all_vertex_handles()) {
    //    std::cout << v->point() << std::endl;
    //    VD::Face
    //    //std::cout << vd.dual().dual(*v) << std::endl;
    //}
    for (VD::Face_iterator fit = vd.faces_begin(); fit!=vd.faces_end();++fit) {
        std::cout << fit->dual()->point() << " #N:" << point2idx[K::Point_2(fit->dual()->point().x(), fit->dual()->point().y())] << std::endl;
    }
    //for (VD::Site_iterator s = vd.sites_begin(); s != vd.sites_end(); s++) {
    //    std::cout << s->cartesian(0) << " " << s->cartesian(1) << std::endl;
    //}
    //auto v = vd.site(0);
    //std::cout << v->point() << std::endl;
    //std::cout << Eigen::MatrixXd(hessian()) << std::endl;
    //std::cout << gradient() << std::endl;
    //hij(4, 0);
    //for (int i = 0; i != this->size; i++)
    //{
    //   for (int j = 0; j != i; j++)
    //        std::cout << "h(" << i << "," << j << ") is " << hij(i, j) << std::endl;
            //if (hij(i, j)!=0)
            //{
            //    std::cout << "edge (" << i << ',' << j << ") is connected.\n";
            //}
            //else
            //{
            //    std::cout << "edge (" << i << ',' << j << ") is not connected.\n";
            //}
    //}
    //std::cout << "hij is" << hij(0, 1);
}

template<int N>
inline K::Point_2 SphericalOMT<N>::sgProj(std::array<double, 3> vertex_position) { // stereographic projection
    //std::cout << vertex_position[0] << vertex_position[1] << vertex_position[2] << std::endl;
    double x = vertex_position[0] / (1 - vertex_position[2]);
	double y = vertex_position[1] / (1 - vertex_position[2]);
	std::cout << "x:" << x << "\t y:" << y << std::endl;
	return K::Point_2(x, y);
}

template<int N>
bool SphericalOMT<N>::constructDiagram(bool verbose) {
    //Create a Regular Triangulation from the points
    this->rt =  RT(wpoints.begin(), wpoints.end());
    if (!this->rt.is_valid())
        return false;
    if (this->rt.number_of_hidden_vertices() != 0) {
        return false;
    }
    //Wrap the triangulation with a Voronoi diagram adaptor. This is necessary to
    //get the Voronoi faces.
    this->vd =  VD(this->rt);
    // construct mapping
    // CGAL returns iterator in random order so we need to find the mapping by querying each point
    //CGAL often returns objects that are either segments or rays. This converts
    //these objects into segments. If the object would have resolved into a ray,
    //that ray is intersected with the bounding box defined above and returned as
    //a segment.
    const auto ConvertToSeg = [&](const CGAL::Object seg_obj, bool outgoing) -> K::Segment_2 {
        //One of these will succeed and one will have a NULL pointer
        auto RAY_LENGTH = 1000;
        const K::Segment_2* dseg = CGAL::object_cast<K::Segment_2>(&seg_obj);
        const K::Ray_2* dray = CGAL::object_cast<K::Ray_2>(&seg_obj);
        if (dseg) { //Okay, we have a segment
            return *dseg;
        }
        else {    //Must be a ray
            const auto& source = dray->source();
            const auto dsx = source.x();
            const auto dsy = source.y();
            const auto& dir = dray->direction();
            auto norm = sqrt((double)(dir.dx() * dir.dx() + dir.dy() * dir.dy()).exact());
            const auto tpoint = K::Point_2(dsx + RAY_LENGTH * dir.dx() / norm, dsy + RAY_LENGTH * dir.dy() / norm);
            if (outgoing)
                return K::Segment_2(
                    dray->source(),
                    tpoint
                );
            else
                return K::Segment_2(
                    tpoint,
                    dray->source()
                );
        }
    };

    for (VD::Face_iterator fit = vd.faces_begin(); fit!=vd.faces_end();++fit) {
        int i = point2idx[K::Point_2(fit->dual()->point().x(), fit->dual()->point().y())];
        //auto t = boost::get<VD::Face_handle>(&lr); // it must be a face
        faces[i] = *fit;
        CGAL::Polygon_2<K> pgon;

        //Edge circulators traverse endlessly around a face. Make a note of the
        //starting point so we know when to quit.
        VD::Face::Ccb_halfedge_circulator ec_start = faces[i]->ccb();

        //Find a bounded edge to start on
        //for(;ec_start->is_unbounded();ec_start++){}

        //Current location of the edge circulator
        VD::Face::Ccb_halfedge_circulator ec = ec_start;

        //In WKT format each polygon must begin and end with the same point
        K::Point_2 first_point;

        do {
            //A half edge circulator representing a ray doesn't carry direction
            //information. To get it, we take the dual of the dual of the half-edge.
            //The dual of a half-edge circulator is the edge of a Delaunay triangle.
            //The dual of the edge of Delaunay triangle is either a segment or a ray.
            const CGAL::Object seg_dual = vd.dual().dual(ec->dual()); // i get myself... it this APi really so stupid or I dont understand how to use it??
            //const CGAL::Object seg_dual = rt.dual(ec->dual());
            //Convert the segment/ray into a segment
            const auto this_seg = ConvertToSeg(seg_dual, ec->has_target());

            pgon.push_back(this_seg.source());
            if (ec == ec_start)
                first_point = this_seg.source();

            //If the segment has no target, it's a ray. This means that the next
            //segment will also be a ray. We need to connect those two rays with a
            //segment. The following accomplishes this.
            if (!ec->has_target()) {
                const CGAL::Object nseg_dual = vd.dual().dual(ec->next()->dual());
                const auto next_seg = ConvertToSeg(nseg_dual, ec->next()->has_target());
                pgon.push_back(next_seg.target());
            }
        } while (++ec != ec_start); //Loop until we get back to the beginning
                                    //Print the polygon as a WKT polygon
        polys[i] = pgon;
        if (verbose) {
            std::cout << "the " << i << "-th polygon:\n";
            printPoly(pgon);
        }
    }
    return true;
}

template<int N>
void SphericalOMT<N>::setHeight(Eigen::VectorXd h){
    for (int i=0; i!=this->size ; i++)
        this->wpoints[i] = K::Weighted_point_2(points[i], h(i));
}


template<int N>
SphericalOMT<N>::SphericalOMT(std::vector<std::array<double, 3>> vertex_positions) {
    //std::cout << "initial construction cgal is stupid:\n";
    //get projection and construct weighted graph
    // construct the mapping 
    this->size = vertex_positions.size();
    for (int i = 0; i != this->size; i++) {
        this->points[i] = (this->sgProj(vertex_positions[i]));
        point2idx.insert(std::make_pair(points[i], i));
    }
    h = Eigen::VectorXd::Zero(size);
    setHeight(h);
    Nu();
    constructDiagram(false);
}

template<int N>
void SphericalOMT<N>::printPoly(CGAL::Polygon_2<K> pgon) {
    std::cout << "POLYGON ((";
    for (auto v = pgon.vertices_begin(); v != pgon.vertices_end(); v++)
        std::cout << v->x() << " " << v->y() << ", ";
    std::cout << pgon.vertices_begin()->x() << " " << pgon.vertices_begin()->y() << "))\n";
}

template<int N>
void SphericalOMT<N>::printPolys() {
    for (auto p : polys)
        printPoly(p);
}

template<int N>
Eigen::SparseMatrix<double> SphericalOMT<N>::hessian()
{
    // iterate over all halfedge
    Eigen::SparseMatrix<double> hess(size, size);
    //for (int i = 0; i!=this->size; )
    for (int i = 0; i != this->size; i++)
    {
        for (int j = 0; j != i; j++) {
            double entry_ij = hij(i, j);
            if (entry_ij) {
                hess.insert(i, j) = entry_ij;    // symmetric
                hess.insert(j, i) = entry_ij;
                hess.coeffRef(i, i) += entry_ij; // set diagonal
                hess.coeffRef(j, j) += entry_ij;
            }
        }
    }
    return hess;
}

template<int N>
double SphericalOMT<N>::hij(int i, int j) {
    //std::cout << " " << i << j << std::endl;
    RT::Vertex_handle vi = faces[i]->dual();
    RT::Vertex_handle vj = faces[j]->dual();
    if (rt.is_edge(vi, vj)) {
        // find the segment that is perpendicular to yi - yj
        // vector p[i] - p[j]
        K::FT pij_x = points[i].x() - points[j].x();
        K::FT pij_y = points[i].y() - points[j].y();
        CGAL::Polygon_2<K> pgon = polys[i]; // the polygon cell
        for (CGAL::Polygon_2<K>::Edge_const_iterator edge = pgon.edges_begin(); edge != pgon.edges_end(); edge++)
        {
            // compute inner product
            const K::Segment_2 seg = (K::Segment_2)(*edge);
            K::FT seg_x = seg.source().x() - seg.target().x();
            K::FT seg_y = seg.source().y() - seg.target().y();
            K::FT inner_prod = seg_x * pij_x + seg_y * pij_y;
            if (abs(inner_prod) < _EPS_OMT)
            {
                double length = sqrt( (double)((pij_x * pij_x+pij_y * pij_y).exact()));
                double denom = lineIntegral(seg);
                return  denom / length;
            }
        }
        //should not reach here, so return a negative value that should not appear for debugging
        return -1;
    }
    else
        return 0;
}

template<int N>
double SphericalOMT<N>::lineIntegral(K::Segment_2 seg) {
    // perform a closed form for computing \int_{x_0}^{x_1} \mu(X, Y) \sqrt(1+k^2) dx where \mu(X,Y) = 4 / (1+X^2+Y^2)^2
    double x0 = (double)seg.source().x().exact(), y0 = (double)seg.source().y().exact(), x1 = (double)seg.target().x().exact(), y1 = (double)seg.target().y().exact();
    // compute coeff for y = kx+b
    double k = (y1 - y0) / (x1 - x0);
    double b = y0 - k * x0;
    // compute coeff that the integrand is in the form of D / (Ax^2+Bx+C)^2
    double A = 1 + k * k;
    double B = 2 * k * b;
    double C = 1 + b * b;
    double D = 4 * sqrt(A);
    double delta = B * B - 4 * A * C;
    double sqrt_delta = sqrt(-delta);
    // compute the close form integral thanks to almighty mathematica
    const auto integral_nominator = [&](double x) -> double {
        return (B + 2 * A * x) / (C + x * (B + A * x)) + 4 * A * atan((B + 2 * A * x) / sqrt_delta) / sqrt_delta;
    };
    double integral =  -(integral_nominator(x1) - integral_nominator(x0)) / delta * D;
    if (x1 > x0)
        return integral; // right direction
    else return -integral;
}

template<int N>
double SphericalOMT<N>::curvatureIntegral(K::Segment_2 seg) {
    // perform a closed form for computing \int_si k_g ds
    double x0 = (double)seg.source().x().exact(), y0 = (double)seg.source().y().exact(), x1 = (double)seg.target().x().exact(), y1 = (double)seg.target().y().exact();
    // compute coeff for y = kx+b
    double k = (y1 - y0) / (x1 - x0);
    double b = y0 - k * x0;
    // compute coeff that the integrand is in the form of D / (Ax^2+Bx+C) = kg(x)*length(dx) where 
    double A = 1 + k * k;
    double B = 2 * k * b;
    double C = 1 + b * b;
    double D = 2 * b;
    double sqrt_delta = sqrt(-B * B + 4 * A * C);
    // compute the close form integral thanks to almighty mathematica again
    const auto integral_nominator = [&](double x) -> double {
        return atan((B + 2 * A * x) / sqrt_delta);
    };
    double integral =  -(integral_nominator(x1) - integral_nominator(x0)) / sqrt_delta * 2 * D;
    return integral;
}


template<int N>
double SphericalOMT<N>::polygonArea(CGAL::Polygon_2<K> pgon) {
    double integral = 0;
    for (CGAL::Polygon_2<K>::Edge_const_iterator edge = pgon.edges_begin(); edge != pgon.edges_end(); edge++)
    {
        // compute inner product
        const K::Segment_2 seg = (K::Segment_2)(*edge);
        K::FT seg_x = seg.source().x() - seg.target().x();
        K::FT seg_y = seg.source().y() - seg.target().y();
        integral += curvatureIntegral(seg);
    }
    return abs(integral); // the area is signed by integral direction but we want the absolute value
}

template<int N>
Eigen::VectorXd SphericalOMT<N>::Wh() {
    Eigen::VectorXd wh(size);
    for (int i = 0; i != size; i++) {
        wh(i) = polygonArea(polys[i]);
    }
    return wh;
}

template<int N>
void SphericalOMT<N>::Nu(){
    nu = Eigen::MatrixXd::Constant(size, 1, 4 * _PI / size);
};

template<int N>
double SphericalOMT<N>::step(double lbda) {
    //std::cout << "hello convex function\n";
    Eigen::VectorXd wh = Wh();
    Eigen::VectorXd grad = wh - nu;
    //std::cout << "sum of wh:" << wh.sum() << std::endl;
    Eigen::SparseMatrix<double> hess = hessian();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(hess);
    Eigen::VectorXd dh = solver.solve(grad);
    setHeight(h - dh * lbda);
    while (!constructDiagram(false)) { // line search
        std::cout << "degenerate at lambda=" << lbda << std::endl;
        lbda = lbda * 0.5;
        setHeight(h - dh * lbda);
    };
    h = h - dh * lbda;
    std::cout << "Wh:\n" << wh << std::endl;
    std::cout << "dh\n" << dh << std::endl;
    std::cout << "h:\n" << h << std::endl;
    std::cout << "norm of grad:" << grad.norm() << std::endl;
    //std::cout << rt.number_of_hidden_vertices() << std::endl;
    return grad.norm();
}