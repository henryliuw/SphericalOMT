#include "SphericalOMT.h"
// implementation start
void SphericalOMT::debug() {
    //std::cout << polygonArea(polys[4]);
    //for (auto v : vd.dual().all_vertex_handles()) {
    //    std::cout << v->point() << std::endl;
    //    VD::Face
    //    //std::cout << vd.dual().dual(*v) << std::endl;
    //}
    //for (VD::Edge_iterator ei = vd.edges_begin(); ei != vd.edges_end(); ei++) {
    //    if (ei->has_source())
    //        std::cout << "source:" << ei->source()->point() << std::endl;
    //    else
    //        std::cout << "no source\n";
    //    if (ei->has_target())
    //        std::cout << "target:" << ei->target()->point() << std::endl;
    //    else
    //        std::cout << "no target\n";
    //    std::cout << "vertex i:" << ei->up()->point() << std::endl;
    //    std::cout << "vertex j:" << ei->down()->point() << std::endl;
    //}
    std::cout << hessian() << std::endl;
    std::cout << hessByEdge() << std::endl;
    //for (VD::Face_iterator fit = vd.faces_begin(); fit!=vd.faces_end();++fit) {
    //    std::cout << fit->dual()->point() << " #N:" << point2idx[K::Point_2(fit->dual()->point().x(), fit->dual()->point().y())] << std::endl;
    //}
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

inline K::Point_2 SphericalOMT::sgProj(std::array<double, 3> vertex_position) { // stereographic projection
                                                                                //std::cout << vertex_position[0] << vertex_position[1] << vertex_position[2] << std::endl;
    double x = vertex_position[0] / (1 - vertex_position[2]);
    double y = vertex_position[1] / (1 - vertex_position[2]);
    //std::cout << "x:" << x << "\t y:" << y << std::endl;
    return K::Point_2(x, y);
}

inline std::array<double, 3> SphericalOMT::isgProj(K::Point_2 point2) {
    double x = (double) point2.x();
    double y = (double) point2.y();
    double A = 1 + x * x + y * y;
    return std::array<double, 3>{ 2*x / A, 2 * y / A, (A - 2) / A  };
};


bool SphericalOMT::constructDiagram(bool verbose) {
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
            auto norm = sqrt((double)(dir.dx() * dir.dx() + dir.dy() * dir.dy()));
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

void SphericalOMT::setHeight(Eigen::VectorXd h){
    for (size_t i=0; i!=this->size ; i++)
        this->wpoints[i] = K::Weighted_point_2(points[i], h(i));
}

SphericalOMT::SphericalOMT(std::vector<std::array<double, 3>> vertex_positions) {
    //std::cout << "initial construction cgal is stupid:\n";
    //get projection and construct weighted graph
    // construct the mapping 
    this->size = vertex_positions.size();
    // initialize vector to known size
    faces = std::vector<VD::Face_handle>{size};
    polys = std::vector<CGAL::Polygon_2<K>>{size};
    points = std::vector<K::Point_2>{size}; 
    wpoints = std::vector<K::Weighted_point_2>{size}; 
    for (size_t i = 0; i != this->size; i++) {
        this->points[i] = (this->sgProj(vertex_positions[i]));
        point2idx.insert(std::make_pair(points[i], i));
    }
    h = Eigen::VectorXd::Zero(size);
    setHeight(h);
    Nu();
    constructDiagram(false);
}


SphericalOMT::SphericalOMT(Eigen::MatrixXd mat) {
    //std::cout << "initial construction cgal is stupid:\n";
    //get projection and construct weighted graph
    // construct the mapping 
    this->size = mat.rows();
    // initialize vector to known size
    faces = std::vector<VD::Face_handle>{size};
    polys = std::vector<CGAL::Polygon_2<K>>{size};
    points = std::vector<K::Point_2>{size}; 
    wpoints = std::vector<K::Weighted_point_2>{size}; 
    for (size_t i = 0; i != this->size; i++) {
        this->points[i] = (this->sgProj(std::array<double, 3>{mat(i, 0), mat(i, 1), mat(i, 2)}));
        point2idx.insert(std::make_pair(points[i], i));
    }
    h = Eigen::VectorXd::Zero(size);
    setHeight(h);
    Nu();
    constructDiagram(false);
}

void SphericalOMT::printPoly(CGAL::Polygon_2<K> pgon) {
    std::cout << "POLYGON ((";
    for (auto v = pgon.vertices_begin(); v != pgon.vertices_end(); v++)
        std::cout << v->x() << " " << v->y() << ", ";
    std::cout << pgon.vertices_begin()->x() << " " << pgon.vertices_begin()->y() << "))\n";
}

void SphericalOMT::printPolys() {
    for (auto p : polys)
        printPoly(p);
}

Eigen::SparseMatrix<double> SphericalOMT::hessByEdge()
{
    Eigen::SparseMatrix<double> hess(size, size);
    std::vector<Eigen::Triplet<double>> entries(size + vd.number_of_halfedges());
    size_t count = 0;
    for (int i = 0; i != size; i++) {
        entries[i] = Eigen::Triplet<double>(i, i, 0);
    }
    for (VD::Edge_iterator ei = vd.edges_begin(); ei != vd.edges_end(); ei++) {
        K::Weighted_point_2 pi = ei->up()->point();
        K::Weighted_point_2 pj = ei->down()->point();
        int i = point2idx[K::Point_2(pi.x(), pi.y())], j = point2idx[K::Point_2(pj.x(), pj.y())];
        double length = sqrt((double)((pi.x() - pj.x()) * (pi.x() - pj.x()) + (pi.y() - pj.y()) * (pi.y() - pj.y())));
        double entry = lineIntegral(ei) / length;
        entries.push_back(Eigen::Triplet<double>(i, j, -entry));
        entries.push_back(Eigen::Triplet<double>(j, i, -entry));
        entries.push_back(Eigen::Triplet<double>(i, i, entry));
        entries.push_back(Eigen::Triplet<double>(j, j, entry));
        //entries[size+count++] = Eigen::Triplet<double>(i, j, entry);
        //entries[size+count++] = Eigen::Triplet<double>(j, i, entry);
        //entries[i] = Eigen::Triplet<double>(i, i, entry + entries[i].value());
        //entries[j] = Eigen::Triplet<double>(j, j, entry + entries[j].value());
    }
    hess.setFromTriplets(entries.begin(), entries.end());
    return hess;
}

Eigen::SparseMatrix<double> SphericalOMT::hessian()
{
    // iterate over all halfedge
    Eigen::SparseMatrix<double> hess(size, size);// ->triplet
    std::vector<Eigen::Triplet<double>> entries;
    entries.resize(size);
    for (int i = 0; i != size; i++) {
        entries[i] = Eigen::Triplet<double>(i, i, 0);
    }
    //for (int i = 0; i!=this->size; )
    for (size_t i = 0; i != this->size; i++)
    {
        for (int j = 0; j != i; j++) {
            double entry_ij = hij(i, j);
            if (entry_ij) {
                //hess.insert(i, j) = entry_ij;    // symmetric
                //hess.insert(j, i) = entry_ij;
                //hess.coeffRef(i, i) += entry_ij; // set diagonal
                //hess.coeffRef(j, j) += entry_ij;
                entries.push_back(Eigen::Triplet<double>(i, j, entry_ij));
                entries.push_back(Eigen::Triplet<double>(j, i, entry_ij));
                //entries.push_back(Eigen::Triplet<double>(i, i, entry_ij));
                //entries.push_back(Eigen::Triplet<double>(j, j, entry_ij));
                entries[i] = Eigen::Triplet<double>(i, i, entry_ij + entries[i].value());
                entries[j] = Eigen::Triplet<double>(j, j, entry_ij + entries[j].value());
            }
        }
    }
    hess.setFromTriplets(entries.begin(), entries.end());
    //Eigen::SparseMatrix<double> hess();
    return hess;
}


double SphericalOMT::hij(int i, int j) {
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
                double length = sqrt( (double)((pij_x * pij_x+pij_y * pij_y)));
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


double SphericalOMT::lineIntegral(K::Segment_2 seg) {
    // perform a closed form for computing \int_{x_0}^{x_1} \mu(X, Y) \sqrt(1+k^2) dx where \mu(X,Y) = 4 / (1+X^2+Y^2)^2
    double x0 = (double)seg.source().x(), y0 = (double)seg.source().y(), x1 = (double)seg.target().x(), y1 = (double)seg.target().y();
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

double SphericalOMT::lineIntegral(VD::Edge_iterator ei) {
    K::Point_2 p1, p2;
    double x0, y0, x1, y1;
    if (ei->has_source() && ei->has_target()){
        p1 = ei->source()->point();
        p2 = ei->target()->point();
        x0 = (double)p1.x(), y0 = (double)p1.y(), x1 = (double)p2.x(), y1 = (double)p2.y();
    }
    else
    {
        double RAY_LENGTH = 1000;
        // cast to ray and found 
        const CGAL::Object seg_dual = vd.dual().dual(ei->dual());
        const K::Ray_2* dray = CGAL::object_cast<K::Ray_2>(&seg_dual);
        const auto& source = dray->source();
        x0 = (double)source.x();
        y0 = (double)source.y();
        const auto& dir = dray->direction();
        auto norm = sqrt((double)(dir.dx() * dir.dx() + dir.dy() * dir.dy()));
        x1 = x0 + RAY_LENGTH * (double)dir.dx() / norm;
        y1 = y0 + RAY_LENGTH * (double)dir.dy() / norm;
    }
    // perform a closed form for computing \int_{x_0}^{x_1} \mu(X, Y) \sqrt(1+k^2) dx where \mu(X,Y) = 4 / (1+X^2+Y^2)^2
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

double SphericalOMT::curvatureIntegral(K::Segment_2 seg) {
    // perform a closed form for computing \int_si k_g ds
    double x0 = (double)seg.source().x(), y0 = (double)seg.source().y(), x1 = (double)seg.target().x(), y1 = (double)seg.target().y();
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


double SphericalOMT::polygonArea(CGAL::Polygon_2<K> pgon) {
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


Eigen::VectorXd SphericalOMT::Wh() {
    Eigen::VectorXd wh(size);
    for (size_t i = 0; i != size; i++) {
        wh(i) = polygonArea(polys[i]);
    }
    return wh;
}


void SphericalOMT::Nu(){
    nu = Eigen::MatrixXd::Constant(size, 1, 4 * _PI / size);
};

void SphericalOMT::setNu(std::string filepath){
    //nu = Eigen::MatrixXd::Constant(size, 1, 4 * _PI / size);
    double area;
    std::ifstream myfile(filepath);
    int i = 0;
    if (myfile.is_open())
    {
        while ( myfile >> area )
            nu(i++, 1) = area * 4 * _PI;
    }
    //std::cout << nu.sum() << std::endl;
};

double SphericalOMT::step(double lbda) {
    //std::cout << "hello convex function\n";
    Eigen::VectorXd wh = Wh();
    Eigen::VectorXd grad = wh - nu;
    //std::cout << "sum of wh:" << wh.sum() << std::endl;
    Eigen::SparseMatrix<double> hess = hessByEdge();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(hess);
    Eigen::VectorXd dh = solver.solve(grad);
    setHeight(h - dh * lbda);
    while (!constructDiagram(false)) { // line search
        //std::cout << "degenerate at lambda=" << lbda << std::endl;
        lbda = lbda * 0.5;
        setHeight(h - dh * lbda);
    };
    h = h - dh * lbda;
    //std::cout << "Wh:\n" << wh << std::endl;
    //std::cout << "dh\n" << dh << std::endl;
    //std::cout << "h:\n" << h << std::endl;
    std::cout << "max error:" << grad.cwiseAbs().maxCoeff() << std::endl;
    //std::cout << rt.number_of_hidden_vertices() << std::endl;
    return grad.cwiseAbs().maxCoeff();
}


Eigen::MatrixXi SphericalOMT::faceMat() {
    // will do the spherical triangulation here 
    Eigen::MatrixXi faces(vd.number_of_vertices(), 3);
    int i = 0;
    for (RT::Face_iterator f = rt.finite_faces_begin(); f != rt.finite_faces_end(); f++) {
        faces(i, 0) = point2idx[K::Point_2(f->vertex(0)->point().x(), f->vertex(0)->point().y())];
        faces(i, 1) = point2idx[K::Point_2(f->vertex(1)->point().x(), f->vertex(1)->point().y())];
        faces(i, 2) = point2idx[K::Point_2(f->vertex(2)->point().x(), f->vertex(2)->point().y())];
        i++;
        //std::cout << f->vertex(0)->point() << " | " << f->vertex(1)->point() << " | " << f->vertex(2)->point() << " | " << std::endl;
    }
    return faces;
}


Eigen::MatrixXd SphericalOMT::mapBack() {
    Eigen::MatrixXd vertices_positions(size, 3);
    for (size_t i = 0; i != size; i++) {
        CGAL::Polygon_2<K> pgon = polys[i];
        double total_x = 0, total_y = 0, total_z = 0;
        // average vertices center
        for (CGAL::Polygon_2<K>::Vertex_iterator vi = pgon.vertices_begin(); vi != pgon.vertices_end(); ++vi) 
        {
            std::array<double, 3>  coordinates = isgProj(*vi);
            total_x += coordinates[0];
            total_y += coordinates[1];
            total_z += coordinates[2];
        }
        double xc = total_x / pgon.size(), yc = total_y / pgon.size(), zc = total_z / pgon.size();
        double length = sqrt(xc * xc + yc * yc + zc * zc);
        // centroid
        //double cx=0, cy=0, cz = 0;
        //double Axy = 0;
        //double Ayz = 0;
        //for (CGAL::Polygon_2<K>::Edge_const_iterator edge = pgon.edges_begin(); edge != pgon.edges_end(); edge++)
        //{
        //    // compute inner product
        //    const K::Segment_2 seg = (K::Segment_2)(*edge);
        //    std::array<double, 3> p1 = isgProj(seg.source());
        //    std::array<double, 3> p2 = isgProj(seg.target());
        //    double axy = (p1[0] * p2[1] - p2[0] * p1[1]);
        //    double ayz = (p1[1] * p2[2] - p2[1] * p1[2]);
        //    cx += (p1[0] + p2[0]) * axy;
        //    cy += (p1[1] + p2[1]) * axy;
        //    cz += (p1[2] + p2[2]) * ayz;
        //    Axy += 0.5 * axy;
        //    Ayz += 0.5 * ayz;
        //}
        //cx /= (6 * Axy);
        //cy /= (6 * Axy);
        //cz /= (6 * Ayz);
        //double length = sqrt(cx * cx + cy * cy + cz * cz);
        vertices_positions(i, 0) = xc / length;
        vertices_positions(i, 1) = yc / length;
        vertices_positions(i, 2) = zc / length;
        //std::cout << i << std::endl;
        //std::cout << "average vertice:\t" << xc / lengthc << " " << yc / lengthc << " " << zc / lengthc << std::endl;
        //std::cout << "centroid\t\t:" << vertices_positions(i, 0) << " " << vertices_positions(i, 1) << " " << vertices_positions(i, 2) << std::endl;
    }
    return vertices_positions;
};

void SphericalOMT::printDistortion() {
    //Eigen::MatrixXd vertices_positions = mapBack();
    //const auto triArea = [&](double v1, CGAL::Vector_3<KS2> v2) -> double { // cross product face area
    //    double cross_x = (double)(v1.y() * v2.z() - v2.y() * v1.z());
    //    double cross_y = (double)(v1.x() * v2.z() - v2.x() * v1.z());
    //    double cross_z = (double)(v1.x() * v2.y() - v2.x() * v1.y());
    //    return sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z) * 0.5;
    //};
    //int i = 0;
    std::cout << "np.array([";
    Eigen::VectorXd cell_areas = Wh();
    for (int i = 0; i != size; i++) {
        std::cout << log10(cell_areas(i) / nu(i)) <<",";
    }
    std::cout << "])" << std::endl;
};


void SphericalOMT::measureDual(Eigen::MatrixXd vertices, Eigen::MatrixXd faces) {
    auto faceArea = [&](int i) {
        int v1 = faces(i, 0), v2 = faces(i, 1), v3 = faces(i, 2);
        double a1 = vertices(v2, 0) - vertices(v1, 0), a2 = vertices(v2, 1) - vertices(v1, 1), a3 = vertices(v2, 2) - vertices(v1, 2);
        double b1 = vertices(v3, 0) - vertices(v1, 0), b2 = vertices(v3, 1) - vertices(v1, 1), b3 = vertices(v3, 2) - vertices(v1, 2);
        double result = sqrt((a2 * b3 - a3 * b2) * (a2 * b3 - a3 * b2) + (a1 * b3 - a3 * b1) * (a1 * b3 - a3 * b1) + (a1 * b2 - a2 * b1) * (a1 * b2 - a2 * b1));
        return result/2;
    };
    Eigen::ArrayXd dual_area(vertices.rows());
    dual_area.setZero();
    for (int i = 0; i != faces.rows(); i++)
    {
        int v1 = faces(i, 0), v2 = faces(i, 1), v3 = faces(i, 2);
        double area = faceArea(i)/3;
        dual_area(v1) += area;
        dual_area(v2) += area;
        dual_area(v3) += area;
    }
    Eigen::ArrayXd area_distortion = dual_area / nu.array();
    std::cout << dual_area.sum() << std::endl;
    std::cout << "local distortion mean:" << area_distortion.mean() << std::endl;
    std::cout << "local distortion std:" << sqrt((area_distortion - area_distortion.mean()).square().sum() / area_distortion.size()) << std::endl;
    std::cout << "relative error:" << (dual_area - nu.array()).abs().sum() / 4 / _PI << std::endl;
    std::ofstream myfile("distortion.txt");
    myfile << area_distortion << std::endl;
    myfile.close();
};

// some utils function
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> readFile(std::string file)
{
    Mesh mesh;
    if (!CGAL::IO::read_OBJ(file, mesh))
    {
        std::cout << "File reading error.";
    }
    Eigen::MatrixXd vertex_positions(mesh.number_of_vertices(), 3); 
    Eigen::MatrixXd faces(mesh.number_of_faces(), 3); // we assume that surfacemesh is a triangular mesh 
    for (auto f : mesh.faces()) {
        //std::cout << "vertices around face:" << f << std::endl;
        int i_v = 0;
        for (auto v : mesh.vertices_around_face(mesh.halfedge(f))){
            //std::cout << v << std::endl;
            faces(f.idx(), i_v) = v.idx();
            i_v++;
        }
    }
    for (auto v : mesh.vertices()) { // get positions
        const auto& point = mesh.point(v);
        vertex_positions(v.idx(), 0) = point.x();
        vertex_positions(v.idx(), 1) = point.y();
        vertex_positions(v.idx(), 2) = point.z();
    }
    return std::make_pair(vertex_positions, faces);
}

std::pair<Eigen::MatrixXd, Eigen::VectorXd> SphereTriangulation(Eigen::MatrixXd vertex_positions) {
    TraitsS2 traits(Point_3(0, 0, 0), 1); // sphere center on (1,1,1), with radius 1
    DToS2 dtos(traits);
    std::unordered_map<Point_3, size_t, HASH_P3> point2idx;
    for (size_t i = 0; i != vertex_positions.rows(); i++) {
        Point_3 p3(vertex_positions(i, 0), vertex_positions(i, 1), vertex_positions(i, 2));
        dtos.insert(p3);
        point2idx.insert(std::make_pair(p3, i));
    }
    Eigen::MatrixXd faces(dtos.number_of_faces(), 3);
    Eigen::VectorXd face_areas(dtos.number_of_faces());
    const auto triArea = [&](CGAL::Vector_3<KS2> v1, CGAL::Vector_3<KS2> v2) -> double { // cross product face area
        double cross_x = (double)(v1.y() * v2.z() - v2.y() * v1.z());
        double cross_y = (double)(v1.x() * v2.z() - v2.x() * v1.z());
        double cross_z = (double)(v1.x() * v2.y() - v2.x() * v1.y());
        return sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z) * 0.5;
    };
    int i = 0;

    for (auto f = dtos.finite_faces_begin(); f != dtos.finite_faces_end(); f++) {
        faces(i, 0) = point2idx[Point_3(f->vertex(0)->point().x(), f->vertex(0)->point().y(), f->vertex(0)->point().z())];
        faces(i, 1) = point2idx[Point_3(f->vertex(1)->point().x(), f->vertex(1)->point().y(), f->vertex(1)->point().z())];
        faces(i, 2) = point2idx[Point_3(f->vertex(2)->point().x(), f->vertex(2)->point().y(), f->vertex(2)->point().z())];
        CGAL::Vector_3<KS2> v1 = f->vertex(1)->point() - f->vertex(0)->point(), v2 = f->vertex(2)->point() - f->vertex(0)->point();
        face_areas(i) = triArea(v1, v2);
        i++;
        //std::cout << f->vertex(0)->point() << " | " << f->vertex(1)->point() << " | " << f->vertex(2)->point() << " | " << std::endl;
    }
    return std::make_pair(faces, face_areas);
}