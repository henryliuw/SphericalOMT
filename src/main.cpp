
#include <vector>
#include <array>
#include <iostream>
#include <stdlib.h>
#include "SphericalOMT.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "imgui.h"


polyscope::SurfaceMesh* psMesh;
int epochs;
double lambda;
int count = 0;
SphericalOMT* omt_ptr;
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> mesh_pair;

void printArea(Eigen::VectorXd areas) {
    std::cout << "np.array([";
    for (int i = 0; i != areas.rows(); i++) {
        std::cout << areas(i) << ",";
    }
    std::cout << "])\n";
}

void functionCallback() { //IMGUi
    ImGui::InputDouble("lambda", &lambda, 0.1, 1);
    ImGui::InputInt("flow epochs", &epochs);
    if (ImGui::Button("Step")) {
        for (int i = 0; i != epochs; i++) {
            omt_ptr->step(lambda);
            count++;
            if (count % 10==0) {
                std::cout << "epoch: " << count << std::endl;
            }
        }
        auto result = omt_ptr->mapBack();
        psMesh->updateVertexPositions(result);
        polyscope::requestRedraw();
        omt_ptr->measureDual(result, mesh_pair.second);
    }
    if (ImGui::Button("print area distortion")) {
        //auto result = SphereTriangulation(omt_ptr->mapBack());
        //polyscope::registerSurfaceMesh("new connectivity", omt_ptr->mapBack(), result.first);
        //polyscope::requestRedraw();
        //printArea(result.second);
        omt_ptr->printDistortion();
    }
    if (ImGui::Button("map new")) {
        auto result = SphereTriangulation(omt_ptr->mapBack());
        polyscope::registerSurfaceMesh("new connectivity", omt_ptr->mapBack(), result.first);
        polyscope::requestRedraw();
        //printArea(result.second);
    }
}

int main(int argc, char** argv) {
    // easy test 
    //std::vector<std::array<double, 3>> vertex_positions = { {0, 0, -1}, {0, 1, 0}, {0.444444, 0.444444, 0.777778}, {0.666667, -0.666667, 0.333333}, {0.833333, 0.555555, -0.1} };
    //SphericalOMT omt(vertex_positions);

    // random mesh test
    //const int n = 10000;
    //Eigen::MatrixXd angle_mat = Eigen::MatrixXd::Random(n, 2);
    //Eigen::MatrixXd sph_mat(n, 3);
    //for (int i = 0; i != n; i++) {
    //    double theta = (angle_mat(i, 0)+1) * _PI, z = angle_mat(i, 1);
    //    sph_mat(i, 0) = sqrt(1-z*z) * cos(theta);
    //    sph_mat(i, 1) = sqrt(1-z*z) * sin(theta);
    //    sph_mat(i, 2) = z;
    //}
    //SphericalOMT omt(sph_mat);

    // bunny test
    mesh_pair = readFile("../input/bunny_conformal.obj");
    SphericalOMT omt(mesh_pair.first);
    omt.setNu("../input/bunny_area.txt");
    // set global ptr
    omt_ptr = &omt;

    //omt.debug();
    // 3D viewer
    polyscope::init();
    polyscope::state::userCallback = functionCallback;
    //polyscope::registerSurfaceMesh("original mesh", sph_mat, SphereTriangulation(sph_mat).first);
    polyscope::registerSurfaceMesh("original mesh", mesh_pair.first,  mesh_pair.second);
    //for (int i = 0; i != 10; i++)
    //    omt_ptr->step(1);
    psMesh = polyscope::registerSurfaceMesh("new mesh", mesh_pair.first,  mesh_pair.second);
    //psMesh = polyscope::registerSurfaceMesh("new mesh", omt.mapBack(), SphereTriangulation(sph_mat).first);
    polyscope::show();

    return 0;
}