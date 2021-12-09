
#include <vector>
#include <array>
#include <iostream>
#include <stdlib.h>
#include "SphericalOMT.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

int main(int argc, char** argv) {

    std::vector<std::array<double, 3>> vertex_positions;
    vertex_positions.push_back(std::array<double, 3> { 0, 0, -1 });
    vertex_positions.push_back(std::array<double, 3> { 0, 1, 0 });
    vertex_positions.push_back(std::array<double, 3> { 0.444444, 0.444444, 0.777778 });
    vertex_positions.push_back(std::array<double, 3> { 0.666667, -0.666667, 0.333333 });
    vertex_positions.push_back(std::array<double, 3> { 0.827586, 0.551724, -0.103448 });

    SphericalOMT<5> omt(vertex_positions);
    for (int i = 0; i != 10; i++){
        std::cout << "------------------step " << i << "------------------\n";
        omt.step(0.5);
    }
    omt.printPolys();
    //omt.debug();
    //polyscope::init();
    //Eigen::MatrixXd meshV{
    //    {0, 0, 0},
    //    {1, 0, 0},
    //    {0, 0, 1},
    //    {1, 1, 1},
    //};
    //Eigen::MatrixXi meshF{
    //    {0, 1, 2},
    //    {1, 2, 3},
    //};
    //Eigen::MatrixXd meshV2{
    //    {7, 0, 0},
    //    {6, 0, -1},
    //    {6, 0, 1},
    //    {6, 1, 1},
    //};
    //Eigen::MatrixXi meshF2{
    //    {0, 1, 2, 3},
    //};
    //polyscope::registerSurfaceMesh("mesh1", meshV, meshF);
    //polyscope::registerSurfaceMesh("mesh2", meshV2, meshF2);
    //polyscope::show();
    

    return 0;
}