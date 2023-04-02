
#include <vector>
#include <array>
#include <iostream>
#include <stdlib.h>
#include "SphericalOMT.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "imgui.h"
#include "happly.h"

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

void readFiles(std::string filename, std::vector<std::array<double, 3>>& vPos, std::vector<std::array<double, 3>>& pPos, std::vector<std::vector<size_t>> & fInd)
{
    happly::PLYData plyIn(filename);

    vPos = plyIn.getVertexPositions();
    fInd = plyIn.getFaceIndices<size_t>();

    std::vector<double> px = plyIn.getElement("vertex").getProperty<double>("px");
    std::vector<double> py = plyIn.getElement("vertex").getProperty<double>("py");
    std::vector<double> pz = plyIn.getElement("vertex").getProperty<double>("pz");
    pPos.resize(vPos.size());
    for (int i = 0; i != pPos.size(); i++)
    {
        pPos[i][0] = px[i];
        pPos[i][1] = py[i]; // translate to view easier
        pPos[i][2] = pz[i];
    }
}

void writeFiles(std::string filenamein, std::string filenameout, Eigen::MatrixXd & pPos)
{
    happly::PLYData plyIn(filenamein);
    happly::PLYData plyOut;
    plyOut.addFaceIndices(plyIn.getFaceIndices<size_t>());
    plyOut.addVertexPositions(plyIn.getVertexPositions());
    plyOut.getElement("vertex").addProperty<unsigned char>("red", plyIn.getElement("vertex").getProperty<unsigned char>("red"));
    plyOut.getElement("vertex").addProperty<unsigned char>("green", plyIn.getElement("vertex").getProperty<unsigned char>("green"));
    plyOut.getElement("vertex").addProperty<unsigned char>("blue", plyIn.getElement("vertex").getProperty<unsigned char>("blue"));
    std::vector<double> px, py, pz;
    for (int i = 0; i != pPos.rows(); i++) {
        px.push_back(pPos(i, 0));
        py.push_back(pPos(i, 1));
        pz.push_back(pPos(i, 2));
    }
    plyOut.getElement("vertex").addProperty<double>("px", px);
    plyOut.getElement("vertex").addProperty<double>("py", py);
    plyOut.getElement("vertex").addProperty<double>("pz", pz);
    plyOut.write(filenameout);
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

    // read images
    std::string filein;
    if (argc>=2)
        filein = argv[1];
    std::vector<std::array<double, 3>> vPos, pPos;
    std::vector<std::vector<size_t>> fInd;
    readFiles(filein, vPos, pPos, fInd);

    // bunny test
    //mesh_pair = readFile("../input/bunny_conformal.obj");
    SphericalOMT omt(pPos);
    omt.setNu(vPos, fInd);
    // set global ptr
    omt_ptr = &omt;
    
    double eps = 1e-3;
    if (argc >= 4)
    {
        std::string temp = argv[3];
        eps = stod(temp);
    }

    // measure time and save
    auto result = omt_ptr->computeAll(eps);
    if (argc >= 3)
    {
        std::string fileout = argv[2];
        writeFiles(filein, fileout, result);
    }
    
    //omt.debug();
    // 3D viewer
    polyscope::init();
    polyscope::state::userCallback = functionCallback;
    //polyscope::registerSurfaceMesh("original mesh", sph_mat, SphereTriangulation(sph_mat).first);
    polyscope::registerSurfaceMesh("original mesh", pPos,  fInd);
    //for (int i = 0; i != 10; i++)
    //    omt_ptr->step(1);
    psMesh = polyscope::registerSurfaceMesh("new mesh", result,  fInd);
    //psMesh = polyscope::registerSurfaceMesh("new mesh", omt.mapBack(), SphereTriangulation(sph_mat).first);
    polyscope::show();
    
    return 0;
}