// SIMPLICIAL COMPLEX

// #include "geometrycentral/surface/manifold_surface_mesh.h"
// #include "geometrycentral/surface/meshio.h"
// #include "geometrycentral/surface/vertex_position_geometry.h"
#include <igl/PI.h>
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/boundary_loop.h>
#include <igl/exact_geodesic.h>
#include <igl/gaussian_curvature.h>
#include <igl/invert_diag.h>
#include <igl/lscm.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/readOBJ.h>
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "dec_util.h"
#include "DECOperators.h"
#include <algorithm>

// using namespace geometrycentral;
// using namespace geometrycentral::surface;
Eigen::MatrixXd meshV;
Eigen::MatrixXi meshF;
std::vector<size_t> edgeIndices;
std::vector<size_t> revEdgeIndices;
// == Geometry-central data
// std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
// std::unique_ptr<VertexPositionGeometry> geometry_uptr;
// so we can more easily pass these to different classes while preserving syntax
// ManifoldSurfaceMesh* mesh;
// VertexPositionGeometry* geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;
std::string MESHNAME;

// Some global variables
DECOperators SCO;
bool isComplexResult = false;
int isPureComplexResult = -1;
double vertexRadius;
double edgeRadius;


std::array<double, 3> BLUE = {0.11, 0.388, 0.89};
glm::vec<3, float> ORANGE_VEC = {1, 0.65, 0};
std::array<double, 3> ORANGE = {1, 0.65, 0};


void flipZ() {
    // Rotate mesh 180 deg about up-axis on startup
    glm::mat4x4 rot = glm::rotate(glm::mat4x4(1.0f), static_cast<float>(igl::PI), glm::vec3(0, 1, 0));
    for (int v=0; v < SCO.nVertices(); v++) {
        auto vec = meshV.row(v);
        glm::vec4 rvec = {vec[0], vec[1], vec[2], 1.0};
        rvec = rot * rvec;
        meshV(v, 0) = rvec[0]; 
        meshV(v, 1) = rvec[1]; 
        meshV(v, 2) = rvec[2]; 
    }
    psMesh->updateVertexPositions(meshV);
}

/*
 * Display the selected simplices.
 * TODO: Use SurfaceVertexCountQuantity* SurfaceMesh::addVertexCountQuantity, etc. instead of SurfaceGraphQuantity for
 * cleaner code
 */
void showSelected() {

    print_vec(psMesh->edgePerm);
    // Show selected vertices. 
    Eigen::MatrixXd vertPos(polyscope::state::subset.vertices.size(), 3);
    std::vector<std::array<size_t, 2>> vertInd;
    int idx = 0;
    for (std::set<size_t>::iterator it = polyscope::state::subset.vertices.begin();
        it != polyscope::state::subset.vertices.end(); ++it) {
        int cur = *it;
        vertPos.row(idx) = meshV.row(cur);
        idx++; 
    }
    polyscope::SurfaceGraphQuantity* showVerts = psMesh->addSurfaceGraphQuantity("selected vertices", vertPos, vertInd);
    showVerts->setEnabled(true);
    showVerts->setRadius(vertexRadius);
    showVerts->setColor(ORANGE_VEC);

    // Show selected edges.
    // std::vector<Vector3> edgePos;
    idx = 0;
    Eigen::MatrixXd edgePos(polyscope::state::subset.edges.size()*2, 3);
    std::vector<std::array<size_t, 2>> edgeInd;
    for (std::set<size_t>::iterator it = polyscope::state::subset.edges.begin();
         it != polyscope::state::subset.edges.end(); ++it) {
        int cur = *it;
        int fir = SCO.E(cur, 0);
        int sec = SCO.E(cur, 1); 
        std::cout<<"fir:"<<fir<<", sec:"<<sec<<std::endl;
        edgePos.row(idx) = meshV.row(fir); 
        edgePos.row(idx+1) = meshV.row(sec);  
        size_t i = idx;
        edgeInd.push_back({i, i + 1});
        idx += 2; 
    }
    polyscope::SurfaceGraphQuantity* showEdges = psMesh->addSurfaceGraphQuantity("selected edges", edgePos, edgeInd);
    showEdges->setEnabled(true);
    showEdges->setRadius(edgeRadius);
    showEdges->setColor(ORANGE_VEC);

//     // Show selected faces.
    std::vector<std::array<double, 3>> faceColors(SCO.nFaces());
    for (size_t i = 0; i < SCO.nFaces(); i++) {
        faceColors[i] = BLUE;
    }
    for (std::set<size_t>::iterator it = polyscope::state::subset.faces.begin();
         it != polyscope::state::subset.faces.end(); ++it) {
        faceColors[*it] = ORANGE;
    }
    polyscope::SurfaceFaceColorQuantity* showFaces = psMesh->addFaceColorQuantity("selected faces", faceColors);
    showFaces->setEnabled(true);
}

void redraw() {
    showSelected();
    polyscope::requestRedraw();
}

std::vector<size_t> getMapping(){
    return psMesh->edgePerm;
}

/*
 * Buttons for the simplicial operators.
 */
void functionCallback() {

    if (ImGui::Button("Reset")) {
        polyscope::state::subset.vertices.clear();
        polyscope::state::subset.edges.clear();
        polyscope::state::subset.faces.clear();
        redraw();
    }

    if (ImGui::Button("isComplex")) {
        isComplexResult = SCO.isComplex(polyscope::state::subset);
    }
    ImGui::SameLine(100);
    ImGui::Text(isComplexResult ? "true" : "false");

    if (ImGui::Button("isPureComplex")) {
        isPureComplexResult = SCO.isPureComplex(polyscope::state::subset);
    }
    ImGui::SameLine(130);
    ImGui::Text("%d", isPureComplexResult);

    if (ImGui::Button("Boundary")) {
        SimplexSubset S = SCO.boundary(polyscope::state::subset);
        polyscope::state::subset = S;
        redraw();
    }

    if (ImGui::Button("Star")) {
        SimplexSubset S = SCO.star(polyscope::state::subset);
        polyscope::state::subset = S;
        redraw();
    }

    if (ImGui::Button("Closure")) {
        SimplexSubset S = SCO.closure(polyscope::state::subset);
        polyscope::state::subset = S;
        redraw();
    }

    if (ImGui::Button("Link")) {
        SimplexSubset S = SCO.link(polyscope::state::subset);
        polyscope::state::subset = S;
        redraw();
    }

    if (ImGui::Button("Face color")) {
        std::vector<std::array<double, 3>> fColor(meshF.rows());
        for (size_t iF = 0; iF < meshF.rows(); iF++) { 
          fColor[iF] = {{polyscope::randomUnit(), polyscope::randomUnit(), polyscope::randomUnit()}};
          std::cout<<"fColor[iF]: "<<fColor[iF][0]<<fColor[iF][1]<<fColor[iF][2]<<std::endl;
        } 
        // Visualize
        polyscope::getSurfaceMesh(MESHNAME)->addFaceColorQuantity("fColor", fColor);
        redraw();
    }
}


int main(int argc, char** argv) {
    // Configure the argument parser
    args::ArgumentParser parser("15-458 HW0");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    // If a mesh name was not given, use default mesh.
    std::string filepath = "../../../input/small_disk.obj";
    if (inputFilename) {
        filepath = args::get(inputFilename);
    }

    MESHNAME = polyscope::guessNiceNameFromPath(filepath);
    igl::readOBJ(filepath, meshV, meshF);
    std::cout<<"meshV:\n"<<meshV<<std::endl;
    std::cout<<"meshF:\n"<<meshF<<std::endl;
    // 
    Eigen::Matrix< size_t, Eigen::Dynamic, Eigen::Dynamic > F = meshF.cast<size_t>();
    Eigen::Matrix< size_t, Eigen::Dynamic, Eigen::Dynamic >  FF = preprocess_matrix(F);
    meshF = FF.cast<int>();
    SCO.initialize(F);
    SCO.vertices_of_diamond(0);
    //
    // Load mesh
    // std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    // mesh = mesh_uptr.release();
    // geometry = geometry_uptr.release();

    // Get indices for element picking
    polyscope::state::facePickIndStart = meshV.rows();
    polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + meshF.rows();
    polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + SCO.nEdges();

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = functionCallback;
    
    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(MESHNAME, meshV, SCO.F);
    // Set edge indices
    typedef std::tuple<size_t, size_t> key_e;
    std::set<key_e> e_set;
    std::vector<size_t> edgeFirst;
    std::vector<size_t> edgeSecond;
    for (int i = 0; i < SCO.nEdges(); ++i){
        edgeFirst.push_back(SCO.E(i, 0));
        edgeSecond.push_back(SCO.E(i, 1));
    }
    int cnt =0;
    for (int i = 0; i < SCO.nFaces(); ++i)
    {
        for(int j = 0; j < 3; j++) {
            size_t vertex_A = SCO.F(i, j);
            size_t vertex_B = SCO.F(i, (j+1) % 3 );
            size_t min = std::min({vertex_A, vertex_B});
            size_t max = std::max({vertex_A, vertex_B});
            key_e cur_key = std::make_tuple(min, max);
            auto search = e_set.find(cur_key);
            if(search != e_set.end()){
                // found
                continue;
            }
            else{
                e_set.insert(cur_key);
                std::cout<<"edge indices: "<<cnt<<", s:"<<min<<", e:"<<max<<", index:"<<SCO.map_e[cur_key]<<std::endl;
                // edgeIndices.push_back(SCO.map_e[cur_key]); 
                edgeIndices.push_back(SCO.map_e[cur_key]); 
                cnt++;
            }
        } 
    }  
    revEdgeIndices.resize(edgeIndices.size());
    for (int i = 0; i < edgeIndices.size(); ++i){
        revEdgeIndices[edgeIndices[i]] = i;
    }
    std::cout<<"edge indices 1: "<<edgeIndices.size()<<std::endl;
    for (int i = 0; i < edgeIndices.size(); ++i) {
        std::cout << edgeIndices[i] << ", ";
    }
    std::cout<<std::endl;
    std::cout<<"edge indices 2: "<<edgeIndices.size()<<std::endl;
    psMesh->setEdgePermutation(edgeIndices, edgeIndices.size());
    print_vec(psMesh->edgePerm);
    // psMesh->buildEdgeInfoGui();

    // psMesh->setEdgePermutation(edgeFirst, edgeSecond);

    // Mesh initialization
    // Add visualization options.
    flipZ();
    double sum = 0;
    for (int i = 0; i < SCO.nEdges(); ++i)
    {
        sum += (meshV.row(SCO.E(i, 0)) - meshV.row(SCO.E(i, 1))).norm();
    }
    sum /= SCO.nEdges();
    double lengthScale = sum;
    vertexRadius = lengthScale * 0.1;
    edgeRadius = lengthScale * 0.05;

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}