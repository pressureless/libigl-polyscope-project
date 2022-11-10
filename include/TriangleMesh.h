//
//  TriangleMesh.h
//  DEC
//
//  Created by pressure on 10/31/22.
//

#ifndef TriangleMesh_h
#define TriangleMesh_h
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <map>
#include "simplex_subset.h"
// using Eigen::Matrix;
// using Eigen::VectorXi;
using Eigen::SparseMatrix;
typedef std::tuple<int, int, int> key_f;
typedef std::tuple<int, int> key_e;
typedef Eigen::Matrix< int, Eigen::Dynamic, 1> Vector;
typedef Eigen::Matrix< int, 1, Eigen::Dynamic> RowVector;
typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic> Matrix;
class TriangleMesh {
public:
    TriangleMesh(Eigen::MatrixXi &T);
    TriangleMesh();
    void initialize(Eigen::MatrixXi &T);
    void create_edges();
    void create_edges_from_faces();
    void create_faces();
    void build_boundary_mat3(); // T -> F, size: |F|x|T|, boundary of tets
    void build_boundary_mat2(); // F -> E, size: |E|x|F|, boundary of triangles
    void build_boundary_mat1(); // E -> V, size: |V|x|E|, boundary of edges
    // mesh API
    // v as input
    std::set<int> get_adjacent_vertices_v(int vindex); 
    std::set<int> get_incident_edges_v(int vindex); 
    std::set<int> get_incident_faces_v(int vindex); 
    // e as input
    std::set<int> get_incident_vertices_e(int eindex); 
    std::set<int> get_incident_faces_e(int eindex); 
    std::set<int> get_diamond_vertices_e(int eindex);
    std::tuple< int, int > get_diamond_vertices_e(int start, int end);
    // f as input
    std::set<int> get_incident_vertices_f(int findex); 
    std::set<int> get_incident_edges_f(int findex); 
    std::set<int> get_adjacent_faces_f(int findex); 
    std::set<int> get_adjacent_faces_f2(int findex);
    //
    int get_opposite_vertex(const RowVector& f, int start, int end);
    // simplicial complex
    Eigen::VectorXi build_vertex_vector(const SimplexSubset& subset) const;
    Eigen::VectorXi build_edge_vector(const SimplexSubset& subset) const;
    Eigen::VectorXi build_face_vector(const SimplexSubset& subset) const;
    Eigen::VectorXi build_vertex_vector(const std::set<int>& vset) const;
    Eigen::VectorXi build_edge_vector(const std::set<int>& eset) const;
    Eigen::VectorXi build_face_vector(const std::set<int>& fset) const;
    SimplexSubset star(const SimplexSubset& subset) const;
    SimplexSubset closure(const SimplexSubset& subset) const;
    SimplexSubset link(const SimplexSubset& subset) const;
    bool is_complex(const SimplexSubset& subset) const;
    int is_pure_complex(const SimplexSubset& subset) const;
    SimplexSubset boundary(const SimplexSubset& subset) const;
    //
    int n_edges() const;
    int n_vertices() const;
    int n_faces() const;
    int n_tets() const;
    int get_edge_index(int i, int j, int &sign);
    int get_face_index(int i, int j, int k, int &sign);
// private:
    int num_v;
    Matrix T;
    Matrix E;
    Matrix F;
    std::map<key_f, int> map_f; // tuple -> face index
    std::map<key_e, int> map_e; // tuple -> edge index
    Eigen::SparseMatrix<int> bm1; // |V|x|E|
    Eigen::SparseMatrix<int> pos_bm1; // |V|x|E|
    Eigen::SparseMatrix<int> bm2; // |E|x|F|
    Eigen::SparseMatrix<int> pos_bm2; // |E|x|F|
    Eigen::SparseMatrix<int> bm3; // |F|x|T|
    Eigen::SparseMatrix<int> pos_bm3; // |F|x|T|
};

#endif /* TriangleMesh_h */
