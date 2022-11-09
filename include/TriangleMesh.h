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
typedef std::tuple<size_t, size_t, size_t> key_f;
typedef std::tuple<size_t, size_t> key_e;
typedef Eigen::Matrix< size_t, Eigen::Dynamic, 1> Vector;
typedef Eigen::Matrix< size_t, Eigen::Dynamic, Eigen::Dynamic> Matrix;
class TriangleMesh {
public:
    TriangleMesh(Eigen::MatrixXi &T);
    TriangleMesh(Matrix &T);
    TriangleMesh();
    void initialize(Eigen::MatrixXi &T);
    void initialize(Matrix &T);
    void create_edges();
    void create_edges_from_faces();
    void create_faces();
    void build_boundary_mat3(); // T -> F, size: |F|x|T|, boundary of tets
    void build_boundary_mat2(); // F -> E, size: |E|x|F|, boundary of triangles
    void build_boundary_mat1(); // E -> V, size: |V|x|E|, boundary of edges
    // mesh API
    // v as input
    std::set<size_t> get_adjacent_vertices_v(size_t vindex); 
    std::set<size_t> get_incident_edges_v(size_t vindex); 
    std::set<size_t> get_incident_faces_v(size_t vindex); 
    // e as input
    std::set<size_t> get_incident_vertices_e(size_t eindex); 
    std::set<size_t> get_incident_faces_e(size_t eindex); 
    std::set<size_t> get_diamond_vertices_e(size_t eindex);
    // f as input
    std::set<size_t> get_incident_vertices_f(size_t findex); 
    std::set<size_t> get_incident_edges_f(size_t findex); 
    std::set<size_t> get_adjacent_faces_f(size_t findex); 
    std::set<size_t> get_adjacent_faces_f2(size_t findex); 
    // simplicial complex
    Eigen::VectorXi build_vertex_vector(const SimplexSubset& subset) const;
    Eigen::VectorXi build_edge_vector(const SimplexSubset& subset) const;
    Eigen::VectorXi build_face_vector(const SimplexSubset& subset) const;
    Eigen::VectorXi build_vertex_vector(const std::set<size_t>& vset) const;
    Eigen::VectorXi build_edge_vector(const std::set<size_t>& eset) const;
    Eigen::VectorXi build_face_vector(const std::set<size_t>& fset) const;
    SimplexSubset star(const SimplexSubset& subset) const;
    SimplexSubset closure(const SimplexSubset& subset) const;
    SimplexSubset link(const SimplexSubset& subset) const;
    bool is_complex(const SimplexSubset& subset) const;
    size_t is_pure_complex(const SimplexSubset& subset) const;
    SimplexSubset boundary(const SimplexSubset& subset) const;
    //
    int n_edges() const;
    int n_vertices() const;
    int n_faces() const;
    int n_tets() const;
    size_t get_edge_index(size_t i, size_t j, int &sign);
    size_t get_face_index(size_t i, size_t j, size_t k, int &sign);
// private:
    size_t num_v;
    Matrix T;
    Matrix E;
    Matrix F;
    std::map<key_f, size_t> map_f; // tuple -> face index
    std::map<key_e, size_t> map_e; // tuple -> edge index
    Eigen::SparseMatrix<int> bm1; // |V|x|E|
    Eigen::SparseMatrix<int> pos_bm1; // |V|x|E|
    Eigen::SparseMatrix<int> bm2; // |E|x|F|
    Eigen::SparseMatrix<int> pos_bm2; // |E|x|F|
    Eigen::SparseMatrix<int> bm3; // |F|x|T|
    Eigen::SparseMatrix<int> pos_bm3; // |F|x|T|
};

#endif /* TriangleMesh_h */
