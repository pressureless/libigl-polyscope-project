//
//  DECOperators.h
//  DEC
//
//  Created by pressure on 10/31/22.
//

#ifndef DECOperators_h
#define DECOperators_h
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
class DECOperators {
public:
    DECOperators(Matrix &T);
    DECOperators();
    void initialize(Matrix &T);
    void create_edges();
    void create_edges_from_faces();
    void create_faces();
    void build_boundary_mat3(); // T -> F, size: |F|x|T|, boundary of tets
    void build_boundary_mat2(); // F -> E, size: |E|x|F|, boundary of triangles
    void build_boundary_mat1(); // E -> V, size: |V|x|E|, boundary of edges
    // mesh API
    std::set<size_t> one_ring_vertices(size_t vindex);
    std::set<size_t> vertices_of_diamond(size_t eindex);
    // simplicial complex
    Eigen::VectorXi buildVertexVector(const SimplexSubset& subset) const;
    Eigen::VectorXi buildEdgeVector(const SimplexSubset& subset) const;
    Eigen::VectorXi buildFaceVector(const SimplexSubset& subset) const;
    SimplexSubset star(const SimplexSubset& subset) const;
    SimplexSubset closure(const SimplexSubset& subset) const;
    SimplexSubset link(const SimplexSubset& subset) const;
    bool isComplex(const SimplexSubset& subset) const;
    size_t isPureComplex(const SimplexSubset& subset) const;
    SimplexSubset boundary(const SimplexSubset& subset) const;
    //
    int nEdges() const;
    int nVertices() const;
    int nFaces() const;
    int nTets() const;
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

#endif /* DECOperators_h */
