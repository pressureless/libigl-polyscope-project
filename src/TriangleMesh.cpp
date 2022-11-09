//
//  TriangleMesh.cpp
//  DEC
//
//  Created by pressure on 10/31/22.
//

#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm> 
#include "dec_util.h"
#include "TriangleMesh.h"
using Eigen::VectorXi;
using Eigen::MatrixXi;
typedef Eigen::Matrix< size_t, Eigen::Dynamic, 1> Vector;
typedef Eigen::Matrix< size_t, 1, Eigen::Dynamic> RowVector;
typedef Eigen::Matrix< size_t, Eigen::Dynamic, Eigen::Dynamic> Matrix;
TriangleMesh::TriangleMesh(){

}
TriangleMesh::TriangleMesh(Eigen::MatrixXi &T){
    Matrix m = T.cast<size_t>();
    this->initialize(m);
}
TriangleMesh::TriangleMesh(Matrix &T){
    this->initialize(T);
}
void TriangleMesh::initialize(Eigen::MatrixXi &T){
    Matrix m = T.cast<size_t>();
    this->initialize(m);
}

void TriangleMesh::initialize(Matrix &T){
    std::cout<<"T cols:"<<T.cols()<<std::endl;
    this->T = preprocess_matrix(T);
    // std::cout<<"T:\n"<<this->T<<std::endl;
    if (T.cols() == 4) {
        // tets, assume each row (tet) already has the positive orientation
        this->T = T;
        Vector maxVal = T.rowwise().maxCoeff();
        this->num_v = maxVal.maxCoeff()+1;
        this->create_edges();
        this->create_faces();
        //
        this->build_boundary_mat1();
        this->build_boundary_mat2();
        this->build_boundary_mat3();
        }
    else if(T.cols() == 3){
        // faces, assume each row (face) already has the positive orientation
        this->F = T;
        std::cout<<"this->F:\n"<<this->F<<std::endl;
        Vector maxVal = this->F.rowwise().maxCoeff();
        this->num_v = maxVal.maxCoeff()+1;
        this->create_edges_from_faces();
        // boundary mat
        this->build_boundary_mat1();
        this->build_boundary_mat2();
    }
    std::cout<<"Total vertices:"<<this->num_v<<", edges:"<<this->E.rows()<<", faces:"<<this->F.rows()<<", tets:"<<this->T.rows()<<std::endl;
}

int TriangleMesh::n_edges() const{
    return this->E.rows();
}
int TriangleMesh::n_vertices() const{
    return this->num_v;
}
int TriangleMesh::n_faces() const{
    return this->F.rows();
}
int TriangleMesh::n_tets() const{
    return this->T.rows();
}

size_t TriangleMesh::get_edge_index(size_t i, size_t j, int &sign){
    if (i < j)
    {
        sign = 1;
        return this->map_e[std::make_tuple(i, j)];
    }
    sign = -1;
    return this->map_e[std::make_tuple(j, i)];
}

size_t TriangleMesh::get_face_index(size_t i, size_t j, size_t k, int &sign){
    RowVector r(3); r << i, j, k;
    RowVector p = permute_rvector(r); // get sorted face
    // find orientation
    key_f key = std::make_tuple(p(0), p(1), p(2));
    auto search = this->map_f.find(key);
    if (search != this->map_f.end()){
        // the input face has the same orientation as the one stored
        sign = 1;
        return this->map_f[key];
    }
    sign = -1;
    return this->map_f[std::make_tuple(p(0), p(2), p(1))];
}

void TriangleMesh::create_edges(){
    Matrix E(6*T.rows(), 2);
    for (int i=0; i<this->T.rows(); i++) {
        Vector v0(2); v0 << this->T(i,0), this->T(i,1);
        E.row(6*i) = sort_vector(v0);
        Vector v1(2); v1 << this->T(i,0), this->T(i,2);
        E.row(6*i+1) = sort_vector(v1);
        Vector v2(2); v2 << this->T(i,0), this->T(i,3);
        E.row(6*i+2) = sort_vector(v2);
        Vector v3(2); v3 << this->T(i,1), this->T(i,2);
        E.row(6*i+3) = sort_vector(v3);
        Vector v4(2); v4 << this->T(i,1), this->T(i,3);
        E.row(6*i+4) = sort_vector(v4);
        Vector v5(2); v5 << this->T(i,2), this->T(i,3);
        E.row(6*i+5) = sort_vector(v5);
    }
    this->E = remove_duplicate_rows(sort_matrix(E));
    std::cout<<"this->E:\n"<<this->E<<std::endl;
    // create the mapping
    for (int i=0; i<this->E.rows(); i++) {
        this->map_e.insert(std::pair<key_e, size_t>(std::make_tuple(this->E(i,0), this->E(i,1)), i));
    }
}

void TriangleMesh::create_edges_from_faces(){
    Matrix E(3*this->F.rows(), 2);
    for (int i=0; i<this->F.rows(); i++) {
        Vector v0(2); v0 << this->F(i,0), this->F(i,1);
        E.row(3*i) = sort_vector(v0);
        Vector v1(2); v1 << this->F(i,0), this->F(i,2);
        E.row(3*i+1) = sort_vector(v1);
        Vector v2(2); v2 << this->F(i,1), this->F(i,2);
        E.row(3*i+2) = sort_vector(v2); 
    }
    // std::cout<<"before this->E:\n"<<E<<std::endl;
    this->E = remove_duplicate_rows(sort_matrix(E));    
    std::cout<<"this->E:\n"<<this->E<<std::endl;
    // create the mapping
    for (int i=0; i<this->E.rows(); i++) {
        this->map_e.insert(std::pair<key_e, size_t>(std::make_tuple(this->E(i,0), this->E(i,1)), i));
    }
}


void TriangleMesh::create_faces(){
    Matrix F(4*T.rows(), 3);
    for (int i=0; i<this->T.rows(); i++) {
        Vector v0(3); v0 << this->T(i,0), this->T(i,1), this->T(i,2);
        F.row(4*i) = sort_vector(v0);
        Vector v1(3); v1 << this->T(i,0), this->T(i,1), this->T(i,3);
        F.row(4*i+1) = sort_vector(v1);
        Vector v2(3); v2 << this->T(i,0), this->T(i,2), this->T(i,3);
        F.row(4*i+2) = sort_vector(v2);
        Vector v3(3); v3 << this->T(i,1), this->T(i,2), this->T(i,3);
        F.row(4*i+3) = sort_vector(v3);
    }
    this->F = remove_duplicate_rows(sort_matrix(F));
    std::cout<<"this->F:\n"<<this->F<<std::endl;
    // create the mapping
    for (int i=0; i<this->F.rows(); i++) {
        this->map_f.insert(std::pair<key_f, size_t>(std::make_tuple(this->F(i,0), this->F(i,1), this->F(i,2)), i));
    }
}

void TriangleMesh::build_boundary_mat3(){
    this->bm3.resize(this->F.rows(), this->T.rows());
    std::vector<Eigen::Triplet<size_t> > tripletList;
    int sign = 1;
    // The faces of a tet +(t0,t1,t2,t3) 
    for(int i=0; i<this->T.rows(); i++){
        // −(t0,t1,t2)
        size_t idx = get_face_index(this->T(i,0), this->T(i,1), this->T(i,2), sign);
        tripletList.push_back(Eigen::Triplet<size_t>(idx, i, -sign));
        // +(t0,t1,t3)
        idx = get_face_index(this->T(i,0), this->T(i,1), this->T(i,3), sign);
        tripletList.push_back(Eigen::Triplet<size_t>(idx, i, sign));
        // −(t0,t2,t3)
        idx = get_face_index(this->T(i,0), this->T(i,2), this->T(i,3), sign);
        tripletList.push_back(Eigen::Triplet<size_t>(idx, i, -sign));
        // +(t1,t2,t3)
        idx = get_face_index(this->T(i,1), this->T(i,2), this->T(i,3), sign);
        tripletList.push_back(Eigen::Triplet<size_t>(idx, i, sign));
        // tripletList.push_back(Eigen::Triplet<size_t>(this->map_f[std::make_tuple(this->T(i,0), this->T(i,1), this->T(i,2))], i, -1));
        // tripletList.push_back(Eigen::Triplet<size_t>(this->map_f[std::make_tuple(this->T(i,0), this->T(i,1), this->T(i,3))], i, 1));
        // tripletList.push_back(Eigen::Triplet<size_t>(this->map_f[std::make_tuple(this->T(i,0), this->T(i,2), this->T(i,3))], i, -1));
        // tripletList.push_back(Eigen::Triplet<size_t>(this->map_f[std::make_tuple(this->T(i,1), this->T(i,2), this->T(i,3))], i, 1));
    }
    this->bm3.setFromTriplets(tripletList.begin(), tripletList.end());
    this->pos_bm3 = this->bm3.cwiseAbs();
    // std::cout<<"this->bm3:\n"<<this->bm3<<std::endl;
    // std::cout<<"this->pos_bm3:\n"<<this->pos_bm3<<std::endl;
}

void TriangleMesh::build_boundary_mat2(){
    this->bm2.resize(this->E.rows(), this->F.rows());
    std::vector<Eigen::Triplet<size_t> > tripletList;
    int sign = 1;
    // The faces of a triangle +(f0,f1,f2)
    for(int i=0; i<this->F.rows(); i++){
        // +(f0,f1)
        size_t idx = get_edge_index(this->F(i,0), this->F(i,1), sign);
        tripletList.push_back(Eigen::Triplet<size_t>(idx, i, sign));
        // −(f0,f2)
        idx = get_edge_index(this->F(i,0), this->F(i,2), sign);
        tripletList.push_back(Eigen::Triplet<size_t>(idx, i, -sign));
        // +(f1, f2)
        idx = get_edge_index(this->F(i,1), this->F(i,2), sign);
        tripletList.push_back(Eigen::Triplet<size_t>(idx, i, sign));
        // tripletList.push_back(Eigen::Triplet<size_t>(this->map_e[std::make_tuple(this->F(i,0), this->F(i,1))], i, 1));
        // tripletList.push_back(Eigen::Triplet<size_t>(this->map_e[std::make_tuple(this->F(i,0), this->F(i,2))], i, -1));
        // tripletList.push_back(Eigen::Triplet<size_t>(this->map_e[std::make_tuple(this->F(i,1), this->F(i,2))], i, 1));
    }
    this->bm2.setFromTriplets(tripletList.begin(), tripletList.end());
    this->pos_bm2 = this->bm2.cwiseAbs();
    // std::cout<<"this->bm2:\n"<<this->bm2<<std::endl;
    // std::cout<<"this->pos_bm2:\n"<<this->pos_bm2<<std::endl;
}

void TriangleMesh::build_boundary_mat1(){
    this->bm1.resize(this->num_v, this->E.rows());
    std::vector<Eigen::Triplet<size_t> > tripletList;
    for(int i=0; i<this->E.rows(); i++){
        tripletList.push_back(Eigen::Triplet<size_t>(this->E(i,0), i, -1));
        tripletList.push_back(Eigen::Triplet<size_t>(this->E(i,1), i, 1));
    }
    this->bm1.setFromTriplets(tripletList.begin(), tripletList.end());
    this->pos_bm1 = this->bm1.cwiseAbs();
    // std::cout<<"this->bm1:\n"<<this->bm1<<std::endl;
    // std::cout<<"this->pos_bm1:\n"<<this->pos_bm1<<std::endl;
}

// v as input
std::set<size_t> TriangleMesh::get_adjacent_vertices_v(size_t vindex){
    // simplex API
    // SimplexSubset verts;
    // verts.addVertex(vindex);
    // SimplexSubset result = link(verts);
    // std::set<size_t> neighbors = result.vertices;
    // print_set(neighbors);
    // return neighbors;
    // matrix API
    SparseMatrix<int> vv = this->pos_bm1 * this->pos_bm1.transpose();
    std::set<size_t> result;
    for (int k=0; k<vv.outerSize(); ++k){
      for (SparseMatrix<int>::InnerIterator it(vv,k); it; ++it) {
        if (it.row() == vindex && it.col() != vindex)
        {
            result.insert(it.col());
        }
      }
    } 
    print_set(result);
    return result;
}

std::set<size_t> TriangleMesh::get_incident_edges_v(size_t vindex){
    std::set<size_t> result;
    for (int k=0; k<this->pos_bm1.outerSize(); ++k){
      for (SparseMatrix<int>::InnerIterator it(this->pos_bm1,k); it; ++it) {
        if (it.row() == vindex)
        {
            result.insert(it.col());
        }
      }
    } 
    print_set(result);
    return result;
}

std::set<size_t> TriangleMesh::get_incident_faces_v(size_t vindex){
    SparseMatrix<int> vf = this->pos_bm1 * this->pos_bm2;
    std::set<size_t> result;
    for (int k=0; k<vf.outerSize(); ++k){
      for (SparseMatrix<int>::InnerIterator it(vf,k); it; ++it) {
        if (it.row() == vindex)
        {
            result.insert(it.col());
        }
      }
    } 
    print_set(result);
    return result;
}
// e as input
std::set<size_t> TriangleMesh::get_incident_vertices_e(size_t eindex){
    std::set<size_t> result;
    for (int k=0; k<this->pos_bm1.outerSize(); ++k){
      for (SparseMatrix<int>::InnerIterator it(this->pos_bm1,k); it; ++it) {
        if (it.col() == eindex)
        {
            result.insert(it.row());
        }
      }
    } 
    print_set(result);
    return result;
}
std::set<size_t> TriangleMesh::get_incident_faces_e(size_t eindex){
    std::set<size_t> result;
    for (int k=0; k<this->pos_bm2.outerSize(); ++k){
      for (SparseMatrix<int>::InnerIterator it(this->pos_bm2,k); it; ++it) {
        if (it.row() == eindex)
        {
            result.insert(it.col());
        }
      }
    } 
    print_set(result);
    return result;
}
//
std::set<size_t> TriangleMesh::get_diamond_vertices_e(size_t eindex){
    SimplexSubset edges;
    edges.addEdge(eindex);
    SimplexSubset result = closure(star(edges));
    std::set<size_t> ver = result.vertices;
    // ver.erase(this->E(eindex, 0));
    // ver.erase(this->E(eindex, 1));
    print_set(ver);
    return ver;
}

// f as input
std::set<size_t> TriangleMesh::get_incident_vertices_f(size_t findex){
    SparseMatrix<int> vf = this->pos_bm1 * this->pos_bm2;
    std::set<size_t> result;
    for (int k=0; k<vf.outerSize(); ++k){
      for (SparseMatrix<int>::InnerIterator it(vf,k); it; ++it) {
        if (it.col() == findex)
        {
            result.insert(it.row());
        }
      }
    } 
    print_set(result);
    return result;
}
std::set<size_t> TriangleMesh::get_incident_edges_f(size_t findex){ 
    std::set<size_t> result;
    for (int k=0; k<this->pos_bm2.outerSize(); ++k){
      for (SparseMatrix<int>::InnerIterator it(this->pos_bm2,k); it; ++it) {
        if (it.col() == findex)
        {
            result.insert(it.row());
        }
      }
    } 
    print_set(result);
    return result;
}
std::set<size_t> TriangleMesh::get_adjacent_faces_f(size_t findex){
    SparseMatrix<int> ff = (this->pos_bm2).transpose()*this->pos_bm2;
    std::set<size_t> result;
    for (int k=0; k<ff.outerSize(); ++k){
      for (SparseMatrix<int>::InnerIterator it(ff,k); it; ++it) {
        if (it.row() == findex && it.col() != findex)
        {
            result.insert(it.col());
        }
      }
    }
    print_set(result);
    return result;
}

std::set<size_t> TriangleMesh::get_adjacent_faces_f2(size_t findex){
    SparseMatrix<int> ff = (this->pos_bm1*this->pos_bm2).transpose()*(this->pos_bm1*this->pos_bm2);
    std::set<size_t> result;
    for (int k=0; k<ff.outerSize(); ++k){
      for (SparseMatrix<int>::InnerIterator it(ff,k); it; ++it) {
        if (it.row() == findex && it.col() != findex)
        {
            result.insert(it.col());
        }
      }
    }
    print_set(result);
    return result;
}

VectorXi TriangleMesh::build_vertex_vector(const SimplexSubset& subset) const{
    VectorXi v = VectorXi::Zero(this->num_v);
    for (size_t idx : subset.vertices)
    {
        v[idx] = 1;
    }
    return v;
}

VectorXi TriangleMesh::build_edge_vector(const SimplexSubset& subset) const{
    VectorXi e = VectorXi::Zero(this->E.rows());
    for (auto idx : subset.edges)
    {
        e[idx] = 1;
    }
    return e;
}

VectorXi TriangleMesh::build_face_vector(const SimplexSubset& subset) const{
    VectorXi f = VectorXi::Zero(this->F.rows());
    for (size_t idx : subset.faces)
    {
        f[idx] = 1;
    }
    return f;
}
Eigen::VectorXi TriangleMesh::build_vertex_vector(const std::set<size_t>& vset) const{
    VectorXi v = VectorXi::Zero(this->num_v);
    for (size_t idx : vset)
    {
        v[idx] = 1;
    }
    return v;
}
Eigen::VectorXi TriangleMesh::build_edge_vector(const std::set<size_t>& eset) const{
    VectorXi e = VectorXi::Zero(this->E.rows());
    for (auto idx : eset)
    {
        e[idx] = 1;
    }
    return e;
}
Eigen::VectorXi TriangleMesh::build_face_vector(const std::set<size_t>& fset) const{
    VectorXi f = VectorXi::Zero(this->F.rows());
    for (size_t idx : fset)
    {
        f[idx] = 1;
    }
    return f;
}

SimplexSubset TriangleMesh::star(const SimplexSubset& subset) const{
    SimplexSubset newSet = subset.deepCopy();
    VectorXi v = this->build_vertex_vector(subset);
    std::cout<<"v:"<<v<<std::endl;
    VectorXi e = this->build_edge_vector(subset);
    // std::cout<<"pos_bm1:"<<pos_bm1<<std::endl;
    // std::cout<<"pos_bm2:"<<pos_bm2<<std::endl;
    SparseMatrix<int> fv = (this->pos_bm1*this->pos_bm2).transpose();
    std::cout<<"fv:"<<fv<<std::endl;
    VectorXi extraE = this->pos_bm1.transpose() * v;
    for (size_t i = 0; i < extraE.size(); ++i)
    {
        if (extraE[i])
        {
            newSet.addEdge(i);
        }
    }
    VectorXi extraF = this->pos_bm2.transpose() * e + fv * v;
    for (size_t i = 0; i < extraF.size(); ++i)
    {
        if (extraF[i])
        {
            newSet.addFace(i);
        }
    }
    std::cout<<"extraF:"<<extraF<<std::endl;
    return newSet;
}
SimplexSubset TriangleMesh::closure(const SimplexSubset& subset) const{
    SimplexSubset newSet = subset.deepCopy();
    VectorXi f = this->build_face_vector(subset);
    VectorXi extraE = this->pos_bm2 * f;
    for (size_t i = 0; i < extraE.size(); ++i)
    {
        if (extraE[i])
        {
            newSet.addEdge(i);
        }
    }
    VectorXi e = this->build_edge_vector(newSet);
    VectorXi extraV = this->pos_bm1 * e;
    for (size_t i = 0; i < extraV.size(); ++i)
    {
        if (extraV[i])
        {
            newSet.addVertex(i);
        }
    }
    return newSet;
}
SimplexSubset TriangleMesh::link(const SimplexSubset& subset) const{
    SimplexSubset newSet = this->closure(this->star(subset));
    newSet.deleteSubset(this->star(this->closure(subset)));
    return newSet;
}
bool TriangleMesh::is_complex(const SimplexSubset& subset) const{
    SimplexSubset newSet = this->closure(subset);
    return newSet.equals(subset);
}
size_t TriangleMesh::is_pure_complex(const SimplexSubset& subset) const{
    size_t degree = -1;
    SparseMatrix<int> fv = (this->pos_bm1*this->pos_bm2).transpose();
    if (this->is_complex(subset))
    {
        if (subset.faces.size() > 0)
        {
            degree = 2;
            // check edges
            for (size_t e : subset.edges)
            {
                bool found = false;
                for (size_t f : subset.faces)
                {
                    if (this->pos_bm2.coeff(e, f) > 0)
                    {
                        found = true;
                        break;
                    }
                }
                if (not found)
                {
                    degree = -1;
                    break;
                }
            }
            // check vertices
            for (size_t v : subset.vertices)
            {
                bool found = false;
                for (size_t f : subset.faces)
                {
                    if (fv.coeff(f, v) > 0)
                    {
                        found = true;
                        break;
                    }
                }
                if (not found)
                {
                    degree = -1;
                    break;
                }
            }
        }
        else if(subset.edges.size() > 0){
            degree = 1;
            // check vertices
            for (size_t v : subset.vertices)
            {
                bool found = false;
                for (size_t e : subset.edges)
                {
                    if (this->pos_bm1.coeff(v, e) > 0)
                    {
                        found = true;
                        break;
                    }
                }
                if (not found)
                {
                    degree = -1;
                    break;
                }
            }
        }
        else if(subset.vertices.size() > 0){
            degree = 0;
        }
    }
    return degree;
} 
SimplexSubset TriangleMesh::boundary(const SimplexSubset& subset) const{
    SimplexSubset newSet;
    SparseMatrix<int> fv = (this->pos_bm1*this->pos_bm2).transpose();
    if (this->is_pure_complex(subset))
    {
        // check edges
        for (size_t e : subset.edges)
        {
            size_t degree = 0;
            for (size_t f : subset.faces)
            {
                if (this->pos_bm2.coeff(e, f) > 0)
                {
                    degree++;
                }
            }
            if (degree == 1)
            {
                newSet.addEdge(e);
            }
        }
        // check vertices
        for (size_t v : subset.vertices)
        {
            size_t degree = 0;
            for (size_t f : subset.faces)
            {
                if (fv.coeff(f, v) > 0)
                {
                    degree++;
                }
            }
            for (size_t e : subset.edges)
            {
                if (this->pos_bm1.coeff(v, e) > 0)
                {
                    degree++;
                }
            }
            if (degree == 1)
            {
                newSet.addVertex(v);
            }
        }
    }
    return this->closure(newSet);
}
    