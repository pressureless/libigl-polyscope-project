//
//  util.h
//  DEC
//
//  Created by pressure on 10/31/22.
//
#pragma once

#include <iostream>
#include <set>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
// using Eigen::Matrix;
// using Eigen::Vector;
#include "TriangleMesh.h"

std::tuple<std::set<int>, std::set<int>, std::set<int>> MeshSets(const TriangleMesh& mesh);
std::tuple<Eigen::SparseMatrix<int>, Eigen::SparseMatrix<int> > BoundaryMatrices(const TriangleMesh& mesh);
std::set<int> vector_to_vertices(const TriangleMesh& mesh, const Eigen::VectorXi& vi);
// std::set<int> vector_to_edges(const TriangleMesh& mesh, const Eigen::VectorXi& ei);
// std::set<int> vector_to_faces(const TriangleMesh& mesh, const Eigen::VectorXi& fi);

