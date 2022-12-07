#include "MeshHelper.h"

std::tuple<std::set<int>, std::set<int>, std::set<int>> MeshSets(const TriangleMesh& mesh){
	return std::tuple<std::set<int>, std::set<int>, std::set<int>>(mesh.Vi, mesh.Ei, mesh.Fi);
}

std::tuple<Eigen::SparseMatrix<int>, Eigen::SparseMatrix<int> > BoundaryMatrices(const TriangleMesh& mesh){
	return std::tuple< Eigen::SparseMatrix<int>, Eigen::SparseMatrix<int> >(mesh.bm1, mesh.bm2);
}

std::set<int> vector_to_vertices(const TriangleMesh& mesh, const Eigen::VectorXi& vi){
	return mesh.vector_to_vertices(vi);
}

// std::set<int> vector_to_edges(const TriangleMesh& mesh, const Eigen::VectorXi& ei){
// 	return mesh.vector_to_vertices(ei);
// }

// std::set<int> vector_to_faces(const TriangleMesh& mesh, const Eigen::VectorXi& fi){
// 	return mesh.vector_to_vertices(fi);
// }