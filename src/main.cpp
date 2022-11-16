//
//  main.cpp
//  DEC
//
//  Created by pressure on 10/31/22.
//

#include <iostream>
#include <Eigen/Dense>
#include "TriangleMesh.h" 
#include <Eigen/Dense>
#include <Eigen/Sparse> 
#include <igl/readOBJ.h>

using Eigen::MatrixXi;
using Eigen::MatrixXd;
/*
acos from trigonometry
get_diamond_vertices_e, edges from triangle_mesh( V, F )

V ∈ ℝ^(n×3)
F ∈ ℤ^(m×3)

normal(i, j, k) = (V_k,* - V_j,*)×(V_i,* - V_j,*)/||(V_i,* - V_j,*)×(V_k,* - V_j,*)|| where i,j,k ∈ ℤ index

dihedral(i,j,k,l) = acos(normal(i,j,k) ⋅ normal(j,i,l)) where i,j,k,l ∈ ℤ index

E = ∑_o dihedral(i,j,k,l)² where i=edges_o,1 ; j =edges_o,2; k, l = get_diamond_vertices_e(i, j) 

*/
#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>
#include "TriangleMesh.h"

struct iheartla {
    double A;

    iheartla(
        const Eigen::Matrix<double, Eigen::Dynamic, 3> & P,
        const Eigen::Matrix<int, Eigen::Dynamic, 3> & Faces)
    {
        const long n = P.rows();
        const long m = Faces.rows();
        assert( P.cols() == 3 );
        assert( Faces.cols() == 3 );
        TriangleMesh triangle_mesh_0(P, Faces);
        double sum_0 = 0;
        for(int f : triangle_mesh_0.Fi){
                // i,j,k = vertices(f)
            std::tuple< int, int, int > tuple = triangle_mesh_0.get_vertices_f(f);
            int i = std::get<0>(tuple);
            int j = std::get<1>(tuple);
            int k = std::get<2>(tuple);
            sum_0 += (1/double(2)) * (((P.row(j).transpose() - P.row(i).transpose())).cross((P.row(k).transpose() - P.row(i).transpose()))).lpNorm<2>();
        }
        // i,j,k = vertices(f)
        A = sum_0;
    }
};

int main(int argc, const char * argv[]) {
    //       2
    //     / | \
    //    0  |  3
    //     \ | /
    //       1
    Eigen::Matrix<double, Eigen::Dynamic, 3> V(4, 3);
    V <<
    0, 1, 0,
    1, 0, 0,
    1, 1, 0,
    2, 1, 0;
    MatrixXi T(2,3);
    T <<
    0,1,2,
    2,1,3; 
    std::cout <<"original:\n"<< T << std::endl;
    iheartla ii(V, T);
    // std::cout<<"ii:"<<ii.E<<std::endl;
    // TriangleMesh dec(V, T);
    // std::tuple< int, int > res = dec.get_diamond_vertices_e(2, 1);
    // std::cout<<"edges:"<<dec.edges<<std::endl;
    std::cout<<"A:\n"<<ii.A<<std::endl;
    // std::cout<<"second_face:"<<std::get<1>(res)<<std::endl;
    // insert code here... 
    return 0;
}
