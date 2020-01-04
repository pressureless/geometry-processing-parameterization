#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/cotmatrix.h> 
#include <igl/massmatrix.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>
#include <iostream>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{ 
	U.resize(V.rows(), 2);
  // Replace with your code
	Eigen::VectorXi bnd; 
	igl::boundary_loop(F, bnd); 

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

    Eigen::MatrixXd UV;
	igl::map_vertices_to_circle(V, bnd, UV);

	Eigen::MatrixXd UU = UV.col(0);
	Eigen::MatrixXd VV = UV.col(1);
    Eigen::MatrixXd U1 = Eigen::MatrixXd::Zero(V.rows(), 1);
    Eigen::MatrixXd U2 = Eigen::MatrixXd::Zero(V.rows(), 1);

	for (int i = 0; i < bnd.size(); ++i)
	{
		U1.row(bnd(i)) = UU.row(i); 
		U2.row(bnd(i)) = VV.row(i); 
	}  

	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(L.rows(), 1); 
	Eigen::MatrixXd Beq = Eigen::MatrixXd::Zero(0,0); 
	igl::min_quad_with_fixed_data<double> mqwf; 
	Eigen::SparseMatrix<double> AeqSparse(0, 0); 

	bool result = igl::min_quad_with_fixed_precompute(L, bnd, AeqSparse, false, mqwf); 
	igl::min_quad_with_fixed_solve(mqwf, B, UU, Beq, U1); 
	igl::min_quad_with_fixed_solve(mqwf, B, VV, Beq, U2);

    // U.col(0) = U1;
    // U.col(1) = U2; 
    //make sure the normal is correct, don't know why
    U.col(1) = U1;
    U.col(0) = U2; 
}

