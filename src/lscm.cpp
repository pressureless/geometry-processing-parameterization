#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/cotmatrix.h> 
#include <igl/massmatrix.h>
#include <igl/boundary_loop.h>
#include <igl/repdiag.h>
#include <igl/eigs.h>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
	int num = V.rows();
	U.resize(num, 2);
	// Replace with your code
	Eigen::SparseMatrix<double> A;
	vector_area_matrix(F, A);
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);
	Eigen::SparseMatrix<double> zeros(L.rows(), L.cols());
	Eigen::SparseMatrix<double> L2;
	igl::repdiag(L, 2, L2);

	Eigen::SparseMatrix<double> Q = L2 - A;

	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	Eigen::SparseMatrix<double> B;
	igl::repdiag(M, 2, B);

    //"Use eigen-decomposition to find an un-biased, non-trivial minimizer"
    bool loop = true;
    int nums_eigs = 1;
    Eigen::MatrixXd sU;
    Eigen::VectorXd sS; 
    Eigen::MatrixXd resultU;
    while(loop){
    	//loop until find a valid one
  	    std::cout<<"loop nums_eigs:"<<nums_eigs<<std::endl;
    	igl::eigs(Q, B, nums_eigs, igl::EIGS_TYPE_SM, sU, sS);
    	for (int i = 0; i < sS.rows(); ++i)
    	{
    		if (std::abs(sS(i)) > 0.0001 )  // non zero
    		{
    			loop = false;
    			resultU = sU.col(i);
  	            std::cout<<"found:"<<sS<<std::endl;
    			break;
    		}
    	}
    	nums_eigs++;
    }  
 
    Eigen::Map<Eigen::MatrixXd> calU(resultU.data(), num, 2);
    //"use singular value decomposition to find a canonical rotation to line the principle axis of U with the 𝑥-axis of the parametric domain."
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(calU.transpose() * calU, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd matU = svd.matrixU();
    Eigen::MatrixXd matV = svd.matrixV(); 
    std::cout<<"matU  :"<<matU<<std::endl;
    std::cout<<"matV  :"<<matV<<std::endl; 
    // U = calU * matU;
    Eigen::MatrixXd res = calU * matU;
    //make sure the normal is correct 
    // U.col(1) = res.col(0);
    // U.col(0) = res.col(1); 
    U = res;
}
