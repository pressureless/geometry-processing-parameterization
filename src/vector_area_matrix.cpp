#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
    // Replace with your code
    int V_size = F.maxCoeff()+1;
    A.resize(V_size*2,V_size*2);

    Eigen::VectorXi bnd; 
	igl::boundary_loop(F, bnd); 

	std::vector<Eigen::Triplet<double>> tripletList; 
	int size = bnd.size();
	for (int i = 0; i < size; ++i)
	{
		int j = (i+1)%size;
		// (u_i,u_j,v_i,v_j)
	    tripletList.push_back(Eigen::Triplet<double>(bnd(i), V_size + bnd(j), 0.25 ));
	    tripletList.push_back(Eigen::Triplet<double>(bnd(j), V_size + bnd(i), -0.25 ));
	    tripletList.push_back(Eigen::Triplet<double>(V_size +bnd(i), bnd(j), -0.25 ));
	    tripletList.push_back(Eigen::Triplet<double>(V_size + bnd(j), bnd(i), 0.25 ));
	    // (u_i,v_i,u_j,v_j)
	    // tripletList.push_back(Eigen::Triplet<double>(2*bnd(i), 2*bnd(j)+1, 0.25 ));
	    // tripletList.push_back(Eigen::Triplet<double>(2*bnd(i) + 1, 2*bnd(j), -0.25 ));
	    // tripletList.push_back(Eigen::Triplet<double>(2*bnd(j), 2*bnd(i) + 1, -0.25 ));
	    // tripletList.push_back(Eigen::Triplet<double>(2*bnd(j)+1, 2*bnd(i), 0.25 ));
	} 
	A.setFromTriplets(tripletList.begin(), tripletList.end());
}

