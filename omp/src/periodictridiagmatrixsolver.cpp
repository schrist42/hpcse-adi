#include <iostream>
#include <assert.h>
#include <algorithm>

#include "periodictridiagmatrixsolver.hpp"
#include "tridiagmatrixsolver.hpp"

#define R(i) (*(result + (i)*inc))


/*
 * inc:     increment in the result matrix
 */
void PeriodicTriDiagMatrixSolver::solve(int n, const TriDiagMatrix& mat, const std::vector<double>& rhs, double *result, unsigned int inc)
{
    assert(mat.size() == rhs.size());
    assert(mat.size() == n);
    
    
//    std::vector<double> l = mat.getL();
//    std::vector<double> m = mat.getM();
//    std::vector<double> u = mat.getU();
    
    double l_first = mat.getL()[0];
    double l_last  = mat.getL()[n-1];
    double m_last  = mat.getM()[n-1];
    double u_2last = mat.getU()[n-2];
    double u_last  = mat.getU()[n-1];
    

    // result = xa + xb*xn, 1 smaller than the original result (compute the last element manually)
    std::vector<double> xa(n-1);
    std::vector<double> xb(n-1);
    
    std::vector<double> rhs2(n-1,0.);
//    rhs2[0]   = -l[0];      // first element of l
//    rhs2[n-2] = -u[n-2];    // second to last element of u, last element of rhs2, since the size is n-1
    rhs2[0]   = -l_first;   // first element of l
    rhs2[n-2] = -u_2last;   // second to last element of u, last element of rhs2, since the size is n-1
    
    
    // ignore the last row in the matrix and the rhs
    TriDiagMatrixSolver::solve(n-1, mat, rhs,  &xa[0], 1);
    TriDiagMatrixSolver::solve(n-1, mat, rhs2, &xb[0], 1);
    
    // compute last element of the results
//    *(result + (n-1)*inc) = (rhs[n-1] - u[n-1]*xa[0] - l[n-1]*xa[n-2]) / (m[n-1] + u[n-1]*xb[0] + l[n-1]*xb[n-2]);
//    *(result + (n-1)*inc) = (rhs[n-1] - u_last*xa[0] - l_last*xa[n-2]) / (m_last + u_last*xb[0] + l_last*xb[n-2]);
      R(n-1) = (rhs[n-1] - u_last*xa[0] - l_last*xa[n-2]) / (m_last + u_last*xb[0] + l_last*xb[n-2]);

    // multiply xb with the last element of result
//    std::transform(xb.begin(), xb.end(), xb.begin(), std::bind1st(std::multiplies<double>(),*(result + (n-1)*inc)));

    for (int i=0; i<=n-2; ++i) {
//        *(result + i*inc) = xa[i] + xb[i] * (*(result + (n-1)*inc));
        R(i) = xa[i] + xb[i] * R(n-1);
    }
}


