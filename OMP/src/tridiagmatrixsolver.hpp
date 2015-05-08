#ifndef TRI_DIAG_MATRIX_SOLVER_HPP
#define TRI_DIAG_MATRIX_SOLVER_HPP

#include "tridiagmatrix.hpp"


/**
 * Solver for a tridiagonal matrix system.
 */
namespace TriDiagMatrixSolver
{
    /**
     * Solve a tridiagonal matrix system using the Thomas algorithm.
     * 
     * @param n         number of elements in the result
     * @param mat       tridiagonal matrix
     * @param rhs       right-hand side of the system
     * @param result    vector for the result, pointer to the first element
     * @param inc       increment for the elements of the result
     */
    void solve(int n, const TriDiagMatrix& mat, const std::vector<double>& rhs, double *result, unsigned int inc);
};


#endif
