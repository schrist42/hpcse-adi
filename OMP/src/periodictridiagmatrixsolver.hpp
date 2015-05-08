#ifndef PERIODIC_TRI_DIAG_MATRIX_SOLVER_HPP
#define PERIODIC_TRI_DIAG_MATRIX_SOLVER_HPP

#include "tridiagmatrix.hpp"


/**
 * Solver for a tridiagonal matrix system.
 */
namespace PeriodicTriDiagMatrixSolver
{
    /*
     * Solve a tridiagonal matrix system.
     * 
     * @param mat       tridiagonal matrix
     * @param rhs       right-hand side of the system
     * @param result    data structure for the result of the solver.
     */
//    void solve(const TriDiagMatrix& mat, const std::vector<double>& rhs, std::vector<double>& result, int start, int stride);
    
    /**
     * Solve a tridiagonal matrix system.
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
