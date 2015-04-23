#ifndef TRI_DIAG_MATRIX_SOLVER
#define TRI_DIAG_MATRIX_SOLVER

#include "tridiagmatrix.hpp"


/**
 * Solver for a tridiagonal matrix system.
 */
namespace TriDiagMatrixSolver
{
    /**
     * Solve a tridiagonal matrix system.
     * 
     * @param mat       tridiagonal matrix
     * @param rhs       right-hand side of the system
     * @param result    data structure for the result of the solver.
     */
    void solve(const TriDiagMatrix& mat, const std::vector<double>& rhs, std::vector<double>& result, int start, int stride);
};


#endif
