#include <iostream>
#include <assert.h>

#include "tridiagmatrixsolver.hpp"


/*
 * inc:     increment in the result matrix
 * n:       number of rows/columns = number of equations
 */

void TriDiagMatrixSolver::solve(int n, const TriDiagMatrix& mat, const std::vector<double>& rhs, double *result, unsigned int inc)
{
    // will ignore some values, if the matrix and rhs are larger than n
    assert(mat.size() >= n);
    assert(rhs.size() >= n);
    

    // NOT by reference -> copy
    std::vector<double> a = mat.getL(); // lower diagonal, indexed 1..n-1
    std::vector<double> b = mat.getM(); // main diagonal, indexed 0..n-1
    std::vector<double> c = mat.getU(); // upper diagonal, indexed 0..n-2
    
    
    n--; // since we start from x0 (not x1) (n now denotes the max index)
    c[0] /= b[0];
    
    std::vector<double> r(rhs);
    r[0] /= b[0];
    
    double tmp;

    for (int i=1; i<n; ++i) {
        tmp = b[i] - a[i]*c[i-1];
        c[i] /= tmp;
        r[i] = (r[i] - a[i]*r[i-1]) / tmp;
    }

    *(result + n*inc) = (r[n] - a[n]*r[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i=n-1; i>=0; --i) {
        *(result + i*inc) = r[i] - c[i] * (*(result + (i+1)*inc));
    }
}


