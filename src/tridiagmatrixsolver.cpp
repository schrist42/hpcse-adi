#include <iostream>
#include <assert.h>

#include "tridiagmatrixsolver.hpp"



/*
 * stride   how long the space between elements is
 */
//void TriDiagMatrixSolver::solve(const TriDiagMatrix& mat, const std::vector<double>& rhs, std::vector<double>& result, int start, int stride)
//{
////    std::cout << "mat.size = " << mat.size() << ", rhs.size = " << rhs.size() << "\n";
//    assert(mat.size() == rhs.size());
//    
//    int n = mat.size();
//    
//    //result.resize(n);
//    
//    std::vector<double> l = mat.getL();
//    std::vector<double> m = mat.getM();
//    std::vector<double> u = mat.getU();
//    
//    
//    n--; // since we start from x0 (not x1) (n now denotes the max index)
//    u[0] /= m[0];
//    
//    std::vector<double> r(rhs);
//    r[0] /= m[0];
//    
//    double temp;

//    for (int i=1; i<n; ++i) {
//        temp = m[i] - l[i]*u[i-1];
//        u[i] /= temp;
//        r[i] = (r[i] - l[i]*r[i-1]) / temp;
//    }

//    result[start + n*stride] = (r[n] - l[n]*r[n-1]) / (m[n] - l[n]*u[n-1]);

//    for (int i=n-1; i>=0; --i) {
//        result[start + i*stride] = r[i] - u[i]*result[start + (i+1)*stride];
//    }
//}


/*
 * inc:     increment in the result matrix
 */

void TriDiagMatrixSolver::solve(const TriDiagMatrix& mat, const std::vector<double>& rhs, double *result, unsigned int inc)
{
//    std::cout << "mat.size = " << mat.size() << ", rhs.size = " << rhs.size() << "\n";
    assert(mat.size() == rhs.size());
    
    int n = mat.size();
    
    //result.resize(n);
    
    std::vector<double> l = mat.getL();
    std::vector<double> m = mat.getM();
    std::vector<double> u = mat.getU();
    
    
    n--; // since we start from x0 (not x1) (n now denotes the max index)
    u[0] /= m[0];
    
    std::vector<double> r(rhs);
    r[0] /= m[0];
    
    double temp;

    for (int i=1; i<n; ++i) {
        temp = m[i] - l[i]*u[i-1];
        u[i] /= temp;
        r[i] = (r[i] - l[i]*r[i-1]) / temp;
    }

    *(result + n*inc) = (r[n] - l[n]*r[n-1]) / (m[n] - l[n]*u[n-1]);

    for (int i=n-1; i>=0; --i) {
        *(result + i*inc) = r[i] - u[i]*(*(result + (i+1)*inc));
    }
}


