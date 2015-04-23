#include <iostream>
#include <iomanip>

#include "tridiagmatrix.hpp"


TriDiagMatrix::TriDiagMatrix()
    : N_(0)
    , l_(0)
    , m_(0)
    , u_(0)
{}


TriDiagMatrix::TriDiagMatrix(int N, double l, double m, double u)
    : N_(N)
    , l_(N,l)
    , m_(N,m)
    , u_(N,u)
{
    // set values that are outside of the matrix to 0
    l_[0] = 0;
    u_[N-1] = 0;
}



std::ostream& operator<< (std::ostream& os, const TriDiagMatrix& matrix)
{
    int n = matrix.size();
    int w = 8;
    
    for (int j=0; j<n; ++j) { // rows
        for (int i=0; i<n; ++i) { // columns
            if (i == j-1) {
                os << std::setw(w) << matrix.l_[j];
            }
            else if (i == j) {
                os << std::setw(w) << matrix.m_[j];
            }
            else if (i == j+1) {
                os << std::setw(w) << matrix.u_[j];
            }
            else {
                os << std::setw(w) << 0;
            }
            os << " ";
        }
        os << "\n";
    }
    return os;
}
