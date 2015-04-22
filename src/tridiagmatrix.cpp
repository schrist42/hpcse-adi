#include <iostream>

#include "tridiagmatrix.hpp"


TriDiagMatrix::TriDiagMatrix()
    : n_(0)
    , l_(0)
    , m_(0)
    , u_(0)
{
}


TriDiagMatrix::TriDiagMatrix(int n, double l, double m, double u)
    : n_(n)
    , l_(n,l)
    , m_(n,m)
    , u_(n,u)
{
    // set values that are outside of the matrix to 0
    l_[0] = 0;
    u_[n-1] = 0;
}

