#ifndef TRI_DIAG_MATRIX
#define TRI_DIAG_MATRIX


#include <vector>
#include <ostream>

class TriDiagMatrix
{
public:
    /**
     * Default constructor
     */
    TriDiagMatrix();
    
    
    
    /**
     * Construct an object of the type TriDiagMatrix
     * 
     * @param n     size of the matrix
     * @param l     value on the lower diagonal
     * @param m     value on the middle diagonal
     * @param u     value on the upper diagonal
     */
    TriDiagMatrix(int N, double l, double m, double u); // specific for the case that the diagonals always have the same values everywhere
    
    int size() const { return N_; }
    std::vector<double> getL() const { return l_; }
    std::vector<double> getM() const { return m_; }
    std::vector<double> getU() const { return u_; }
    
    
    friend std::ostream& operator<<(std::ostream& os, const TriDiagMatrix& matrix);

private:
    /**
     * size of the matrix
     */
    const int N_;
    
    
    /**
     * Diagonals of the matrix.
     * 
     * They all have the same length (n_), but have 0 values where they are out of the matrix.
     * The indices determine the row.
     * 
     * | m0 u0 0  |
     * | l1 m1 u1 |
     * | 0  l2 m2 |
     * 
     * i.e.: l_[0] and u_[2] are 0.
     */
    std::vector<double> l_; // lower diagonal
    std::vector<double> m_; // middle diagonal
    std::vector<double> u_; // upper diagonal
};








#endif
