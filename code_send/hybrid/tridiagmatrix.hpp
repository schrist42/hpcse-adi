#ifndef TRI_DIAG_MATRIX_HPP
#define TRI_DIAG_MATRIX_HPP


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
     * @param N     size of the matrix
     * @param l     value on the lower diagonal
     * @param m     value on the middle diagonal
     * @param u     value on the upper diagonal
     */
    TriDiagMatrix(int N, double l, double m, double u); // specific for the case that the diagonals always have the same values everywhere
    
    /**
     * Return the size of the matrix in one dimension.
     * 
     * @return      the matrix size in one dimension
     */
    int size() const { return N_; }
    
    
    /**
     * @defgroup getdiagonals Return the diagonals.
     * 
     * @{
     */
    std::vector<double> getL() const { return l_; } /**< Return the lower diagonal */
    std::vector<double> getM() const { return m_; } /**< Return the middle diagonal */
    std::vector<double> getU() const { return u_; } /**< Return the upper diagonal */
    /** @} */
    
    
    /**
     * Print the matrix to the stream os.
     * 
     * @param os        stream to print the matrix to
     * @param matrix    matrix
     */
    friend std::ostream& operator<<(std::ostream& os, const TriDiagMatrix& matrix);

private:
    /**
     * size of the matrix
     */
    const int N_;
    
    
    /**
     * @defgroup Diagonals      Diagonals of the tridiagonal matrix.
     * 
     * The diagonals all have the same length (n_), but have 0 values where they are outside of the matrix.
     * The indices determine the row.
     * 
     * | m0 u0 0  |
     * | l1 m1 u1 |
     * | 0  l2 m2 |
     * 
     * i.e.: l_[0] and u_[2] are 0.
     * 
     * @{
     */
    std::vector<double> l_; /**< Lower diagonal */
    const std::vector<double> m_; /**< Middle diagonal */
    std::vector<double> u_; /**< Upper diagonal */
    /** @} */
};


#endif // TRI_DIAG_MATRIX_HPP
