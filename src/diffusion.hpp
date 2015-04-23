#include <vector>

#include "tridiagmatrix.hpp"


class Diffusion
{
public:
    Diffusion(int N, double L, double dt, double Du, double Dv, double F, double k, int nSteps);
    
    /**
     * Run the simulation.
     */
    void run();
    
    /**
     * Print the fields to the specified file.
     * 
     * @param uName     filename of the file to print u to
     * @param vName     filename of the file to print v to
     */
    void print_fields(std::string uName, std::string vName);
    
private:
    
    /**
     * Perform one simulation step.
     */
    void step();



    ///////////////////////////// PARAMETERS ///////////////////////////////////
    /**
     * Grid size.
     * 
     * Number of cells in each direction. The grid is square.
     */
    const int N_;
    /**
     * Total number of grid cells.
     */
    int Ntot_;
    
    /**
     * Size of the domain.
     * 
     * Length of the computation domain in each direction. The grid is square.
     */
//    const double L_;
    
    const double dx_;
    
    /**
     * Time step.
     */
    const double dt_;
    
    
    /**
     * Number of simulation steps
     */
    const int nSteps_;
    
        
        
    ///////////////////////////// FIELDS ///////////////////////////////////////
    /**
     * Quantities to diffuse.
     * 
     * The quantities are two different chemical species.
     */
    std::vector<double> u_;
    std::vector<double> v_;
    
    /**
     * Diffusion coefficients.
     */
    const double Du_;
    const double Dv_;
    /**
     * Model parameters.
     */
    const double F_;
    const double k_;
    
        
    /**
     * Tridiagonal matrices.
     * 
     * These matrices don't depend on the time step or the diffusing quantities and are thus constant.
     */
     
    TriDiagMatrix matU1_; /**< Matrix for the first step of u. */
    TriDiagMatrix matU2_; /**< Matrix for the second step of u. */
    
    TriDiagMatrix matV1_;
    TriDiagMatrix matV2_;
    
    
    /////////////////////// HELPER FUNCTIONS ///////////////////////////////////
    void initialize_fields();
};




