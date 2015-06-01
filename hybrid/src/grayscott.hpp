#ifndef GRAYSCOTT_HPP
#define GRAYSCOTT_HPP


#include <vector>

#include "tridiagmatrix.hpp"
#include "world.hpp"




//struct world_info
//{
//    int size;
//    int dims_x;
//    int dims_y;
//    
//    int left_proc;
//    int right_proc;
//    int top_proc;
//    int bottom_proc;
//    
//    int rank;
//    int cart_rank;
//    int coord_x;
//    int coord_y;
//    
//} world;




class GrayScott
{
public:
    /**
     * Construct a new simulation object.
     * 
     * @param N         number of grid cells in one dimension
     * @param L         length of the domain in one dimension
     * @param dt        time step
     * @param Du        diffusion coefficient for u
     * @param Dv        diffusion coefficient for v
     * @param F         model parameter
     * @param k         model parameter
     * @param nSteps    number of steps in the simulation
     */
    GrayScott(int N, double rmin, double rmax, double dt, double Du, double Dv, double F, double k, int nRep, int nSteps, std::string pngname, world_info w, bool localtranspose, unsigned int nthreads);
    
    /**
     * Destructor
     */
    ~GrayScott();
    
    /**
     * Run the simulation.
     */
    void run();
    
    /**
     * Benchmark the simulation.
     * 
     * Run multiple times and compute error.
     */
    void benchmark();
    
    /**
     * Perform one simulation step.
     * 
     * Public because of the visualization
     */
    void step();
    
    
    
    /**
     * Return the size of the system in one dimension.
     * 
     * @return  number of grid cells in one dimension.
     */
    int size() const { return N_; }
    
    /**
     * Get the field U of the simulation.
     * 
     * @return  the field U
     */
    std::vector<double> getU() const { return u_; }
    
    /**
     * Get the current step in the simulation.
     * 
     * @return  the current step.
     */
    int getCurrStep() const { return currStep_; }
    
    /**
     * Get the time step for the simulation.
     * 
     * @return  the timestep
     */
    double getDt() const { return dt_; }
    
    /**
     * Get the time that has passed.
     * 
     * @return  time
     */
    double getTime() const { return (double)currStep_*dt_; }
    

    
private:




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
     * Number of repetitions for benchmark
     */
    const int nRep_;
    /**
     * Number of simulation steps
     */
    const int nSteps_;
    /**
     * Current simulation step
     */
    int currStep_;
    
        
        
    ///////////////////////////// FIELDS ///////////////////////////////////////
    /**
     * Quantities to diffuse.
     * 
     * The quantities are two different chemical species.
     */
    std::vector<double> u_;
    std::vector<double> v_;
    
    /**
     * GrayScott coefficients.
     */
    const double Du_;
    const double Dv_;
    
    const double uCoeff;
    const double vCoeff;
    
    /**
     * Model parameters.
     */
    const double F_;
    const double k_;
    
        
    /**
     * @defgroup matrices   Tridiagonal matrices.
     * 
     * These matrices don't depend on the time step or the diffusing quantities and are thus constant.
     * 
     * @{
     */
     
    TriDiagMatrix matU1_; /**< Matrix for the first step of u. */
//    TriDiagMatrix matU2_; /**< Matrix for the second step of u. */
    
    TriDiagMatrix matV1_;
//    TriDiagMatrix matV2_;
    
    /** @} */
    
    
    /////////////////////// HELPER FUNCTIONS ///////////////////////////////////
    /**
     * Initialize the fields.
     * 
     * Initialize the chemical quantities U and V for the simulation.
     */
    void initialize_fields();
    
    /**
     * Print the fields to the specified file.
     * 
     * @param uName     filename of the file to print u to
     * @param vName     filename of the file to print v to
     */
    void save_fields();
    
    /**
     * Directory to save the output to
     */
    std::string dirPath_;
    
    std::string pngName_;
    
    void save_png();
    
    
    // added for mpi version
    
    world_info world;
    
    
    const double rmin_, rmax_;
    
    
    double xmin_loc, ymin_loc;
    double xmax_loc, ymax_loc;
    
    int NN_loc, Nx_loc, Ny_loc, Nb_loc;
    int NN_glo, Nx_glo, Ny_glo;
    
    MPI_Datatype bottom_boundary, top_boundary, block_resized_send, block_resized_recv;
    
    MPI_Comm cart_comm;
    
    
    // locally tranpose blocks or use datatype for transpose
    bool localtranspose_;
    
    
    /**
     * Number of threads used
     */
    const unsigned int nthreads_;
};



#endif

