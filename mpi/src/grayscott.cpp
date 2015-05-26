#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <ctime>
#include <boost/filesystem.hpp>
#include <chrono>
#include <mpi.h>

#include "grayscott.hpp"
#include "tridiagmatrixsolver.hpp"
#include "periodictridiagmatrixsolver.hpp"
#include "lodepng.h" // NOT loadpng


// takes into account the ghost rows (no need to consider them in coordinates)
#define U(x,y) u_[((x)+1)*Ny_loc + (y)]
#define V(x,y) v_[((x)+1)*Ny_loc + (y)]

#define TAG 0


GrayScott::GrayScott(int N, double rmin, double rmax, double dt, double Du, double Dv, double F, double k, int nSteps, std::string pngname, world_info w)
    : N_(N)
    , Ntot_(N*N)
//    , L_(L)
    , dx_((double) (rmax-rmin) / (double) (N-1))
    , dt_(dt)//(dx_*dx_ / (2.*std::max(Du,Dv)))
    , nSteps_(nSteps)
    , currStep_(0)
    , Du_(Du)
    , Dv_(Dv)
    , F_(F)
    , k_(k)
    , matU1_(N, -Du*dt/(2.*dx_*dx_), 1.+Du*dt/(dx_*dx_), -Du*dt/(2.*dx_*dx_))
    , matU2_(N, -Du*dt/(2.*dx_*dx_), 1.+Du*dt/(dx_*dx_), -Du*dt/(2.*dx_*dx_)) // equal to matU1_, since we have a square grid (dx==dy)
    , matV1_(N, -Dv*dt/(2.*dx_*dx_), 1.+Dv*dt/(dx_*dx_), -Dv*dt/(2.*dx_*dx_))
    , matV2_(N, -Dv*dt/(2.*dx_*dx_), 1.+Dv*dt/(dx_*dx_), -Dv*dt/(2.*dx_*dx_))
    , pngName_(pngname)
    , world(w)
    , rmin_(rmin)
    , rmax_(rmax)
{
    if (world.rank == 0) {
        // create directory to save output to
        time_t rawtime;
        struct tm * timeinfo;
        char buffer[80];
        time (&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(buffer,80,"%d-%m-%Y_%H-%M-%S",timeinfo);
        std::string timeString(buffer);
	
        dirPath_ = "data/" + timeString + "/";
	
        boost::filesystem::path dir(dirPath_);
        boost::filesystem::create_directory(dir);
    }
    
    
    
    
    // global grid
    Nx_glo = N;
    Ny_glo = N;
    NN_glo = Nx_glo * Ny_glo;
    
    // local grid
    Nx_loc = Nx_glo / world.dims_x;
    Ny_loc = Ny_glo / world.dims_y;
    NN_loc = Nx_loc * Ny_loc;
    
    
    
    // build periodic process geometry with cartesian communicator
    int periods[2] = {false, false};
    int dims[2] = {world.dims_x, world.dims_y};
    
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &cart_comm);
    
    MPI_Comm_rank(cart_comm, &world.cart_rank);
    
    MPI_Cart_shift(cart_comm, 1, 1, &world.top_proc, &world.bottom_proc);
    
    int coords[2];
    MPI_Cart_coords(cart_comm, world.cart_rank, 2, coords);
    
    world.coord_x = coords[0];
    world.coord_y = coords[1];
    
    
    
    // datatypes
    
    // build contiguous (rows) vectors for boundaries.
    // each process has multiple rows in the grid
    MPI_Type_contiguous(Ny_loc, MPI_DOUBLE, &bottom_boundary);
    MPI_Type_commit(&bottom_boundary);
    
    MPI_Type_contiguous(Ny_loc, MPI_DOUBLE, &top_boundary);
    MPI_Type_commit(&top_boundary);
    
    
    // build datatype for transpose
    int sizes[2]    = {Nx_loc, Ny_loc}; // size of global array
    int subsizes[2] = {Nx_loc, Nx_loc}; // size of sub-region (square)
    int starts[2]   = {0,0};            // let's say we're looking at region "0", which begins at index [0,0]
    
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &block);
//    MPI_Type_vector(Ny_loc, Nx_loc, Ny_loc, MPI_INT, &block);
    MPI_Type_commit(&block);

    // resize -> make contiguous
    MPI_Type_create_resized(block, 0, Nx_loc*sizeof(double), &block_resized);
//    MPI_Type_free(&block);
    MPI_Type_commit(&block_resized);
    
    
    
    // sub-domain boundaries
    xmin_loc = rmin + world.coord_x * Nx_loc * dx_;
    xmax_loc = xmin_loc + (Nx_loc - 1) * dx_;
    ymin_loc = rmin + world.coord_y * Ny_loc * dx_;
    ymax_loc = ymin_loc + (Ny_loc - 1) * dx_;
    
    
    
    initialize_fields();
    
    MPI_Barrier(MPI_COMM_WORLD);
}


GrayScott::~GrayScott()
{
    if (world.rank == 0) {
        boost::filesystem::path dir(dirPath_);
        if (boost::filesystem::exists(dir) && boost::filesystem::is_empty(dir)) {
            boost::filesystem::remove(dir);
        }
    }
    
    MPI_Type_free(&bottom_boundary);
    MPI_Type_free(&top_boundary);
    MPI_Type_free(&block);
    MPI_Type_free(&block_resized);
    
    MPI_Comm_free(&cart_comm);
}


void GrayScott::run()
{
    //save_fields();
	// timer init and start
//	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
//	start = std::chrono::high_resolution_clock::now();

    double start = MPI_Wtime();
    
    for (int i=0; i<nSteps_; ++i) {
//        if (world.rank == 0) std::cout << "starting step " << i << "\n";
        step();
//        if (world.rank == 0)
//            std::cout << "Step " << i << " done!\n";
    }
    
    
    // timer end and output of time
//	end = std::chrono::high_resolution_clock::now();
//	double elapsed = std::chrono::duration<double>(end-start).count();    
	
	double end = MPI_Wtime();
    double elapsed = end-start;
    
//    save_png();
    
    if (world.rank == 0) {
        std::string sep = "_";
        
        std::cout
        << "performance: " << '\t'
        << world.size << '\t'
        << elapsed << '\t'
        << N_ << '\t'
        << "\n*********\n"
        << std::endl;
    }
    
    save_png();
        

//    if (world.rank == 0) {
//        std::cout << "\n";
//        std::cout << "exec time: " << '\t' << elapsed << std::endl;
//        std::cout << "\n";
//        
//        save_png();
//    }
}

#define UHALF(x,y) uHalf[((x)+1)*Ny_loc + (y)]
#define VHALF(x,y) vHalf[((x)+1)*Ny_loc + (y)]

void GrayScott::step()
{
    // update step
    if (world.rank == 0) {
        ++currStep_;
    }
    
    
    MPI_Request request[8];
    MPI_Status status[8];
    
    // exchange boundaries along y-direction
    if (world.coord_x % 2 == 0) { // first send top, then botton
        
        MPI_Isend(&u_[(Nx_loc)*(Ny_loc)],   1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[0]);
        MPI_Irecv(&u_[(Nx_loc+1)*(Ny_loc)], 1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[1]);
        MPI_Isend(&u_[(Ny_loc)],            1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[2]);
        MPI_Irecv(&u_[0],                   1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[3]);
        
        MPI_Isend(&v_[(Nx_loc)*(Ny_loc)],   1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[4]);
        MPI_Irecv(&v_[(Nx_loc+1)*(Ny_loc)], 1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[5]);
        MPI_Isend(&v_[(Ny_loc)],            1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[6]);
        MPI_Irecv(&v_[0],                   1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[7]);
    }
    else { // first send botton, then top
        
        MPI_Irecv(&u_[0],                   1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[0]);
        MPI_Isend(&u_[(Ny_loc)],            1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[1]);
        MPI_Irecv(&u_[(Nx_loc+1)*(Ny_loc)], 1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[2]);
        MPI_Isend(&u_[(Nx_loc)*(Ny_loc)],   1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[3]);
        
        MPI_Irecv(&v_[0],                   1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[4]);
        MPI_Isend(&v_[(Ny_loc)],            1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[5]);
        MPI_Irecv(&v_[(Nx_loc+1)*(Ny_loc)], 1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[6]);
        MPI_Isend(&v_[(Nx_loc)*(Ny_loc)],   1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[7]);
    }
    
    
    
    
    // u and v at the half step
    std::vector<double> uHalf(u_.size());
    std::vector<double> vHalf(v_.size());
    
    // right hand sides for u and for v
    std::vector<double> uRhs(N_);
    std::vector<double> vRhs(N_);
    
    double uCoeff = Du_*dt_/(2.*dx_*dx_);
    double vCoeff = Dv_*dt_/(2.*dx_*dx_);
    
    
    /****************** DIFFUSION (ADI) ***************************************/
    
    // perform the first half-step
    // loop over all rows
    
    
    // parallelize outer loop (y-direction) with mpi, inner loop with openmp (TODO)
    
    // inner grid points
    for (int i=1; i<Nx_loc-1; ++i) {
        // create right-hand side of the systems
        for (int j=0; j<N_; ++j) {
            uRhs[j] = U(i,j) + uCoeff * (U(i+1,j) - 2.*U(i,j) + U(i-1,j));
            vRhs[j] = V(i,j) + vCoeff * (V(i+1,j) - 2.*V(i,j) + V(i-1,j));
        }
        
        TriDiagMatrixSolver::solve(N_, matU1_, uRhs, &UHALF(i,0), 1);
        TriDiagMatrixSolver::solve(N_, matV1_, vRhs, &VHALF(i,0), 1);
    }
    
    
    
    // wait for boundaries to arrive
    MPI_Waitall(8,request,status);
    
    
    // update local boundaries

    if (world.rank == 0) {
        // i=0 local and global
        for (int j=0; j<N_; ++j) {
            uRhs[j] = U(0,j) + uCoeff * (U(1,j) - U(0,j));
            vRhs[j] = V(0,j) + vCoeff * (V(1,j) - V(0,j));
        }
        TriDiagMatrixSolver::solve(N_, matU1_, uRhs, &UHALF(0,0), 1);
        TriDiagMatrixSolver::solve(N_, matV1_, vRhs, &VHALF(0,0), 1);
    }
    else {
        // i=0 local
        for (int j=0; j<N_; ++j) {
            uRhs[j] = U(0,j) + uCoeff * (U(0+1,j) - 2.*U(0,j) + U(0-1,j));
            vRhs[j] = V(0,j) + vCoeff * (V(0+1,j) - 2.*V(0,j) + V(0-1,j));
        }
        
        TriDiagMatrixSolver::solve(N_, matU1_, uRhs, &UHALF(0,0), 1);
        TriDiagMatrixSolver::solve(N_, matV1_, vRhs, &VHALF(0,0), 1);
    }
    
    
    if (world.rank == world.size-1) {
        // i=Nx_loc-1 local and i=N_-1 global
        for (int j=0; j<N_; ++j) {
            uRhs[j] = U(Nx_loc-1,j) + uCoeff * (- U(Nx_loc-1,j) + U(Nx_loc-2,j));
            vRhs[j] = V(Nx_loc-1,j) + vCoeff * (- V(Nx_loc-1,j) + V(Nx_loc-2,j));
        }
        TriDiagMatrixSolver::solve(N_, matU1_, uRhs, &UHALF(Nx_loc-1,0), 1);
        TriDiagMatrixSolver::solve(N_, matV1_, vRhs, &VHALF(Nx_loc-1,0), 1);
    }
    else {
        // i=Nx_loc-1 local and global
        for (int j=0; j<N_; ++j) {
            uRhs[j] = U(Nx_loc-1,j) + uCoeff * (U(Nx_loc-1+1,j) - 2.*U(Nx_loc-1,j) + U(Nx_loc-1-1,j));
            vRhs[j] = V(Nx_loc-1,j) + vCoeff * (V(Nx_loc-1+1,j) - 2.*V(Nx_loc-1,j) + V(Nx_loc-1-1,j));
        }
        
        TriDiagMatrixSolver::solve(N_, matU1_, uRhs, &UHALF(Nx_loc-1,0), 1);
        TriDiagMatrixSolver::solve(N_, matV1_, vRhs, &VHALF(Nx_loc-1,0), 1);
    }
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    // transpose matrix
    
    // transpose global blocks (send from uHalf to u_)
    // start at Ny_loc, because we ignore the ghost cells
    MPI_Alltoall(&uHalf[Ny_loc], 1, block_resized, &u_[Ny_loc], 1, block_resized, MPI_COMM_WORLD);
    MPI_Alltoall(&vHalf[Ny_loc], 1, block_resized, &v_[Ny_loc], 1, block_resized, MPI_COMM_WORLD);
    
    // locally transpose blocks
    // loop over blocks TODO parallelize block loop with openmp
    int Nb_loc = Ny_loc/Nx_loc;
    int ind1, ind2;
//    int tmp2;
    for (int b=0; b<Nb_loc; ++b) {
        for (int i=0; i<Nx_loc; ++i) {
            for (int j=0; j<i; ++j) {
                ind1 = (i+1)*Ny_loc + j + b*Nx_loc; // regular index + offset of block
                ind2 = (j+1)*Ny_loc + i + b*Nx_loc; // switch i and j
                
                std::swap(u_[ind1], u_[ind2]);
                std::swap(v_[ind1], v_[ind2]);
            }
        }
    }
    
    
    
    
    
    // exchange new boundaries
    if (world.coord_x % 2 == 0) { // first send top, then botton
        
        MPI_Isend(&u_[(Nx_loc)*(Ny_loc)],   1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[0]);
        MPI_Irecv(&u_[(Nx_loc+1)*(Ny_loc)], 1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[1]);
        MPI_Isend(&u_[(Ny_loc)],            1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[2]);
        MPI_Irecv(&u_[0],                   1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[3]);
        
        MPI_Isend(&v_[(Nx_loc)*(Ny_loc)],   1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[4]);
        MPI_Irecv(&v_[(Nx_loc+1)*(Ny_loc)], 1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[5]);
        MPI_Isend(&v_[(Ny_loc)],            1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[6]);
        MPI_Irecv(&v_[0],                   1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[7]);
    }
    else { // first send botton, then top
        
        MPI_Irecv(&u_[0],                   1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[0]);
        MPI_Isend(&u_[(Ny_loc)],            1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[1]);
        MPI_Irecv(&u_[(Nx_loc+1)*(Ny_loc)], 1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[2]);
        MPI_Isend(&u_[(Nx_loc)*(Ny_loc)],   1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[3]);
        
        MPI_Irecv(&v_[0],                   1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[4]);
        MPI_Isend(&v_[(Ny_loc)],            1, top_boundary,    world.top_proc,    TAG, cart_comm, &request[5]);
        MPI_Irecv(&v_[(Nx_loc+1)*(Ny_loc)], 1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[6]);
        MPI_Isend(&v_[(Nx_loc)*(Ny_loc)],   1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[7]);
    }

    

    // inner grid points
    for (int i=1; i<Nx_loc-1; ++i) {
        // create right-hand side of the systems
        for (int j=0; j<N_; ++j) {
            uRhs[j] = U(i,j) + uCoeff * (U(i+1,j) - 2.*U(i,j) + U(i-1,j));
            vRhs[j] = V(i,j) + vCoeff * (V(i+1,j) - 2.*V(i,j) + V(i-1,j));
        }
        
        TriDiagMatrixSolver::solve(N_, matU1_, uRhs, &UHALF(i,0), 1);
        TriDiagMatrixSolver::solve(N_, matV1_, vRhs, &VHALF(i,0), 1);
    }
    
    
    // wait for boundaries to arrive
    MPI_Waitall(8,request,status);
    
    
    // update local boundaries

    // top
    if (world.rank == 0) {
        // i=0 local and global
        for (int j=0; j<N_; ++j) {
            uRhs[j] = U(0,j) + uCoeff * (U(1,j) - U(0,j));
            vRhs[j] = V(0,j) + vCoeff * (V(1,j) - V(0,j));
        }
        TriDiagMatrixSolver::solve(N_, matU1_, uRhs, &UHALF(0,0), 1);
        TriDiagMatrixSolver::solve(N_, matV1_, vRhs, &VHALF(0,0), 1);
    }
    else {
        // i=0 local, but not globally
        for (int j=0; j<N_; ++j) {
            uRhs[j] = U(0,j) + uCoeff * (U(0+1,j) - 2.*U(0,j) + U(0-1,j));
            vRhs[j] = V(0,j) + vCoeff * (V(0+1,j) - 2.*V(0,j) + V(0-1,j));
        }
        
        TriDiagMatrixSolver::solve(N_, matU1_, uRhs, &UHALF(0,0), 1);
        TriDiagMatrixSolver::solve(N_, matV1_, vRhs, &VHALF(0,0), 1);
    }
    
    // bottom
    if (world.rank == world.size-1) {
        // i=Nx_loc-1 local and i=N_-1 global
        for (int j=0; j<N_; ++j) {
            uRhs[j] = U(Nx_loc-1,j) + uCoeff * (- U(Nx_loc-1,j) + U(Nx_loc-2,j));
            vRhs[j] = V(Nx_loc-1,j) + vCoeff * (- V(Nx_loc-1,j) + V(Nx_loc-2,j));
        }
        TriDiagMatrixSolver::solve(N_, matU1_, uRhs, &UHALF(Nx_loc-1,0), 1);
        TriDiagMatrixSolver::solve(N_, matV1_, vRhs, &VHALF(Nx_loc-1,0), 1);
    }
    else {
        // i=Nx_loc-1 local
        for (int j=0; j<N_; ++j) {
            uRhs[j] = U(Nx_loc-1,j) + uCoeff * (U(Nx_loc-1+1,j) - 2.*U(Nx_loc-1,j) + U(Nx_loc-1-1,j));
            vRhs[j] = V(Nx_loc-1,j) + vCoeff * (V(Nx_loc-1+1,j) - 2.*V(Nx_loc-1,j) + V(Nx_loc-1-1,j));
        }
        
        TriDiagMatrixSolver::solve(N_, matU1_, uRhs, &UHALF(Nx_loc-1,0), 1);
        TriDiagMatrixSolver::solve(N_, matV1_, vRhs, &VHALF(Nx_loc-1,0), 1);
    }
    
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    // transpose back
    
    // transpose global blocks (send from uHalf to u_)
    // start at Ny_loc, because we ignore the ghost cells
    MPI_Alltoall(&uHalf[Ny_loc], 1, block_resized, &u_[Ny_loc], 1, block_resized, MPI_COMM_WORLD);
    MPI_Alltoall(&vHalf[Ny_loc], 1, block_resized, &v_[Ny_loc], 1, block_resized, MPI_COMM_WORLD);
    
    // locally transpose blocks
    // loop over blocks TODO parallelize block loop with openmp
//    int Nb_loc = Ny_loc/Nx_loc;
//    int ind1, ind2;
//    int tmp2;
    for (int b=0; b<Nb_loc; ++b) {
        for (int i=0; i<Nx_loc; ++i) {
            for (int j=0; j<i; ++j) {
                ind1 = (i+1)*Ny_loc + j + b*Nx_loc; // regular index + offset of block
                ind2 = (j+1)*Ny_loc + i + b*Nx_loc; // switch i and j
                
                std::swap(u_[ind1], u_[ind2]);
                std::swap(v_[ind1], v_[ind2]);
            }
        }
    }
    
    
    
    /****************** REACTION **********************************************/

    double uind, vind;
    for (int j=0; j<Ny_loc; ++j) {
        for (int i=0; i<Nx_loc; ++i) {
//            const int ind = (i+1)*Ny_loc + j;
//            uind = u_[ind];
//            vind = v_[ind];
//            u_[ind] += dt_ * ( -uind*vind*vind + F_*(1.-uind) );
//            v_[ind] += dt_ * ( uind*vind*vind - (F_+k_)*vind );
            uind = U(i,j);
            vind = V(i,j);
            U(i,j) += dt_ * ( -uind*vind*vind + F_*(1.-uind) );
            V(i,j) += dt_ * ( uind*vind*vind - (F_+k_)*vind );
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
}



void GrayScott::initialize_fields()
{
    // domain [-1,1]²
    
    // initialize with all U
    u_.clear();
    u_.resize((Nx_loc+2) * (Ny_loc), 0.);
    v_.clear();
    v_.resize((Nx_loc+2) * (Ny_loc), 0.);
    
    
    
    std::mt19937 rng(world.rank+1); // different rng for each processor -> problem?
    std::normal_distribution<double> dist(0,1);
    
    double chi;
    double x, y;
    
    for (int i=0; i<Nx_loc; ++i) {
        for (int j=0; j<Ny_loc; ++j) {
            x = -1. + (double)i*dx_ + xmin_loc;
            y = -1. + (double)j*dx_ + ymin_loc;
            
            chi = 0.;
            if (x>=-0.2 && x<=0.2 && y>=-0.2 && y<=0.2) {
                chi = 1.;
            }
            
            U(i,j) = (1.-chi) + chi*(0.5 + dist(rng)/100.);
            V(i,j) = chi * (0.25 + dist(rng)/100.);
        }
    }
}


void GrayScott::save_fields()
{	
    char stepString[10];
    std::sprintf(stepString, "%05d_", currStep_);
    
    std::string fname(dirPath_ + stepString + "u.dat");
    std::ofstream uOut(fname);
    
    fname = dirPath_ + stepString + "v.dat";
    std::ofstream vOut(fname);
    
    
    int w = 12;
    
    for (int j=0; j<N_; ++j) {
        for (int i=0; i<N_; ++i) {
            uOut << std::setw(w) << U(i,j) << " ";
            vOut << std::setw(w) << V(i,j) << " ";
        }
        uOut << "\n";
        vOut << "\n";
    }
    
    uOut.close();
    vOut.close();
}




namespace png {
double red( double gray ) {
    if (gray < 0.6)
        return 637.5 * std::max(gray-0.2, 0.);
    else
        return 255;
}
double green( double gray ) {
    if (gray < 0.6)
        return 637.5 * std::max(gray-0.2, 0.);
    else
        return -637.5 * (gray-1.);
}
double blue( double gray ) {
    if (gray < 0.6)
        return -637.5 * (gray-0.6);
    else
        return 0.;
}
}



#define UA(x,y) uAll[(x)*N_ + (y)]
void GrayScott::save_png()
{
    // first send all data to proc 0 -> gather
//    if (world.rank == 0) {
        std::vector<double> uAll(N_*N_);
//    }
    
    MPI_Gather(&u_[Ny_loc], Nx_loc*Ny_loc, MPI_DOUBLE, &uAll[0], Nx_loc*Ny_loc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
    if (world.rank == 0) {
    
        char parameters[50];
        std::sprintf(parameters, "_k_%1.03f_F_%1.03f_step_%06d", k_, F_, currStep_);
        std::string fname = "images/" + pngName_ + parameters + ".png";
        
        char filename[1024];
        strncpy(filename, fname.c_str(), sizeof(filename));
        filename[sizeof(filename) - 1] = 0;
        
        // 4 values per pixel: RGBA
        std::vector<unsigned char> data(4*N_*N_);
        
        int idx = 0;
        for (int i=0; i<N_; ++i) {
            for (int j=0; j<N_; ++j) {
                double color = UA(i,j);
                
                data[idx] = png::red(color); ++idx;
                data[idx] = png::green(color); ++idx;
                data[idx] = png::blue(color); ++idx;
                data[idx] = 0xFFu; ++idx;
            }
        }
        
        lodepng::encode(filename, data, N_, N_);
    }
}










