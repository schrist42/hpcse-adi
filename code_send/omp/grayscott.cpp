#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <ctime>
#include <boost/filesystem.hpp>
#include <chrono>
#include <omp.h>

#include "grayscott.hpp"
#include "tridiagmatrixsolver.hpp"
#include "lodepng.h" // NOT loadpng


#define U(x,y) u_[(x) + (y)*N_]
#define V(x,y) v_[(x) + (y)*N_]


GrayScott::GrayScott(int N, double L, double dt, double Du, double Dv, double F, double k, int nSteps, std::string pngname, unsigned int nthreads)
    : N_(N)
    , Ntot_(N*N)
    , dx_((double) L / (double) N)
    , dt_(dt)
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
    , nthreads_(nthreads)
{
    initialize_fields();
}


GrayScott::~GrayScott()
{}


void GrayScott::run()
{
	// timer init and start
	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	start = std::chrono::high_resolution_clock::now();
    
    for (int i=0; i<nSteps_; ++i) {
        step();

        if (!(i%10000)) {
            std::cout << "Step " << i+1 << " done!\n";
        }
    }
    // timer end and output of time
	end = std::chrono::high_resolution_clock::now();
	double elapsed = std::chrono::duration<double>(end-start).count();  
	  
    #pragma omp parallel num_threads(nthreads_)
    #pragma omp single
    std::cout << "# omp threads = " << omp_get_num_threads() << std::endl;    
    std::cout << "\n";
    std::cout << "exec time: " << '\t' << elapsed << std::endl;
    std::cout << "\n";
    std::cout << "size: " << '\t' << N_ << std::endl;
    std::cout << "\n";
    
    save_png();
}


void GrayScott::benchmark()
{
    double elapsed;
    double avg_time = 0.;
    double avg_squared_time = 0.;
    
    // timer initialization
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    
    
     for (int i=0; i<nSteps_; ++i) {
	    // timer start
	    start = std::chrono::high_resolution_clock::now();
        
        step();
        
      	// timer end
	    end = std::chrono::high_resolution_clock::now();
	    elapsed = std::chrono::duration<double>(end-start).count();    
	    
	    avg_time += elapsed;
        avg_squared_time += elapsed*elapsed;
    }
    
    avg_time /= (double)nSteps_;
    avg_squared_time /= (double)nSteps_;
    
    double err = std::sqrt((avg_squared_time - avg_time*avg_time)/(double)nSteps_);
    
    
    #pragma omp parallel num_threads(nthreads_)
    #pragma omp single
    std::cout << omp_get_num_threads() << "\t";
    std::cout << avg_time << "\t";
    std::cout << N_ << "\t";
    std::cout << err << "\n";
}


#define UHALF(x,y) uHalf[(x) + (y)*N_]
#define VHALF(x,y) vHalf[(x) + (y)*N_]

void GrayScott::step()
{
    // update step
    ++currStep_;
    
    // u and v at the half step
    std::vector<double> uHalf(Ntot_);
    std::vector<double> vHalf(Ntot_);
    
    // right hand sides for u and for v
    std::vector<double> uRhs(N_);
    std::vector<double> vRhs(N_);
    
    double uCoeff = Du_*dt_/(2.*dx_*dx_);
    double vCoeff = Dv_*dt_/(2.*dx_*dx_);
    
    
    /****************** DIFFUSION (ADI) ***************************************/
    
    // perform the first half-step
    // loop over all rows
    
    // j=0
    for (int i=0; i<N_; ++i) {
        uRhs[i] = U(i,0) + uCoeff * (U(i,1) - U(i,0));
        vRhs[i] = V(i,0) + vCoeff * (V(i,1) - V(i,0));
    }
    TriDiagMatrixSolver::solve(N_, matU1_, uRhs, &UHALF(0,0), 1);
    TriDiagMatrixSolver::solve(N_, matV1_, vRhs, &VHALF(0,0), 1);
    
    // inner grid points

    #pragma omp parallel num_threads(nthreads_)
    { // PARALLEL REGION BEGIN
    // right hand sides for u and for v
    std::vector<double> puRhs(N_);
    std::vector<double> pvRhs(N_);

    #pragma omp for
    for (int j=1; j<N_-1; ++j) {
        // create right-hand side of the systems
        for (int i=0; i<N_; ++i) {
            puRhs[i] = U(i,j) + uCoeff * (U(i,j+1) - 2.*U(i,j) + U(i,j-1));
            pvRhs[i] = V(i,j) + vCoeff * (V(i,j+1) - 2.*V(i,j) + V(i,j-1));
        }
        
        TriDiagMatrixSolver::solve(N_, matU1_, puRhs, &UHALF(0,j), 1);
        TriDiagMatrixSolver::solve(N_, matV1_, pvRhs, &VHALF(0,j), 1);
    }
    } // PARALLEL REGION END
    // j=N_-1
    for (int i=0; i<N_; ++i) {
        uRhs[i] = U(i,N_-1) + uCoeff * (- U(i,N_-1) + U(i,N_-2));
        vRhs[i] = V(i,N_-1) + vCoeff * (- V(i,N_-1) + V(i,N_-2));
    }
    TriDiagMatrixSolver::solve(N_, matU1_, uRhs, &UHALF(0,N_-1), 1);
    TriDiagMatrixSolver::solve(N_, matV1_, vRhs, &VHALF(0,N_-1), 1);
    
    // perform the second half-step
    // loop over all columns
    
    // i=0
    for (int j=0; j<N_; ++j) {
        uRhs[j] = U(0,j) + uCoeff * (U(1,j) - U(0,j));
        vRhs[j] = V(0,j) + vCoeff * (V(1,j) - V(0,j));
    }
    TriDiagMatrixSolver::solve(N_, matU2_, uRhs, &U(0,0), N_);
    TriDiagMatrixSolver::solve(N_, matV2_, vRhs, &V(0,0), N_);
    
    // inner grid points

    #pragma omp parallel num_threads(nthreads_)
    {  // PARALLEL REGION BEGIN
    // right hand sides for u and for v
    std::vector<double> puRhs(N_);
    std::vector<double> pvRhs(N_);

    #pragma omp for
    for (int i=1; i<N_-1; ++i) {
        // create right-hand side of the systems

        for (int j=0; j<N_; ++j) {
            puRhs[j] = UHALF(i,j) + uCoeff * (UHALF(i+1,j) - 2.*UHALF(i,j) + UHALF(i-1,j));
            pvRhs[j] = VHALF(i,j) + vCoeff * (VHALF(i+1,j) - 2.*VHALF(i,j) + VHALF(i-1,j));
        }
        
        TriDiagMatrixSolver::solve(N_, matU2_, puRhs, &U(i,0), N_);
        TriDiagMatrixSolver::solve(N_, matV2_, pvRhs, &V(i,0), N_);
    }
    } // PARALLEL REGION END
    // i=N_-1
    for (int j=0; j<N_; ++j) {
        uRhs[j] = U(N_-1,j) + uCoeff * (- U(N_-1,j) + U(N_-2,j));
        vRhs[j] = V(N_-1,j) + vCoeff * (- V(N_-1,j) + V(N_-2,j));
    }
    TriDiagMatrixSolver::solve(N_, matU2_, uRhs, &U(N_-1,0), N_);
    TriDiagMatrixSolver::solve(N_, matV2_, vRhs, &V(N_-1,0), N_);
    
    
    
    /****************** REACTION **********************************************/
    #pragma omp parallel num_threads(nthreads_)
    { // PARALLEL REGION BEGIN
    double tmp;
    #pragma omp for collapse(2)
    for (int i=0; i<N_; ++i) {
        for (int j=0; j<N_; ++j) {
            tmp = U(i,j);
            U(i,j) += dt_ * ( -tmp*V(i,j)*V(i,j) + F_*(1.-tmp) );
            V(i,j) += dt_ * ( tmp*V(i,j)*V(i,j) - (F_+k_)*V(i,j) );
        }
    }
    } // PARALLEL REGION END
}



void GrayScott::initialize_fields()
{
    // domain [-1,1]²
    
    // initialize with all U
    u_.clear();
    u_.resize(Ntot_);
    v_.clear();
    v_.resize(Ntot_);
    
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(0,1);
    
    double chi;
    double x, y;

    for (int i=0; i<N_; ++i) {
        for (int j=0; j<N_; ++j) {
            x = -1 + (double)i*dx_;
            y = -1 + (double)j*dx_;
            
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


void GrayScott::save_png()
{
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
            double color = U(i,j);
            
            data[idx] = png::red(color); ++idx;
            data[idx] = png::green(color); ++idx;
            data[idx] = png::blue(color); ++idx;
            data[idx] = 0xFFu; ++idx;
        }
    }
    
    lodepng::encode(filename, data, N_, N_);
}










