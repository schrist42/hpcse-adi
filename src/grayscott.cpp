#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <ctime>
#include <boost/filesystem.hpp>

#include "grayscott.hpp"
#include "tridiagmatrixsolver.hpp"


//#define U(x,y) u_[(x) + (y)*N_]
//#define V(x,y) v_[(x) + (y)*N_]

// periodic boundary conditions
#define U(x,y) u_[((x)+N_)%N_ + (((y)+N_)%N_)*N_]
#define V(x,y) v_[((x)+N_)%N_ + (((y)+N_)%N_)*N_]


GrayScott::GrayScott(int N, double L, double dt, double Du, double Dv, double F, double k, int nSteps)
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
{
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
	if (boost::filesystem::create_directory(dir)) {
//		std::cout << "Success" << "\n";
	}
    
    initialize_fields();
    
    print_fields();
    
}


void GrayScott::run()
{
    for (int i=0; i<nSteps_; ++i) {
        step();
        
        if (i == 0) {
            print_fields();
        }
    }
    
    
    print_fields();
}



//#define UHALF(x,y) uHalf[(x) + (y)*N_]
//#define VHALF(x,y) vHalf[(x) + (y)*N_]

// periodic boundary conditions
#define UHALF(x,y) uHalf[((x)+N_)%N_ + (((y)+N_)%N_)*N_]
#define VHALF(x,y) vHalf[((x)+N_)%N_ + (((y)+N_)%N_)*N_]

void GrayScott::step()
{
//    std::cout << "step\n";
    // update step
    ++currStep_;
    
    // u and v at the half step
    std::vector<double> uHalf(Ntot_);
    std::vector<double> vHalf(Ntot_);
    
    std::vector<double> uRhs(N_);
    std::vector<double> vRhs(N_);
    
    // perform the first half-step
    // loop over all rows
    for (int j=0; j<N_; ++j) {
        // create right-hand side of the systems
        for (int i=0; i<N_; ++i) {
            uRhs[j] = U(i,j) + Du_*dt_/(2.*dx_*dx_) * ((j+1<N_ ? U(i,j+1) : 0) - 2.*U(i,j) + (j-1>=0 ? U(i,j-1) : 0)) + dt_/2. * (-U(i,j)*V(i,j)*V(i,j) + F_*(1-U(i,j)));
            vRhs[j] = V(i,j) + Dv_*dt_/(2.*dx_*dx_) * ((j+1<N_ ? V(i,j+1) : 0) - 2.*V(i,j) + (j-1>=0 ? V(i,j-1) : 0)) + dt_/2. * (U(i,j)*V(i,j)*V(i,j) - (F_+k_)*V(i,j));
        }
        
        TriDiagMatrixSolver::solve(matU1_, uRhs, &uHalf[j], 1);
        TriDiagMatrixSolver::solve(matV1_, vRhs, &vHalf[j], 1);
    }
    
    
    
    // MPI TRANSPOSE MATRIX (ALL-TO-ALL) ?
    
    
    
    // perform the second half-step
    // loop over all columns
    for (int i=0; i<N_; ++i) {
        // create right-hand side of the systems
        for (int j=0; j<N_; ++j) {
            uRhs[i] = UHALF(i,j) + Du_*dt_/(2.*dx_*dx_) * ((i+1<N_ ? UHALF(i+1,j) : 0) - 2.*UHALF(i,j) + (i-1>=0 ? UHALF(i-1,j) : 0)) + dt_/2. * (-UHALF(i,j)*VHALF(i,j)*VHALF(i,j) + F_*(1-UHALF(i,j)));
            vRhs[i] = VHALF(i,j) + Dv_*dt_/(2.*dx_*dx_) * ((i+1<N_ ? VHALF(i+1,j) : 0) - 2.*VHALF(i,j) + (i-1>=0 ? VHALF(i-1,j) : 0)) + dt_/2. * (UHALF(i,j)*VHALF(i,j)*VHALF(i,j) - (F_+k_)*VHALF(i,j));
        }
        
        TriDiagMatrixSolver::solve(matU2_, uRhs, &u_[i], N_);
        TriDiagMatrixSolver::solve(matV2_, vRhs, &v_[i], N_);
    }
}



void GrayScott::initialize_fields()
{
    // initialize with all U
    u_.clear();
    u_.resize(Ntot_,1.);
    v_.clear();
    v_.resize(Ntot_,0.);
    
    
    
    // Source: https://github.com/derekrb/gray-scott/blob/master/main.py
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0,1);

    int perturbNum = 1;
    double perturbU = 0.5;
    double perturbV = 0.25;
    double perturbMag = 0.5; // how much of the domain has perturbations? (in each dimension)
    
    // Apply the specified number of randomly-sized perturbations
    for (int i=0; i<perturbNum; ++i) {

        int xStart = dist(rng) * N_ * 0.9;
        int yStart = dist(rng) * N_ * 0.9;
        int xEnd = xStart + (dist(rng)*N_ + N_) * perturbMag;
        int yEnd = yStart + (dist(rng)*N_ + N_) * perturbMag;
        
        // Apply perturbations
        for (int x=xStart; x<xEnd; ++x) {
            for (int y=yStart; y<yEnd; ++y) {

                // Constant perturbation
                U(x,y) = perturbU;
                V(x,y) = perturbV;

                // Random perturbation
                double perturb = 0.01 * dist(rng);
                U(x,y) -= perturb;
                V(x,y) += perturb;
            }
        }
   }
    
//    for (int y=0; y<N_/3; ++y) {
//        for (int x=0; x<N_/3; ++x) {
//            U(x,y) = 2.;
//            V(x,y) = 4.;
//        }
//    }
}


void GrayScott::print_fields()
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







