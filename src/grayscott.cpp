#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <ctime>
#include <boost/filesystem.hpp>

#include "grayscott.hpp"
#include "tridiagmatrixsolver.hpp"
#include "periodictridiagmatrixsolver.hpp"


//#define U(x,y) u_[(x) + (y)*N_]
//#define V(x,y) v_[(x) + (y)*N_]

// periodic boundary conditions
#define U(x,y) u_[((x)+N_)%N_ + (((y)+N_)%N_)*N_]
#define V(x,y) v_[((x)+N_)%N_ + (((y)+N_)%N_)*N_]


GrayScott::GrayScott(int N, double L, double dt, double Du, double Dv, double F, double k, int nSteps)
    : N_(N)
    , Ntot_(N*N)
    , L_(L)
    , dx_((double) L / (double) N)
    , dt_(dx_*dx_ / (2.*std::max(Du,Dv)))
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
    std::cout << "dt = " << dt_ << "\n";
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
    
    initialize_fields();
}


GrayScott::~GrayScott()
{
    boost::filesystem::path dir(dirPath_);
    if (boost::filesystem::exists(dir) && boost::filesystem::is_empty(dir)) {
        boost::filesystem::remove(dir);
    }
}


void GrayScott::run()
{
    save_fields();
    
    for (int i=0; i<nSteps_; ++i) {
        step();
        
        if (i == 0) {
            save_fields();
        }
    }
    
    
    save_fields();
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
            uRhs[i] = U(i,j);// + Du_*dt_/(2.*dx_*dx_) * (U(i,j+1) - 2.*U(i,j) + U(i,j-1)) + dt_/2. * (-U(i,j)*V(i,j)*V(i,j) + F_*(1-U(i,j)));
            vRhs[i] = V(i,j);// + Dv_*dt_/(2.*dx_*dx_) * (V(i,j+1) - 2.*V(i,j) + V(i,j-1)) + dt_/2. * (U(i,j)*V(i,j)*V(i,j) - (F_+k_)*V(i,j));
        }
        
        PeriodicTriDiagMatrixSolver::solve(N_, matU1_, uRhs, &uHalf[j], 1);
        PeriodicTriDiagMatrixSolver::solve(N_, matV1_, vRhs, &vHalf[j], 1);
        
        for (int i=0; i<N_; ++i) {
            UHALF(i,j) += Du_*dt_/(2.*dx_*dx_) * (U(i,j+1) - 2.*U(i,j) + U(i,j-1)) + dt_/2. * (-U(i,j)*V(i,j)*V(i,j) + F_*(1-U(i,j)));
            VHALF(i,j) += Dv_*dt_/(2.*dx_*dx_) * (V(i,j+1) - 2.*V(i,j) + V(i,j-1)) + dt_/2. * (U(i,j)*V(i,j)*V(i,j) - (F_+k_)*V(i,j));
        }
    }
    
    
    
    
    // MPI TRANSPOSE MATRIX (ALL-TO-ALL) ?
    
    
    
    // perform the second half-step
    // loop over all columns
    for (int i=0; i<N_; ++i) {
        // create right-hand side of the systems
        for (int j=0; j<N_; ++j) {
            uRhs[j] = UHALF(i,j);// + Du_*dt_/(2.*dx_*dx_) * (UHALF(i+1,j) - 2.*UHALF(i,j) + UHALF(i-1,j)) + dt_/2. * (-UHALF(i,j)*VHALF(i,j)*VHALF(i,j) + F_*(1-UHALF(i,j)));
            vRhs[j] = VHALF(i,j);// + Dv_*dt_/(2.*dx_*dx_) * (VHALF(i+1,j) - 2.*VHALF(i,j) + VHALF(i-1,j)) + dt_/2. * (UHALF(i,j)*VHALF(i,j)*VHALF(i,j) - (F_+k_)*VHALF(i,j));
        }
        
        PeriodicTriDiagMatrixSolver::solve(N_, matU2_, uRhs, &u_[i], N_);
        PeriodicTriDiagMatrixSolver::solve(N_, matV2_, vRhs, &v_[i], N_);
        
        for (int j=0; j<N_; ++j) {
            U(i,j) += Du_*dt_/(2.*dx_*dx_) * (UHALF(i+1,j) - 2.*UHALF(i,j) + UHALF(i-1,j)) + dt_/2. * (-UHALF(i,j)*VHALF(i,j)*VHALF(i,j) + F_*(1-UHALF(i,j)));
            V(i,j) += Dv_*dt_/(2.*dx_*dx_) * (VHALF(i+1,j) - 2.*VHALF(i,j) + VHALF(i-1,j)) + dt_/2. * (UHALF(i,j)*VHALF(i,j)*VHALF(i,j) - (F_+k_)*VHALF(i,j));
        }
    }
    
//    std::cout << "step done\n";
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
    
    
//    
//    // Source: https://github.com/derekrb/gray-scott/blob/master/main.py
//    std::mt19937 rng(42);
//    std::uniform_real_distribution<double> dist(0,1);

//    int perturbNum = 1;
//    double perturbU = 0.5;
//    double perturbV = 0.25;
//    double perturbMag = 0.5; // how much of the domain has perturbations? (in each dimension)
//    
//    // Apply the specified number of randomly-sized perturbations
//    for (int i=0; i<perturbNum; ++i) {

//        int xStart = dist(rng) * N_ * 0.9;
//        int yStart = dist(rng) * N_ * 0.9;
//        int xEnd = xStart + (dist(rng)*N_ + N_) * perturbMag;
//        int yEnd = yStart + (dist(rng)*N_ + N_) * perturbMag;
//        
//        // Apply perturbations
//        for (int x=xStart; x<xEnd; ++x) {
//            for (int y=yStart; y<yEnd; ++y) {

//                // Constant perturbation
//                U(x,y) = perturbU;
//                V(x,y) = perturbV;

//                // Random perturbation
//                double perturb = 0.01 * dist(rng);
//                U(x,y) -= perturb;
//                V(x,y) += perturb;
//            }
//        }
//   }
    
//    for (int y=0; y<N_/3; ++y) {
//        for (int x=0; x<N_/3; ++x) {
//            U(x,y) = 2.;
//            V(x,y) = 4.;
//        }
//    }
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







