#include <iostream>
#include <fstream>

#include "diffusion.hpp"
#include "tridiagmatrixsolver.hpp"


#define U(x,y) u_[(x) + (y)*N_]
#define V(x,y) v_[(x) + (y)*N_]


Diffusion::Diffusion(int N, double L, double dt, double Du, double Dv, double F, double k, int nSteps)
    : N_(N)
    , Ntot_(N*N)
//    , L_(L)
    , dx_((double) L / (double) N)
    , dt_(dt)
    , nSteps_(nSteps)
    , Du_(Du)
    , Dv_(Dv)
    , F_(F)
    , k_(k)
    , matU1_(N, -Du*dt/(2.*dx_*dx_), 1.+Du*dt/(dx_*dx_), -Du*dt/(2.*dx_*dx_))
    , matU2_(N, -Du*dt/(2.*dx_*dx_), 1.+Du*dt/(dx_*dx_), -Du*dt/(2.*dx_*dx_))
    , matV1_(N, -Dv*dt/(2.*dx_*dx_), 1.+Dv*dt/(dx_*dx_), -Dv*dt/(2.*dx_*dx_))
    , matV2_(N, -Dv*dt/(2.*dx_*dx_), 1.+Dv*dt/(dx_*dx_), -Dv*dt/(2.*dx_*dx_))
{
//    u_.resize(Ntot_,0.);
//    v_.resize(Ntot_,0.);

    
    initialize_fields();
    
    std::ofstream uOut("data/uInit.dat");
    std::ofstream vOut("data/vInit.dat");
    
    for (int j=0; j<N_; ++j) {
        for (int i=0; i<N_; ++i) {
            uOut << U(i,j) << " ";
            vOut << V(i,j) << " ";
        }
        uOut << "\n";
        vOut << "\n";
    }
    
    uOut.close();
    vOut.close();
    
}


void Diffusion::run()
{
    for (int i=0; i<nSteps_; ++i) {
        step();
    }
    
    
    std::ofstream uOut("data/u.dat");
    std::ofstream vOut("data/v.dat");
    
    for (int j=0; j<N_; ++j) {
        for (int i=0; i<N_; ++i) {
            uOut << U(i,j) << " ";
            vOut << V(i,j) << " ";
        }
        uOut << "\n";
        vOut << "\n";
    }
    
    uOut.close();
    vOut.close();
}



#define UHALF(x,y) uHalf[(x) + (y)*N_]
#define VHALF(x,y) vHalf[(x) + (y)*N_]

void Diffusion::step()
{
    // temporary u and v for the result
    std::vector<double> uHalf(Ntot_,0.);
    std::vector<double> vHalf(Ntot_,0.);
    
    std::vector<double> uRhs(N_,0.);
    std::vector<double> vRhs(N_,0.);
    
    // perform the first half-step
    // loop over all rows
    for (int j=0; j<N_; ++j) {
        // create rhs
        for (int i=0; i<N_; ++i) {
            uRhs[j] = U(i,j) + Du_*dt_/(2.*dx_*dx_) * ((j+1<N_ ? U(i,j+1) : 0) - 2.*U(i,j) + (j-1>=0 ? U(i,j-1) : 0)) + dt_/2. * (-U(i,j)*V(i,j)*V(i,j) + F_*(1-U(i,j)));
            vRhs[i] = V(i,j) + Dv_*dt_/(2.*dx_*dx_) * ((j+1<N_ ? V(i,j+1) : 0) - 2.*V(i,j) + (j-1>=0 ? V(i,j-1) : 0)) + dt_/2. * (U(i,j)*V(i,j)*V(i,j) - (F_+k_)*V(i,j));
        }
        
        TriDiagMatrixSolver::solve(matU1_, uRhs, uHalf);
        TriDiagMatrixSolver::solve(matV1_, vRhs, vHalf);
    }
    
    // perform the second half-step
    // loop over all columns
    for (int i=0; i<N_; ++i) {
        // create rhs
        for (int j=0; j<N_; ++j) {
            uRhs[i] = UHALF(i,j) + Du_*dt_/(2.*dx_*dx_) * ((i+1<N_ ? UHALF(i+1,j) : 0) - 2.*UHALF(i,j) + (i-1>=0 ? UHALF(i-1,j) : 0)) + dt_/2. * (-UHALF(i,j)*VHALF(i,j)*VHALF(i,j) + F_*(1-UHALF(i,j)));
            vRhs[i] = VHALF(i,j) + Dv_*dt_/(2.*dx_*dx_) * ((i+1<N_ ? VHALF(i+1,j) : 0) - 2.*VHALF(i,j) + (i-1>=0 ? VHALF(i-1,j) : 0)) + dt_/2. * (UHALF(i,j)*VHALF(i,j)*VHALF(i,j) - (F_+k_)*VHALF(i,j));
        }
        
        TriDiagMatrixSolver::solve(matU1_, uRhs, u_);
        TriDiagMatrixSolver::solve(matV1_, vRhs, v_);
    }
}



void Diffusion::initialize_fields()
{
    u_.clear();
    u_.resize(Ntot_,1.);
    v_.clear();
    v_.resize(Ntot_,0.);
    
    for (int y=0; y<N_/10; ++y) {
        for (int x=0; x<N_/10; ++x) {
            U(x,y) = 0.5;
            V(x,y) = 0.25;
        }
    }
}







