#include <numeric>

#include "Solver.h"
#include "DMFT.h"
#include "DF.h"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;
typedef grid_object<complex_type,fmatsubara_grid> GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, real_type tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}


int main()
{
    INFO("Hi! Doing Falicov-Kimball. ");
    real_type U = 5.0;
    real_type mu = U/2.0;
    real_type e_d = 0.0;
    real_type beta = 4.0;
    real_type T=1.0/beta;
    real_type t = 1.0; 
    size_t maxit = 1000;
    real_type mix = 0.5;

    size_t n_freq = 256;
    size_t n_b_freq = 2;

    static const size_t KPOINTS=16;
    static const size_t D=2;

    //Log.setDebugging(true);
    fmatsubara_grid gridF(-n_freq, n_freq, beta);
    bmatsubara_grid gridB(1, n_b_freq, beta);

    GF Delta(gridF);
    std::function<complex_type(complex_type)> f1;
    f1 = [t](complex_type w) -> complex_type {return t*2.0*D*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    real_type diff=1.0;
    kmesh kGrid(KPOINTS);
    kmesh_patch qGrid(kGrid);
    //std::array<kmesh_patch,2> qGrids( {{ qGrid, qGrid }}) ; 
    CubicDMFTSC<D> SC(Solver, kmesh(KPOINTS), t);
    
    for (int i=0; i<maxit && diff>1e-8; ++i) {
        INFO("Iteration " << i);
        Solver.run();
        auto Delta_new = SC();
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);
        Delta = Delta_new*mix + (1.0-mix)*Solver.Delta;
        Solver.Delta = Delta; 
        }

    DFLadderCubic<D> SC_DF(Solver, gridF, SC._kGrid, t);
    SC_DF._n_GD_iter = 1;
    SC_DF._GDmix = 0.0;
    auto data_out = SC_DF.calculateLatticeData(gridB);
    auto LatticeSusc = std::get<0>(data_out);
  
    /* Checking the charge susc as opposed to a pure DMFT value. */ 
    auto iW = gridB[0];
    GF Gw = Solver.gw;
    GF Sigma = Solver.Sigma;
    GF Gw_shift = Solver.gw.shift(iW);
    GF Sigma_shift = Solver.Sigma.shift(iW);
    auto GLat = SC_DF.GLat;
    auto ChiLat0 = Diagrams::getBubble(GLat,std::make_tuple(iW,kGrid[0],kGrid[0]));

    if (!is_equal(Gw_shift(FMatsubara(gridF.size()-2,beta)), Gw(FMatsubara(gridF.size()-1,beta)))) {
        ERROR("Error in shifted gw");
        return EXIT_FAILURE;
        }; 
    
    auto tmp2 = SC_DF.GLat.shift(iW,kGrid[0],kGrid[0]);
    if (!is_equal(GLat(FMatsubara(gridF.w_max_,beta),0.0,0.0), tmp2(FMatsubara(gridF.w_max_-1,beta),0.0,0.0)) || !is_equal(GLat(FMatsubara(gridF.w_max_+3,beta),0.0,0.0), tmp2(FMatsubara(gridF.w_max_+2,beta),0.0,0.0))) {
        ERROR("Error in shifted glat");
        return EXIT_FAILURE;
        }; 

    GF Chi0DMFT = -T*(Gw-Gw_shift)/(complex_type(iW)+Sigma-Sigma_shift);
    if (!is_equal(ChiLat0.diff(Chi0DMFT),0.0,1e-5)) { 
        ERROR("Bare lattice susc doesn't correspond to the analytical value at q=0");
        return EXIT_FAILURE;
        }
    
    GF IrrVertexDMFT(gridF);
    IrrVertexDMFT = (-1.0/T)*(Sigma-Sigma_shift)/(Gw-Gw_shift);
    
    auto DMFTVertex2 = 1.0/(1.0/IrrVertexDMFT - ChiLat0);
    auto SuscDMFT = ChiLat0 + ChiLat0*DMFTVertex2*ChiLat0;
    auto SuscDMFT2 = T*(-1.0)*(Gw-Gw_shift)/iW;
    INFO("DMFT Susc diff: " << SuscDMFT2.diff(SuscDMFT));
    if (!is_equal(SuscDMFT2.diff(SuscDMFT),0,1e-5)) { 
        ERROR("Full lattice susc doesn't correspond to the analytical value at q=0");
        return EXIT_FAILURE;
        };
    auto ChiVal = SuscDMFT.sum();
    DEBUG(ChiVal);
    DEBUG(LatticeSusc(gridB[0],0,0));
    if (!is_equal(LatticeSusc(gridB[0],0,0), ChiVal)) { 
        ERROR("Full lattice susc doesn't correspond to the analytical value at q=0");
        return EXIT_FAILURE;
    }
   
    return EXIT_SUCCESS;
}
