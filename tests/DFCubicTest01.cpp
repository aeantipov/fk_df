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
    real_type U = 4.0;
    real_type mu = U/2.0;
    real_type e_d = 0.0;
    real_type beta = 3.0;
    //real_type T=1.0/beta;
    real_type t = 1.0; 
    real_type mix = 1.0;

    size_t n_freq = 32;
    size_t n_b_freq = 2;

    static const size_t KPOINTS=16;
    static const size_t D=1;

    //Log.setDebugging(true);
    fmatsubara_grid gridF(-n_freq, n_freq, beta);
    bmatsubara_grid gridB(1, n_b_freq, beta);

    GF Delta(gridF);
    std::function<complex_type(complex_type)> f1;
    f1 = [t](complex_type w) -> complex_type {return t*2.0*D*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    Solver.w_0 = 0.5;
    Solver.w_1 = 0.5;
    real_type diff=1.0;
    kmesh kGrid(KPOINTS);
    kmesh_patch qGrid(kGrid);
    //std::array<kmesh_patch,2> qGrids( {{ qGrid, qGrid }}) ; 
    CubicDMFTSC<D> SC_DMFT(Solver, kmesh(KPOINTS), t);
    DFLadderCubic<D> SC_DF(Solver, gridF, SC_DMFT._kGrid, t);
    
    bool calc_DMFT = true;
    size_t i_dmft = 0; 
    size_t i_df = 0;
    real_type DFCutoff=1e-7;

    size_t NDMFTRuns = 100;
    size_t NDFRuns = 1;

    for (; i_dmft<=NDMFTRuns-calc_DMFT && i_df<=NDFRuns && diff>1e-8+(1-calc_DMFT)*(DFCutoff-1e-8); (calc_DMFT)?i_dmft++:i_df++) {
        INFO("Iteration " << i_dmft+i_df <<". Mixing = " << mix);
        Solver.run(false);
        if (calc_DMFT) {  
            Delta = SC_DMFT();
            }
        else { 
            Delta = SC_DF();
             }
        auto Delta_new = Delta*mix+(1.0-mix)*Solver.Delta;
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);

        Solver.Delta = Delta_new; 
        if (diff<=1e-8 && calc_DMFT) { 
            diff = 1.0; calc_DMFT = false; 
            }; // now continue with DF 
        }

    auto gloc = SC_DF.GLatLoc;
    DEBUG(gloc[n_freq]);
    if (!is_equal(gloc[n_freq],-2.285624025387e-01*I,1e-4)) return EXIT_FAILURE;
    DEBUG(gloc[n_freq+1]);
    if (!is_equal(gloc[n_freq+1],-2.095296954079e-01*I,1e-4)) return EXIT_FAILURE;
       
    std::array<real_type, D> q;
    q.fill(PI);

    auto ChiPiVal = SC_DF.getStaticLatticeSusceptibility<real_type>(q,fmatsubara_grid(-32,32,beta));
    INFO(ChiPiVal);
    ChiPiVal = SC_DF.getStaticLatticeSusceptibility<real_type>(q,fmatsubara_grid(-33,33,beta));
    INFO(ChiPiVal);
    ChiPiVal = SC_DF.getStaticLatticeSusceptibility<real_type>(q,fmatsubara_grid(-1024,1024,beta));
    INFO(ChiPiVal);
    if (!is_equal(ChiPiVal,0.923208,1e-2)) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
