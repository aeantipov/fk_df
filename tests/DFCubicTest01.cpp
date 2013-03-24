#include <numeric>

#include "Solver.h"
#include "SelfConsistency.h"
#include "DF.h"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;
typedef GFWrap GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}


int main()
{
    INFO("Hi! Doing Falicov-Kimball. ");
    RealType U = 4.0;
    RealType mu = U/2.0;
    RealType e_d = 0.0;
    RealType beta = 3.0;
    //RealType T=1.0/beta;
    RealType t = 1.0; 
    RealType mix = 1.0;

    size_t n_freq = 32;
    size_t n_b_freq = 2;

    static const size_t KPOINTS=16;
    static const size_t D=1;

    //Log.setDebugging(true);
    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    BMatsubaraGrid gridB(1, n_b_freq, beta);

    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*2.0*D*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    Solver.w_0 = 0.5;
    Solver.w_1 = 0.5;
    RealType diff=1.0;
    KMesh kGrid(KPOINTS);
    KMeshPatch qGrid(kGrid);
    //std::array<KMeshPatch,2> qGrids( {{ qGrid, qGrid }}) ; 
    CubicDMFTSC<D> SC_DMFT(Solver, t, KMesh(KPOINTS));
    DFLadder<D> SC_DF(Solver, gridF, SC_DMFT._kGrid, t);
    
    bool calc_DMFT = true;
    size_t i_dmft = 0; 
    size_t i_df = 0;
    RealType DFCutoff=1e-7;

    size_t NDMFTRuns = 100;
    size_t NDFRuns = 100;

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
    if (!is_equal(gloc[n_freq],-2.156531283632e-01*I,1e-4)) return EXIT_FAILURE;
    DEBUG(gloc[n_freq+1]);
    if (!is_equal(gloc[n_freq+1],-2.071054102717e-01*I,1e-4)) return EXIT_FAILURE;
       
    std::array<RealType, D> q;
    q.fill(PI);
    auto ChiPiVal = SC_DF.getStaticLatticeSusceptibility<RealType>(q,FMatsubaraGrid(-512,512,beta));
    INFO(ChiPiVal);

    if (!is_equal(ChiPiVal,0.8715,1e-3)) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
