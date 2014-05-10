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
    real_type mu = U/2.0+1.0;
    real_type e_d = 0.0;
    real_type beta = 1.0;
    real_type T=1.0/beta;
    real_type t = 1.0; 
    size_t maxit = 1000;
    real_type mix = 0.5;

    size_t n_freq = 64;

    static const size_t KPOINTS=16;
    static const size_t D=2;

    //Log.setDebugging(true);
    fmatsubara_grid gridF(-n_freq, n_freq, beta);
    fmatsubara_grid gridF2(-n_freq-1, n_freq+1, beta);
    fmatsubara_grid gridF3(-n_freq*2, n_freq*2, beta);
    fmatsubara_grid gridF4(-n_freq*4, n_freq*4, beta);

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
    SC_DF._n_GD_iter = 0;
    SC_DF._GDmix = 0.0;
    SC_DF();

    std::array<real_type,2> q_pi ={{PI,PI}};
    INFO("CC susc with a basic grid");
    auto t1 = SC_DF.getStaticLatticeSusceptibility(q_pi,gridF);
    INFO("CC susc with a +1 point grid");
    auto t2 = SC_DF.getStaticLatticeSusceptibility(q_pi,gridF2);
    INFO("CC susc with a *2 points grid");
    auto t3 = SC_DF.getStaticLatticeSusceptibility(q_pi,gridF3);
    //INFO("CC susc with a *4 points grid");
    //auto t4 = SC_DF.getStaticLatticeSusceptibility(q_pi,gridF4);

    //std::vector<GF> bubbles = { -T*gw*gw, BubbleqPI, Bubbleq0 };
    auto BubbleqPI = SC.getBubblePI(0.0); 
    GF BubbleqPI2(gridF2);
    BubbleqPI2.copy_interpolate(BubbleqPI);
    /*DEBUG(BubbleqPI(FMatsubara(gridF.w_max_,beta)));
    DEBUG(BubbleqPI._f(FMatsubara(gridF.w_max_,beta)));
    DEBUG(-T/std::pow(FMatsubara(gridF.w_max_,beta),2));
    DEBUG(BubbleqPI);
    DEBUG(BubbleqPI2);
    */
    auto t1_dmft = getStaticLatticeDMFTSusceptibility(Solver,BubbleqPI,gridF);
    auto t2_dmft = getStaticLatticeDMFTSusceptibility(Solver,BubbleqPI,gridF2);
    auto t3_dmft = getStaticLatticeDMFTSusceptibility(Solver,BubbleqPI,gridF3);
 
    auto t1_dmft_skel = getStaticLatticeDMFTSkeletonSusceptibility(Solver,BubbleqPI,gridF)[0];
    auto t2_dmft_skel = getStaticLatticeDMFTSkeletonSusceptibility(Solver,BubbleqPI,gridF2)[0];
    auto t3_dmft_skel = getStaticLatticeDMFTSkeletonSusceptibility(Solver,BubbleqPI,gridF3)[0];
    if (!is_equal(t1,t1_dmft) || !is_equal(t1,t1_dmft_skel)) return EXIT_FAILURE; 
    INFO("PASSED");
    if (!is_equal(t2,t2_dmft) || !is_equal(t2,t2_dmft_skel)) return EXIT_FAILURE; 
    INFO("PASSED");
    if (!is_equal(t3,t3_dmft) || !is_equal(t3,t3_dmft_skel)) return EXIT_FAILURE; 
    INFO("PASSED");
    return EXIT_SUCCESS;
}
