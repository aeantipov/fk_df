#include <numeric>

#include "Solver.h"
#include "DMFT.h"
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
    RealType U = 5.0;
    RealType mu = U/2.0+1.0;
    RealType e_d = 0.0;
    RealType beta = 1.0;
    RealType T=1.0/beta;
    RealType t = 1.0; 
    size_t maxit = 1000;
    RealType mix = 0.5;

    size_t n_freq = 64;

    static const size_t KPOINTS=16;
    static const size_t D=2;

    //Log.setDebugging(true);
    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    FMatsubaraGrid gridF2(-n_freq-1, n_freq+1, beta);
    FMatsubaraGrid gridF3(-n_freq*2, n_freq*2, beta);
    FMatsubaraGrid gridF4(-n_freq*4, n_freq*4, beta);

    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*2.0*D*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    RealType diff=1.0;
    KMesh kGrid(KPOINTS);
    KMeshPatch qGrid(kGrid);
    //std::array<KMeshPatch,2> qGrids( {{ qGrid, qGrid }}) ; 
    CubicDMFTSC<D> SC(Solver, KMesh(KPOINTS), t);
    
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

    std::array<RealType,2> q_pi ={{PI,PI}};
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
    BubbleqPI2.copyInterpolate(BubbleqPI);
    /*DEBUG(BubbleqPI(FMatsubara(gridF._w_max,beta)));
    DEBUG(BubbleqPI._f(FMatsubara(gridF._w_max,beta)));
    DEBUG(-T/std::pow(FMatsubara(gridF._w_max,beta),2));
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
