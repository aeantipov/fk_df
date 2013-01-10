#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"
#include "Solver.h"
#include "SelfConsistency.h"

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
    RealType mu = 2.0;
    RealType e_d = 0.0;
    RealType beta = 10;
    RealType T = 1.0/beta;
    RealType t = 1.0; 
    size_t maxit = 1000;
    RealType mix = 0.5;
    
    static const size_t KPOINTS=16;
    static const size_t D=2;

    size_t n_freq = 128;
    size_t n_b_freq = 15;
    Log.setDebugging(true);
    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    BMatsubaraGrid gridB(-n_b_freq, n_b_freq+1, beta);
    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*2.0*D*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    RealType diff=1.0;
    CubicDMFTSC<FKImpuritySolver, D, KPOINTS> SC(Solver, t);
    DEBUG(SC._kGrid);

    if (!is_equal((SC._ek + SC._ek.shift(PI,PI)).sum(),0.0)) { ERROR("Glat(k+pi)!=conj(Glat(k))"); return EXIT_FAILURE;};

    for (int i=0; i<maxit && diff>1e-8; ++i) {
        INFO("Iteration " << i);
        Solver.run();
        auto Delta_new = SC();
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);
        Delta = Delta_new*mix + (1.0-mix)*Solver.Delta;
        Solver.Delta = Delta; 
        }

    FMatsubaraGrid gridF_half(0,n_freq*2,beta);
    GF Delta_half(gridF_half); 
    Delta_half = Solver.Delta;
    GF gw_half(gridF_half); 
    gw_half = Solver.gw;
    GF sigma_half(gridF_half); 
    sigma_half = Solver.Sigma;

    sigma_half.savetxt("Sigma.dat");
    gw_half.savetxt("Gw.dat");
    Delta_half.savetxt("Delta.dat");

    bool success = false;

    auto glat = SC.getGLat(gridF);

    GF glat_cut0(gridF_half);
    GF glat_cutPI(gridF_half);
    GF glat_cut4(gridF_half);
    glat_cut0.fill([&](ComplexType w){return glat(w,0.0,0.0);});
    glat_cutPI.fill([&](ComplexType w){return glat(w,PI,PI);});
    glat_cut4.fill([&](ComplexType w){return glat(w,PI/4,PI/4);});
    glat_cut0.savetxt("GlatK0,0.dat");
    glat_cutPI.savetxt("GlatKpi,pi.dat");
    glat_cut4.savetxt("GlatKpi|4,pi|4.dat");

    auto glat_shift1 = glat.shift(0.0,PI,PI);
    DEBUG((glat.conj()+glat_shift1).sum());

    GF glat_sum(gridF);
    for (auto w:gridF.getVals()) {
        glat_sum[size_t(w)] = glat[size_t(w)].sum()/RealType(__power<KPOINTS,D>::value); 
    }

    success = is_equal((glat_sum - Solver.gw).sum(),0.0);
    if (!success) {ERROR("Gw out of solver != sum_k glat."); return EXIT_FAILURE; };


    GridObject<RealType,BMatsubaraGrid> chi0_q0(gridB);
    GridObject<RealType,BMatsubaraGrid> chi0_q0_2(gridB);
    GridObject<RealType,BMatsubaraGrid> chi0_qPI(gridB);
    GridObject<RealType,BMatsubaraGrid> chi0_qPI_2(gridB);
    GF iw_gf(gridF); 
    iw_gf.fill([](ComplexType w){return w;});
    decltype(chi0_q0_2)::PointFunctionType chi0_q0_f = [&](BMatsubaraGrid::point in)->RealType { 
            auto g_shift = Solver.gw.shift(in);
            auto sigma_shift = Solver.Sigma.shift(in);
            if (is_equal(ComplexType(in),0.0)) return (-T)*std::real((glat*glat).sum()/RealType(__power<KPOINTS,D>::value));
            return -T*std::real(((Solver.gw - g_shift)/(ComplexType(in)+Solver.Sigma - sigma_shift)).sum());
        };
    decltype(chi0_qPI)::PointFunctionType chi0_qPI_f = [&](BMatsubaraGrid::point in)->RealType { 
            auto g_shift = Solver.gw.shift(in);
            auto sigma_shift = Solver.Sigma.shift(in);
            //if (is_equal(ComplexType(in),0.0)) return (-T)*std::real((glat*glat).sum()/RealType(__power<KPOINTS,D>::value));
            return -T*std::real(((Solver.gw + g_shift)/(2*iw_gf + ComplexType(in)+ 2*mu - Solver.Sigma - sigma_shift)).sum());
        };

    chi0_q0_2.fill(chi0_q0_f);
    chi0_qPI_2.fill(chi0_qPI_f);
    for (auto iW : gridB.getVals()) {
        INFO("iW = " << iW);
        auto glat_shift = glat.shift(iW,0.0,0.0);
        auto glat_shift_pi = glat.shift(iW,PI,PI);
        chi0_q0[size_t(iW)] = -T*std::real((glat*glat_shift).sum()/RealType(__power<KPOINTS,D>::value));
        chi0_qPI[size_t(iW)] = -T*std::real((glat*glat_shift_pi).sum()/RealType(__power<KPOINTS,D>::value));
        };

    DEBUG(chi0_q0);
    DEBUG(chi0_q0_2);
    DEBUG(chi0_qPI);
    DEBUG(chi0_qPI_2);
    
    success = is_equal((chi0_q0.diff(chi0_q0_2)),0.0,1e-4);
    if (!success) {ERROR("Numerical q=0 susc != Analytical q=0 susc."); return EXIT_FAILURE; };

    success = is_equal((chi0_qPI.diff(chi0_qPI_2)),0.0,1e-3);
    if (!success) {ERROR("Numerical q=PI susc != Analytical q=PI susc."); return EXIT_FAILURE; };
    
    return EXIT_SUCCESS;
}
