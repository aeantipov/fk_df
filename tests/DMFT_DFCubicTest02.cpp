#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"
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
    RealType mu = 2.2;
    RealType e_d = 0.0;
    RealType beta = 10;
    RealType T=1.0/beta;
    RealType t = 1.0; 
    size_t maxit = 1000;
    RealType mix = 0.5;

    size_t n_freq = 256;
    size_t n_b_freq = 15;

    static const size_t KPOINTS=16;
    static const size_t D=2;

    Log.setDebugging(true);
    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    BMatsubaraGrid gridB(-n_b_freq, n_b_freq+1, beta);

    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*2.0*D*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    RealType diff=1.0;
    KMesh kGrid(KPOINTS);
    KMeshPatch qGrid(kGrid);
    std::array<KMeshPatch,2> qGrids( {{ qGrid, qGrid }}) ; 
    CubicDMFTSC<FKImpuritySolver, D> SC(Solver, t, KMesh(KPOINTS));
    
    for (int i=0; i<maxit && diff>1e-8; ++i) {
        INFO("Iteration " << i);
        Solver.run();
        auto Delta_new = SC();
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);
        Delta = Delta_new*mix + (1.0-mix)*Solver.Delta;
        Solver.Delta = Delta; 
        }

    DFLadder<FKImpuritySolver, 2> SC_DF(Solver, gridF, kGrid, gridB, t);
    SC_DF._n_GD_iter = 0;
    SC_DF();

    FMatsubaraGrid gridF_half(0,n_freq*2,beta);
    GF Delta_half(gridF_half); 
    Delta_half = Delta;
    GF gw_half(gridF_half); 
    gw_half = Solver.gw;
    GF sigma_half(gridF_half); 
    sigma_half = Solver.Sigma;

    sigma_half.savetxt("Sigma.dat");
    gw_half.savetxt("Gw.dat");
    Delta_half.savetxt("Delta.dat");

    bool success = false;

    std::vector<ComplexType> right_vals = {{ 1.901098273610857259e-02-2.130828360852165537e-01*I,2.161465786590606453e-02 -2.485358842730046314e-01*I }}; 

    for (size_t t=0; t<right_vals.size(); ++t) { 
        success = (is_equal(gw_half[t],right_vals[t],1e-3));
        //if (!success) { ERROR(gw_half[t] << " != " << right_vals[t]); return EXIT_FAILURE; };
        };

    auto glat = SC.getGLat(gridF);
    GridObject<RealType,BMatsubaraGrid> chiDMFT0_q0(gridB);
    GridObject<RealType,BMatsubaraGrid> chiDMFT0_q0_2(gridB);
    GridObject<RealType,BMatsubaraGrid> chiDMFT0_qPI(gridB);
    GridObject<RealType,BMatsubaraGrid> chiDMFT0_qPI_2(gridB);

    GridObject<ComplexType,BMatsubaraGrid> chiDF0_q0(gridB);
    GridObject<ComplexType,BMatsubaraGrid> chiDF1_q0(gridB);
    GridObject<ComplexType,BMatsubaraGrid> chiDF0_qPI(gridB);
    GridObject<ComplexType,BMatsubaraGrid> chiDF1_qPI(gridB);

    GF iw_gf(gridF); 
    iw_gf.fill([](ComplexType w){return w;});
    
    decltype(chiDMFT0_q0_2)::PointFunctionType chiDMFT0_q0_f = [&](BMatsubaraGrid::point in)->RealType { 
            auto g_shift = Solver.gw.shift(in);
            auto sigma_shift = Solver.Sigma.shift(in);
            if (is_equal(ComplexType(in),0.0)) return (-T)*std::real((glat*glat).sum()/RealType(__power<KPOINTS,D>::value));
            return -T*std::real(((Solver.gw - g_shift)/(ComplexType(in)+Solver.Sigma - sigma_shift)).sum());
        };
    decltype(chiDMFT0_qPI)::PointFunctionType chiDMFT0_qPI_f = [&](BMatsubaraGrid::point in)->RealType { 
            auto g_shift = Solver.gw.shift(in);
            auto sigma_shift = Solver.Sigma.shift(in);
            //if (is_equal(ComplexType(in),0.0)) return (-T)*std::real((glat*glat).sum()/RealType(__power<KPOINTS,D>::value));
            return -T*std::real(((Solver.gw + g_shift)/(2*iw_gf + ComplexType(in)+ 2*mu - Solver.Sigma - sigma_shift)).sum());
        };
    
    chiDMFT0_q0_2.fill(chiDMFT0_q0_f);
    chiDMFT0_qPI_2.fill(chiDMFT0_qPI_f);

    
    for (auto iW : gridB.getPoints()) {
        INFO("iW = " << iW);
        auto glat_shift = glat.shift(iW,0.0,0.0);
        auto glat_shift_pi = glat.shift(iW,PI,PI);
        chiDMFT0_q0[size_t(iW)] = -T*std::real((glat*glat_shift).sum()/RealType(__power<KPOINTS,D>::value));
        chiDMFT0_qPI[size_t(iW)] = -T*std::real((glat*glat_shift_pi).sum()/RealType(__power<KPOINTS,D>::value));
        auto chiD0 = Diagrams::getBubble(SC_DF.GD,std::make_tuple(iW,kGrid[0],kGrid[0]));
        auto chiDPI = Diagrams::getBubble(SC_DF.GD,std::make_tuple(iW,kGrid[KPOINTS/2],kGrid[KPOINTS/2]));
        chiDF0_q0[size_t(iW)] = chiD0.sum();
        chiDF1_q0[size_t(iW)] = chiDMFT0_q0[size_t(iW)]+T*(Solver.gw*Solver.gw.shift(iW)).sum(); 
        chiDF0_qPI[size_t(iW)] = chiDPI.sum();
        chiDF1_qPI[size_t(iW)] = chiDMFT0_qPI[size_t(iW)]+T*(Solver.gw*Solver.gw.shift(iW)).sum(); 
        };

    INFO("Checking DMFT susceptibilities.");
    DEBUG(chiDMFT0_q0);
    DEBUG(chiDMFT0_q0_2);
    if (!is_equal(chiDMFT0_q0.diff(chiDMFT0_q0_2),0,1e-5)) { ERROR("Susceptibilities don't match with analytic value"); return EXIT_FAILURE; }
    DEBUG(chiDMFT0_qPI);
    DEBUG(chiDMFT0_qPI_2);
    if (!is_equal(chiDMFT0_qPI.diff(chiDMFT0_qPI_2),0,1e-5)) { ERROR("Susceptibilities don't match with analytic value"); return EXIT_FAILURE; }
    
    DEBUG("----------------");
    DEBUG(chiDF0_q0);
    DEBUG(chiDF1_q0);
    DEBUG(chiDF1_q0.diff(chiDF0_q0));
    if (!is_equal(chiDF0_q0.diff(chiDF1_q0),0,1e-5)) { ERROR("Susceptibilities don't match with DF value"); return EXIT_FAILURE; }
    DEBUG(chiDF0_qPI);
    DEBUG(chiDF1_qPI);
    DEBUG(chiDF1_qPI-chiDF0_qPI);
    
    DEBUG(chiDF1_qPI.diff(chiDF0_qPI));
    if (!is_equal(chiDF0_qPI.diff(chiDF1_qPI),0,1e-5)) { ERROR("Susceptibilities don't match with DF value"); return EXIT_FAILURE; }
    
    return EXIT_SUCCESS;
}
