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
    real_type mu = 2.2;
    real_type e_d = 0.0;
    real_type beta = 10;
    real_type T=1.0/beta;
    real_type t = 1.0; 
    size_t maxit = 1000;
    real_type mix = 0.5;

    size_t n_freq = 256;
    size_t n_b_freq = 15;

    static const size_t KPOINTS=16;
    static const size_t D=2;

    //Log.setDebugging(true);
    fmatsubara_grid gridF(-n_freq, n_freq, beta);
    bmatsubara_grid gridB(-n_b_freq, n_b_freq+1, beta);

    GF Delta(gridF);
    std::function<complex_type(complex_type)> f1;
    f1 = [t](complex_type w) -> complex_type {return t*2.0*D*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    real_type diff=1.0;
    kmesh kGrid(KPOINTS);
    kmesh_patch qGrid(kGrid);
    std::array<kmesh_patch,2> qGrids( {{ qGrid, qGrid }}) ; 
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

    DFLadderCubic<2> SC_DF(Solver, gridF, kGrid, t);
    SC_DF._n_GD_iter = 0;
    SC_DF();

    fmatsubara_grid grid_half(0,n_freq*2,beta);
    GF Delta_half(grid_half); 
    Delta_half.copy_interpolate(Solver.Delta);
    GF gw_half(grid_half); 
    gw_half.copy_interpolate(Solver.gw);
    GF sigma_half(grid_half); 
    sigma_half.copy_interpolate(Solver.Sigma);

    sigma_half.savetxt("Sigma.dat");
    gw_half.savetxt("Gw.dat");
    Delta_half.savetxt("Delta.dat");

    bool success = false;

    std::vector<complex_type> right_vals = {{ 1.901098273610857259e-02-2.130828360852165537e-01*I,2.161465786590606453e-02 -2.485358842730046314e-01*I }}; 

    for (size_t t=0; t<right_vals.size(); ++t) { 
        success = (is_equal(gw_half[t],right_vals[t],1e-3));
        //if (!success) { ERROR(gw_half[t] << " != " << right_vals[t]); return EXIT_FAILURE; };
        };

    auto glat = SC.getGLat(gridF);
    grid_object<real_type,bmatsubara_grid> chiDMFT0_q0(gridB);
    grid_object<real_type,bmatsubara_grid> chiDMFT0_q0_2(gridB);
    grid_object<real_type,bmatsubara_grid> chiDMFT0_qPI(gridB);
    grid_object<real_type,bmatsubara_grid> chiDMFT0_qPI_2(gridB);

    grid_object<complex_type,bmatsubara_grid> chiDF0_q0(gridB);
    grid_object<complex_type,bmatsubara_grid> chiDF1_q0(gridB);
    grid_object<complex_type,bmatsubara_grid> chiDF0_qPI(gridB);
    grid_object<complex_type,bmatsubara_grid> chiDF1_qPI(gridB);

    GF iw_gf(gridF); 
    iw_gf.fill(typename GF::function_type([](complex_type w){return w;}));
    
    grid_object<real_type,bmatsubara_grid> ::point_function_type chiDMFT0_q0_f = [&](bmatsubara_grid::point in)->real_type { 
            auto g_shift = Solver.gw.shift(in.val_);
            auto sigma_shift = Solver.Sigma.shift(in.val_);
            if (is_equal(complex_type(in),0.0)) return (-T)*std::real((glat*glat).sum()/real_type(std::pow(KPOINTS,D)));
            return -T*std::real(((Solver.gw - g_shift)/(complex_type(in)+Solver.Sigma - sigma_shift)).sum());
        };
    grid_object<real_type,bmatsubara_grid>::point_function_type chiDMFT0_qPI_f = [&](bmatsubara_grid::point in)->real_type { 
            auto g_shift = Solver.gw.shift(in.val_);
            auto sigma_shift = Solver.Sigma.shift(in.val_);
            //if (is_equal(complex_type(in),0.0)) return (-T)*std::real((glat*glat).sum()/real_type(std::pow(KPOINTS,D)));
            return -T*std::real(((Solver.gw + g_shift)/(2*iw_gf + complex_type(in)+ 2*mu - Solver.Sigma - sigma_shift)).sum());
        };
    
    chiDMFT0_q0_2.fill(chiDMFT0_q0_f);
    chiDMFT0_qPI_2.fill(chiDMFT0_qPI_f);

    
    for (auto iW : gridB.points()) {
        INFO("iW = " << iW);
        auto glat_shift = glat.shift(iW,0.0,0.0);
        auto glat_shift_pi = glat.shift(iW,PI,PI);
        chiDMFT0_q0[size_t(iW)] = -T*std::real((glat*glat_shift).sum()/real_type(std::pow(KPOINTS,D)));
        chiDMFT0_qPI[size_t(iW)] = -T*std::real((glat*glat_shift_pi).sum()/real_type(std::pow(KPOINTS,D)));
        auto chiD0 = Diagrams::getBubble(SC_DF.GD,std::make_tuple(iW.val_,kGrid[0],kGrid[0]));
        auto chiDPI = Diagrams::getBubble(SC_DF.GD,std::make_tuple(iW.val_,kGrid[KPOINTS/2],kGrid[KPOINTS/2]));
        chiDF0_q0[size_t(iW)] = chiD0.sum();
        chiDF1_q0[size_t(iW)] = chiDMFT0_q0[size_t(iW)]+T*(Solver.gw*Solver.gw.shift(iW.val_)).sum();
        chiDF0_qPI[size_t(iW)] = chiDPI.sum();
        chiDF1_qPI[size_t(iW)] = chiDMFT0_qPI[size_t(iW)]+T*(Solver.gw*Solver.gw.shift(iW.val_)).sum();
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
