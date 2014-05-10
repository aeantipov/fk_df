#include <numeric>

#include "Solver.h"
#include "DMFT.h"

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
    real_type mu = 2.0;
    real_type e_d = 0.0;
    real_type beta = 10;
    real_type T = 1.0/beta;
    real_type t = 1.0; 
    size_t maxit = 1000;
    real_type mix = 0.5;
    
    static const size_t KPOINTS=16;
    static const size_t D=2;

    size_t n_freq = 1024;
    size_t n_b_freq = 15;
    //Log.setDebugging(true);
    fmatsubara_grid gridF(-n_freq, n_freq, beta);
    bmatsubara_grid gridB(-n_b_freq, n_b_freq+1, beta);
    GF Delta(gridF);
    std::function<complex_type(complex_type)> f1;
    f1 = [t](complex_type w) -> complex_type {return t*2.0*D*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    real_type diff=1.0;
    CubicDMFTSC<D> SC(Solver, kmesh(KPOINTS), t);
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

    const auto gw = Solver.gw;
    const auto Sigma = Solver.Sigma;

    auto iW = BMatsubara(5,beta);
    GF Vertex4_1(gridF);
    Vertex4_1.fill(typename GF::point_function_type([&](fmatsubara_grid::point w){return Solver.getBVertex4(iW,w);}));
    auto irrDMFTV4 = Diagrams::BS(Diagrams::getBubble(gw,iW), Vertex4_1, false);
    auto gw_shift = Solver.gw.shift(iW);
    auto Sigma_shift = Solver.Sigma.shift(iW);
    auto irrDMFT2 = (Sigma - Sigma_shift)/(gw-gw_shift)*(-1.0)/T;
    auto irrDMFT3 = (-1.0)/T*Solver.w_0*Solver.w_1*std::pow(Solver.U,2)/( (1.+gw*Sigma)*(1.+gw*(Sigma-U))*(1.+gw_shift*(Sigma_shift))*(1.+gw_shift*(Sigma_shift - U))+ Solver.w_0*Solver.w_1*U*U*gw*gw_shift);

    if (!is_equal(irrDMFTV4.diff(irrDMFT3),0,1e-4)) { ERROR("DMFT vertex is incorrect. "); return EXIT_FAILURE; };

    //DEBUG(irrDMFTV4);
    //DEBUG(irrDMFT3);
    //DEBUG(irrDMFT3 - irrDMFTV4);
    DEBUG(irrDMFTV4.diff(irrDMFT3));
    auto glat = SC.getGLat(gridF);

    GF glat_cut0(grid_half);
    GF glat_cutPI(grid_half);
    GF glat_cut4(grid_half);
    glat_cut0.fill(typename GF::function_type([&](complex_type w){return glat(w,0.0,0.0);}));
    glat_cutPI.fill(typename GF::function_type([&](complex_type w){return glat(w,PI,PI);}));
    glat_cut4.fill(typename GF::function_type([&](complex_type w){return glat(w,PI/4,PI/4);}));
    glat_cut0.savetxt("GlatK0,0.dat");
    glat_cutPI.savetxt("GlatKpi,pi.dat");
    glat_cut4.savetxt("GlatKpi|4,pi|4.dat");

    auto glat_shift1 = glat.shift(0.0,PI,PI);
    DEBUG((glat.conj()+glat_shift1).sum());

    GF glat_sum(gridF);
    for (auto w:gridF.points()) {
        glat_sum[size_t(w)] = glat[size_t(w)].sum()/real_type(std::pow(KPOINTS,D));
    }

    success = is_equal((glat_sum - Solver.gw).sum(),0.0);
    if (!success) {ERROR("Gw out of solver != sum_k glat."); return EXIT_FAILURE; };


    grid_object<real_type,bmatsubara_grid> chi0_q0_vals(gridB), chi0_q0_dmft_vals(gridB), chi0_qPI_vals(gridB), chi0_qPI_dmft_vals(gridB);
    grid_object<real_type,bmatsubara_grid> chi_q0_vals(gridB), chi_qPI_vals(gridB), chi_q0_dmft_vals(gridB), chi_qPI_dmft_vals(gridB);
    GF iw_gf(gridF); 
    iw_gf.fill(typename GF::function_type([](complex_type w){return w;}));
    grid_object<real_type,bmatsubara_grid>::point_function_type chi0_q0_vals_f = [&](bmatsubara_grid::point in)->real_type { 
            auto g_shift = Solver.gw.shift(in.val_);
            auto sigma_shift = Solver.Sigma.shift(in.val_);
            if (is_equal(complex_type(in),0.0)) return (-T)*std::real((glat*glat).sum()/real_type(std::pow(KPOINTS,D)));
            return -T*std::real(((Solver.gw - g_shift)/(complex_type(in)+Solver.Sigma - sigma_shift)).sum());
        };
    grid_object<real_type,bmatsubara_grid>::point_function_type chi0_qPI_vals_f = [&](bmatsubara_grid::point in)->real_type { 
            auto g_shift = Solver.gw.shift(in.val_);
            auto sigma_shift = Solver.Sigma.shift(in.val_);
            auto chi_w = (Solver.gw + g_shift)/(2.*iw_gf + complex_type(in)+ 2.*mu - Solver.Sigma - sigma_shift); 
            auto w = gridF.find_nearest(-in.val_ / 2.); 
            if (std::abs(in.val_ + w.val_*2.)<1e-8) 
                chi_w[w.index_] = std::real((glat[w.index_]*glat[w.index_]*(-1.)).sum())/real_type(std::pow(KPOINTS,D)); // Patch a resonance -2iw = i\Omega
            auto a = std::real(chi_w.sum());
            return -T*a;
        };

    chi0_q0_dmft_vals.fill(chi0_q0_vals_f);
    chi0_qPI_dmft_vals.fill(chi0_qPI_vals_f);

    for (auto iW : gridB.points()) {
        INFO("iW = " << iW);
        size_t iwn = size_t(iW);
        GF Vertex4(gridF);
        Vertex4.fill(typename GF::point_function_type([&](fmatsubara_grid::point w){return Solver.getBVertex4(iW,w);}));
        auto gw_bubble = Diagrams::getBubble(gw, iW.val_);
        auto Chi0q0 = Diagrams::getBubble(glat, iW.val_, 0.0, 0.0);
        auto Chiq0 = Diagrams::getSusc<GF>(Chi0q0, Diagrams::BS(Chi0q0 - gw_bubble, Vertex4 , true));
        auto Chi0qPI = Diagrams::getBubble(glat, iW.val_, PI, PI);
        auto ChiqPI = Diagrams::getSusc<GF>(Chi0qPI, Diagrams::BS(Chi0qPI - gw_bubble, Vertex4, true)); 
        chi_q0_vals[iwn] = std::real(Chiq0.sum());
        chi0_q0_vals[iwn] = std::real(Chi0q0.sum());
        chi0_qPI_vals[iwn] = std::real(Chi0qPI.sum());
        chi_qPI_vals[iwn] = std::real(ChiqPI.sum());
        auto chiq0_dmft = -T/complex_type(iW.val_)*(gw-gw.shift(iW.val_)).sum();
        chi_q0_dmft_vals[size_t(iW)] = std::real(chiq0_dmft); 
        if (is_equal(complex_type(iW),0.0)) { 
            INFO("Static val = " << chi_q0_vals[size_t(iW)])
            chi_q0_dmft_vals[size_t(iW)] = chi_q0_vals[size_t(iW)]; 
            };
        };

    INFO("Chi0[q=0]     = " << chi0_q0_vals);
    INFO("Chi0DMFT[q=0] = " << chi0_q0_dmft_vals);
    INFO("Chi0, q=0 diff = " << chi0_q0_vals.diff(chi0_q0_dmft_vals));
    
    INFO("Chi0[q=PI]     = " << chi0_qPI_vals);
    INFO("Chi0DMFT[q=PI] = " << chi0_qPI_dmft_vals);
    INFO("Chi0, q=pi diff = " << chi0_qPI_vals.diff(chi0_qPI_dmft_vals));


    INFO("Chi[q=0]     = " << chi_q0_vals);
    INFO("Chi0DMFT[q=0] = " << chi_q0_dmft_vals);
    INFO("Full Chi, q=0 diff = " << chi_q0_vals.diff(chi_q0_dmft_vals));
    
    success = is_equal((chi0_q0_vals.diff(chi0_q0_dmft_vals)),0.0,1e-4);
    if (!success) {ERROR("Numerical q=0 susc != Analytical q=0 susc."); return EXIT_FAILURE; };

    success = is_equal((chi0_qPI_vals.diff(chi0_qPI_dmft_vals)),0.0,3e-4);
    if (!success) {ERROR("Numerical q=PI susc != Analytical q=PI susc."); return EXIT_FAILURE; };

    success = is_equal((chi_q0_vals.diff(chi_q0_dmft_vals)),0.0,1e-4);
    if (!success) {ERROR("Numerical q=0 full susc != Analytical q=0 DMFT susc."); return EXIT_FAILURE; };

    return EXIT_SUCCESS;
}
