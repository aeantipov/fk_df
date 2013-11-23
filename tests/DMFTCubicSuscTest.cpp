#include <numeric>

#include "Solver.h"
#include "DMFT.h"

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

    size_t n_freq = 1024;
    size_t n_b_freq = 15;
    //Log.setDebugging(true);
    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    BMatsubaraGrid gridB(-n_b_freq, n_b_freq+1, beta);
    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*2.0*D*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    RealType diff=1.0;
    CubicDMFTSC<D> SC(Solver, t, KMesh(KPOINTS));
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

    const auto gw = Solver.gw;
    const auto Sigma = Solver.Sigma;

    auto iW = BMatsubara(5,beta);
    GF Vertex4_1(gridF);
    Vertex4_1.fill(typename GF::PointFunctionType([&](FMatsubaraGrid::point w){return Solver.getBVertex4(iW,w);}));
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

    GF glat_cut0(gridF_half);
    GF glat_cutPI(gridF_half);
    GF glat_cut4(gridF_half);
    glat_cut0.fill(typename GF::FunctionType([&](ComplexType w){return glat(w,0.0,0.0);}));
    glat_cutPI.fill(typename GF::FunctionType([&](ComplexType w){return glat(w,PI,PI);}));
    glat_cut4.fill(typename GF::FunctionType([&](ComplexType w){return glat(w,PI/4,PI/4);}));
    glat_cut0.savetxt("GlatK0,0.dat");
    glat_cutPI.savetxt("GlatKpi,pi.dat");
    glat_cut4.savetxt("GlatKpi|4,pi|4.dat");

    auto glat_shift1 = glat.shift(0.0,PI,PI);
    DEBUG((glat.conj()+glat_shift1).sum());

    GF glat_sum(gridF);
    for (auto w:gridF.getPoints()) {
        glat_sum[size_t(w)] = glat[size_t(w)].sum()/RealType(__power<KPOINTS,D>::value); 
    }

    success = is_equal((glat_sum - Solver.gw).sum(),0.0);
    if (!success) {ERROR("Gw out of solver != sum_k glat."); return EXIT_FAILURE; };


    GridObject<RealType,BMatsubaraGrid> chi0_q0_vals(gridB), chi0_q0_dmft_vals(gridB), chi0_qPI_vals(gridB), chi0_qPI_dmft_vals(gridB);
    GridObject<RealType,BMatsubaraGrid> chi_q0_vals(gridB), chi_qPI_vals(gridB), chi_q0_dmft_vals(gridB), chi_qPI_dmft_vals(gridB);
    GF iw_gf(gridF); 
    iw_gf.fill(typename GF::FunctionType([](ComplexType w){return w;}));
    decltype(chi0_q0_dmft_vals)::PointFunctionType chi0_q0_vals_f = [&](BMatsubaraGrid::point in)->RealType { 
            auto g_shift = Solver.gw.shift(in);
            auto sigma_shift = Solver.Sigma.shift(in);
            if (is_equal(ComplexType(in),0.0)) return (-T)*std::real((glat*glat).sum()/RealType(__power<KPOINTS,D>::value));
            return -T*std::real(((Solver.gw - g_shift)/(ComplexType(in)+Solver.Sigma - sigma_shift)).sum());
        };
    decltype(chi0_qPI_vals)::PointFunctionType chi0_qPI_vals_f = [&](BMatsubaraGrid::point in)->RealType { 
            auto g_shift = Solver.gw.shift(in);
            auto sigma_shift = Solver.Sigma.shift(in);
            auto chi_w = (Solver.gw + g_shift)/(2.*iw_gf + ComplexType(in)+ 2.*mu - Solver.Sigma - sigma_shift); 
            auto a = std::real(chi_w.sum());
            if (std::isnan(a)) { // Hit a resonance -2iw = i\Omega
                auto w = gridF.findClosest(-in._val / 2.);
                chi_w[w._index] = std::real((glat[w._index]*glat[w._index]*(-1.)).sum())/RealType(__power<KPOINTS,D>::value);
                a = std::real(chi_w.sum());
                };
            return -T*a;
        };

    chi0_q0_dmft_vals.fill(chi0_q0_vals_f);
    chi0_qPI_dmft_vals.fill(chi0_qPI_vals_f);

    for (auto iW : gridB.getPoints()) {
        INFO("iW = " << iW);
        size_t iwn = size_t(iW);
        GF Vertex4(gridF);
        Vertex4.fill(typename GF::PointFunctionType([&](FMatsubaraGrid::point w){return Solver.getBVertex4(iW,w);}));
        auto gw_bubble = Diagrams::getBubble(gw, iW);
        auto Chi0q0 = Diagrams::getBubble(glat, iW, 0.0, 0.0);    
        auto Chiq0 = Diagrams::getSusc<GF>(Chi0q0, Diagrams::BS(Chi0q0 - gw_bubble, Vertex4 , true));
        auto Chi0qPI = Diagrams::getBubble(glat, iW, PI, PI);    
        auto ChiqPI = Diagrams::getSusc<GF>(Chi0qPI, Diagrams::BS(Chi0qPI - gw_bubble, Vertex4, true)); 
        chi_q0_vals[iwn] = std::real(Chiq0.sum());
        chi0_q0_vals[iwn] = std::real(Chi0q0.sum());
        chi0_qPI_vals[iwn] = std::real(Chi0qPI.sum());
        chi_qPI_vals[iwn] = std::real(ChiqPI.sum());
        auto chiq0_dmft = -T/ComplexType(iW)*(gw-gw.shift(iW)).sum();
        chi_q0_dmft_vals[size_t(iW)] = std::real(chiq0_dmft); 
        if (is_equal(ComplexType(iW),0.0)) { 
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
