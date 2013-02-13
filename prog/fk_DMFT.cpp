#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"
#include "GFWrap.h"
#include "Solver.h"
#include "SelfConsistency.h"

#include "FKOptionParserDMFT.h"

#include <iostream>
#include <ctime>
#include <array>
#include <unordered_map>
#include <csignal>

using namespace FK;

RealType beta;
size_t D=0;
size_t extra_ops;
 
typedef GFWrap GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}

bool INTERRUPT = false;
void sighandler(int signal)
{
    static size_t count = 0;
    count++;
    INFO("Caught INTERRUPT, signal " << signal <<" " << count << " times. ")
    INTERRUPT = true;
    if (count >= 3) { INFO("Force exiting"); exit(signal); }
}

template <class SCType> void calcStats(const SCType& SC, const FMatsubaraGrid& gridF)
{
    const auto &Solver = SC._S;
    auto gw = Solver.gw;
    auto Sigma = Solver.Sigma;
    auto w_0 = Solver.w_0;
    auto w_1 = Solver.w_1;
    auto K0 = Solver.K0;
    auto K1 = Solver.K1;
    auto U = Solver.U;
    //auto mu = Solver.mu;
    RealType beta = Solver.beta;
    RealType T=1.0/beta;

    INFO("Calculating additional statistics.");
    INFO("Static susceptibility");

    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> Vertex4_2(std::forward_as_tuple(gridF,gridF)); 
    decltype(Vertex4_2)::PointFunctionType VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return Solver.getVertex4(0.0, w1,w2);};
    Vertex4_2.fill(VertexF2);
    auto V4 = Vertex4_2.getData().getAsMatrix();

    auto Bubbleq0 = SC.getBubble0(0.0);
    auto BubbleqPI = SC.getBubblePI(0.0); 
    std::vector<std::string> names = {"local", "pi", "zero"};
    std::vector<GF> bubbles = { -T*gw*gw, BubbleqPI, Bubbleq0 };
    std::map<std::string,ComplexType> chi_cc_vals, chi_ff_vals, chi_cf_vals;

    std::map<std::string,ComplexType> chi_cc_vertex_vals, chi_cc_bare_susc_vals;

    auto Lambda = 1.0 + gw*(2.0*Sigma - U); 
    auto Lambda2 = (1. + gw*Sigma)*(1.+gw*(Sigma-U));

    for (size_t i=0; i<bubbles.size(); ++i) { 

        auto bubble = bubbles[i];
        auto dual_bubble = bubbles[i]+T*gw*gw;
        auto dual_bubble_matrix = dual_bubble.getData().getAsDiagonalMatrix();
        auto FullVertex = Diagrams::BS(dual_bubble_matrix, V4, true);
        ComplexType susc = 0.0;
        ComplexType bare_susc = 0.0;

        /** Skeleton expansion. */
        auto nu = gw*(-1.0/gw/gw - T/bubble);
        auto d1 = Lambda*gw*nu + Lambda2;
        auto ugamma_down = 1.0 - (w_0*w_1*U*U*gw*gw*gw*nu/(Lambda2 * (d1))).sum();
        auto ugamma = (w_0*w_1*U*U*gw*gw/(d1)).sum() / (ugamma_down); 

        auto chi_cf = ugamma/U;
        auto chi_ff = w_0*w_1/T/ugamma_down;
        auto chi_cc = -T*((Lambda - ugamma)*gw*gw/d1).sum();
        
        /** Vertex expansion. */
        for (auto w1: gridF.getPoints()) { 
            bare_susc+=bubble(w1);
            for (auto w2: gridF.getPoints()) {
                susc+=bubble(w1)*FullVertex(size_t(w1),size_t(w2))*bubble(w2); 
                }
            };
        susc+=bare_susc;
        chi_cc_vals[names[i]] = chi_cc;
        chi_cc_vertex_vals[names[i]] = susc;
        chi_cc_vals[names[i]] = chi_cc;
        chi_cf_vals[names[i]] = chi_cf;
        chi_ff_vals[names[i]] = chi_ff;
        chi_cc_bare_susc_vals[names[i]] = bare_susc;
        INFO2("Static susc " << names[i] <<" (bs) = " << susc);
        INFO2("Static susc " << names[i] <<" (exact) = " << chi_cc);

        __num_format<ComplexType>(susc).savetxt("StaticChiCC_" + names[i] + ".dat");
        __num_format<ComplexType>(bare_susc).savetxt("StaticChi0CC_" + names[i] + ".dat");
        __num_format<ComplexType>(chi_cc).savetxt("StaticChiCC_" + names[i] + "_skeleton.dat");
        __num_format<ComplexType>(chi_cf).savetxt("StaticChiCF_" + names[i] + "_skeleton.dat");
        __num_format<ComplexType>(chi_ff).savetxt("StaticChiFF_" + names[i] + "_skeleton.dat");
        };

    if (extra_ops>=2) { 
    INFO("Dynamic susceptibility");
    size_t n_b_freq = std::max(std::min(Solver.w_grid._w_max/2,int(2*beta)),10);
    BMatsubaraGrid gridB(-n_b_freq, n_b_freq+1, beta);
    GridObject<RealType,BMatsubaraGrid> chi0_q0_vals(gridB), chi0_qPI_vals(gridB);
    GridObject<RealType,BMatsubaraGrid> chi_q0_vals(gridB), chi_qPI_vals(gridB), chi_q0_dmft_vals(gridB), chi_qPI_dmft_vals(gridB);

    for (auto iW : gridB.getPoints()) {
        if (INTERRUPT) exit(0);
        INFO("iW = " << iW);
        size_t iWn = size_t(iW);
        GF Vertex4(gridF);
        Vertex4.fill(typename GF::PointFunctionType([&](FMatsubaraGrid::point w){return Solver.getBVertex4(iW,w);}));
        auto gw_bubble = Diagrams::getBubble(gw, iW);
        auto Chi0q0 = SC.getBubble0(iW);
        auto Chiq0 = Diagrams::getSusc<GFWrap>(Chi0q0, Diagrams::BS(Chi0q0 - gw_bubble, Vertex4 , true));
        auto Chi0qPI = SC.getBubblePI(iW);
        auto ChiqPI = Diagrams::getSusc<GFWrap>(Chi0qPI, Diagrams::BS(Chi0qPI - gw_bubble, Vertex4, true)); 
        chi_q0_vals[iWn] = std::real(Chiq0.sum());
        chi0_q0_vals[iWn] = std::real(Chi0q0.sum());
        chi0_qPI_vals[iWn] = std::real(Chi0qPI.sum());
        chi_qPI_vals[iWn] = std::real(ChiqPI.sum());
        auto chiq0_dmft = -T/ComplexType(iW)*(gw-gw.shift(iW)).sum();
        chi_q0_dmft_vals[size_t(iW)] = std::real(chiq0_dmft); 

        if (is_equal(ComplexType(iW),0.0)) { 
            INFO("Static val = " << chi_q0_vals[iWn])
            chi_q0_dmft_vals[iWn] = chi_q0_vals[iWn];
            };
        };

    INFO("Chi0[q=0]     = " << chi0_q0_vals);
    INFO("Chi0[q=PI]     = " << chi0_qPI_vals);

    INFO("Chi[q=0]     = " << chi_q0_vals);
    INFO("ChiDMFT[q=0] = " << chi_q0_dmft_vals);
    INFO("Full Chi, q=0 diff = " << chi_q0_vals.diff(chi_q0_dmft_vals));

    chi0_q0_vals.savetxt("DynamicChi0q0.dat");
    chi0_qPI_vals.savetxt("DynamicChi0qPI.dat");
    chi_q0_vals.savetxt("DynamicChiq0.dat");
    chi_q0_dmft_vals.savetxt("DynamicChiq0_DMFT.dat");
    chi_qPI_vals.savetxt("DynamicChiqPI.dat");
    };
}

int main(int argc, char *argv[])
{
  // Catch CTRL-C
  std::signal(SIGABRT, &sighandler);
  std::signal(SIGTERM, &sighandler);
  std::signal(SIGINT , &sighandler);

  FKOptionParserDMFT opt;
   try {
        opt.parse(&argv[1], argc-1); // Skip argv[0].
        INFO("Hi! Doing Falicov-Kimball. ");
        std::cout << "FK. Parameters " << std::endl;
        std::cout << "beta                 : " << opt.beta << std::endl;
        std::cout << "U                    : " << opt.U    << std::endl;
        std::cout << "t                    : " << opt.t    << std::endl;
        std::cout << "mu                   : " << opt.mu   << std::endl;
        std::cout << "e_d                  : " << opt.e_d << std::endl;
        std::cout << "Selfconsistency      : " << opt.sc_type << std::endl;
        std::cout << "Number Of Matsubaras : " << opt.n_freq << std::endl;
        std::cout << "Max number of iterations : " << opt.n_iter << std::endl;
    } catch (const optparse::unrecognized_option& e) {
        std::cout << "unrecognized option: " << e.what() << std::endl;
        return 1;
    } catch (const optparse::invalid_value& e) {
        std::cout << "invalid value: " << e.what() << std::endl;
        return 1;
    }
    
    RealType U = opt.U;
    RealType mu = opt.mu;
    RealType e_d = opt.e_d;
    beta = opt.beta;
    RealType t = opt.t; 
    size_t n_freq = opt.n_freq;
    size_t maxit = opt.n_iter;
    RealType mix = opt.mix;
    auto sc_switch = opt.sc_index;
    extra_ops = opt.extra_ops;
    size_t kpoints = opt.kpts;

    KMesh kgrid(kpoints);
    
    Log.setDebugging(true);

    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1, f2;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};

    try { GF Delta2("Delta_full.dat", beta); Delta = Delta2;} 
    catch (std::exception &e) { Delta.fill(f1); };

    Delta.savetxt("Delta_0.dat");
    
    FKImpuritySolver Solver(U,mu,e_d,Delta);

    std::unique_ptr<SelfConsistency> SC_ptr;

    typedef FKOptionParserDMFT::SC enumSC;
    switch (sc_switch) {
        case enumSC::Bethe:       SC_ptr.reset(new BetheSC(Solver, t)); break;
        case enumSC::DMFTCubic1d: SC_ptr.reset(new CubicDMFTSC<1>(Solver, t, kgrid)); D=1; break;
        case enumSC::DMFTCubic2d: SC_ptr.reset(new CubicDMFTSC<2>(Solver, t, kgrid)); D=2; break;
        case enumSC::DMFTCubic3d: SC_ptr.reset(new CubicDMFTSC<3>(Solver, t, kgrid)); D=3; break;
        case enumSC::DMFTCubic4d: SC_ptr.reset(new CubicDMFTSC<4>(Solver, t, kgrid)); D=4; break;
        case enumSC::DMFTCubicInfd: SC_ptr.reset(new CubicInfDMFTSC(Solver,t,RealGrid(-6.0*t,6.0*t,1024))); break;
        default:                  ERROR("No self-consistency provided. Exiting..."); exit(1); 
    };
    auto &SC = *SC_ptr;

    RealType diff=1.0;
    for (int i=0; i<maxit && diff>1e-8 &&!INTERRUPT; ++i) {
        INFO("Iteration " << i <<". Mixing = " << mix);
        if (diff/mix>1e-3) Solver.run(true);
        else Solver.run(false);
        Delta = SC();
        auto Delta_new = Delta*mix+(1.0-mix)*Solver.Delta;
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);
        Solver.Delta = Delta_new;
        }
   
    FMatsubaraGrid gridF_half(0, std::max(n_freq*3,size_t(100)), beta);
    GF Delta_half(gridF_half); Delta_half = Delta;
    GF gw_half(gridF_half); gw_half = Solver.gw;
    GF sigma_half(gridF_half); sigma_half = Solver.Sigma;
    sigma_half.savetxt("Sigma.dat");
    gw_half.savetxt("Gw.dat");
    Delta_half.savetxt("Delta.dat");
    int __msize = std::max(n_freq*5,size_t(1024));
    FMatsubaraGrid gridF_large(-__msize, __msize, beta);
    GF Delta_large(gridF_large); Delta_large = Delta;
    Delta_large.savetxt("Delta_full.dat");

    // Everything important is finished by here
    if (extra_ops) {
        switch (sc_switch) {
            case enumSC::DMFTCubic1d: 
                calcStats(*(static_cast<CubicDMFTSC<1>*> (SC_ptr.get())), gridF); 
                break;
            case enumSC::DMFTCubic2d: 
                calcStats(*(static_cast<CubicDMFTSC<2>*> (SC_ptr.get())), gridF); 
                break;
            case enumSC::DMFTCubic3d: 
                calcStats(*(static_cast<CubicDMFTSC<3>*> (SC_ptr.get())), gridF); 
                break;
            case enumSC::DMFTCubic4d: 
                calcStats(*(static_cast<CubicDMFTSC<4>*> (SC_ptr.get())), gridF); 
                break;
            case enumSC::DMFTCubicInfd: 
                calcStats(*(static_cast<CubicInfDMFTSC*> (SC_ptr.get())), gridF); 
                break;
 
            default: break;
            }; 
        };
    SC_ptr.release(); 
}

