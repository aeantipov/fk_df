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
#ifdef K8
    static const size_t KPOINTS = 8;
#elif K16
    static const size_t KPOINTS = 16;
#elif K32
    static const size_t KPOINTS = 32;
#elif K64
    static const size_t KPOINTS = 64;
#else 
    static const size_t KPOINTS = 16;
#endif

RealType beta;
size_t D=0;
size_t extra_ops;
 
typedef GFWrap GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}

bool interrupt = false;
void sighandler(int signal)
{
    INFO("Caught interrupt, signal " << signal <<". Exiting...")
    interrupt = true;
}

template <class SCType> void calcStats(const SCType& SC, const FMatsubaraGrid& gridF)
{
    const auto &Solver = SC._S;
    RealType beta = Solver.beta;
    RealType T=1.0/beta;
    INFO("Calculating additional statistics.");
    INFO("Static susceptibility");

    auto Bubbleq0 = SC.getBubble0(0.0);
    auto BubbleqPI = SC.getBubblePI(0.0); 
    auto StaticV4 = SC.getStaticLatticeDMFTVertex4();

    auto FullVertexqPI = StaticV4.getData().getAsMatrix();
    auto FullVertexq0 = StaticV4.getData().getAsMatrix();
    auto Chiq0 = Bubbleq0.getData().getAsDiagonalMatrix();
    auto ChiqPI = BubbleqPI.getData().getAsDiagonalMatrix();

    FullVertexq0 = Diagrams::BS(Chiq0, FullVertexq0, true);
    FullVertexqPI = Diagrams::BS(ChiqPI, FullVertexqPI, true);

    GF susc0(gridF), suscPI(gridF);
    for (auto w1: gridF.getVals()) { 
        susc0[size_t(w1)]+=Bubbleq0(w1);
        suscPI[size_t(w1)]+=BubbleqPI(w1);
        for (auto w2: gridF.getVals()) {
            susc0[size_t(w1)]+=Bubbleq0(w1)*FullVertexq0(size_t(w1),size_t(w2))*Bubbleq0(w2); 
            suscPI[size_t(w1)]+=BubbleqPI(w1)*FullVertexqPI(size_t(w1),size_t(w2))*BubbleqPI(w2); 
            }
        };
    INFO("T = " << T);
    RealType susc0_val = std::real(susc0.sum());
    RealType suscPI_val = std::real(suscPI.sum());
    INFO("Static q=0 susc = " << susc0_val);
    INFO("Static q=PI susc = " << suscPI_val);
    __num_format<RealType>(susc0_val).savetxt("StaticChiq0.dat");
    __num_format<RealType>(suscPI_val).savetxt("StaticChiqPI.dat");
    __num_format<RealType>(std::real(Bubbleq0.sum())).savetxt("StaticChi0q0.dat");
    __num_format<RealType>(std::real(BubbleqPI.sum())).savetxt("StaticChi0qPI.dat");

    if (extra_ops>=2) { 
    INFO("Dynamic susceptibility");
    size_t n_b_freq = std::min(Solver.w_grid._w_max/2,int(2*beta));
    BMatsubaraGrid gridB(-n_b_freq, n_b_freq+1, beta);
    GridObject<RealType,BMatsubaraGrid> chi0_q0_vals(gridB), chi0_qPI_vals(gridB);
    GridObject<RealType,BMatsubaraGrid> chi_q0_vals(gridB), chi_qPI_vals(gridB), chi_q0_dmft_vals(gridB), chi_qPI_dmft_vals(gridB);
    auto gw = SC._S.gw;
    auto Sigma = SC._S.Sigma;

    for (auto iW : gridB.getVals()) {
        if (interrupt) exit(0);
        INFO("iW = " << iW);
        size_t iwn = size_t(iW);
        auto Chi0q0 = SC.getBubble0(iW); 
        auto Chiq0 = Diagrams::getSusc<GFWrap>(Chi0q0, Diagrams::BS(Chi0q0,SC.getLatticeDMFTVertex4(iW), true));
        auto Chi0qPI = SC.getBubblePI(iW);
        auto ChiqPI = Diagrams::getSusc<GFWrap>(Chi0qPI, Diagrams::BS(Chi0qPI,SC.getLatticeDMFTVertex4(iW), true)); 
        chi_q0_vals[iwn] = std::real(Chiq0.sum());
        chi0_q0_vals[iwn] = std::real(Chi0q0.sum());
        chi0_qPI_vals[iwn] = std::real(Chi0qPI.sum());
        chi_qPI_vals[iwn] = std::real(ChiqPI.sum());
        auto chiq0_dmft = -T/ComplexType(iW)*(gw-gw.shift(iW)).sum();
        chi_q0_dmft_vals[size_t(iW)] = std::real(chiq0_dmft); 
        if (is_equal(ComplexType(iW),0.0)) { 
            INFO("Static val = " << chi_q0_vals[size_t(iW)])
        //    chi_q0_vals[size_t(iW)] += susc0_val;
        //    chi_qPI_vals[size_t(iW)] += suscPI_val;
            chi_q0_dmft_vals[size_t(iW)] = chi_q0_vals[size_t(iW)]; 
            //chi_qPI_dmft_vals[size_t(iW)] += suscPI_val;
            };
        };

    INFO("Chi0[q=0]     = " << chi0_q0_vals);
    INFO("Chi0[q=PI]     = " << chi0_qPI_vals);

    INFO("Chi[q=0]     = " << chi_q0_vals);
    INFO("Chi0DMFT[q=0] = " << chi_q0_dmft_vals);
    INFO("Full Chi, q=0 diff = " << chi_q0_vals.diff(chi_q0_dmft_vals));

    chi0_q0_vals.savetxt("DynamicChi0q0.dat");
    chi0_qPI_vals.savetxt("DynamicChi0qPI.dat");
    chi_q0_vals.savetxt("DynamicChiq0.dat");
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
    
    Log.setDebugging(true);

    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    FMatsubaraGrid gridF_half(0, n_freq*2, beta);
    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1, f2;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};

    try { Delta.loadtxt("Delta_full.dat"); } 
    catch (std::exception &e) { Delta.fill(f1); };

    Delta.savetxt("Delta_0.dat");
    
    FKImpuritySolver Solver(U,mu,e_d,Delta);

    std::unique_ptr<SelfConsistency<FKImpuritySolver>> SC_ptr;

    typedef FKOptionParserDMFT::SC enumSC;
    switch (sc_switch) {
        case enumSC::Bethe:       SC_ptr.reset(new BetheSC<FKImpuritySolver>(Solver, t)); break;
        case enumSC::DMFTCubic1d: SC_ptr.reset(new CubicDMFTSC<FKImpuritySolver,1, KPOINTS>(Solver, t)); D=1; break;
        case enumSC::DMFTCubic2d: SC_ptr.reset(new CubicDMFTSC<FKImpuritySolver,2, KPOINTS>(Solver, t)); D=2; break;
        case enumSC::DMFTCubic3d: SC_ptr.reset(new CubicDMFTSC<FKImpuritySolver,3, KPOINTS>(Solver, t)); D=3; break;
        case enumSC::DMFTCubic4d: SC_ptr.reset(new CubicDMFTSC<FKImpuritySolver,4, KPOINTS>(Solver, t)); D=4; break;
        case enumSC::DMFTCubicInfd: SC_ptr.reset(new CubicInfDMFTSC<FKImpuritySolver>(Solver,t,RealGrid(-6.0*t,6.0*t,1024))); break;
        default:                  ERROR("No self-consistency provided. Exiting..."); exit(1); 
    };
    auto &SC = *SC_ptr;

    RealType diff=1.0;
    for (int i=0; i<maxit && diff>1e-8 &&!interrupt; ++i) {
        INFO("Iteration " << i <<". Mixing = " << mix);
        if (diff/mix>1e-3) Solver.run(true);
        else Solver.run(false);
        Delta = SC();
        auto Delta_new = Delta*mix+(1.0-mix)*Solver.Delta;
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);
        Solver.Delta = Delta_new;
        }
   
    GF Delta_half(gridF_half); Delta_half = Delta;
    GF gw_half(gridF_half); gw_half = Solver.gw;
    GF sigma_half(gridF_half); sigma_half = Solver.Sigma;
    sigma_half.savetxt("Sigma.dat");
    gw_half.savetxt("Gw.dat");
    Delta_half.savetxt("Delta.dat");
    Solver.Delta.savetxt("Delta_full.dat");

    // Everything important is finished by here
    if (extra_ops) {
        switch (sc_switch) {
            case enumSC::DMFTCubic1d: 
                calcStats(*(static_cast<CubicDMFTSC<FKImpuritySolver,1, KPOINTS>*> (SC_ptr.get())), gridF); 
                break;
            case enumSC::DMFTCubic2d: 
                calcStats(*(static_cast<CubicDMFTSC<FKImpuritySolver,2, KPOINTS>*> (SC_ptr.get())), gridF); 
                break;
            case enumSC::DMFTCubic3d: 
                calcStats(*(static_cast<CubicDMFTSC<FKImpuritySolver,3, KPOINTS>*> (SC_ptr.get())), gridF); 
                break;
            case enumSC::DMFTCubic4d: 
                calcStats(*(static_cast<CubicDMFTSC<FKImpuritySolver,4, KPOINTS>*> (SC_ptr.get())), gridF); 
                break;
            case enumSC::DMFTCubicInfd: 
                calcStats(*(static_cast<CubicInfDMFTSC<FKImpuritySolver>*> (SC_ptr.get())), gridF); 
                break;
 
            default: break;
            }; 
        };
    SC_ptr.release(); 
}

