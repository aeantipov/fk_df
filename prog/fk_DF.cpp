#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"
#include "GFWrap.h"
#include "Solver.h"
#include "SelfConsistency.h"
#include "DF.h"

#include "OptionParser.h"

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
size_t D = 0;

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

template <class SCType> void getGwBubble(const SCType& SC, const FMatsubaraGrid& gridF)
{
    assert(D);
    INFO("Calculating additional statistics.");
    const auto &Solver = SC._S;
    RealType beta = Solver.beta;
    RealType T=1.0/beta;
    size_t n_b_freq = 40; // std::max(30,int(beta));
    BMatsubaraGrid gridB(-n_b_freq, n_b_freq+1, beta);
    auto glat = SC.getGLat();
    GF iw_gf(gridF); 
    iw_gf.fill([](ComplexType w){return w;});
    GridObject<RealType,BMatsubaraGrid> chi0_q0(gridB), chi0_qPI(gridB);
    for (auto iW : gridB.getVals()) {
        INFO("iW = " << iW);
        std::array<KMesh::point,SCType::NDim> a0, api;
        a0.fill(SC._kGrid[0]);
        api.fill(SC._kGrid[KPOINTS/2]);
        auto args_0 = std::tuple_cat(std::forward_as_tuple(iW),a0);
        auto args_pi = std::tuple_cat(std::forward_as_tuple(iW),api);
        auto glat_shift = glat.shift(args_0);
        auto glat_shift_pi = glat.shift(args_pi);
        chi0_q0[size_t(iW)] = -T*std::real((glat*glat_shift).sum()/RealType(pow(KPOINTS,D)));
        chi0_qPI[size_t(iW)] = -T*std::real((glat*glat_shift_pi).sum()/RealType(pow(KPOINTS,D)));
        };

    chi0_q0.savetxt("Chi0q0.dat");
    chi0_qPI.savetxt("Chi0qPI.dat");

    RealType chi0_q0_0 = T*chi0_q0.sum();
    RealType chi0_qPI_0 = T*chi0_qPI.sum();
    INFO("Chi0(q=0) sum  = " << chi0_q0_0);
    INFO("Chi0(q=pi) sum = " << chi0_qPI_0);
}


int main(int argc, char *argv[])
{
  // Catch CTRL-C
  std::signal(SIGABRT, &sighandler);
  std::signal(SIGTERM, &sighandler);
  std::signal(SIGINT , &sighandler);

  FKOptionParser opt;
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
    size_t n_dual_freq = opt.n_dual_freq;
    size_t maxit = opt.n_iter;
    RealType mix = opt.mix;
    auto sc_switch = opt.sc_index;
    bool extra_ops = opt.extra_ops;

    KMesh kGrid(KPOINTS);

    Log.setDebugging(true);

    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    FMatsubaraGrid gridF_half(0, 2*n_freq, beta);

    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1, f2;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};
    try { Delta.loadtxt("Delta_full.dat"); } 
    catch (std::exception &e) { Delta.fill(f1); };
    Delta.savetxt("Delta_0.dat");
    
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    
    std::unique_ptr<SelfConsistency<FKImpuritySolver>> SC_DF_ptr, SC_DMFT_ptr;
    typedef FKOptionParser::SC enumSC;
    KMeshPatch qGrid(kGrid);
    switch (sc_switch) {
        case enumSC::DFCubic1d: 
            SC_DMFT_ptr.reset(new CubicDMFTSC<FKImpuritySolver,1, KPOINTS>(Solver, t));
            SC_DF_ptr.reset(new DFLadder<FKImpuritySolver, 1, KPOINTS>(Solver, gridF, BMatsubaraGrid(0,n_dual_freq, beta), {{qGrid}}, t)); 
            D=1; break;
        case enumSC::DFCubic2d: 
            SC_DMFT_ptr.reset(new CubicDMFTSC<FKImpuritySolver,2, KPOINTS>(Solver, t));
            SC_DF_ptr.reset(new DFLadder<FKImpuritySolver, 2, KPOINTS>(Solver, gridF, BMatsubaraGrid(0,n_dual_freq, beta), {{qGrid,qGrid}}, t)); 
            D=2; break;
        case enumSC::DFCubic3d: 
            SC_DMFT_ptr.reset(new CubicDMFTSC<FKImpuritySolver,3, KPOINTS>(Solver, t));
            SC_DF_ptr.reset(new DFLadder<FKImpuritySolver, 3, KPOINTS>(Solver, gridF, BMatsubaraGrid(0,n_dual_freq, beta), {{qGrid,qGrid,qGrid}}, t)); 
            D=3; break;
        case enumSC::DFCubic4d: 
            SC_DMFT_ptr.reset(new CubicDMFTSC<FKImpuritySolver,4, KPOINTS>(Solver, t));
            SC_DF_ptr.reset(new DFLadder<FKImpuritySolver, 4, KPOINTS>(Solver, gridF, BMatsubaraGrid(0,n_dual_freq, beta), {{qGrid,qGrid,qGrid,qGrid}}, t)); 
            D=4; break;
        default:
            ERROR("Couldn't find the desired SC type. Exiting.");
            exit(0);
    };
    auto &SC_DMFT = *SC_DMFT_ptr;
    auto &SC_DF   = *SC_DF_ptr;
  
    RealType diff=1.0;
    bool calc_DMFT = true;
    for (int i=0; i<maxit && diff>1e-8 &&!interrupt; ++i) {
        INFO("Iteration " << i <<". Mixing = " << mix);
        if (diff/mix>1e-3) Solver.run(true);
        else Solver.run(false);
        if (calc_DMFT) {  
            Delta = SC_DMFT();
            }
        else { 
            Delta = SC_DF();
//            break; 
             }
        auto Delta_new = Delta*mix+(1.0-mix)*Solver.Delta;
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);
        Solver.Delta = Delta_new;
        if (diff<=1e-8 && calc_DMFT) { diff = 1.0; calc_DMFT = false; }; // now continue with DF 
        }
   
    GF Delta_half(gridF_half); Delta_half = Delta;
    GF gw_half(gridF_half); gw_half = Solver.gw;
    GF sigma_half(gridF_half); sigma_half = Solver.Sigma;
    sigma_half.savetxt("Sigma.dat");
    gw_half.savetxt("Gw.dat");
    Delta_half.savetxt("Delta.dat");
    Solver.Delta.savetxt("Delta_full.dat");

    if (extra_ops && D) {
        switch (sc_switch) {
            case enumSC::DFCubic1d: 
                getGwBubble(*(static_cast<DFLadder<FKImpuritySolver,1, KPOINTS>*> (SC_DF_ptr.get())), gridF); 
                break;
            case enumSC::DFCubic2d: 
                getGwBubble(*(static_cast<DFLadder<FKImpuritySolver,2, KPOINTS>*> (SC_DF_ptr.get())), gridF); 
                break;
            case enumSC::DFCubic3d: 
                getGwBubble(*(static_cast<DFLadder<FKImpuritySolver,3, KPOINTS>*> (SC_DF_ptr.get())), gridF); 
                break;
            case enumSC::DFCubic4d: 
                getGwBubble(*(static_cast<DFLadder<FKImpuritySolver,4, KPOINTS>*> (SC_DF_ptr.get())), gridF); 
                break;
            default: break;
            }; 
        };


    SC_DMFT_ptr.release(); 
    SC_DF_ptr.release();
}

/*
    DEBUG(Delta(FMatsubara(gridF._w_max-1, beta)));
    DEBUG(Delta(FMatsubara(gridF._w_max, beta)));
    DEBUG(Solver.gw(FMatsubara(gridF._w_max-1, beta)));
    DEBUG(Solver.gw(FMatsubara(gridF._w_max, beta)));
    DEBUG(Solver.Sigma(FMatsubara(gridF._w_max-1, beta)));
    DEBUG(Solver.Sigma(FMatsubara(gridF._w_max, beta)));
*/
   // exit(0);
 
