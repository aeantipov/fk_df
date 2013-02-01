#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"
#include "GFWrap.h"
#include "Solver.h"
#include "SelfConsistency.h"
#include "DF.h"

#include "FKOptionParserDF.h"

#include <iostream>
#include <ctime>
#include <array>
#include <unordered_map>
#include <csignal>

using namespace FK;

RealType beta;
size_t D = 0;
size_t extraops=0;

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

template <class SCType> void getExtraData(SCType& SC, const FMatsubaraGrid& gridF)
{
    constexpr size_t D = SCType::NDim;
    INFO("Calculating additional statistics.");
    const auto &Solver = SC._S;
    RealType beta = Solver.beta;
    RealType T=1.0/beta;
    SC.GLatLoc.savetxt("gloc.dat");

    if (extraops>=1) {
    }

    if (extraops >= 2) {
        size_t n_b_freq = gridF._w_max/2; // std::max(gridF._w_max/2,int(beta));
        BMatsubaraGrid gridB(-n_b_freq, n_b_freq+1, beta);

        auto glat = SC.getGLat();
        auto gloc = SC.GLatLoc;
        size_t ksize = SC._kGrid.getSize();
        
        GF iw_gf(gridF); 
        iw_gf.fill([](ComplexType w){return w;});
        GridObject<ComplexType,BMatsubaraGrid> chi_q0(gridB), chi_qPI(gridB);


        KMeshPatch grid0PI(SC._kGrid, {{0, SC._kGrid.getSize()/2}} );
        std::array<KMeshPatch, D> qgrids = __repeater<KMeshPatch,D>::get_array(grid0PI); 
        auto data = SC.calculateLatticeData(gridB, qgrids); // Heavy operation

        auto LatticeSusc = std::get<0>(data);

        std::array<KMesh::point,SCType::NDim> q_0, q_PI;
        q_0.fill(SC._kGrid[0]);
        q_PI.fill(SC._kGrid[ksize/2]);

        for (auto iW : gridB.getPoints()) {
            auto args_0 = std::tuple_cat(std::forward_as_tuple(iW),q_0);
            auto args_pi = std::tuple_cat(std::forward_as_tuple(iW),q_PI);
            chi_q0[size_t(iW)] = LatticeSusc(args_0);
            chi_qPI[size_t(iW)] = LatticeSusc(args_pi);
            };

        chi_q0.savetxt("Chiq0.dat");
        chi_qPI.savetxt("ChiqPI.dat");
        auto chi_q0_0 = -T*chi_q0.sum();
        auto chi_qPI_0 = -T*chi_qPI.sum();
        INFO("Chi0(q=0) sum  = " << chi_q0_0);
        INFO("Chi0(q=pi) sum = " << chi_qPI_0);
    };

}


int main(int argc, char *argv[])
{
  // Catch CTRL-C
  std::signal(SIGABRT, &sighandler);
  std::signal(SIGTERM, &sighandler);
  std::signal(SIGINT , &sighandler);

  FKOptionParserDF opt;
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
        std::cout << "Max number of DMFT iterations : " << opt.NDMFTRuns << std::endl;
        std::cout << "Max number of DF   iterations : " << opt.NDFRuns << std::endl;
        INFO("Extra options flag: " << opt.extraops);
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
    size_t ksize = opt.kpts;
    size_t n_dual_freq = opt.n_dual_freq;
    RealType mix = opt.mix;
    auto sc_switch = opt.sc_index;
    extraops = opt.extraops;
    size_t NDMFTRuns = opt.NDMFTRuns;
    size_t NDFRuns = opt.NDFRuns;

    KMesh kGrid(ksize);

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
    typedef FKOptionParserDF::SC enumSC;
    KMeshPatch qGrid(kGrid);
    switch (sc_switch) {
        case enumSC::DFCubic1d: 
  //          SC_DMFT_ptr.reset(new CubicDMFTSC<FKImpuritySolver,1, ksize>(Solver, t));
  //          SC_DF_ptr.reset(new DFLadder<FKImpuritySolver, 1, ksize>(Solver, gridF, BMatsubaraGrid(-n_dual_freq+1,n_dual_freq, beta), t)); 
            D=1; break;
        case enumSC::DFCubic2d: 
            SC_DMFT_ptr.reset(new CubicDMFTSC<FKImpuritySolver,2>(Solver, t, kGrid));
            typedef DFLadder<FKImpuritySolver, 2> DFSCType;
            SC_DF_ptr.reset(new DFSCType(Solver, gridF, kGrid, BMatsubaraGrid(-n_dual_freq+1,n_dual_freq, beta), t)); 
            static_cast<DFSCType*> (SC_DF_ptr.get())->_n_GD_iter = opt.DFNumberOfSelfConsistentIterations;
            static_cast<DFSCType*> (SC_DF_ptr.get())->_GDmix = opt.DFSCMixing;
            static_cast<DFSCType*> (SC_DF_ptr.get())->_eval_BS_SC = opt.DFEvaluateBSSelfConsistent;
            static_cast<DFSCType*> (SC_DF_ptr.get())->_n_BS_iter = opt.DFNumberOfBSIterations;
            static_cast<DFSCType*> (SC_DF_ptr.get())->_BSmix = opt.DFBSMixing;
            D=2; break;
        case enumSC::DFCubic3d: 
//            SC_DMFT_ptr.reset(new CubicDMFTSC<FKImpuritySolver,3, ksize>(Solver, t));
//            SC_DF_ptr.reset(new DFLadder<FKImpuritySolver, 3, ksize>(Solver, gridF, BMatsubaraGrid(-n_dual_freq+1,n_dual_freq, beta), t)); 
            D=3; break;
        case enumSC::DFCubic4d: 
//            SC_DMFT_ptr.reset(new CubicDMFTSC<FKImpuritySolver,4, ksize>(Solver, t));
//            SC_DF_ptr.reset(new DFLadder<FKImpuritySolver, 4, ksize>(Solver, gridF, BMatsubaraGrid(-n_dual_freq+1,n_dual_freq, beta), t)); 
            D=4; break;
        default:
            ERROR("Couldn't find the desired SC type. Exiting.");
            exit(0);
    };
    auto &SC_DMFT = *SC_DMFT_ptr;
    auto &SC_DF   = *SC_DF_ptr;
  
    RealType diff=1.0;
    bool calc_DMFT = true;

    size_t i_dmft = 0; 
    size_t i_df = 0;

    for (; i_dmft<NDMFTRuns && i_df<=NDFRuns && diff>1e-8 &&!interrupt; (calc_DMFT)?i_dmft++:i_df++) {
        INFO("Iteration " << i_dmft+i_df <<". Mixing = " << mix);
        if (diff/mix>1e-3) Solver.run(true);
        else Solver.run(false);
        if (calc_DMFT) {  
            Delta = SC_DMFT();
            }
        else { 
            Delta = SC_DF();
             }
        auto Delta_new = Delta*mix+(1.0-mix)*Solver.Delta;
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);
        Solver.Delta = Delta_new;
        if (diff<=1e-8 && calc_DMFT) { 
            Solver.Sigma.savetxt("SigmaDMFT.dat"); 
            Solver.Delta.savetxt("DeltaDMFT.dat");
            Solver.gw.savetxt("GwDMFT.dat");
            diff = 1.0; calc_DMFT = false; mix = 1.0; }; // now continue with DF 
        }
   
    GF Delta_half(gridF_half); Delta_half = Delta;
    GF gw_half(gridF_half); gw_half = Solver.gw;
    GF sigma_half(gridF_half); sigma_half = Solver.Sigma;
    sigma_half.savetxt("Sigma.dat");
    gw_half.savetxt("Gw.dat");
    Delta_half.savetxt("Delta.dat");
    Solver.Delta.savetxt("Delta_full.dat");

    if (extraops>0) {
        switch (sc_switch) {
            case enumSC::DFCubic1d: 
                getExtraData(*(static_cast<DFLadder<FKImpuritySolver,1>*> (SC_DF_ptr.get())), gridF); 
                break;
            case enumSC::DFCubic2d: 
                getExtraData(*(static_cast<DFLadder<FKImpuritySolver,2>*> (SC_DF_ptr.get())), gridF); 
                break;
            case enumSC::DFCubic3d: 
                getExtraData(*(static_cast<DFLadder<FKImpuritySolver,3>*> (SC_DF_ptr.get())), gridF); 
                break;
            case enumSC::DFCubic4d: 
                getExtraData(*(static_cast<DFLadder<FKImpuritySolver,4>*> (SC_DF_ptr.get())), gridF); 
                break;
            default: break;
            }; 
        };


    SC_DMFT_ptr.release(); 
    SC_DF_ptr.release();
}

