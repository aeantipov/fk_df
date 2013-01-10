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

//typedef GridObject<ComplexType,FMatsubaraGrid> GF;
typedef GFWrap GF;
typedef GridObject<ComplexType,KMesh,FMatsubaraGrid> GF1d;
typedef GridObject<ComplexType,KMesh,KMesh,FMatsubaraGrid> GF2d;
typedef GridObject<ComplexType,KMesh,KMesh,KMesh,FMatsubaraGrid> GF3d;

bool interrupt = false;

void sighandler(int signal)
{
    INFO("Caught interrupt, signal " << signal <<". Exiting...")
    interrupt = true;
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
    RealType beta = opt.beta;
    RealType t = opt.t; 
    size_t n_freq = opt.n_freq;
    size_t maxit = opt.n_iter;
    RealType mix = opt.mix;
    auto sc_index = opt.sc_index;
    
    static const size_t KPOINTS = 64;

    Log.setDebugging(true);

    FMatsubaraGrid grid(-n_freq, n_freq, beta);
    FMatsubaraGrid grid_half(0, n_freq*2, beta);
    GF Delta(grid);
    std::function<ComplexType(ComplexType)> f1, f2;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};

    try { Delta.loadtxt("Delta_full.dat"); } 
    catch (std::exception &e) { Delta.fill(f1); };

    Delta.savetxt("Delta_0.dat");
    
    FKImpuritySolver Solver(U,mu,e_d,Delta);

    std::unique_ptr<SelfConsistency<FKImpuritySolver>> SC_ptr;
    typedef FKOptionParser::SC enumSC;
    switch (sc_index) {
        case enumSC::Bethe:       SC_ptr.reset(new BetheSC<FKImpuritySolver>(Solver, t)); break;
        case enumSC::DMFTCubic1d: SC_ptr.reset(new CubicDMFTSC<FKImpuritySolver,1, KPOINTS>(Solver, t)); break;
        case enumSC::DMFTCubic2d: SC_ptr.reset(new CubicDMFTSC<FKImpuritySolver,2, KPOINTS>(Solver, t)); break;
        case enumSC::DMFTCubic3d: SC_ptr.reset(new CubicDMFTSC<FKImpuritySolver,2, KPOINTS>(Solver, t)); break;
        case enumSC::DMFTCubic4d: SC_ptr.reset(new CubicDMFTSC<FKImpuritySolver,2, KPOINTS>(Solver, t)); break;
        case enumSC::DMFTCubicInfd: SC_ptr.reset(new CubicInfDMFTSC<FKImpuritySolver>(Solver,t,RealGrid(-6.0*t,6.0*t,1024)));
        default:                  SC_ptr.reset(new BetheSC<FKImpuritySolver>(Solver, t)); break;
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
   
    GF Delta_half(grid_half); Delta_half = Delta;
    GF gw_half(grid_half); gw_half = Solver.gw;
    GF sigma_half(grid_half); sigma_half = Solver.Sigma;
    sigma_half.savetxt("Sigma.dat");
    gw_half.savetxt("Gw.dat");
    Delta_half.savetxt("Delta.dat");
    Solver.Delta.savetxt("Delta_full.dat");
   
    SC_ptr.release(); 
}

