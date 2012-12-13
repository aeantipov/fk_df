#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"
#include "GFWrap.h"
#include "Solver.h"
#include "SelfConsistency.h"

#include "OptionParser.h"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;

//typedef GridObject<ComplexType,FMatsubaraGrid> GF;
typedef GFWrap GF;
typedef GridObject<ComplexType,KMesh,FMatsubaraGrid> GF1d;
typedef GridObject<ComplexType,KMesh,KMesh,FMatsubaraGrid> GF2d;
typedef GridObject<ComplexType,KMesh,KMesh,KMesh,FMatsubaraGrid> GF3d;

int main(int argc, char *argv[])
{

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


    Log.setDebugging(true);
    FMatsubaraGrid grid(-n_freq, n_freq, beta);
    GF Delta(grid);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};
    Delta = f1;
    //DEBUG(Delta(FMatsubara(11,beta)));
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    RealType diff=1.0;
    CubicDMFTSC<2> SC(t,32);
    //BetheSC SC(t);
    //CubicInfDMFTSC SC(t,RealGrid(-5,5,1024));

    for (int i=0; i<maxit && diff>1e-8; ++i) {
        INFO("Iteration " << i);
        Solver.run();
        //DEBUG(Solver.gw);
        //exit(0);
        auto G1 = SC(Solver.gw, Solver.Delta);
        auto diffG = Solver.Delta - G1;
        diff = std::real(grid.integrate(diffG.conj()*diffG));
        INFO("diff = " << diff);
        Solver.Delta = G1;
        }
    Solver.Delta.savetxt("Delta.dat");
}
