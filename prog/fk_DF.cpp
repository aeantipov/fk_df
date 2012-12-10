#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"
#include "Solver.h"
#include "SelfConsistency.h"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;

typedef GridObject<ComplexType,FMatsubaraGrid> GF;
typedef GridObject<ComplexType,KMesh,FMatsubaraGrid> GF1d;
typedef GridObject<ComplexType,KMesh,KMesh,FMatsubaraGrid> GF2d;
typedef GridObject<ComplexType,KMesh,KMesh,KMesh,FMatsubaraGrid> GF3d;

int main()
{
    INFO("Hi! Doing Falicov-Kimball. ");
    RealType U = 4.0;
    RealType mu = U/2;
    RealType e_d = 0.0;
    RealType beta = 10;
    RealType t = 0.5; 
    size_t maxit = 100;

    size_t n_freq = 10;
    Log.setDebugging(true);
    FMatsubaraGrid grid(-n_freq, n_freq, beta);
    GF Delta(grid);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};
    Delta = f1;
    //DEBUG(Delta);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    RealType diff=1.0;
    CubicDMFTSC<3> SC(t,32);
    DEBUG(SC.dispersion(0.0,PI,PI/2.0));
    
    for (int i=0; i<maxit && diff>1e-8; ++i) {
        INFO("Iteration " << i);
        Solver.run();
        auto G1 = Solver.gw*(t*t);
        auto diffG = Solver.Delta - G1;
        diff = std::real(grid.integrate(diffG.conj()*diffG));
        INFO("diff = " << diff);
        Solver.Delta = G1;
        }
    auto G1 = SC.glat(Solver.gw, Solver.Delta, 0.0, PI, PI/2.0);
    auto G1_k = std::bind(&CubicDMFTSC<3>::glat<RealType,RealType,RealType>, SC, Solver.gw, Solver.Delta, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    auto G1_k1 = std::bind(&CubicDMFTSC<3>::glat<RealType,RealType,RealType>, SC, Solver.gw, Solver.Delta, 0.0, PI, std::placeholders::_1);
    auto G1_k2 = std::bind(&CubicDMFTSC<3>::glat<RealType,RealType,RealType>, SC, Solver.gw, Solver.Delta, 0.0, std::placeholders::_1, std::placeholders::_2);
    //DEBUG(G1_k(1.0,2.0,3.0));
    //DEBUG(G1_k1(1.0));
    //DEBUG(G1);
    auto ff1 = [&G1_k1](RealType k){return G1_k1(k);};
    auto ff2 = [&G1_k2](RealType k1, RealType k2){return G1_k2(k1,k2);};
    DEBUG(SC._kgrid.integrate(ff1)); 
    RecursiveGridIntegrator<1,decltype(ff1)> t1;
    auto t1_int = t1.integrate(SC._kgrid,ff1);
    DEBUG(t1_int);
    RecursiveGridIntegrator<2,decltype(ff2)> t2;
    exit(0);
 //std::function<ComplexType(ComplexType,ComplexType)> f1 = std::bind(g4, std::placeholders::_1, std::placeholders::_2);
    //std::function<ComplexType(ComplexType,ComplexType)> f2 = std::bind(&FKImpuritySolver::getVertex4, Solver, std::placeholders::_1, std::placeholders::_2);
    //GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> g44(std::make_tuple(grid,grid));
    //g44.fill(f2);
    //DEBUG(g44);


}
