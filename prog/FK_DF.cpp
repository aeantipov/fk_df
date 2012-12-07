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
    RealType mu = 0.5;
    RealType e_d = 0.01;
    RealType beta = 100;
    RealType t = 0.5; 
    size_t maxit = 100;

    size_t n_freq = 1000;
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
    
    exit(0);
    for (int i=0; i<maxit && diff>1e-8; ++i) {
        INFO("Iteration " << i);
        Solver.run();
        auto G1 = Solver.gw*(t*t);
        auto diffG = Solver.Delta - G1;
        diff = std::real(grid.integrate(diffG.conj()*diffG));
        INFO("diff = " << diff);
        Solver.Delta = G1;
        }
 //std::function<ComplexType(ComplexType,ComplexType)> f1 = std::bind(g4, std::placeholders::_1, std::placeholders::_2);
    std::function<ComplexType(ComplexType,ComplexType)> f2 = std::bind(&FKImpuritySolver::getVertex4, Solver, std::placeholders::_1, std::placeholders::_2);
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> g44(std::make_tuple(grid,grid));
    //g44.fill(f2);
    //DEBUG(g44);


}
