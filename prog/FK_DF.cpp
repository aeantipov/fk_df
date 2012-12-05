#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"
#include "Solver.h"

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
    RealType U = 1.0;
    RealType mu = 0.5;
    RealType e_d = 0.0;
    RealType beta = 10;
    RealType t = 0.5; 

    size_t n_freq = 10;
    Log.setDebugging(true);
    FMatsubaraGrid grid(-n_freq, n_freq, beta);
    GF Delta(grid);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};
    Delta = f1;
    DEBUG(Delta);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    for (int i=0; i<10; ++i) {
        Solver.run();
        auto G1 = Solver.gw*(t*t);
        RealType diff=0.0;
        auto diffG = Solver.Delta - G1;
        diffG.getData().conj();
        Solver.Delta = G1;
        }
 //std::function<ComplexType(ComplexType,ComplexType)> f1 = std::bind(g4, std::placeholders::_1, std::placeholders::_2);
    std::function<ComplexType(ComplexType,ComplexType)> f2 = std::bind(&FKImpuritySolver::getVertex4, Solver, std::placeholders::_1, std::placeholders::_2);
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> g44(std::make_tuple(grid,grid));
    g44.fill(f2);
    //DEBUG(g44);


}
