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

typedef GFWrap GF;
typedef GridObject<ComplexType,KMesh,FMatsubaraGrid> GF1d;
typedef GridObject<ComplexType,KMesh,KMesh,FMatsubaraGrid> GF2d;
typedef GridObject<ComplexType,KMesh,KMesh,KMesh,FMatsubaraGrid> GF3d;

int main()
{
    RealType beta = 10;
    size_t n_freq = 10;
    Log.setDebugging(true);
    FMatsubaraGrid grid(-n_freq, n_freq, beta);
    GF Delta(grid);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [](ComplexType w) -> ComplexType {return 1.0/w;};
    Delta = f1;
    RealType U = 1.0;
    RealType mu = 0.5;
    RealType e_d = 0.0;
    
    FKImpuritySolver Solver(U,mu,e_d,Delta);
       //DEBUG("Delta = " << Delta);
    Solver.run();
 //std::function<ComplexType(ComplexType,ComplexType)> f1 = std::bind(g4, std::placeholders::_1, std::placeholders::_2);
    std::function<ComplexType(ComplexType,ComplexType)> f2 = std::bind(&FKImpuritySolver::getVertex4, Solver, std::placeholders::_1, std::placeholders::_2);
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> g44(std::make_tuple(grid,grid));
    g44.fill(f2);
    DEBUG(g44);


}
