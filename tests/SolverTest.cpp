#include <numeric>

#include "Solver.h"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;

//typedef GridObject<ComplexType,FMatsubaraGrid> GF;
typedef GridObject<ComplexType,FMatsubaraGrid> GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}

int main()
{
    RealType beta = 4;
    size_t n_freq = 100;
    //Log.setDebugging(true);
    FMatsubaraGrid grid(-n_freq, n_freq, beta);
    BMatsubaraGrid gridB(0,100,beta);
    GF Delta(grid);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [](ComplexType w) -> ComplexType {return 1.0/w;};
    Delta.fill(f1);
    RealType U = 5.0;
    RealType mu = U/2.0;
    RealType e_d = 0.0;
    
    FKImpuritySolver Solver(U,mu,e_d,Delta);
       //DEBUG("Delta = " << Delta);
    Solver.run();
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid>::PointFunctionType f2 = 
        [Solver](typename FMatsubaraGrid::point w1, typename FMatsubaraGrid::point w2) { return Solver.getFVertex4<FMatsubaraGrid::point, FMatsubaraGrid::point>(w1,w2); };
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> g44(std::make_tuple(grid,grid));
    g44.fill(f2);

    GridObject<ComplexType,FMatsubaraGrid> VertexG(grid);
    auto iW = gridB[80];
    typename GridObject<ComplexType,FMatsubaraGrid>::PointFunctionType fV = 
        [Solver,iW](typename FMatsubaraGrid::point w1) {
        return Solver.getBVertex4<BMatsubaraGrid::point, FMatsubaraGrid::point>(iW, w1);};
    typename GridObject<ComplexType,FMatsubaraGrid>::FunctionType fV_2 =
        [Solver,iW](ComplexType w1) { 
            return Solver.getBVertex4<ComplexType, ComplexType>(iW, w1); };
    VertexG.fill(fV);
    VertexG._f = fV_2;

    if (!is_equal(VertexG(FMatsubara(grid.getSize()-1,beta)),-beta*U*U/4.0,1e-1)) return EXIT_FAILURE;;
    if (!is_equal(VertexG(FMatsubara(grid.getSize(),beta)),-beta*U*U/4.0,1e-1)) return EXIT_FAILURE;;
    if (!is_equal(VertexG(FMatsubara(grid.getSize()+1,beta)),-beta*U*U/4.0,1e-1)) return EXIT_FAILURE;;

    iW = gridB[0];
    fV = [Solver,iW](typename FMatsubaraGrid::point w1) { return Solver.getBVertex4<BMatsubaraGrid::point, FMatsubaraGrid::point>(iW, w1);};
    VertexG.fill(fV);

    auto &Gw = Solver.gw;
    auto Gw_shift = Solver.gw.shift(gridB[3]);
    DEBUG(Gw_shift(FMatsubara(grid.getSize()-1,beta)));
    DEBUG(Gw._f(FMatsubara(grid.getSize(),beta)));

}
