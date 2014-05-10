#include <numeric>

#include "Solver.h"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;

//typedef grid_object<complex_type,fmatsubara_grid> GF;
typedef grid_object<complex_type,fmatsubara_grid> GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, real_type tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}

int main()
{
    real_type beta = 4;
    size_t n_freq = 100;
    //Log.setDebugging(true);
    fmatsubara_grid grid(-n_freq, n_freq, beta);
    bmatsubara_grid gridB(0,100,beta);
    GF Delta(grid);
    std::function<complex_type(complex_type)> f1;
    f1 = [](complex_type w) -> complex_type {return 1.0/w;};
    Delta.fill(f1);
    real_type U = 5.0;
    real_type mu = U/2.0;
    real_type e_d = 0.0;
    
    FKImpuritySolver Solver(U,mu,e_d,Delta);
       //DEBUG("Delta = " << Delta);
    Solver.run();
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid>::point_function_type f2 = 
        [Solver](typename fmatsubara_grid::point w1, typename fmatsubara_grid::point w2) { return Solver.getFVertex4<fmatsubara_grid::point, fmatsubara_grid::point>(w1,w2); };
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> g44(std::make_tuple(grid,grid));
    g44.fill(f2);

    grid_object<complex_type,fmatsubara_grid> VertexG(grid);
    auto iW = gridB[80];
    typename grid_object<complex_type,fmatsubara_grid>::point_function_type fV = 
        [Solver,iW](typename fmatsubara_grid::point w1) {
        return Solver.getBVertex4<bmatsubara_grid::point, fmatsubara_grid::point>(iW, w1);};
    typename grid_object<complex_type,fmatsubara_grid>::function_type fV_2 =
        [Solver,iW](complex_type w1) { 
            return Solver.getBVertex4<complex_type, complex_type>(iW, w1); };
    VertexG.fill(fV);
    VertexG.tail_ = fV_2;

    if (!is_equal(VertexG(FMatsubara(grid.size()-1,beta)),-beta*U*U/4.0,1e-1)) return EXIT_FAILURE;;
    if (!is_equal(VertexG(FMatsubara(grid.size(),beta)),-beta*U*U/4.0,1e-1)) return EXIT_FAILURE;;
    if (!is_equal(VertexG(FMatsubara(grid.size()+1,beta)),-beta*U*U/4.0,1e-1)) return EXIT_FAILURE;;

    iW = gridB[0];
    fV = [Solver,iW](typename fmatsubara_grid::point w1) { return Solver.getBVertex4<bmatsubara_grid::point, fmatsubara_grid::point>(iW, w1);};
    VertexG.fill(fV);

    auto &Gw = Solver.gw;
    auto Gw_shift = Solver.gw.shift(gridB[3].val_);
    DEBUG(Gw_shift(FMatsubara(grid.size()-1,beta)));
    DEBUG(Gw.tail_(FMatsubara(grid.size(),beta)));

}
