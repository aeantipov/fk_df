#include <numeric>

#include "Solver.h"
#include "DMFT.h"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;
typedef GFWrap GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}


int main()
{
    INFO("Hi! Doing Falicov-Kimball. ");
    RealType U = 4.0;
    RealType mu = 2.0;
    RealType e_d = 0.0;
    RealType beta = 3;
    RealType t = 1.0; 
    size_t maxit = 1000;
    RealType mix = 1.0;

    size_t n_freq = 16;
    //Log.setDebugging(true);
    FMatsubaraGrid grid(-n_freq, n_freq, beta);
    GF Delta(grid);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    RealType diff=1.0;
    CubicDMFTSC<2> SC(Solver, t, KMesh(16));

    for (int i=0; i<maxit && diff>1e-8; ++i) {
        INFO("Iteration " << i);
        Solver.run();
        auto Delta_new = SC();
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);
        Solver.Delta = Delta_new*mix + (1.0-mix)*Solver.Delta;
        DEBUG(Solver.Delta);
        }

    FMatsubaraGrid grid_half(0,n_freq*2,beta);
    GF Delta_half(grid_half); 
    Delta_half = Solver.Delta;
    GF gw_half(grid_half); 
    gw_half = Solver.gw;
    GF sigma_half(grid_half); 
    sigma_half = Solver.Sigma;

    sigma_half.savetxt("Sigma.dat");
    gw_half.savetxt("Gw.dat");
    Delta_half.savetxt("Delta.dat");

    bool success = false;

    std::vector<ComplexType> right_vals = {{ -2.499740280930e-01*I, -2.026008712834e-01*I }}; 

    for (size_t t=0; t<right_vals.size(); ++t) { 
        success = (is_equal(gw_half[t],right_vals[t],1e-3));
        if (!success) { ERROR(gw_half[t] << " != " << right_vals[t]); return EXIT_FAILURE; };
        };
    
    return EXIT_SUCCESS;
}
