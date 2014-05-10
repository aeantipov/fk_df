#include <numeric>

#include "Solver.h"
#include "DMFT.h"
#include "DF.h"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;
typedef grid_object<complex_type,fmatsubara_grid> GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, real_type tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}


int main()
{
    INFO("Hi! Doing Falicov-Kimball. ");
    real_type U = 4.0;
    real_type mu = 2.2;
    real_type e_d = 0.1;
    real_type beta = 10;
    real_type T=1.0/beta; T*=1;
    real_type t = 1.0; 
    size_t maxit = 1000;
    real_type mix = 0.5;

    size_t n_freq = 256;
    size_t n_b_freq = 15;

    static const size_t KPOINTS=32;
    static const size_t D=2;

    //Log.setDebugging(true);
    fmatsubara_grid gridF(-n_freq, n_freq, beta);
    bmatsubara_grid gridB(-n_b_freq, n_b_freq+1, beta);

    GF Delta(gridF);
    std::function<complex_type(complex_type)> f1;
    f1 = [t](complex_type w) -> complex_type {return t*2.0*D*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    real_type diff=1.0;
    kmesh kGrid(KPOINTS);
    kmesh_patch qGrid(kGrid);
    std::array<kmesh_patch,2> qGrids( {{ qGrid, qGrid }}) ; 
    DFLadderCubic<2> SC(Solver, gridF, kGrid, t);
    SC._n_GD_iter = 0;
    
    for (int i=0; i<maxit && diff>1e-8; ++i) {
        INFO("Iteration " << i);
        Solver.run();
        auto Delta_new = SC();
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);
        Delta = Delta_new*mix + (1.0-mix)*Solver.Delta;
        Solver.Delta = Delta; 
        }


    fmatsubara_grid grid_half(0,n_freq*2,beta);
    GF Delta_half(grid_half); 
    Delta_half.copy_interpolate(Solver.Delta);
    GF gw_half(grid_half); 
    gw_half.copy_interpolate(Solver.gw);
    GF sigma_half(grid_half); 
    sigma_half.copy_interpolate(Solver.Sigma);

    sigma_half.savetxt("Sigma.dat");
    gw_half.savetxt("Gw.dat");
    Delta_half.savetxt("Delta.dat");

    bool success = false;

    std::vector<complex_type> right_vals = {{ 1.901098273610857259e-02-2.130828360852165537e-01*I,2.161465786590606453e-02 -2.485358842730046314e-01*I }}; 

    for (size_t t=0; t<right_vals.size(); ++t) { 
        success = (is_equal(gw_half[t],right_vals[t],1e-3));
        if (!success) { ERROR(gw_half[t] << " != " << right_vals[t]); return EXIT_FAILURE; };
        };
    
    return EXIT_SUCCESS;
}
