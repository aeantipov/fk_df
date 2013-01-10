#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"
#include "Solver.h"
#include "SelfConsistency.h"
#include "DF.h"

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
    RealType mu = 2.2;
    RealType e_d = 0.1;
    RealType beta = 10;
    RealType T=1.0/beta; T*=1;
    RealType t = 1.0; 
    size_t maxit = 1000;
    RealType mix = 0.5;

    size_t n_freq = 256;
    size_t n_b_freq = 15;

    static const size_t KPOINTS=32;
    static const size_t D=2;

    Log.setDebugging(true);
    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    BMatsubaraGrid gridB(-n_b_freq, n_b_freq+1, beta);

    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*2.0*D*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    RealType diff=1.0;
    KMesh kGrid(KPOINTS);
    KMeshPatch qGrid(kGrid);
    std::array<KMeshPatch,2> qGrids( {{ qGrid, qGrid }}) ; 
    DFLadder<FKImpuritySolver, 2, KPOINTS> SC(Solver, gridF, gridB, qGrids, t);
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

    FMatsubaraGrid gridF_half(0,n_freq*2,beta);
    GF Delta_half(gridF_half); 
    Delta_half = Delta;
    GF gw_half(gridF_half); 
    gw_half = Solver.gw;
    GF sigma_half(gridF_half); 
    sigma_half = Solver.Sigma;

    sigma_half.savetxt("Sigma.dat");
    gw_half.savetxt("Gw.dat");
    Delta_half.savetxt("Delta.dat");

    bool success = false;

    std::vector<ComplexType> right_vals = {{ 1.901098273610857259e-02-2.130828360852165537e-01*I,2.161465786590606453e-02 -2.485358842730046314e-01*I }}; 

    for (size_t t=0; t<right_vals.size(); ++t) { 
        success = (is_equal(gw_half[t],right_vals[t],1e-3));
        if (!success) { ERROR(gw_half[t] << " != " << right_vals[t]); return EXIT_FAILURE; };
        };
    
    return EXIT_SUCCESS;
}
