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
    RealType U = 5.0;
    RealType mu = U/2.0;
    RealType e_d = 0.0;
    RealType beta = 4.0;
    RealType T=1.0/beta;
    RealType t = 1.0; 
    size_t maxit = 1000;
    RealType mix = 0.5;

    size_t n_freq = 256;
    size_t n_b_freq = 2;

    static const size_t KPOINTS=16;
    static const size_t D=2;
    Log.setDebugging(true);

    MatrixType<ComplexType> Chi0(5,5), V4(5,5);
    Chi0.setIdentity(); Chi0*=0.1;
    V4.setIdentity();

    DEBUG(V4);
    DEBUG(Chi0);

    if ((Diagrams::BS(Chi0,V4,true) - Diagrams::BS(Chi0,V4,true,true,50)).norm() > 1e-8) return EXIT_FAILURE;

    V4.setOnes();
    V4.diagonal().setZero();
    DEBUG(V4);
    if (((Diagrams::BS(Chi0,V4,true) - Diagrams::BS(Chi0,V4,true,true,50)).norm()) > 1e-8) return EXIT_FAILURE;

/*

    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    BMatsubaraGrid gridB(1, n_b_freq, beta);

    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*2.0*D*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    RealType diff=1.0;
    KMesh kGrid(KPOINTS);
    KMeshPatch qGrid(kGrid);
    //std::array<KMeshPatch,2> qGrids( {{ qGrid, qGrid }}) ; 
    CubicDMFTSC<D> SC(Solver, t, KMesh(KPOINTS));
    
    for (int i=0; i<maxit && diff>1e-8; ++i) {
        INFO("Iteration " << i);
        Solver.run();
        auto Delta_new = SC();
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);
        Delta = Delta_new*mix + (1.0-mix)*Solver.Delta;
        Solver.Delta = Delta; 
        }

    DFLadder<D> SC_DF(Solver, gridF, SC._kGrid, gridB, t);
*/
 
    return EXIT_SUCCESS;
}
