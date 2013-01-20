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
    RealType mu = U/2;
    RealType e_d = 0.0;
    RealType beta = 10;
    RealType T = 1.0/beta;
    RealType t = 1.0; 
    size_t maxit = 1000;
    RealType mix = 0.5;
    
    size_t n_freq = 512;
    size_t n_b_freq = 15;
    Log.setDebugging(true);
    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    BMatsubaraGrid gridB(-n_b_freq, n_b_freq+1, beta);
    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*2.0*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    RealType diff=1.0;
    CubicInfDMFTSC<FKImpuritySolver> SC(Solver, t, RealGrid(-6.0*t,6.0*t,1024));

    for (int i=0; i<maxit && diff>1e-8; ++i) {
        INFO("Iteration " << i);
        if (i<3) Solver.run(true);
        else Solver.run(false);
        auto Delta_new = SC();
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);
        Delta = Delta_new*mix + (1.0-mix)*Solver.Delta;
        Solver.Delta = Delta; 
        }

    bool success = false;

    INFO("Static susc");
    auto Bubbleq0 = SC.getBubble0(0.0);
    auto BubbleqPI = SC.getBubblePI(0.0); 
    auto StaticV4 = SC.getStaticLatticeDMFTVertex4();

    auto FullVertexqPI = StaticV4.getData().getAsMatrix();
    auto FullVertexq0 = StaticV4.getData().getAsMatrix();
    auto Chiq0 = Bubbleq0.getData().getAsDiagonalMatrix();
    auto ChiqPI = BubbleqPI.getData().getAsDiagonalMatrix();

    FullVertexq0 = Diagrams::BS(Chiq0, FullVertexq0, true);
    FullVertexqPI = Diagrams::BS(ChiqPI, FullVertexqPI, true);

    GF susc0(gridF), suscPI(gridF);
    for (auto w1: gridF.getVals()) { 
        susc0[size_t(w1)]+=Bubbleq0(w1);
        suscPI[size_t(w1)]+=BubbleqPI(w1);
        for (auto w2: gridF.getVals()) {
            susc0[size_t(w1)]+=Bubbleq0(w1)*FullVertexq0(size_t(w1),size_t(w2))*Bubbleq0(w2); 
            suscPI[size_t(w1)]+=BubbleqPI(w1)*FullVertexqPI(size_t(w1),size_t(w2))*BubbleqPI(w2); 
            }
        };
    INFO("T = " << T);
    INFO("Static bare q=0 susc = " << Bubbleq0.sum());
    INFO("Static bare q=PI susc = " << BubbleqPI.sum());
    INFO("Static q=0 susc = " << susc0.sum());
    INFO("Static q=PI susc = " << suscPI.sum());
    success = true;
    return EXIT_SUCCESS;
}
