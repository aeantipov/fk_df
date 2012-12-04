#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"

#include <iostream>
#include <ctime>
#include <array>
#include<Eigen/Sparse>

using namespace FK;

typedef GridObject<ComplexType,FMatsubaraGrid> GF;
typedef GridObject<ComplexType,KMesh,FMatsubaraGrid> GF1d;
typedef GridObject<ComplexType,KMesh,KMesh,FMatsubaraGrid> GF2d;
typedef GridObject<ComplexType,KMesh,KMesh,KMesh,FMatsubaraGrid> GF3d;

template <class GFType>
struct FK_ImpuritySolver
{
    const FMatsubaraGrid w_grid;
    const RealType beta;
    GFType& Delta;
    GFType gw; 
    GFType K0;
    GFType K1;
    GFType Sigma;
    RealType U;
    RealType mu;
    RealType e_d;
    ComplexType w_0;
    ComplexType w_1;
    FK_ImpuritySolver(GF& Delta);
    void run();
};

template <class GFType>
FK_ImpuritySolver<GFType>::FK_ImpuritySolver(GF& Delta): 
    w_grid(Delta.getGrid()), 
    beta(w_grid._beta), 
    Delta(Delta), 
    gw(GF(w_grid)), K0(GF(w_grid)), K1(GF(w_grid)), Sigma(GF(w_grid))
{
};

template <class GFType>
void FK_ImpuritySolver<GFType>::run()
{
    INFO("Running FK Solver, beta = " << beta << ", U = " << U << ", mu = " << mu << ", e_d = " << e_d);
    std::function<ComplexType(ComplexType)> K0f, K1f;
    K0f = [this](ComplexType w){return 1.0/(w+mu-Delta(w));};
    K1f = [this](ComplexType w){return 1.0/(w+mu-Delta(w)-U);};
    K0 = K0f;
    K1 = K1f;
    auto Zf1 = [this](ComplexType w, RealType mu1){return 1.0+I*(Delta(w)-mu1)/w;};
    std::function<ComplexType(RealType)> Zfprod = [this,Zf1](RealType mu1){return w_grid.prod(std::bind(Zf1, std::placeholders::_1, mu1));};

    ComplexType Z0 = Zfprod(mu);
    ComplexType Z1 = Zfprod(mu-U)*std::exp(beta*(mu-e_d-U/2));
    ComplexType Z=Z0+Z1;

    DEBUG("Z0 = " << Z0 << ", Z1 = " << Z1);
    w_0 = Z0/Z;
    w_1 = Z1/Z;
    
    INFO("w_0 = " << w_0 << "; w_1 = " << w_1 );
    
    
    //ComplexType Z0 = 
}

int main()
{
    RealType beta = 10;
    size_t n_freq = 5;
    Log.setDebugging(true);
    FMatsubaraGrid grid(-n_freq, n_freq, beta);
    GF Delta(grid);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [](ComplexType w) -> ComplexType {return 1.0/w;};
    Delta = f1;
    
    FK_ImpuritySolver<GF> Solver(Delta);
    Solver.U = 1.0;
    Solver.mu = 0.5;
    Solver.e_d = 0.0;

    INFO("Delta = " << Delta);
    Solver.run();
}
