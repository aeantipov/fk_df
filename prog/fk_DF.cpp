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

typedef GridObject<ComplexType,FMatsubaraGrid> GF;
typedef GridObject<ComplexType,KMesh,FMatsubaraGrid> GF1d;
typedef GridObject<ComplexType,KMesh,KMesh,FMatsubaraGrid> GF2d;
typedef GridObject<ComplexType,KMesh,KMesh,KMesh,FMatsubaraGrid> GF3d;

int f44(int a, int b, int c, int d) { return a+b+c+d; };

template <typename...> struct map1;
template<class T1, typename ...OtherArgs>  
struct map1<T1,OtherArgs...> { 
    std::function< int (OtherArgs...)> operator() (
    std::function<int(T1,OtherArgs...)> &in, const T1& t1) { 
        std::function<int(OtherArgs...)> out = [&](const OtherArgs&... other){return in(t1,other...);};
        return out;
        };
};

std::function<int(int,int,int,int)> f44_1 = f44;
std::function<int(int,int,int)> f24;

int main()
{
    
    INFO("Hi! Doing Falicov-Kimball. ");
    RealType U = 4.0;
    RealType mu = U/2;
    RealType e_d = 0.0;
    RealType beta = 10;
    RealType t = 0.5; 
    size_t maxit = 100;

    size_t n_freq = 1024;
    Log.setDebugging(true);
    FMatsubaraGrid grid(-n_freq, n_freq, beta);
    GF Delta(grid);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};
    Delta = f1;
    //DEBUG(Delta);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    RealType diff=1.0;
    CubicDMFTSC<2> SC(t,32);
    //BetheSC SC(t);
    //CubicInfDMFTSC SC(t,RealGrid(-5,5,1024));

    for (int i=0; i<maxit && diff>1e-8; ++i) {
        INFO("Iteration " << i);
        Solver.run();
        auto G1 = SC(Solver.gw, Solver.Delta);
        auto diffG = Solver.Delta - G1;
        diff = std::real(grid.integrate(diffG.conj()*diffG));
        INFO("diff = " << diff);
        Solver.Delta = G1;
        }
    Solver.Delta.savetxt("Delta.dat");
}
