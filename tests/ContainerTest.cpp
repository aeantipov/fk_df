#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;

int main()
{
    Log.setDebugging(true);
    std::cout << "Hi!" << std::endl;

    typedef GridObject<ComplexType,FMatsubaraGrid> GF;
    FMatsubaraGrid n1(-100,100,10);
    FMatsubaraGrid n2(0,32,20);

    std::function<ComplexType(ComplexType)> F2=[](ComplexType x) {return x;};
    ComplexType out = n1.integrate<std::function<ComplexType(ComplexType)> > (F2);
    INFO(out);
    DEBUG(n2);

    auto a1 = std::make_tuple(n1,n2);

    std::array<size_t,2> Ar1 {{1,3}};
    std::array<size_t,3> Ar2 {{1,1,3}};
    Container<2,ComplexType> B(Ar1);
    Container<2,ComplexType> C(Ar1);
    Container<3,ComplexType> D(Ar2);
    Container<3,ComplexType> E(D);
    B[0][2]=3.0;
    C[0][2]=-2.0;
    DEBUG(B+C*2);
    DEBUG(B[0]+C[0]);
    DEBUG((B+C)[0]);
    DEBUG((B[0]*3.0));
    DEBUG((B-C)[0]);
    D[0][0][0]=-1.0;
    D[0][0][1]=1.0;
    E*=(-1);
    DEBUG(D);
    DEBUG(E);
    DEBUG(D+E);
    DEBUG(D*5+2.0);
    DEBUG(D[0]+B);

    return EXIT_SUCCESS;
}
