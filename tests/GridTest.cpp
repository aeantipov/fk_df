#include <numeric>

#include "FKCommon.h"
#include "Grid.h"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;

int main()
{
    Log.setDebugging(true);

    FMatsubaraGrid n1(-100,100,10);
    FMatsubaraGrid n2(0,32,20);

    std::function<ComplexType(ComplexType)> F2=[](ComplexType x) {return x;};
    ComplexType out = n1.integrate<std::function<ComplexType(ComplexType)> > (F2);
    if (std::abs(out)>1e-8) return EXIT_FAILURE;
    INFO(out);
    DEBUG(n2);

    auto a1 = std::make_tuple(n1,n2);
    
    std::vector<int> x(n2.getSize());
    for (int i=0; i<x.size(); ++i) x[i]=i*i;

    INFO(std::get<1>(n2.find(FMatsubara(10,20))));
    if (std::get<1>(n2.find(FMatsubara(10,20)))!=10) return EXIT_FAILURE;
    INFO(std::get<1>(n2.find(FMatsubara(20,20))));
    if (std::get<1>(n2.find(FMatsubara(20,20)))!=20) return EXIT_FAILURE;
    INFO(std::get<1>(n2.find(FMatsubara(32,20))));
    if (std::get<1>(n2.find(FMatsubara(32,20)))!=32) return EXIT_FAILURE;
    INFO(std::get<1>(n2.find(FMatsubara(35,20))));

    INFO(n2.getValue(x,FMatsubara(10,20)));
    if (n2.getValue(x,FMatsubara(10,20))!=100) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
