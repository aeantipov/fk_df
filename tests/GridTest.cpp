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
    INFO(out);
    DEBUG(n2);

    auto a1 = std::make_tuple(n1,n2);

    return EXIT_SUCCESS;
}
