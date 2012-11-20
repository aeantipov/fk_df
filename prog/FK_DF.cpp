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
    //GF D1(n2);

    std::array<size_t,2> Ar1 {{10,20}};
    Container<2,ComplexType> B(Ar1);
/*
    //INFO(D1);
   // INFO(D1[4]);
    INFO_NONEWLINE("Count: " << b1_c << ". "); std::cin >> b1; b1_c++;
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> D2(std::make_tuple(n1,n1));
    INFO_NONEWLINE("Count: " << b1_c << ". "); std::cin >> b1; b1_c++;
    DEBUG(D2.getGrid().getSize());
    DEBUG(D2[1].getGrid().getSize());
    D2[0].set(F2);
//    INFO(D2[0]);
//    INFO(D2[5]);
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid,FMatsubaraGrid> D3(std::make_tuple(n1,n1,n2));
    INFO_NONEWLINE("Count: " << b1_c << ". "); std::cin >> b1; b1_c++;
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid,FMatsubaraGrid,FMatsubaraGrid> D4(std::make_tuple(n1,n1,n2,n2));
    INFO_NONEWLINE("Count: " << b1_c << ". "); std::cin >> b1; b1_c++;
    //DEBUG(D2.getGrid().getSize());
*/
}
