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

int main()
{
    Log.setDebugging(true);
    std::cout << "Hi!" << std::endl;
    typedef GridObject<ComplexType,FMatsubaraGrid> GF;
    FMatsubaraGrid n1(0,2,10);
    FMatsubaraGrid n2(0,5,20);
    auto a1 = std::make_tuple(n1,n2);

    GF D1(n2);
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> D2(std::make_tuple(n1,n2));
    
    D2.getData()[0][1]=4.0;
    D2.getData()[1][2]=3.1;
    DEBUG(D2);

    auto C1 = D2.getData();
    DEBUG(C1[0]);
    //decltype (C1[0]) x(std::make_tuple(3));
    const std::array<size_t, 1> arr{{5}};
    Container<1,ComplexType> C2(C1[0]);
    DEBUG(C2[1]);
    decltype (C1[0])& C2_2 = decltype(C1[0])(arr);
    DEBUG(n2.getValue(C2, FMatsubara(1,20)));
    DEBUG(n1.getValue(C1, FMatsubara(1,10)));
    auto C22 = n1.getValue(C1, FMatsubara(0,10));
    DEBUG(C22+C2-C22*2.0);
    DEBUG(n2.getValue(n1.getValue(C1, FMatsubara(1,10)), FMatsubara(1,20)));

    DEBUG(D2(FMatsubara(1,10),FMatsubara(2,20)));
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
