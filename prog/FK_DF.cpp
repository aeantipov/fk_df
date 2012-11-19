#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "GridObject.h"

#include <iostream>
#include <ctime>

using namespace FK;

int main()
{
    Log.setDebugging(true);
    std::cout << "Hi!" << std::endl;
    typedef GridObject<ComplexType,FMatsubaraGrid> GF;
    //Grid<RealType> a1;
    //std::function<RealType(const int &)> F1=[&](const int & n) {return PI/10*(2*n+1);};
    //std::cout<<F1(-100) << std::endl;
    //Grid<RealType> n1(-100,100,F1);
    //INFO(n1);

    
    int b1, b1_c=0;
//    INFO_NONEWLINE("Count: " << b1_c << ". "); std::cin >> b1; b1_c++;
   // Container<1,ComplexType> a_(std::make_tuple(100));
//    INFO_NONEWLINE("Count: " << b1_c << ". "); std::cin >> b1; b1_c++;
   // Container<2,ComplexType> b_(std::make_tuple(100,100));
//    INFO_NONEWLINE("Count: " << b1_c << ". "); std::cin >> b1; b1_c++;
   // Container<3,ComplexType> c_(std::make_tuple(100,100,100));
//    INFO_NONEWLINE("Count: " << b1_c << ". "); std::cin >> b1; b1_c++;
    //Container<4,ComplexType> d_(std::make_tuple(100,100,32,32));
    //Container<6,ComplexType> e_(std::make_tuple(100,100,32,32,2,2));
//    INFO_NONEWLINE("Count: " << b1_c << ". "); std::cin >> b1; b1_c++;
    FMatsubaraGrid n1(-100,100,10);
    //FMatsubaraGrid n1(-3,3,10);
    FMatsubaraGrid n2(0,32,20);
    std::function<ComplexType(ComplexType)> F2=[](ComplexType x) {return x;};
    ComplexType out = n1.integrate<std::function<ComplexType(ComplexType)> > (F2);
    INFO(out);
/*
    INFO(std::get<1>(n2.findSmaller(FMatsubara(5,7))));

    auto a1 = std::make_tuple(n1,n2);
    std::cout << "type check: " << std::is_same<FMatsubaraGrid, std::tuple_element<0,std::tuple<FMatsubaraGrid,FMatsubaraGrid> >::type >::value << std::endl;
    static_assert(std::is_same<FMatsubaraGrid, std::tuple_element<0,std::tuple<FMatsubaraGrid,FMatsubaraGrid> >::type >::value,"!!!");

    INFO_NONEWLINE("Count: " << b1_c << ". "); std::cin >> b1; b1_c++;
    GF D1(n2);
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
