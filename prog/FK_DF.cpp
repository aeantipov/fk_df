
#include "FKCommon.h"
#include "Grid.h"
#include "GridObject.h"

#include <iostream>

using namespace FK;

int main()
{
    Log.setDebugging(true);
    std::cout << "Hi!" << std::endl;
    typedef GridObject<1,ComplexType,FMatsubaraGrid> GF;
    //Grid<RealType> a1;
    //std::function<RealType(const int &)> F1=[&](const int & n) {return PI/10*(2*n+1);};
    //std::cout<<F1(-100) << std::endl;
    //Grid<RealType> n1(-100,100,F1);
    //INFO(n1);
    FMatsubaraGrid n2(-10,10,10);
    FMatsubaraGrid n1(-3,3,20);
    std::shared_ptr<FMatsubaraGrid> a = std::make_shared<FMatsubaraGrid>(n2);
    //INFO(n2);
    //INFO(n2[2]);
    std::function<ComplexType(ComplexType)> F2=[](ComplexType x) {return x*x;};
    //INFO(F2(3));
    ComplexType out = n2.integrate<std::function<ComplexType(ComplexType)> > (F2);
    INFO(out);
    
    INFO(std::get<1>(n2.findSmaller(FMatsubara(5,7))));

    auto a1 = std::make_tuple(n1,n2);
    std::cout << "type check: " << std::is_same<FMatsubaraGrid, std::tuple_element<0,std::tuple<FMatsubaraGrid,FMatsubaraGrid> >::type >::value << std::endl;
    static_assert(std::is_same<FMatsubaraGrid, std::tuple_element<0,std::tuple<FMatsubaraGrid,FMatsubaraGrid> >::type >::value,"!!!");

   //GF D1(n2);
   // INFO(D1);
   // INFO(D1[4]);
    GridObject<2,ComplexType,FMatsubaraGrid,FMatsubaraGrid> D2(std::make_tuple(n1,n2));
    DEBUG(D2.getGrid().getSize());
    DEBUG(D2[0].getGrid().getSize());
    //D2[0].set(F2);
}
