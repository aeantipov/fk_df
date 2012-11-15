
#include "FKCommon.h"
#include "Grid.h"
#include "GfWrap.h"

#include <iostream>

using namespace FK;

int main()
{
    Log.setDebugging(true);
    std::cout << "Hi!" << std::endl;
    //Grid1d<RealType> a1;
    //std::function<RealType(const int &)> F1=[&](const int & n) {return PI/10*(2*n+1);};
    //std::cout<<F1(-100) << std::endl;
    //Grid1d<RealType> n1(-100,100,F1);
    //INFO(n1);
    FMatsubaraGrid1d n2(-10,10,10);
    //INFO(n2);
    //INFO(n2[2]);
    std::function<ComplexType(ComplexType)> F2=[](ComplexType x) {return x;};
    //INFO(F2(3));
    ComplexType out = n2.integrate<std::function<ComplexType(ComplexType)> > (F2);
    INFO(out);
}
