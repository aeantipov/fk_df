
#include "FKCommon.h"
#include "Grid.h"
#include "GridObject.h"

#include <iostream>

using namespace FK;

int main()
{
    Log.setDebugging(true);
    std::cout << "Hi!" << std::endl;
    typedef Grid1dObject<1,ComplexType,FMatsubaraGrid1d> GF;
    //Grid1d<RealType> a1;
    //std::function<RealType(const int &)> F1=[&](const int & n) {return PI/10*(2*n+1);};
    //std::cout<<F1(-100) << std::endl;
    //Grid1d<RealType> n1(-100,100,F1);
    //INFO(n1);
    FMatsubaraGrid1d n2(-10,10,10);
    std::shared_ptr<FMatsubaraGrid1d> a = std::make_shared<FMatsubaraGrid1d>(n2);
    //INFO(n2);
    //INFO(n2[2]);
    std::function<ComplexType(ComplexType)> F2=[](ComplexType x) {return x*x;};
    //INFO(F2(3));
    ComplexType out = n2.integrate<std::function<ComplexType(ComplexType)> > (F2);
    INFO(out);
    
    INFO(std::get<1>(n2.findSmaller(FMatsubara(5,7))));

    auto a1 = std::make_tuple(n2,n2);
    std::cout << "type check: " << std::is_same<FMatsubaraGrid1d, std::tuple_element<0,std::tuple<FMatsubaraGrid1d,FMatsubaraGrid1d> >::type >::value << std::endl;
    static_assert(std::is_same<FMatsubaraGrid1d, std::tuple_element<0,std::tuple<FMatsubaraGrid1d,FMatsubaraGrid1d> >::type >::value,"!!!");

    GF D1(n2);
   // INFO(D1);
   // INFO(D1[4]);
    Grid1dObject<2,ComplexType,FMatsubaraGrid1d,FMatsubaraGrid1d> D2(a1);
    //D2[0].set(F2);
    //FMatsubaraGrid1d n3;
    Eigen::VectorXf x (Eigen::VectorXf::Constant(2,3)); 
    DEBUG(x);
    std::tuple<int,double,int> t = std::make_tuple(1,2.3,2);
    std::tuple<double,int> t2;
    std::tuple<int> t2_2;
    std::tuple_cat(t2_2, t2) = t;
    DEBUG(std::get<0>(t2_2));
    DEBUG(std::get<0>(t));
}
