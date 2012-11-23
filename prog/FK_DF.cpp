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
/*    
    GF D1(n2);
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> D2(std::make_tuple(n1,n1));

    MatrixType<ComplexType> d3(3,3);
    MatrixType<ComplexType>::InnerIterator d3_it2(d3,0);
    //VectorType<ComplexType> d3(3);
    auto d3_it = index_begin<ComplexType>(d3);
    auto d3_it_end = index_end<ComplexType>(d3);
    //std::for_each(d3_it,d3_it_end,[](ComplexType &x){std::cout << x << " " << std::endl;});
  */  
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
