#ifndef ___FK_SOLVER_H___
#define ___FK_SOLVER_H___

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"
#include "GFWrap.h"

namespace FK { 

struct FKImpuritySolver
{
private:
    ComplexType _v_mult;
public:
    //typedef GridObject<ComplexType,FMatsubaraGrid> GFType;
    typedef GFWrap GFType;
    const RealType U;
    const RealType mu;
    const RealType e_d;
    const FMatsubaraGrid w_grid;
    const FMatsubaraGrid half_grid;
    const RealType beta;
    GFType Delta;
    GFType gw; 
    GFType K0;
    GFType K1;
    GFType Sigma;
    RealType w_0;
    RealType w_1;

    friend struct ImpurityVertex4;

    FKImpuritySolver(RealType U, RealType mu, RealType e_d, GFType& Delta);
    void run(bool calc_weight = true);
    //ComplexType getVertex4(ComplexType w1, ComplexType w2) const {return _v_mult*K0(w1)*K0(w2)*K1(w1)*K1(w2);};
    template <typename Arg1, typename Arg2, 
              typename std::enable_if<std::is_convertible<Arg1, ComplexType>::value, int>::type=0, 
              typename std::enable_if<std::is_convertible<Arg2, ComplexType>::value, int>::type=0> 
        ComplexType getVertex4(Arg1, Arg2) const;
    //ComplexType getVertex4<FMatsubaraGrid::point, FMatsubaraGrid::point>(FMatsubaraGrid::point, FMatsubaraGrid::point) const;
 //{return _v_mult*K0(wF)*K0(ComplexType(wF)+ComplexType(WB))*K1(wF)*K1(ComplexType(wF)+ComplexType(WB));};
};

struct ImpurityVertex4 { 
private:
    FKImpuritySolver& _parent;
    const ComplexType mult = _parent.U*_parent.U*_parent.w_0*_parent.w_1*_parent.beta;
public:
    ImpurityVertex4(FKImpuritySolver& in):_parent(in){};
    //decltype(_parent.gw[0]) operator()(ComplexType w1, ComplexType w2) const
    ComplexType operator()(ComplexType w1, ComplexType w2) const /// -> ComplexType // decltype(_parent.gw[0])
    {
        return mult*_parent.K0(w1)*_parent.K0(w2)*_parent.K1(w1)*_parent.K1(w2);
    }
};

} // end of namespace FK
#endif // endif :: #ifndef ___FK_SOLVER_H___
