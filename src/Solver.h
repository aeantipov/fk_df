#ifndef ___FK_SOLVER_H___
#define ___FK_SOLVER_H___

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"

namespace FK { 

struct FKImpuritySolver
{
private:
    ComplexType _v_mult;
public:
    typedef GridObject<ComplexType,FMatsubaraGrid> GFType;
    const RealType U;
    const RealType mu;
    const RealType e_d;
    const FMatsubaraGrid w_grid;
    const RealType beta;
    GFType& Delta;
    GFType gw; 
    GFType K0;
    GFType K1;
    GFType Sigma;
    ComplexType w_0;
    ComplexType w_1;

    friend struct ImpurityVertex4;

    FKImpuritySolver(RealType U, RealType mu, RealType e_d, GFType& Delta);
    void run();
    ComplexType getVertex4(ComplexType w1, ComplexType w2) const {return _v_mult*K0(w1)*K0(w2)*K1(w1)*K1(w2);};
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
