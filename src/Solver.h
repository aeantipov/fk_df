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
        inline ComplexType getFVertex4(Arg1 w1, Arg2 w2) const 
        { return _v_mult*K0(w1)*K0(w2)*K1(w1)*K1(w2)/(gw(w1)*gw(w1)*gw(w2)*gw(w2)); };
    template <typename Arg1, typename Arg2, 
              typename std::enable_if<std::is_convertible<Arg1, ComplexType>::value, int>::type=0, 
              typename std::enable_if<std::is_convertible<Arg2, ComplexType>::value, int>::type=0> 
        inline ComplexType getBVertex4(Arg1 O1, Arg2 w1) const { 
        auto w2 = w_grid.shift(w1,O1);
        return -this->getFVertex4(w1,w2);
        };
     template <typename Arg1, typename Arg2, typename Arg3, 
              typename std::enable_if<std::is_convertible<Arg1, ComplexType>::value, int>::type=0, 
              typename std::enable_if<std::is_convertible<Arg2, ComplexType>::value, int>::type=0, 
              typename std::enable_if<std::is_convertible<Arg3, ComplexType>::value, int>::type=0> 
        inline ComplexType getVertex4(Arg1 O1, Arg2 w1, Arg3 w2) const { 
            return ((std::abs(ComplexType(O1))<std::numeric_limits<RealType>::epsilon())?getFVertex4(w1,w2):0.0) + ((w1==w2)?getBVertex4(O1,w1):0.0);
        
        };
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> getBubble() const;
};

} // end of namespace FK
#endif // endif :: #ifndef ___FK_SOLVER_H___
