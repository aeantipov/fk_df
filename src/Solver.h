#ifndef ___FK_SOLVER_H___
#define ___FK_SOLVER_H___

#include "gftools/matsubara_grid.hpp"
#include "gftools/container.hpp"
#include "gftools/grid_object.hpp"

namespace FK { 

using namespace gftools;

struct FKImpuritySolver
{
private:
    complex_type _v_mult;
public:
    typedef grid_object<complex_type,fmatsubara_grid> GFType;
    const real_type U;
    const real_type mu;
    const real_type e_d;
    const fmatsubara_grid w_grid;
    const fmatsubara_grid half_grid;
    const real_type beta;
    GFType Delta;
    GFType gw; 
    GFType K0;
    GFType K1;
    GFType Lambda;
    GFType Sigma;
    real_type w_0;
    real_type w_1;

    friend struct ImpurityVertex4;

    FKImpuritySolver(real_type U, real_type mu, real_type e_d, GFType& Delta);
    void run(bool calc_weight = true);
    //complex_type getVertex4(complex_type w1, complex_type w2) const {return _v_mult*K0(w1)*K0(w2)*K1(w1)*K1(w2);};
    template <typename Arg1, typename Arg2, 
              typename std::enable_if<std::is_convertible<Arg1, complex_type>::value, int>::type=0,
              typename std::enable_if<std::is_convertible<Arg2, complex_type>::value, int>::type=0>
        inline complex_type getFVertex4(Arg1 w1, Arg2 w2) const
        { return _v_mult*K0(w1)*K0(w2)*K1(w1)*K1(w2)/(gw(w1)*gw(w1)*gw(w2)*gw(w2)); };
    template <typename Arg1, typename Arg2, 
              typename std::enable_if<std::is_convertible<Arg1, complex_type>::value, int>::type=0,
              typename std::enable_if<std::is_convertible<Arg2, complex_type>::value, int>::type=0>
        inline complex_type getBVertex4(Arg1 O1, Arg2 w1) const {
        try { 
            auto w2 = w_grid.shift(w1,complex_type(O1));
            return -this->getFVertex4(w1,w2);
            }
        catch (typename gftools::fmatsubara_grid::ex_wrong_index)
            {
            complex_type w2 = w_grid.shift(complex_type(w1),complex_type(O1));
            return -this->getFVertex4(w1,w2);
            }
        };
     template <typename Arg1, typename Arg2, typename Arg3, 
              typename std::enable_if<std::is_convertible<Arg1, complex_type>::value, int>::type=0,
              typename std::enable_if<std::is_convertible<Arg2, complex_type>::value, int>::type=0,
              typename std::enable_if<std::is_convertible<Arg3, complex_type>::value, int>::type=0>
        inline complex_type getVertex4(Arg1 O1, Arg2 w1, Arg3 w2) const {
            return ((std::abs(complex_type(O1))<std::numeric_limits<real_type>::epsilon())?getFVertex4(w1,w2):0.0) + ((w1==w2)?getBVertex4(O1,w1):0.0);
        
        };
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> getBubble() const;

    GFType getLambda() const;
};

} // end of namespace FK
#endif // endif :: #ifndef ___FK_SOLVER_H___
