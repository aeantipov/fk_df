#include "DMFT.h"
#include "Diagrams.h"
#include <Eigen/LU>

namespace FK { 

//
// DMFT
//

grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> DMFTBase::getBubblePI() const
{
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> out(std::forward_as_tuple(this->_S.w_grid, this->_S.w_grid));
    real_type T = 1.0/_S.w_grid._beta;
    GFType iwn(this->_S.w_grid); iwn.fill(typename GFType::function_type([](complex_type w){return w;}));
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> ::point_function_type f = [&](fmatsubara_grid::point w1, fmatsubara_grid::point w2) { 
        return -T*(_S.gw(w1)+_S.gw(w2))/(complex_type(w1)+complex_type(w2)+2.0*_S.mu-_S.Sigma(w1)-_S.Sigma(w2));
    };
    out.fill(f);
    return out;
}

//
// Bethe SC
//

BetheSC::BetheSC(const FKImpuritySolver &S, real_type t):DMFTBase(S),_t(t)
{
}


typename BetheSC::GFType BetheSC::operator()()
{
    INFO("Using DMFT self-consistency on a Bethe lattice in infinite dimensions");
    return (this->_S.gw)*(_t*_t);
}

//
// CubicDMFT
//


// Most important line for compilation.
//template struct LatticeDMFTSC<1>;
//template struct LatticeDMFTSC<2>;
//template struct LatticeDMFTSC<3>;
//template struct LatticeDMFTSC<4>;

//
// CubicInfDMFTSC
//


CubicInfDMFTSC :: CubicInfDMFTSC(const FKImpuritySolver &S, real_type t, const real_grid& realgrid):
    DMFTBase(S),
    _t(t),
    _realgrid(realgrid),
    _nominator(ComplW(realgrid))
{
    std::function<complex_type(real_type)> f1 = [=](real_type w){return std::exp(-1.0*w*w/t/t);};
    _nominator.fill(f1);
}



typename CubicInfDMFTSC::GFType CubicInfDMFTSC::operator()() 
{

    INFO("Using DMFT self-consistency on a hybercubic lattice in infinite dimensions");
    CubicInfDMFTSC::GFType gl(this->_S.w_grid); 
    CubicInfDMFTSC::GFType out(this->_S.w_grid); 
    ComplW denominator(_realgrid);
    for (auto iw : this->_S.w_grid.points()) { 
        //std::function<complex_type(real_grid::point)> f2 = [&](real_grid::point w){return 1.0/this->_S.gw(iw)+this->_S.Delta(iw)-real_type(w);};
        std::function<complex_type(real_grid::point)> f2 = [&](real_grid::point w){return complex_type(iw)+this->_S.mu-this->_S.Sigma(iw)-real_type(w);};
        denominator.fill(f2);
        auto tmp=_nominator/denominator;
        gl.get(iw) = _realgrid.integrate(_nominator/denominator)/_t/sqrt(PI);
        out.get(iw) = -1.0/gl(iw)+this->_S.mu-this->_S.Sigma(iw)+complex_type(iw);
    }; 
    return out;
}

grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> CubicInfDMFTSC::getBubble0() const
{
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> out(std::forward_as_tuple(_S.w_grid,_S.w_grid));
    auto T = 1.0/_S.w_grid._beta;
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid>::point_function_type f = [&](fmatsubara_grid::point w1, fmatsubara_grid::point w2)->complex_type {
        if (w1 == w2) return 2.0*T/_t/_t*(1.0-(complex_type(w1)+_S.mu-_S.Sigma(w1))*_S.gw(w1));
        else return -T*(_S.gw(w1) - _S.gw(w2))/( complex_type(w2) - complex_type(w1) + _S.Sigma(w1) - _S.Sigma(w2));
    };
    out.fill(f);
    return out;
} 

} // end of namespace FK
