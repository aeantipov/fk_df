#include "DMFT.h"
#include "Diagrams.h"
#include <Eigen/LU>

namespace FK { 

//
// DMFT
//

GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> DMFTBase::getBubblePI() const
{
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> out(std::forward_as_tuple(this->_S.w_grid, this->_S.w_grid));
    RealType T = 1.0/_S.w_grid._beta;
    GFType iwn(this->_S.w_grid); iwn.fill(typename GFType::FunctionType([](ComplexType w){return w;}));
    decltype(out)::PointFunctionType f = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2) { 
        return -T*(_S.gw(w1)+_S.gw(w2))/(ComplexType(w1)+ComplexType(w2)+2.0*_S.mu-_S.Sigma(w1)-_S.Sigma(w2));
    };
    out.fill(f);
    return out;
}

//
// Bethe SC
//

BetheSC::BetheSC(const FKImpuritySolver &S, RealType t):DMFTBase(S),_t(t)
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


CubicInfDMFTSC :: CubicInfDMFTSC(const FKImpuritySolver &S, RealType t, const RealGrid& realgrid):
    DMFTBase(S),
    _t(t),
    _realgrid(realgrid),
    _nominator(ComplW(realgrid))
{
    std::function<ComplexType(RealType)> f1 = [=](RealType w){return std::exp(-1.0*w*w/t/t);};
    _nominator.fill(f1);
}



typename CubicInfDMFTSC::GFType CubicInfDMFTSC::operator()() 
{

    INFO("Using DMFT self-consistency on a hybercubic lattice in infinite dimensions");
    CubicInfDMFTSC::GFType gl(this->_S.w_grid); 
    CubicInfDMFTSC::GFType out(this->_S.w_grid); 
    ComplW denominator(_realgrid);
    for (auto iw : this->_S.w_grid.getPoints()) { 
        //std::function<ComplexType(RealGrid::point)> f2 = [&](RealGrid::point w){return 1.0/this->_S.gw(iw)+this->_S.Delta(iw)-RealType(w);};
        std::function<ComplexType(RealGrid::point)> f2 = [&](RealGrid::point w){return ComplexType(iw)+this->_S.mu-this->_S.Sigma(iw)-RealType(w);};
        denominator.fill(f2);
        auto tmp=_nominator/denominator;
        gl.get(iw) = _realgrid.integrate(_nominator/denominator)/_t/sqrt(PI);
        out.get(iw) = -1.0/gl(iw)+this->_S.mu-this->_S.Sigma(iw)+ComplexType(iw);
    }; 
    return out;
}

GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> CubicInfDMFTSC::getBubble0() const
{
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> out(std::forward_as_tuple(_S.w_grid,_S.w_grid));
    auto T = 1.0/_S.w_grid._beta;
    decltype(out)::PointFunctionType f = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2)->ComplexType {
        if (w1 == w2) return 2.0*T/_t/_t*(1.0-(ComplexType(w1)+_S.mu-_S.Sigma(w1))*_S.gw(w1));
        else return -T*(_S.gw(w1) - _S.gw(w2))/( ComplexType(w2) - ComplexType(w1) + _S.Sigma(w1) - _S.Sigma(w2));
    };
    out.fill(f);
    return out;
} 

} // end of namespace FK
