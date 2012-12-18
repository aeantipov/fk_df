#ifndef ___FK_SELFCONSISTENCY_HPP___
#define ___FK_SELFCONSISTENCY_HPP___

namespace FK {

//
// Bethe SC
//

template <class Solver>
inline BetheSC<Solver>::BetheSC(const Solver &S, RealType t):SelfConsistency<Solver>(S),_t(t)
{
}

template <class Solver>
inline typename BetheSC<Solver>::GFType BetheSC<Solver>::operator()() const
{
    return (this->_S.gw)*(_t*_t);
}

//
// CubicDMFT
//

template <class Solver, size_t D>
inline CubicDMFTSC<Solver, D>::CubicDMFTSC ( const Solver &S, RealType t, size_t npoints):
    SelfConsistency<Solver>(S),
    _t(t),
    _npoints(npoints), 
    _kgrid(KMesh(_ksize)),
    _gloc(this->_S.w_grid)
{
    __fill_ek<D,_ksize>::template fill<index_iterator<ComplexType,EkStorage>>(index_begin<ComplexType, EkStorage>(_ek_vals), _t, _kgrid);
}

template <size_t M, size_t ksize> 
template <class IteratorType, typename ...ArgTypes> 
inline void __fill_ek<M, ksize>::fill(IteratorType in, RealType t, const KMesh& grid, ArgTypes ... otherpos)
{
    size_t mult = __power<ksize,M-1>::value;
    auto move_it = in;
    for (size_t kx = 0; kx<ksize; ++kx) { 
        __fill_ek<M-1, ksize>::template fill<IteratorType, ArgTypes..., RealType> (move_it, t, grid, otherpos..., grid[kx]);
        move_it += mult;
        }
}

template <size_t ksize> 
template <class IteratorType, typename ...ArgTypes> 
inline void __fill_ek<0, ksize>::fill (IteratorType in, RealType t, const KMesh& grid, ArgTypes ... otherpos)
{
    *in = ek(t,otherpos...);
}

template <size_t ksize> 
template <typename ArgType1, typename ...ArgTypes> 
inline RealType __fill_ek<0, ksize>::ek(RealType t, ArgType1 kpoint1, ArgTypes... kpoints) 
{
    static_assert(std::is_convertible<ArgType1, RealType>::value,"Wrong kpoint");
    assert (kpoint1>=0 && kpoint1 < 2*PI);
    return -2.0*t*cos(kpoint1)+ek(t, kpoints...);
}
 
template <size_t ksize> 
template <typename ArgType1> 
inline RealType __fill_ek<0, ksize>::ek(RealType t, ArgType1 kpoint1)
{
    static_assert(std::is_convertible<ArgType1, RealType>::value,"Wrong kpoint");
    assert (kpoint1>=0 && kpoint1 < 2*PI);
    return -2*t*cos(kpoint1);
}
 
template <class Solver, size_t D>
template <typename ...ArgTypes>
inline RealType CubicDMFTSC<Solver, D>::dispersion(const ArgTypes&... kpoints)
{
    static_assert(sizeof...(ArgTypes) == D, "Number of points mismatch!" );
    return __fill_ek<0,_ksize>::ek(_t, kpoints...);
}

template <class Solver, size_t D>
template <typename ...ArgTypes>
inline typename CubicDMFTSC<Solver, D>::GFType CubicDMFTSC<Solver,D>::glat(ArgTypes... kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D,"!");
    GFType out = 1.0/(1.0/this->_S.gw+this->_S.Delta-__fill_ek<0,_ksize>::ek(_t, kpoints...));
    return out;

}

template <class Solver, size_t D>
inline typename CubicDMFTSC<Solver,D>::GFType CubicDMFTSC<Solver,D>::operator()()
{
    GFType out(this->_S.w_grid); 
    out=0.0; 
    for (auto w : _gloc.getGrid().getVals()) {
        EkStorage e1 = (1.0/(1.0/this->_S.gw(w)+this->_S.Delta(w)-_ek_vals)); 
        _gloc.get(w) = e1.sum()/RealType(_ksize*_ksize);
        out.get(w) = -1.0/_gloc(w)+this->_S.mu-this->_S.Sigma(w)+ComplexType(w);
    }
    return out;
}

//
// CubicInfDMFTSC
//

template <class Solver>
inline CubicInfDMFTSC<Solver> :: CubicInfDMFTSC(const Solver &S, RealType t, const RealGrid& realgrid):
    SelfConsistency<Solver>(S),
    _t(t),
    _realgrid(realgrid),
    _nominator(ComplW(realgrid))
{
    std::function<ComplexType(RealType)> f1 = [=](RealType w){return std::exp(-1.0*w*w/t/t);};
    _nominator = f1;
}


template <class Solver>
inline typename CubicInfDMFTSC<Solver>::GFType CubicInfDMFTSC<Solver>::operator()() const 
{
    CubicInfDMFTSC::GFType gl(this->_S.w_grid); 
    CubicInfDMFTSC::GFType out(this->_S.w_grid); 
    ComplW denominator(_realgrid);
    for (auto iw : this->_S.w_grid.getVals()) { 
        //std::function<ComplexType(RealGrid::point)> f2 = [&](RealGrid::point w){return 1.0/this->_S.gw(iw)+this->_S.Delta(iw)-RealType(w);};
        std::function<ComplexType(RealGrid::point)> f2 = [&](RealGrid::point w){return ComplexType(iw)+this->_S.mu-this->_S.Sigma(iw)-RealType(w);};
        denominator.fill(f2);
        auto tmp=_nominator/denominator;
        gl.get(iw) = _realgrid.integrate(_nominator/denominator)/_t/sqrt(PI);
        out.get(iw) = -1.0/gl(iw)+this->_S.mu-this->_S.Sigma(iw)+ComplexType(iw);
    }; 
    return out;
}



} // end of namespace FK
#endif // endif :: #ifndef ___FK_SELFCONSISTENCY_HPP___
