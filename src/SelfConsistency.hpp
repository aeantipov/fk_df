#ifndef ___FK_SELFCONSISTENCY_HPP___
#define ___FK_SELFCONSISTENCY_HPP___

namespace FK {

template <size_t D>
inline CubicDMFTSC<D>::CubicDMFTSC ( RealType t, size_t npoints):
//    SelfConsistency(gw),
    _t(t),
    _npoints(npoints), 
    _kgrid(KMesh(npoints))
//    _wgrid(gw.getGrid()),
//    _Delta_old(Delta_old)
{
}

template <size_t D>
template <typename ArgType1, typename ...ArgTypes> 
inline RealType CubicDMFTSC<D>::ek(const ArgType1& kpoint1, const ArgTypes&... kpoints) 
{
    static_assert(std::is_convertible<ArgType1, RealType>::value,"Wrong kpoint");
    assert (kpoint1>=0 && kpoint1 < 2*PI);
    return -2*_t*cos(kpoint1)+this->ek(kpoints...);
}
 
template <size_t D>
template <typename ArgType1> 
inline RealType CubicDMFTSC<D>::ek(const ArgType1& kpoint1)
{
    static_assert(std::is_convertible<ArgType1, RealType>::value,"Wrong kpoint");
    assert (kpoint1>=0 && kpoint1 < 2*PI);
    return -2*_t*cos(kpoint1);
}
 

/*
template <size_t D>
template <typename ArgType1, typename ...ArgTypes> 
inline auto CubicDMFTSC<D>::getSC(const ArgType1& kpoint1, const ArgTypes&... kpoints, const GFType& gw, const GFType& Delta) const 
    -> decltype(this->getSC(kpoints...,gw,Delta)) 
{
}

template <size_t D>
template <> 
inline typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::getSC<1>(const GFType& gw, const GFType &Delta)
{
}
*/

template <size_t D>
template <typename ...ArgTypes>
inline RealType CubicDMFTSC<D>::dispersion(const ArgTypes&... kpoints)
{
    static_assert(sizeof...(ArgTypes) == D, "Number of points mismatch!" );
    return this->ek(kpoints...);
}

template <size_t D>
template <typename ...ArgTypes>
inline typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::glat(const GFType& gw, const GFType& Delta, const ArgTypes&... kpoints)
{
    assert(Delta.getGrid().getVals() == gw.getGrid().getVals());
    GFType out(Delta.getGrid()); 
    RealType ek_ = this->ek(kpoints...);
    std::function<ComplexType(ComplexType)> f1 = [&gw,&Delta,&ek_](ComplexType w){return 1.0/(1.0/gw(w)+Delta(w)-ek_);};
    out.fill(f1);
    return out;
}

template <size_t D>
inline typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::operator()(const GFType& gw, const GFType& Delta)
{
    assert(Delta.getGrid().getVals() == gw.getGrid().getVals());
    //return k_grid.integrate(
    return Delta;
}

} // end of namespace FK
#endif // endif :: #ifndef ___FK_SELFCONSISTENCY_HPP___
