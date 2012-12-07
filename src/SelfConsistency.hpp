#ifndef ___FK_SELFCONSISTENCY_HPP___
#define ___FK_SELFCONSISTENCY_HPP___

namespace FK {

template <size_t D>
inline CubicDMFTSC<D>::CubicDMFTSC ( RealType t, size_t npoints) : 
    _t(t),
    _npoints(npoints), 
    _kgrid(KMesh(npoints)) 
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
 
template <size_t D>
template <typename ...ArgTypes>
inline RealType CubicDMFTSC<D>::dispersion(const ArgTypes&... kpoints)
{
    static_assert(sizeof...(ArgTypes) == D, "Number of points mismatch!" );
    return this->ek(kpoints...);
}

template <size_t D>
template <typename ...ArgTypes>
inline ComplexType CubicDMFTSC<D>::glat(const CubicDMFTSC<D>::GFType &gw, const CubicDMFTSC<D>::GFType &Delta, const ArgTypes&... kpoints)
{
}

template <size_t D>
inline typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::operator()(const CubicDMFTSC<D>::GFType &gw)
{
    GFType Delta(gw.getGrid());
    Delta*=0.0;
    return Delta;
}

} // end of namespace FK
#endif // endif :: #ifndef ___FK_SELFCONSISTENCY_HPP___
