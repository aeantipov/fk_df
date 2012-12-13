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
inline RealType CubicDMFTSC<D>::ek(RealType t, ArgType1 kpoint1, ArgTypes... kpoints) 
{
    static_assert(std::is_convertible<ArgType1, RealType>::value,"Wrong kpoint");
    assert (kpoint1>=0 && kpoint1 < 2*PI);
    return -2.0*t*cos(kpoint1)+ek(t, kpoints...);
}
 
template <size_t D>
template <typename ArgType1> 
inline RealType CubicDMFTSC<D>::ek(RealType t, ArgType1 kpoint1)
{
    static_assert(std::is_convertible<ArgType1, RealType>::value,"Wrong kpoint");
    assert (kpoint1>=0 && kpoint1 < 2*PI);
    return -2*t*cos(kpoint1);
}
 
template <size_t D>
template <typename ...ArgTypes>
inline RealType CubicDMFTSC<D>::dispersion(const ArgTypes&... kpoints)
{
    static_assert(sizeof...(ArgTypes) == D, "Number of points mismatch!" );
    return ek(_t, kpoints...);
}

template <size_t D>
template <typename ...ArgTypes>
inline typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::glat(const GFType& gw, const GFType& Delta, ArgTypes... kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D,"!");
    assert(Delta.getGrid().getVals() == gw.getGrid().getVals());
    GFType out = 1.0/(1.0/gw+Delta-ek(_t, kpoints...));
    //DEBUG(out);
    //exit(0);
    //GFType out = 1.0/(1.0/gw+Delta-ek(_t, kpoints...)); 
    //DEBUG("!");
    //std::function<ComplexType(FMatsubaraGrid::point)> f1 = [&](typename FMatsubaraGrid::point w){return glat_val(gw(w),Delta(w),kpoints...);};
    //std::function<ComplexType(FMatsubaraGrid::point)> f1 = [&](typename FMatsubaraGrid::point w){return 1.0/gw(w)+Delta(w)-1.0;};
    //out.fill(f1);
    //DEBUG(out);
    return out;

}

template <size_t D>
template <typename ...ArgTypes>
inline ComplexType CubicDMFTSC<D>::glat_val(const ComplexType &gw, const ComplexType &Delta, ArgTypes... kpoints) const
{
    return (1.0/(1.0/gw+Delta-ek(_t, kpoints...)));
}

/* Tried to do it better. In the end here it is acceptible. */
//1d
template <>
inline typename CubicDMFTSC<1>::GFType CubicDMFTSC<1>::operator()(const GFType& gw, const GFType& Delta) const
{
    assert(Delta.getGrid().getVals() == gw.getGrid().getVals());
    RecursiveGridIntegrator<KMesh, ArgFunType> t;
    ArgFunType f = [&](RealType k1){return glat(gw, Delta, k1);};
    return t.integrate(_kgrid, f);
}


//2d
template <>
inline typename CubicDMFTSC<2>::GFType CubicDMFTSC<2>::operator()(const GFType& gw, const GFType& Delta) const
{
    assert(Delta.getGrid().getVals() == gw.getGrid().getVals());
/*
    RecursiveGridIntegrator<KMesh, ArgFunType> t;
    ArgFunType f = [&](RealType k1, RealType k2){return glat(gw, Delta, k1,k2);};
    return t.integrate(_kgrid, f);
*/

    size_t ksize = _kgrid._points;
    GridObject<ComplexType, KMesh, KMesh> e_k(std::make_tuple(_kgrid,_kgrid));
    GFType out(Delta); out=0.0; 
    for (auto w : out.getGrid().getVals()) {
        std::function<ComplexType(RealType,RealType)> f1 = [&](RealType k1, RealType k2){return glat_val(gw(w), Delta(w), k1, k2);};
        e_k.fill(f1);
        out.get(w) = e_k.sum()/RealType(ksize*ksize);
    }
    return out;

/*
    GFType out(Delta); out=0.0; 
    for ( RealType kx : _kgrid.getVals() ) 
        for ( RealType ky : _kgrid.getVals() ){
            GFType temp(Delta.getGrid());
            temp=1.0;
            temp/=gw;
            temp/=(temp+Delta-ek(_t,kx,ky));
            temp*=gw;
            out+=temp/_kgrid._points/_kgrid._points;
                };
    return out;
*/
}

//3d
template <>
inline typename CubicDMFTSC<3>::GFType CubicDMFTSC<3>::operator()(const GFType& gw, const GFType& Delta) const
{
    assert(Delta.getGrid().getVals() == gw.getGrid().getVals());
    RecursiveGridIntegrator<KMesh, ArgFunType> t;
    ArgFunType f = [&](RealType k1, RealType k2, RealType k3){return glat(gw, Delta, k1, k2, k3);};
    return t.integrate(_kgrid, f);
}

//4d
template <>
inline typename CubicDMFTSC<4>::GFType CubicDMFTSC<4>::operator()(const GFType& gw, const GFType& Delta) const
{
    assert(Delta.getGrid().getVals() == gw.getGrid().getVals());

    RecursiveGridIntegrator<KMesh, ArgFunType> t;
    ArgFunType f = [&](RealType k1, RealType k2, RealType k3, RealType k4){return glat(gw, Delta, k1, k2, k3, k4);};
    return t.integrate(_kgrid, f);

/*
    size_t ksize = _kgrid._points;
    GridObject<ComplexType, KMesh, KMesh, KMesh, KMesh> e_k(std::make_tuple(_kgrid,_kgrid,_kgrid,_kgrid));
    GFType out(Delta); out=0.0; 
    for (auto w : out.getGrid().getVals()) {
        std::function<ComplexType(RealType,RealType,RealType,RealType)> f1 = [&](RealType k1, RealType k2, RealType k3, RealType k4){return glat_val(gw(w), Delta(w), k1, k2, k3, k4);};
        e_k.fill(f1);
        out(w) = e_k.sum()/RealType(ksize*ksize*ksize*ksize);
    }
    return out;
*/
}



} // end of namespace FK
#endif // endif :: #ifndef ___FK_SELFCONSISTENCY_HPP___
