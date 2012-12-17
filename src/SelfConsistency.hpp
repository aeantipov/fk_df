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
    _kgrid(KMesh(npoints))
{
//    fill_ek<D-1>::fill(_ek_vals,0);
}

/*
template <class Solver, size_t D>
template <size_t M, typename ...ArgTypes> 
inline void CubicDMFTSC<Solver, D>::fill_ek<M, ArgTypes...>::fill(typename CubicDMFTSC<Solver, D>::EkStorage&in, size_t pos, ArgTypes ... otherpos)
{
    size_t mult = __power<_ksize,M>::value;
    for (size_t kx = 0; kx<_ksize; ++kx) fill_ek<M-1, ArgTypes...>::fill(in, pos+kx*mult, otherpos..., kx);
}


template <class Solver, size_t D>
template <typename ...ArgTypes> 
inline void CubicDMFTSC<Solver, D>::fill_ek<1, ArgTypes...>::fill(typename CubicDMFTSC<Solver, D>::EkStorage&in, size_t pos, ArgTypes ... otherpos)
{
    //size_t mult = __power<_ksize,M>::value;
    //for (size_t kx = 0; kx<_ksize; ++kx) fill_ek<D-1>(pos+kx*mult, otherpos..., _kgrid[kx])
}
*/


template <class Solver, size_t D>
template <typename ArgType1, typename ...ArgTypes> 
inline RealType CubicDMFTSC<Solver, D>::ek(RealType t, ArgType1 kpoint1, ArgTypes... kpoints) 
{
    static_assert(std::is_convertible<ArgType1, RealType>::value,"Wrong kpoint");
    assert (kpoint1>=0 && kpoint1 < 2*PI);
    return -2.0*t*cos(kpoint1)+ek(t, kpoints...);
}
 
template <class Solver, size_t D>
template <typename ArgType1> 
inline RealType CubicDMFTSC<Solver, D>::ek(RealType t, ArgType1 kpoint1)
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
    return ek(_t, kpoints...);
}

template <class Solver, size_t D>
template <typename ...ArgTypes>
inline typename CubicDMFTSC<Solver, D>::GFType CubicDMFTSC<Solver,D>::glat(const GFType& gw, const GFType& Delta, ArgTypes... kpoints) const
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

template <class Solver, size_t D>
template <typename ...ArgTypes>
inline ComplexType CubicDMFTSC<Solver,D>::glat_val(const ComplexType &gw, const ComplexType &Delta, ArgTypes... kpoints) const
{
    return (1.0/(1.0/gw+Delta-ek(_t, kpoints...)));
}

/* Tried to do it better. In the end here it is acceptible. */
//1d
template <>
inline typename CubicDMFTSC<FKImpuritySolver,1>::GFType CubicDMFTSC<FKImpuritySolver,1>::operator()() const
{
    //assert(Delta.getGrid().getVals() == gw.getGrid().getVals());
    RecursiveGridIntegrator<KMesh, ArgFunType> t;
    ArgFunType f = [&](RealType k1){return glat(this->_S.gw, this->_S.Delta, k1);};
    assert(1);
    return t.integrate(_kgrid, f);
}


//2d
template <>
inline typename CubicDMFTSC<FKImpuritySolver,2>::GFType CubicDMFTSC<FKImpuritySolver, 2>::operator()() const
{
//    assert(Delta.getGrid().getVals() == gw.getGrid().getVals());


    /*RecursiveGridIntegrator<KMesh, ArgFunType> t;
    ArgFunType f = [&](RealType k1, RealType k2){return glat(this->_S.gw, this->_S.Delta, k1,k2);};
    return t.integrate(_kgrid, f);
*/

    GFType out(this->_S.w_grid); 
    out=0.0; 
    // No lapack
/*    size_t ksize = _kgrid._points;
    GridObject<ComplexType, KMesh, KMesh> e_k(std::make_tuple(_kgrid,_kgrid));
    for (auto w : out.getGrid().getVals()) {
        std::function<ComplexType(RealType,RealType)> f1 = [&](RealType k1, RealType k2){return glat_val(this->_S.gw(w), this->_S.Delta(w), k1, k2);};
        e_k.fill(f1);
        out.get(w) = e_k.sum()/RealType(ksize*ksize);
        out.get(w) = -1.0/out(w)+this->_S.mu-this->_S.Sigma(w)+ComplexType(w);
    }
  */  

    // Eigen
    //size_t ksize = _kgrid._points;
    Eigen::Array<ComplexType,_ksize,_ksize,Eigen::AutoAlign> ek_matrix;
    Eigen::Array<RealType,_ksize,1> v;
    v.setLinSpaced(0,2*PI);
    for (auto w : out.getGrid().getVals()) {
        std::function<ComplexType(RealType,RealType)> f1 = [&](RealType k1, RealType k2){return glat_val(this->_S.gw(w), this->_S.Delta(w), k1, k2);};
        //out.get(w) = e_k.sum()/RealType(ksize*ksize);
        //out.get(w) = -1.0/out(w)+this->_S.mu-this->_S.Sigma(w)+ComplexType(w);
    }

    return out;
/*
    GFType out(this->_S.Delta); out=0.0; 
    for ( RealType kx : _kgrid.getVals() ) 
        for ( RealType ky : _kgrid.getVals() ){
            GFType temp(this->_S.Delta.getGrid());
            temp=1.0;
            temp/=gw;
            temp/=(temp+this->_S.Delta-ek(_t,kx,ky));
            temp*=gw;
            out+=temp/_kgrid._points/_kgrid._points;
                };
    return out;
*/
}

//3d
template <>
inline typename CubicDMFTSC<FKImpuritySolver,3>::GFType CubicDMFTSC<FKImpuritySolver,3>::operator()() const
{
    RecursiveGridIntegrator<KMesh, ArgFunType> t;
    ArgFunType f = [&](RealType k1, RealType k2, RealType k3){return glat(this->_S.gw, this->_S.Delta, k1, k2, k3);};
    assert(1);
    return t.integrate(_kgrid, f);
}

//4d
template <>
inline typename CubicDMFTSC<FKImpuritySolver,4>::GFType CubicDMFTSC<FKImpuritySolver,4>::operator()() const
{
    //assert(Delta.getGrid().getVals() == gw.getGrid().getVals());

    RecursiveGridIntegrator<KMesh, ArgFunType> t;
    ArgFunType f = [&](RealType k1, RealType k2, RealType k3, RealType k4){return glat(this->_S.gw, this->_S.Delta, k1, k2, k3, k4);};
    assert(1);
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
