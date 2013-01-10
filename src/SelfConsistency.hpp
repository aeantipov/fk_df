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
inline typename BetheSC<Solver>::GFType BetheSC<Solver>::operator()()
{
    INFO("Using DMFT self-consistency on a Bethe lattice in infinite dimensions");
    return (this->_S.gw)*(_t*_t);
}

//
// CubicTraits
//

template <size_t M, size_t ksize> 
template <class IteratorType, typename ...ArgTypes> 
inline void CubicTraits<M, ksize>::fill(IteratorType in, RealType t, const KMesh& grid, ArgTypes ... otherpos)
{
    assert(grid.getSize() == ksize);
    size_t mult = __power<ksize,M-1>::value;
    auto move_it = in;
    for (size_t kx = 0; kx<ksize; ++kx) { 
        CubicTraits<M-1, ksize>::template fill<IteratorType, ArgTypes..., RealType> (move_it, t, grid, otherpos..., grid[kx]);
        move_it += mult;
        }
}

template <size_t M, size_t ksize> 
template <class ContainerType, typename ...ArgTypes> 
inline void CubicTraits<M, ksize>::fillContainer(ContainerType &in, RealType t, const KMesh& grid, ArgTypes ... otherpos)
{
    assert(grid.getSize() == ksize);
    for (size_t kx = 0; kx<ksize; ++kx) { 
        CubicTraits<M-1, ksize>::fillContainer(in[kx], t, grid, otherpos..., grid[kx]);
        }
}



template <size_t ksize> 
template <class IteratorType, typename ...ArgTypes> 
inline void CubicTraits<0, ksize>::fill (IteratorType in, RealType t, const KMesh& grid, ArgTypes ... otherpos)
{
    *in = ek(t,otherpos...);
}

template <size_t ksize> 
template <class DataType, typename ...ArgTypes> 
inline void CubicTraits<0, ksize>::fillContainer(DataType &in, RealType t, const KMesh& grid, ArgTypes ... otherpos)
{
    assert(grid.getSize() == ksize);
    in = ek(t,otherpos...);
}


template <size_t ksize> 
template <typename ArgType1, typename ...ArgTypes> 
inline RealType CubicTraits<0, ksize>::ek(RealType t, ArgType1 kpoint1, ArgTypes... kpoints) 
{
    static_assert(std::is_convertible<ArgType1, RealType>::value,"Wrong kpoint");
    assert (kpoint1>=0 && kpoint1 < 2*PI);
    return -2.0*t*cos(kpoint1)+ek(t, kpoints...);
}
 
template <size_t ksize> 
template <typename ArgType1> 
inline RealType CubicTraits<0, ksize>::ek(RealType t, ArgType1 kpoint1)
{
    static_assert(std::is_convertible<ArgType1, RealType>::value,"Wrong kpoint");
    assert (kpoint1>=0 && kpoint1 < 2*PI);
    return -2*t*cos(kpoint1);
}

//
// CubicDMFT
//

template <class Solver, size_t D, size_t ksize>
inline CubicDMFTSC<Solver,D,ksize>::CubicDMFTSC ( const Solver &S, RealType t):
    SelfConsistency<Solver>(S),
    _t(t),
    _kGrid(KMesh(ksize)),
    _ek(CubicTraits<D,ksize>::getTuples(_kGrid)),
    _gloc(this->_S.w_grid)
{
    //CubicTraits<D,ksize>::template fill<index_iterator<ComplexType,EkStorage>>(index_begin<ComplexType, EkStorage>(_ek_vals), _t, _kGrid);
    CubicTraits<D,ksize>::template fillContainer<Container<D,ComplexType>>(_ek.getData(), _t, _kGrid);
    //CubicTraits<D,ksize>::template fillContainer<EkStorage>(_ek, _t, _kGrid);
    _ek._f = CubicTraits<D,ksize>::template get_dispersion<typename EkStorage::FunctionType> (t); 
}


template <class Solver, size_t D, size_t ksize>
template <typename ...ArgTypes>
inline RealType CubicDMFTSC<Solver,D,ksize>::dispersion(ArgTypes... kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D, "Number of points mismatch!" );
    return CubicTraits<0,ksize>::ek(_t, kpoints...);
}

template <class Solver, size_t D, size_t ksize>
template <typename ...ArgTypes>
inline RealType CubicDMFTSC<Solver,D,ksize>::dispersion(const std::tuple<ArgTypes...>& kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D, "Number of points mismatch!" );
    typename EkStorage::PointFunctionType f1 = [&](ArgTypes... kpoints)->RealType{return dispersion(kpoints...);};
    auto f2 =__fun_traits<typename EkStorage::PointFunctionType>::getTupleF(f1); 
    return std::real(f2(kpoints));
}


template <class Solver, size_t D, size_t ksize>
template <typename ...ArgTypes>
inline typename CubicDMFTSC<Solver,D,ksize>::GFType CubicDMFTSC<Solver,D,ksize>::glat(ArgTypes... kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D,"!");
    auto e = dispersion<ArgTypes...>(kpoints...);
    GFType out = 1.0/(1.0/_S.gw+_S.Delta-e);
    return out;

}

template <class Solver, size_t D, size_t ksize>
template <typename MPoint, typename ...ArgTypes> 
ComplexType CubicDMFTSC<Solver,D,ksize>::glat_val(MPoint w, ArgTypes... kpoints) const
{
    return 1.0/(1.0/_S.gw(w)+_S.Delta(w)-dispersion(kpoints...));
}

template <class Solver, size_t D, size_t ksize>
typename CubicDMFTSC<Solver,D,ksize>::GKType CubicDMFTSC<Solver,D,ksize>::getGLat(const FMatsubaraGrid& fGrid) const
{
    std::array<KMesh,D> kgrids;
    kgrids.fill(_kGrid);
    GKType out(std::tuple_cat(std::forward_as_tuple(fGrid),kgrids));
    //for (auto iw : fGrid.getVals()) {
    auto f1 = [&](const typename GKType::PointTupleType& in)->ComplexType {
        FMatsubaraGrid::point w = std::get<0>(in);
        auto ktuple = __tuple_tail(in);
        typename ArgBackGenerator<D,KMesh::point,__caller,RealType>::type t1; //( {{in, &dispersion}});
        return 1.0/(1.0/_S.gw(w)+_S.Delta(w)-_ek(ktuple));
    };
    typename GKType::PointFunctionType f2 = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(f1);
    out.fill(f2);

    auto glatdmft_f = [&](const typename GKType::ArgTupleType &in)->ComplexType{
        ComplexType w = std::get<0>(in);
        auto ktuple = __tuple_tail(in);
        return (_S.mu - _S.Sigma._f(w)-_ek(ktuple))/std::abs(w*w)+1.0/w;
        };
    out._f = __fun_traits<typename GKType::FunctionType>::getFromTupleF(glatdmft_f);

    return out;
}

template <class Solver, size_t D, size_t ksize>
inline typename CubicDMFTSC<Solver,D,ksize>::GFType CubicDMFTSC<Solver,D,ksize>::operator()()
{
    INFO("Using DMFT self-consistency on a cubic lattice in " << D << " dimensions on a lattice of " << ksize << "^" << D << " atoms.");
    GFType out(this->_S.w_grid); 
    out=0.0; 
    for (auto w : _gloc.getGrid().getVals()) {
        EkStorage e1 = (1.0/(1.0/_S.gw(w)+_S.Delta(w)-_ek)); 
        _gloc.get(w) = e1.sum()/RealType(__power<ksize,D>::value);
        out.get(w) = -1.0/_gloc(w)+_S.mu-_S.Sigma(w)+ComplexType(w);
    }
    //out._f = std::bind([&](ComplexType w){return _t*_t*2*RealType(D)/w;}, std::placeholders::_1);
    out._f = std::bind([&](ComplexType w)->ComplexType{return _t*_t*2*RealType(D)*((_S.mu-_S.w_1*_S.U)/std::abs(w*w) + 1.0/w);}, std::placeholders::_1);
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
    _nominator.fill(f1);
}


template <class Solver>
inline typename CubicInfDMFTSC<Solver>::GFType CubicInfDMFTSC<Solver>::operator()() 
{

    INFO("Using DMFT self-consistency on a hybercubic lattice in infinite dimensions");
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
