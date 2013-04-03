#include "SelfConsistency.h"
#include "Diagrams.h"
#include <Eigen/LU>

namespace FK { 

//
// SelfConsistency
//

GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> SelfConsistency::getBubblePI() const
{
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> out(std::forward_as_tuple(this->_S.w_grid, this->_S.w_grid));
    RealType T = 1.0/_S.w_grid._beta;
    GFType iwn(this->_S.w_grid); iwn.fill([](ComplexType w){return w;});
    decltype(out)::PointFunctionType f = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2) { 
        return -T*(_S.gw(w1)+_S.gw(w2))/(ComplexType(w1)+ComplexType(w2)+2.0*_S.mu-_S.Sigma(w1)-_S.Sigma(w2));
    };
    out.fill(f);
    return out;
}

//
// Bethe SC
//

BetheSC::BetheSC(const FKImpuritySolver &S, RealType t):SelfConsistency(S),_t(t)
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

template <size_t D>
inline CubicDMFTSC<D>::CubicDMFTSC ( const FKImpuritySolver &S, RealType t, KMesh kGrid):
    SelfConsistency(S),
    _t(t),
    _kGrid(kGrid),
    _ek(__repeater<KMesh,D>::get_tuple(_kGrid)),
    _gloc(this->_S.w_grid)
{
    //CubicTraits<D>::template fill<index_iterator<ComplexType,EkStorage>>(index_begin<ComplexType, EkStorage>(_ek_vals), _t, _kGrid);
    CubicTraits<D>::template fillContainer<Container<D,ComplexType>>(_ek.getData(), _t, _kGrid);
    //CubicTraits<D>::template fillContainer<EkStorage>(_ek, _t, _kGrid);
    _ek._f = CubicTraits<D>::template get_dispersion<typename EkStorage::FunctionType> (t); 
}

template <size_t D>
typename CubicDMFTSC<D>::GKType CubicDMFTSC<D>::getGLat(const FMatsubaraGrid& fGrid) const
{
    std::array<KMesh,D> kgrids;
    kgrids.fill(_kGrid);
    GKType out(std::tuple_cat(std::forward_as_tuple(fGrid),kgrids));
    //for (auto iw : fGrid.getPoints()) {
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
        //return (_S.mu - _S.w_1*_S.U-_ek(ktuple))/std::abs(w*w)+1.0/w;
        };
    out._f = __fun_traits<typename GKType::FunctionType>::getFromTupleF(glatdmft_f);

    return out;
}

template <size_t D>
inline typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::operator()()
{
    INFO("Using DMFT self-consistency on a cubic lattice in " << D << " dimensions on a lattice of " << _kGrid.getSize() << "^" << D << " atoms.");
    GFType out(this->_S.w_grid); 
    out=0.0; 
    size_t ksize = _kGrid.getSize();
    RealType knorm = pow(ksize,D);
    for (auto w : _gloc.getGrid().getPoints()) {
        EkStorage e1 = (1.0/(1.0/_S.gw(w)+_S.Delta(w)-_ek)); 
        _gloc.get(w) = e1.sum()/knorm;
        out.get(w) = -1.0/_gloc(w)+_S.mu-_S.Sigma(w)+ComplexType(w);
    }
    //out._f = std::bind([&](ComplexType w){return _t*_t*2*RealType(D)/w;}, std::placeholders::_1);
    out._f = std::bind([&](ComplexType w)->ComplexType{return _t*_t*2*RealType(D)*((_S.mu-_S.w_1*_S.U)/std::abs(w*w) + 1.0/w);}, std::placeholders::_1);
    return out;
}


// Most important line for compilation.
template struct CubicDMFTSC<1>;
template struct CubicDMFTSC<2>;
template struct CubicDMFTSC<3>;
template struct CubicDMFTSC<4>;

//
// CubicInfDMFTSC
//


CubicInfDMFTSC :: CubicInfDMFTSC(const FKImpuritySolver &S, RealType t, const RealGrid& realgrid):
    SelfConsistency(S),
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
