#ifndef ___FK_SELFCONSISTENCY_HPP___
#define ___FK_SELFCONSISTENCY_HPP___

#include "Diagrams.h"

#include <Eigen/LU>

namespace FK {

template <class Solver>
template <typename MPoint>
inline typename SelfConsistency<Solver>::GFType SelfConsistency<Solver>::getLatticeDMFTVertex4(MPoint in) const
{
    GFType Vertex4(_S.w_grid);
    typename GFType::PointFunctionType VertexFillf = 
        std::bind(&FKImpuritySolver::getBVertex4<MPoint, FMatsubaraGrid::point>, std::cref(_S), in, std::placeholders::_1); 
    typename GFType::FunctionType Vertexf = 
        std::bind(&FKImpuritySolver::getBVertex4<ComplexType, ComplexType>, std::cref(_S), ComplexType(in), std::placeholders::_1); 
    Vertex4.fill(VertexFillf);
    Vertex4._f = Vertexf;

    GFType Chi0 = Diagrams::getBubble(_S.gw,in);
    auto Vertex_out = Diagrams::BS(Chi0, Vertex4, false);
    return Vertex_out; 
}

template <class Solver>
inline GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> SelfConsistency<Solver>::getStaticLatticeDMFTVertex4() const
{
    INFO_NONEWLINE("\tObtaining vertex from Solver and static 2-freq susceptibility...");
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> Vertex4_out(std::forward_as_tuple(_S.w_grid,_S.w_grid)); 
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> Chi0(std::forward_as_tuple(_S.w_grid,_S.w_grid)); 
    decltype(Vertex4_out)::PointFunctionType VertexF = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return _S.getFVertex4(w1,w2);};
    decltype(Chi0)::PointFunctionType Chi0F = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return _S.gw(w1)*_S.gw(w2)/_S.w_grid._beta;};
    Vertex4_out.fill(VertexF);
    Chi0.fill(Chi0F);

    auto Vertex4 = Vertex4_out.getData().getAsMatrix();
    auto Chi0Matrix = Chi0.getData().getAsMatrix();

    INFO("done");

    Vertex4 = Diagrams::BS(Chi0Matrix, Vertex4, false);
    
    INFO_NONEWLINE("\tFilling output vertex...");
    Vertex4_out.getData() = Vertex4;
    INFO("done.");
    return Vertex4_out;
}

template <class Solver>
template <typename MPoint>
inline typename SelfConsistency<Solver>::GFType SelfConsistency<Solver>::getBubblePI(MPoint in) const
{
    GFType out(this->_S.w_grid);
    GFType gw_shift(_S.gw), Sigma_shift(_S.Sigma);
    RealType T = 1.0/_S.w_grid._beta;
    if (std::abs(ComplexType(in))<PI*T) {
         gw_shift = _S.gw.shift(in);
         Sigma_shift = _S.Sigma.shift(in);
        };
    GFType iwn(this->_S.w_grid); iwn.fill([](ComplexType w){return w;});
    out = -T*(_S.gw+gw_shift)/(2*iwn+ComplexType(in)+2.0*_S.mu-_S.Sigma-Sigma_shift);
    return out;
}

template <class Solver>
inline GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> SelfConsistency<Solver>::getBubblePI() const
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


/*
template <class Solver>
template <typename MPoint>
inline typename CubicInfDMFTSC<Solver>::GFType CubicInfDMFTSC<Solver>::getBubble0(MPoint in)
{
    GFType out(this->_S.w_grid);
    return out;
}
*/

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

template <size_t M> 
template <class IteratorType, typename ...ArgTypes> 
inline void CubicTraits<M>::fill(IteratorType in, RealType t, const KMesh& grid, ArgTypes ... otherpos)
{
    //assert(grid.getSize() == ksize);
    size_t ksize = grid.getSize();
    size_t mult = pow(ksize,M-1); 
    auto move_it = in;
    for (size_t kx = 0; kx<ksize; ++kx) { 
        CubicTraits<M-1>::template fill<IteratorType, ArgTypes..., RealType> (move_it, t, grid, otherpos..., grid[kx]);
        move_it += mult;
        }
}

template <size_t M> 
template <class ContainerType, typename ...ArgTypes> 
inline void CubicTraits<M>::fillContainer(ContainerType &in, RealType t, const KMesh& grid, ArgTypes ... otherpos)
{
    assert(grid.getSize() == in.getSize());
    size_t ksize = grid.getSize();
    for (size_t kx = 0; kx<ksize; ++kx) { 
        CubicTraits<M-1>::fillContainer(in[kx], t, grid, otherpos..., grid[kx]);
        }
}



template <class IteratorType, typename ...ArgTypes> 
inline void CubicTraits<0>::fill (IteratorType in, RealType t, const KMesh& grid, ArgTypes ... otherpos)
{
    *in = ek(t,otherpos...);
}

template <class DataType, typename ...ArgTypes> 
inline void CubicTraits<0>::fillContainer(DataType &in, RealType t, const KMesh& grid, ArgTypes ... otherpos)
{
    assert(grid.getSize() == in.getSize());
    in = ek(t,otherpos...);
}


template <typename ArgType1, typename ...ArgTypes> 
inline RealType CubicTraits<0>::ek(RealType t, ArgType1 kpoint1, ArgTypes... kpoints) 
{
    static_assert(std::is_convertible<ArgType1, RealType>::value,"Wrong kpoint");
    assert (kpoint1>=0 && kpoint1 < 2*PI);
    return -2.0*t*cos(kpoint1)+ek(t, kpoints...);
}
 
template <typename ArgType1> 
inline RealType CubicTraits<0>::ek(RealType t, ArgType1 kpoint1)
{
    static_assert(std::is_convertible<ArgType1, RealType>::value,"Wrong kpoint");
    assert (kpoint1>=0 && kpoint1 < 2*PI);
    return -2*t*cos(kpoint1);
}

//
// CubicDMFT
//

template <class Solver, size_t D>
inline CubicDMFTSC<Solver,D>::CubicDMFTSC ( const Solver &S, RealType t, KMesh kGrid):
    SelfConsistency<Solver>(S),
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


template <class Solver, size_t D>
template <typename ...ArgTypes>
inline RealType CubicDMFTSC<Solver,D>::dispersion(ArgTypes... kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D, "Number of points mismatch!" );
    return CubicTraits<0>::ek(_t, kpoints...);
}

template <class Solver, size_t D>
template <typename ...ArgTypes>
inline RealType CubicDMFTSC<Solver,D>::dispersion(const std::tuple<ArgTypes...>& kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D, "Number of points mismatch!" );
    typename EkStorage::PointFunctionType f1 = [&](ArgTypes... kpoints)->RealType{return dispersion(kpoints...);};
    auto f2 =__fun_traits<typename EkStorage::PointFunctionType>::getTupleF(f1); 
    return std::real(f2(kpoints));
}


template <class Solver, size_t D>
template <typename ...ArgTypes>
inline typename CubicDMFTSC<Solver,D>::GFType CubicDMFTSC<Solver,D>::glat(ArgTypes... kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D,"!");
    auto e = dispersion<ArgTypes...>(kpoints...);
    GFType out = 1.0/(1.0/_S.gw+_S.Delta-e);
    return out;

}

template <class Solver, size_t D>
template <typename MPoint, typename ...ArgTypes> 
ComplexType CubicDMFTSC<Solver,D>::glat_val(MPoint w, ArgTypes... kpoints) const
{
    return 1.0/(1.0/_S.gw(w)+_S.Delta(w)-dispersion(kpoints...));
}

template <class Solver, size_t D>
typename CubicDMFTSC<Solver,D>::GKType CubicDMFTSC<Solver,D>::getGLat(const FMatsubaraGrid& fGrid) const
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
        //return (_S.mu - _S.w_1*_S.U-_ek(ktuple))/std::abs(w*w)+1.0/w;
        };
    out._f = __fun_traits<typename GKType::FunctionType>::getFromTupleF(glatdmft_f);

    return out;
}

template <class Solver, size_t D>
inline typename CubicDMFTSC<Solver,D>::GFType CubicDMFTSC<Solver,D>::operator()()
{
    INFO("Using DMFT self-consistency on a cubic lattice in " << D << " dimensions on a lattice of " << _kGrid.getSize() << "^" << D << " atoms.");
    GFType out(this->_S.w_grid); 
    out=0.0; 
    size_t ksize = _kGrid.getSize();
    RealType knorm = pow(ksize,D);
    for (auto w : _gloc.getGrid().getVals()) {
        EkStorage e1 = (1.0/(1.0/_S.gw(w)+_S.Delta(w)-_ek)); 
        _gloc.get(w) = e1.sum()/knorm;
        out.get(w) = -1.0/_gloc(w)+_S.mu-_S.Sigma(w)+ComplexType(w);
    }
    //out._f = std::bind([&](ComplexType w){return _t*_t*2*RealType(D)/w;}, std::placeholders::_1);
    out._f = std::bind([&](ComplexType w)->ComplexType{return _t*_t*2*RealType(D)*((_S.mu-_S.w_1*_S.U)/std::abs(w*w) + 1.0/w);}, std::placeholders::_1);
    return out;
}

template <class Solver, size_t D>
template <typename MPoint>
inline typename CubicDMFTSC<Solver,D>::GFType CubicDMFTSC<Solver,D>::getBubble0(MPoint in) const
{
    std::array<KMesh::point,D> q;
    auto q1 = _kGrid.findClosest(0.0);
    q.fill(q1);
    auto args = std::tuple_cat(std::forward_as_tuple(in),q);
    return Diagrams::getBubble(this->getGLat(_S.w_grid),args);
}

template <class Solver, size_t D>
template <typename MPoint>
inline typename CubicDMFTSC<Solver,D>::GFType CubicDMFTSC<Solver,D>::getBubblePI(MPoint in) const
{
    std::array<KMesh::point,D> q;
    auto q1 = _kGrid.findClosest(PI);
    DEBUG(q1);
    DEBUG(_kGrid[size_t(q1)+1]);
    q.fill(q1);
    auto args = std::tuple_cat(std::forward_as_tuple(in),q);
    return Diagrams::getBubble(this->getGLat(_S.w_grid),args);
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

template <class Solver>
template <typename MPoint>
inline typename CubicInfDMFTSC<Solver>::GFType CubicInfDMFTSC<Solver>::getBubble0(MPoint in) const
{
    auto T = 1.0/_S.w_grid._beta;
    GFType iwn(_S.w_grid);
    iwn.fill([](ComplexType w){return w;});
    return 2.0*T/_t/_t*(1.0-(iwn+_S.mu-_S.Sigma)*_S.gw);
} 

template <class Solver>
inline GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> CubicInfDMFTSC<Solver>::getBubble0() const
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
#endif // endif :: #ifndef ___FK_SELFCONSISTENCY_HPP___
