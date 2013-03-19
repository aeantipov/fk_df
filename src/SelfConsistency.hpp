#ifndef ___FK_SELFCONSISTENCY_HPP___
#define ___FK_SELFCONSISTENCY_HPP___

#include "Diagrams.h"

namespace FK {

template <typename MPoint>
inline typename SelfConsistency::GFType SelfConsistency::getBubblePI(MPoint in) const
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

//
// CubicTraits
//

template <size_t M> 
template <typename FunctionType, typename ...ArgTypes> 
inline FunctionType CubicTraits<M>::get_dispersion(RealType t)
{
    return CubicTraits<M-1>::template get_dispersion<FunctionType,ArgTypes...,RealType>(t); 
}

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

template <size_t M>
inline std::vector<BZPoint<M>> CubicTraits<M>::getAllBZPoints(const KMesh& kGrid)
{
    size_t ksize = kGrid.getSize();
    size_t totalqpts = size_t(pow(ksize,M));
    std::vector<BZPoint<M>> out(totalqpts);

    std::array<KMesh::point, M> q;
    for (size_t nq=0; nq<totalqpts; ++nq) { // iterate over all kpoints
        size_t offset = 0;
        for (size_t i=0; i<M; ++i) { 
            q[M-1-i]=kGrid[(nq-offset)/(int(pow(ksize,i)))%ksize]; 
            offset+=(int(pow(ksize,i)))*size_t(q[M-1-i]); 
            };
        out[nq]=q;
        }
    return out;
}

template <size_t M>
inline std::map<std::array<KMesh::point, M>, size_t> CubicTraits<M>::getUniqueBZPoints(const KMesh& kGrid)
{
    auto all_pts = getAllBZPoints(kGrid);
    auto totalqpts = all_pts.size();
    std::map<std::array<KMesh::point, M>, size_t> unique_pts;
    for (size_t nq=0; nq<totalqpts; ++nq) {
        auto q = all_pts[nq];
      //  INFO_NONEWLINE("Considering: " << nq << "/" << totalqpts << " ");
      //  printKPoint(q);
        BZPoint<M> q_unique = q;
        // Flip all pi+x to pi-x
        for (size_t i=0; i<M; ++i) {
            if (RealType(q[i])>PI) q_unique[i]=kGrid.findClosest(2.0*PI-RealType(q_unique[i]));
            }
        // Order x,y,z. Ensures x<=y<=z
        std::sort(q_unique.begin(), q_unique.end());

    //    INFO_NONEWLINE("Got "); printKPoint(q_unique);
        if (unique_pts.find(q_unique)==unique_pts.end()) 
            unique_pts[q_unique]=1;
        else if (q_unique != q)  
            unique_pts[q_unique]++;
        };
    size_t count = 0;
    for (auto it = unique_pts.begin(); it!=unique_pts.end(); it++) { /*printKPoint(it->first); DEBUG(" : " << it->second);*/ count+=it->second; };
    //DEBUG(totalqpts << " == " << count);
    assert(totalqpts == count);
    return unique_pts;
} 

//
// CubicDMFT
//


template <size_t D>
template <typename ...ArgTypes>
inline RealType CubicDMFTSC<D>::dispersion(ArgTypes... kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D, "Number of points mismatch!" );
    return CubicTraits<0>::ek(_t, kpoints...);
}

template <size_t D>
template <typename ...ArgTypes>
inline RealType CubicDMFTSC<D>::dispersion(const std::tuple<ArgTypes...>& kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D, "Number of points mismatch!" );
    std::function<RealType(ArgTypes...)> f1 = [&](ArgTypes... kpoints)->RealType{return dispersion(kpoints...);};
    auto f2 =__fun_traits<decltype(f1)>::getTupleF(f1); 
    return std::real(f2(kpoints));
}


template <size_t D>
template <typename ...ArgTypes>
inline typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::glat(ArgTypes... kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D,"!");
    auto e = dispersion<ArgTypes...>(kpoints...);
    GFType out = 1.0/(1.0/_S.gw+_S.Delta-e);
    return out;

}

template <size_t D>
template <typename MPoint, typename ...ArgTypes> 
ComplexType CubicDMFTSC<D>::glat_analytic(MPoint w, ArgTypes... kpoints) const
{
    return 1.0/(1.0/_S.gw(w)+_S.Delta(w)-dispersion(kpoints...));
}

template <size_t D>
template <typename MPoint, typename ...ArgTypes> 
ComplexType CubicDMFTSC<D>::glat_analytic(std::tuple<MPoint,ArgTypes...> in) const
{
    auto w = std::get<0>(in);
    auto kpts = __tuple_tail(in);
    return 1.0/(1.0/_S.gw(w)+_S.Delta(w)-dispersion(kpts));
}

template <size_t D>
template <typename MPoint>
inline typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::getBubble0(MPoint in) const
{
    std::array<KMesh::point,D> q;
    auto q1 = _kGrid.findClosest(0.0);
    q.fill(q1);
    auto args = std::tuple_cat(std::forward_as_tuple(in),q);
    return Diagrams::getBubble(this->getGLat(_S.w_grid),args);
}

template <size_t D>
template <typename MPoint>
inline typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::getBubblePI(MPoint in) const
{
    std::array<KMesh::point,D> q;
    auto q1 = _kGrid.findClosest(PI);
    q.fill(q1);
    auto args = std::tuple_cat(std::forward_as_tuple(in),q);
    return Diagrams::getBubble(this->getGLat(_S.w_grid),args);
}

//
// CubicInfDMFTSC
//


template <typename MPoint>
inline typename CubicInfDMFTSC::GFType CubicInfDMFTSC::getBubble0(MPoint in) const
{
    auto T = 1.0/_S.w_grid._beta;
    GFType iwn(_S.w_grid);
    iwn.fill([](ComplexType w){return w;});
    return 2.0*T/_t/_t*(1.0-(iwn+_S.mu-_S.Sigma)*_S.gw);
} 



} // end of namespace FK
#endif // endif :: #ifndef ___FK_SELFCONSISTENCY_HPP___
