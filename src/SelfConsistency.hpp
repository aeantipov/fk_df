#ifndef ___FK_SELFCONSISTENCY_HPP___
#define ___FK_SELFCONSISTENCY_HPP___

#include "Diagrams.h"

namespace FK {

template <typename MPoint>
typename SelfConsistency::GFType SelfConsistency::getBubblePI(MPoint in) const
{
    GFType out(this->_S.w_grid);
    GFType gw_shift(_S.gw), Sigma_shift(_S.Sigma);
    RealType T = 1.0/_S.w_grid._beta;
    if (std::abs(ComplexType(in))<PI*T) {
         gw_shift = _S.gw.shift(in);
         Sigma_shift = _S.Sigma.shift(in);
        };
    GFType iwn(this->_S.w_grid); iwn.fill(typename GFType::FunctionType([](ComplexType w){return w;}));
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
std::vector<BZPoint<M>> CubicTraits<M>::getAllBZPoints(const KMesh& kGrid)
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
inline std::array<RealType,M> CubicTraits<M>::findSymmetricBZPoint(const std::array<RealType,M>& in)
{
    std::array<RealType,M> out;
    for (size_t i=0; i<M; ++i) {
            in[i] = std::fmod(in[i],2.0*PI);
            if (RealType(in[i])>PI) out[i]=2.0*PI-in[i];
            }
        // Order x,y,z. Ensures x<=y<=z
        std::sort(out.begin(), out.end());
    return out;
}

template <size_t M>
inline BZPoint<M> CubicTraits<M>::findSymmetricBZPoint(const BZPoint<M>& in, const KMesh& kGrid)
{
    BZPoint<M> out(in);
    // Flip all pi+x to pi-x
    for (size_t i=0; i<M; ++i) {
        if (RealType(in[i])>PI) out[i]=kGrid.findClosest(2.0*PI-RealType(in[i]));
        }
    // Order x,y,z. Ensures x<=y<=z
    std::sort(out.begin(), out.end());
    return out;
}

template <size_t M>
std::map<BZPoint<M>, std::vector<BZPoint<M>>> CubicTraits<M>::getUniqueBZPoints(const KMesh& kGrid)
{
    auto all_pts = getAllBZPoints(kGrid);
    auto totalqpts = all_pts.size();
    std::map<std::array<KMesh::point, M>, std::vector<BZPoint<M>>> unique_pts;
    for (size_t nq=0; nq<totalqpts; ++nq) {
        auto q = all_pts[nq];
//        INFO_NONEWLINE("Considering: " << nq << "/" << totalqpts << " " << q);
        BZPoint<M> q_unique = findSymmetricBZPoint(q, kGrid);

        if (unique_pts.find(q_unique)==unique_pts.end()) {
            unique_pts[q_unique]=std::vector<BZPoint<M>>();
            unique_pts[q_unique].push_back(q_unique);
            }
        else if (q_unique != q)  
            unique_pts[q_unique].push_back(q);
        };
    size_t count = 0;
    for (auto it = unique_pts.begin(); it!=unique_pts.end(); it++) { 
//        DEBUG(it->first << " : " << it->second.size()); 
        count+=it->second.size(); 
        };
//    DEBUG(totalqpts << " == " << count);
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
    return _ek._f(kpoints...);
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
typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::glat(ArgTypes... kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D,"!");
    auto e = dispersion<ArgTypes...>(kpoints...);
    GFType out = 1.0/(1.0/_S.gw+_S.Delta-e);
    return out;

}

template <size_t D>
template <typename MPoint, typename KPoint> 
inline typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::getBubble(MPoint in, std::array<KPoint,D> q) const
{
    auto args = std::tuple_cat(std::forward_as_tuple(in),q);
    auto out = Diagrams::getBubble(this->getGLat(_S.w_grid),args);
    auto T = 1.0/_S.beta; auto w_1 = _S.w_1; auto mu = _S.mu; auto U = _S.U;
    out._f = std::bind([T,in,w_1,mu,U](ComplexType w){return -T/(w*(w+ComplexType(in)))*(1.0 - (1.0/w + 1.0/((w+ComplexType(in))))*(mu-w_1*U));},std::placeholders::_1); 
    return out;
}

template <size_t D>
template <typename MPoint>
inline typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::getBubble0(MPoint in) const
{
    std::array<KMesh::point,D> q;
    auto q1 = _kGrid.findClosest(0.0);
    q.fill(q1);
    return getBubble(in,q);
}

template <size_t D>
template <typename MPoint>
inline typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::getBubblePI(MPoint in) const
{
    std::array<KMesh::point,D> q;
    auto q1 = _kGrid.findClosest(PI);
    q.fill(q1);
    return getBubble(in,q);
}

template <size_t D>
CubicDMFTSC<D>::CubicDMFTSC ( const FKImpuritySolver &S, RealType t, KMesh kGrid):
    SelfConsistency(S),
    _t(t),
    _kGrid(kGrid),
    _ek(__repeater<KMesh,D>::get_tuple(_kGrid)),
    _gloc(this->_S.w_grid)
{
    _ek._f = CubicTraits<D>::template get_dispersion<typename EkStorage::FunctionType> (t); 
    _ek.fill(_ek._f);
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
typename CubicDMFTSC<D>::GFType CubicDMFTSC<D>::operator()()
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


//
// CubicInfDMFTSC
//


template <typename MPoint>
inline typename CubicInfDMFTSC::GFType CubicInfDMFTSC::getBubble0(MPoint in) const
{
    auto T = 1.0/_S.w_grid._beta;
    GFType iwn(_S.w_grid);
    iwn.fill(typename GFType::FunctionType([](ComplexType w){return w;}));
    return 2.0*T/_t/_t*(1.0-(iwn+_S.mu-_S.Sigma)*_S.gw);
} 

//
// Extra tools
//

template <typename SolverType>
std::vector<RealType> getStaticLatticeDMFTSusceptibility(const SolverType& Solver, const std::vector<GFWrap>& Bubbles_in, const FMatsubaraGrid& gridF)
{
    GFWrap gw(gridF), Sigma(gridF), K0(gridF), K1(gridF);
    gw.copyInterpolate(Solver.gw);
    Sigma.copyInterpolate(Solver.Sigma);
    RealType beta = Solver.beta;
    RealType T=1.0/beta;

    GFWrap bubble(gridF);

    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> Vertex4_2(std::forward_as_tuple(gridF,gridF)); 
    decltype(Vertex4_2)::PointFunctionType VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return Solver.getVertex4(0.0, w1,w2);};
    Vertex4_2.fill(VertexF2);
    auto V4 = Vertex4_2.getData().getAsMatrix();

    std::vector<RealType> out;

    for (const auto& bubble_in : Bubbles_in) {

        bubble.copyInterpolate(bubble_in);

        auto dual_bubble = bubble+T*gw*gw;
        auto dual_bubble_matrix = dual_bubble.getData().getAsDiagonalMatrix();
        auto FullVertex = Diagrams::BS(dual_bubble_matrix, V4, true);
        ComplexType susc = 0.0;
        ComplexType bare_susc = 0.0;

        /** Vertex expansion. */
        for (auto w1: gridF.getPoints()) { 
            bare_susc+=bubble(w1);
            for (auto w2: gridF.getPoints()) {
                susc+=bubble(w1)*FullVertex(size_t(w1),size_t(w2))*bubble(w2); 
                }
            };
        susc+=bare_susc;
        if (std::abs(std::imag(susc))>1e-5) throw 1; 
        out.push_back(std::real(susc));
    };
    return out;
}

template <typename SolverType>
RealType getStaticLatticeDMFTSusceptibility(const SolverType& Solver, const GFWrap& Bubble_in, const FMatsubaraGrid& gridF)
{
    std::vector<GFWrap> bubbles;
    bubbles.push_back(Bubble_in);
    return getStaticLatticeDMFTSusceptibility(Solver,bubbles,gridF)[0];
}


template <typename SolverType>
std::vector<std::array<RealType,3>> getStaticLatticeDMFTSkeletonSusceptibility(const SolverType& Solver, const std::vector<GFWrap>& Bubbles_in, const FMatsubaraGrid& gridF)
{
    GFWrap gw(gridF), Sigma(gridF), K0(gridF), K1(gridF);
    gw.copyInterpolate(Solver.gw);
    Sigma.copyInterpolate(Solver.Sigma);
    K0.copyInterpolate(Solver.K0);
    K1.copyInterpolate(Solver.K1);
    auto U = Solver.U;
    auto w_0 = Solver.w_0;
    auto w_1 = Solver.w_1;
    //auto mu = Solver.mu;
    RealType beta = Solver.beta;
    RealType T=1.0/beta;
    GFWrap bubble(gridF);

    auto Lambda = 1.0 + gw*(2.0*Sigma - U); 
    auto Lambda2 = (1. + gw*Sigma)*(1.+gw*(Sigma-U));

    std::vector<std::array<RealType,3>> out;

    for (const auto& bubble_in : Bubbles_in) {

        bubble.copyInterpolate(bubble_in);
        /** Skeleton expansion. */
        auto nu = gw*(-1.0/gw/gw - T/bubble);
        auto d1 = Lambda*gw*nu + Lambda2;
        auto ugamma_down = 1.0 - (w_0*w_1*U*U*gw*gw*gw*nu/(Lambda2 * (d1))).sum();
        auto ugamma = (w_0*w_1*U*U*gw*gw/(d1)).sum() / (ugamma_down); 

        auto chi_cc = -T*((Lambda - ugamma)*gw*gw/d1).sum();
        auto chi_cf = ugamma/U;
        auto chi_ff = w_0*w_1/T/ugamma_down;
        
        out.push_back({{std::real(chi_cc), std::real(chi_cf), std::real(chi_ff)}});
    };
    return out;
}

template <typename SolverType>
std::array<RealType,3> getStaticLatticeDMFTSkeletonSusceptibility(const SolverType& Solver, const GFWrap& Bubble_in, const FMatsubaraGrid& gridF)
{
    std::vector<GFWrap> bubbles;
    bubbles.push_back(Bubble_in);
    return getStaticLatticeDMFTSkeletonSusceptibility(Solver,bubbles,gridF)[0];
}


} // end of namespace FK
#endif // endif :: #ifndef ___FK_SELFCONSISTENCY_HPP___
