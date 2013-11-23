#ifndef ___FK_DMFT_HPP___
#define ___FK_DMFT_HPP___

#include "Diagrams.h"

namespace FK {

template <typename MPoint>
typename DMFTBase::GFType DMFTBase::getBubblePI(MPoint in) const
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
// LatticeDMFT
//


template <typename LatticeT, size_t D>
template <typename ...ArgTypes>
inline RealType LatticeDMFTSC<LatticeT,D>::dispersion(ArgTypes... kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D, "Number of points mismatch!" );
    return _ek._f(kpoints...);
}

template <typename LatticeT, size_t D>
template <typename ...ArgTypes>
inline RealType LatticeDMFTSC<LatticeT,D>::dispersion(const std::tuple<ArgTypes...>& kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D, "Number of points mismatch!" );
    std::function<RealType(ArgTypes...)> f1 = [&](ArgTypes... kpoints)->RealType{return dispersion(kpoints...);};
    auto f2 =__fun_traits<decltype(f1)>::getTupleF(f1); 
    return std::real(f2(kpoints));
}


template <typename LatticeT, size_t D>
template <typename ...ArgTypes>
typename LatticeDMFTSC<LatticeT,D>::GFType LatticeDMFTSC<LatticeT,D>::glat(ArgTypes... kpoints) const
{
    static_assert(sizeof...(ArgTypes) == D,"!");
    auto e = dispersion<ArgTypes...>(kpoints...);
    GFType out = 1.0/(1.0/_S.gw+_S.Delta-e);
    return out;

}

template <typename LatticeT, size_t D>
template <typename MPoint, typename KPoint> 
inline typename LatticeDMFTSC<LatticeT,D>::GFType LatticeDMFTSC<LatticeT,D>::getBubble(MPoint in, std::array<KPoint,D> q) const
{
    auto args = std::tuple_cat(std::forward_as_tuple(in),q);
    auto out = Diagrams::getBubble(this->getGLat(_S.w_grid),args);
    auto T = 1.0/_S.beta; auto w_1 = _S.w_1; auto mu = _S.mu; auto U = _S.U;
    out._f = std::bind([T,in,w_1,mu,U](ComplexType w){return -T/(w*(w+ComplexType(in)))*(1.0 - (1.0/w + 1.0/((w+ComplexType(in))))*(mu-w_1*U));},std::placeholders::_1); 
    return out;
}

template <typename LatticeT, size_t D>
template <typename MPoint>
inline typename LatticeDMFTSC<LatticeT,D>::GFType LatticeDMFTSC<LatticeT,D>::getBubble0(MPoint in) const
{
    std::array<KMesh::point,D> q;
    auto q1 = _kGrid.findClosest(0.0);
    q.fill(q1);
    return getBubble(in,q);
}

template <typename LatticeT, size_t D>
template <typename MPoint>
inline typename LatticeDMFTSC<LatticeT,D>::GFType LatticeDMFTSC<LatticeT,D>::getBubblePI(MPoint in) const
{
    std::array<KMesh::point,D> q;
    auto q1 = _kGrid.findClosest(PI);
    q.fill(q1);
    return getBubble(in,q);
}

template <typename LatticeT, size_t D>
LatticeDMFTSC<LatticeT,D>::LatticeDMFTSC ( const FKImpuritySolver &S, KMesh kGrid, RealType t):
    DMFTBase(S),
    _t(t),
    _kGrid(kGrid),
    _ek(__repeater<KMesh,D>::get_tuple(_kGrid)),
    _gloc(this->_S.w_grid)
{
    _ek._f = lattice_traits::template get_dispersion<typename EkStorage::FunctionType> (t); 
    _ek.fill(_ek._f);
}

template <typename LatticeT, size_t D>
typename LatticeDMFTSC<LatticeT,D>::GKType LatticeDMFTSC<LatticeT,D>::getGLat(const FMatsubaraGrid& fGrid) const
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

template <typename LatticeT, size_t D>
typename LatticeDMFTSC<LatticeT,D>::GFType LatticeDMFTSC<LatticeT,D>::operator()()
{
    INFO("Using DMFT self-consistency on a cubic lattice in " << D << " dimensions on a lattice of " << _kGrid.getSize() << "^" << D << " atoms.");
    GFType out(this->_S.w_grid); 
    out=0.0; 
    size_t ksize = _kGrid.getSize();
    RealType knorm = pow(ksize,D);
    for (auto w : _gloc.getGrid().getPoints()) {
        EkStorage e1 = (1.0/(1.0/_S.gw(w)+_S.Delta(w)-_ek)); 
        _gloc.get(w) = e1.sum()/knorm;
        //out.get(w) = -1.0/_gloc(w)+_S.mu-_S.Sigma(w)+ComplexType(w);
        out.get(w) = -1.0/_gloc(w)+_S.Delta(w) + 1.0/_S.gw(w); 
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
std::vector<RealType> getStaticLatticeDMFTSusceptibility(const SolverType& Solver, const std::vector<typename SolverType::GFType>& Bubbles_in, const FMatsubaraGrid& gridF)
{
    typename SolverType::GFType gw(gridF), Sigma(gridF), K0(gridF), K1(gridF);
    gw.copyInterpolate(Solver.gw);
    Sigma.copyInterpolate(Solver.Sigma);
    RealType beta = Solver.beta;
    RealType T=1.0/beta;

    typename SolverType::GFType bubble(gridF);

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
RealType getStaticLatticeDMFTSusceptibility(const SolverType& Solver, const typename SolverType::GFType& Bubble_in, const FMatsubaraGrid& gridF)
{
    std::vector<typename SolverType::GFType> bubbles;
    bubbles.push_back(Bubble_in);
    return getStaticLatticeDMFTSusceptibility(Solver,bubbles,gridF)[0];
}


template <typename SolverType>
std::vector<std::array<RealType,3>> getStaticLatticeDMFTSkeletonSusceptibility(const SolverType& Solver, const std::vector<typename SolverType::GFType>& Bubbles_in, const FMatsubaraGrid& gridF)
{
    typename SolverType::GFType gw(gridF), Sigma(gridF), K0(gridF), K1(gridF);
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
    typename SolverType::GFType bubble(gridF);

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
std::array<RealType,3> getStaticLatticeDMFTSkeletonSusceptibility(const SolverType& Solver, const typename SolverType::GFType& Bubble_in, const FMatsubaraGrid& gridF)
{
    std::vector<typename SolverType::GFType> bubbles;
    bubbles.push_back(Bubble_in);
    return getStaticLatticeDMFTSkeletonSusceptibility(Solver,bubbles,gridF)[0];
}


} // end of namespace FK
#endif // endif :: #ifndef ___FK_DMFT_HPP___
