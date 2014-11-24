#ifndef ___FK_DMFT_HPP___
#define ___FK_DMFT_HPP___

#include "Diagrams.h"

namespace FK {

template <typename MPoint>
typename DMFTBase::GFType DMFTBase::getBubblePI(MPoint in) const
{
    GFType out(this->_S.w_grid);
    GFType gw_shift(_S.gw), Sigma_shift(_S.Sigma);
    real_type T = 1.0/_S.w_grid.beta();
    if (std::abs(complex_type(in))<PI*T) {
         gw_shift = _S.gw.shift(in);
         Sigma_shift = _S.Sigma.shift(in);
        };
    GFType iwn(this->_S.w_grid); iwn.fill(typename GFType::function_type([](complex_type w){return w;}));
    out = -T*(_S.gw+gw_shift)/(2*iwn+complex_type(in)+2.0*_S.mu-_S.Sigma-Sigma_shift);
    return out;
}

//
// LatticeDMFT
//

template <typename LatticeT>
template <typename ...LatticeParams> 
LatticeDMFTSC<LatticeT>::LatticeDMFTSC(const FKImpuritySolver &S, kmesh kGrid, LatticeParams ... lattice_p):
    DMFTBase(S),
    lattice(lattice_traits(lattice_p...)),
    _kGrid(kGrid),
    _ek(tuple_tools::repeater<kmesh,_D>::get_tuple(_kGrid)),
    _gloc(this->_S.w_grid)
{
    //_ek.tail_ = lattice.template get_dispersion<typename EkStorage::function_type> (lattice_p...);
    _ek.tail_ = lattice.get_dispersion();//<typename EkStorage::function_type> (lattice_p...);
    _ek.fill(_ek.tail_);
}

template <typename LatticeT>
template <typename ...ArgTypes>
inline real_type LatticeDMFTSC<LatticeT>::dispersion(ArgTypes... kpoints) const
{
    static_assert(sizeof...(ArgTypes) == _D, "Number of points mismatch!" );
    return _ek.tail_(kpoints...);
}

template <typename LatticeT>
template <typename ...ArgTypes>
inline real_type LatticeDMFTSC<LatticeT>::dispersion(const std::tuple<ArgTypes...>& kpoints) const
{
    static_assert(sizeof...(ArgTypes) == _D, "Number of points mismatch!" );
    std::function<real_type(ArgTypes...)> f1 = [&](ArgTypes... kpoints)->real_type{return dispersion(kpoints...);};
    auto f2 =tools::tuple_f(f1);
    return std::real(f2(kpoints));
}


template <typename LatticeT>
template <typename ...ArgTypes>
typename LatticeDMFTSC<LatticeT>::GFType LatticeDMFTSC<LatticeT>::glat(ArgTypes... kpoints) const
{
    static_assert(sizeof...(ArgTypes) == _D,"!");
    auto e = dispersion<ArgTypes...>(kpoints...);
    GFType out = 1.0/(1.0/_S.gw+_S.Delta-e);
    return out;

}

template <typename LatticeT>
template <typename MPoint, typename KPoint> 
inline typename LatticeDMFTSC<LatticeT>::GFType LatticeDMFTSC<LatticeT>::getBubble(MPoint in, std::array<KPoint,_D> q) const
{
    auto args = std::tuple_cat(std::forward_as_tuple(in),q);
    auto out = Diagrams::getBubble(this->getGLat(_S.w_grid),args);
    auto T = 1.0/_S.beta; auto w_1 = _S.w_1; auto mu = _S.mu; auto U = _S.U;
    out.tail_ = std::bind([T,in,w_1,mu,U](complex_type w){return -T/(w*(w+complex_type(in)))*(1.0 - (1.0/w + 1.0/((w+complex_type(in))))*(mu-w_1*U));},std::placeholders::_1);
    return out;
}

template <typename LatticeT>
template <typename MPoint>
inline typename LatticeDMFTSC<LatticeT>::GFType LatticeDMFTSC<LatticeT>::getBubble0(MPoint in) const
{
	std::array<kmesh::point,_D> q = tuple_tools::repeater<typename kmesh::point, _D>::get_array(_kGrid.find_nearest(0.0));
    return getBubble(in,q);
}

template <typename LatticeT>
template <typename MPoint>
inline typename LatticeDMFTSC<LatticeT>::GFType LatticeDMFTSC<LatticeT>::getBubblePI(MPoint in) const
{
    std::array<kmesh::point,_D> q = tuple_tools::repeater<typename kmesh::point, _D>::get_array(_kGrid.find_nearest(PI));
    return getBubble(in,q);
}

template <typename LatticeT>
typename LatticeDMFTSC<LatticeT>::GKType LatticeDMFTSC<LatticeT>::getGLat(const fmatsubara_grid& fGrid) const
{
    std::array<kmesh,_D> kgrids;
    kgrids.fill(_kGrid);
    GKType out(std::tuple_cat(std::forward_as_tuple(fGrid),kgrids));
    //for (auto iw : fGrid.points()) {
    std::function<complex_type(typename GKType::point_tuple in)> f1 =
    	[&](typename GKType::point_tuple in)->complex_type {
        fmatsubara_grid::point w = std::get<0>(in);
        auto ktuple = tuple_tools::tuple_tail(in);
        return 1.0/(1.0/_S.gw(w.value())+_S.Delta(w.value())-_ek(ktuple));
    };
    typename GKType::point_function_type f2 = tools::extract_tuple_f(f1);
    out.fill(f2);

    std::function<complex_type(typename GKType::arg_tuple in)>glatdmft_f =
    	[&](typename GKType::arg_tuple in)->complex_type{
        complex_type w = std::get<0>(in);
        auto ktuple = tuple_tools::tuple_tail(in);
        return (_S.mu - _S.Sigma.tail_eval(w)-_ek(ktuple))/std::abs(w*w)+1.0/w;
        //return (_S.mu - _S.w_1*_S.U-_ek(ktuple))/std::abs(w*w)+1.0/w;
        };
    out.tail_ = tools::extract_tuple_f(glatdmft_f);

    return out;
}

template <typename LatticeT>
typename LatticeDMFTSC<LatticeT>::GFType LatticeDMFTSC<LatticeT>::operator()()
{
    INFO("Using DMFT self-consistency on a cubic lattice in " << _D << " dimensions on a lattice of " << _kGrid.size() << "^" << _D << " atoms.");
    GFType out(this->_S.w_grid); 
    out=0.0; 
    size_t ksize = _kGrid.size();
    real_type knorm = pow(ksize,_D);
    for (auto w : _gloc.grid().points()) {
        EkStorage e1 = (1.0/(1.0/_S.gw(w)+_S.Delta(w)-_ek)); 
        _gloc.get(w) = e1.sum()/knorm;
        out.get(w) = -1.0/_gloc(w)+_S.Delta(w) + 1.0/_S.gw(w); 
    }
    out.tail_ = std::bind([&](complex_type w)->complex_type{return lattice.disp_square_sum()*((_S.mu-_S.w_1*_S.U)/std::abs(w*w) + 1.0/w);}, std::placeholders::_1);
    return out;
}


//
// CubicInfDMFTSC
//


template <typename MPoint>
inline typename CubicInfDMFTSC::GFType CubicInfDMFTSC::getBubble0(MPoint in) const
{
    auto T = 1.0/_S.w_grid.beta();
    GFType iwn(_S.w_grid);
    iwn.fill(typename GFType::function_type([](complex_type w){return w;}));
    return 2.0*T/_t/_t*(1.0-(iwn+_S.mu-_S.Sigma)*_S.gw);
} 

//
// Extra tools
//

template <typename SolverType>
std::vector<real_type> getStaticLatticeDMFTSusceptibility(const SolverType& Solver, const std::vector<typename SolverType::GFType>& Bubbles_in, const fmatsubara_grid& gridF)
{
    typename SolverType::GFType gw(gridF), Sigma(gridF), K0(gridF), K1(gridF);
    gw.copy_interpolate(Solver.gw);
    Sigma.copy_interpolate(Solver.Sigma);
    real_type beta = Solver.beta;
    real_type T=1.0/beta;

    typename SolverType::GFType bubble(gridF);

    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> Vertex4_2(std::forward_as_tuple(gridF,gridF)); 
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid>::point_function_type VertexF2 = [&](fmatsubara_grid::point w1, fmatsubara_grid::point w2){return Solver.getVertex4(0.0, w1,w2);};
    Vertex4_2.fill(VertexF2);
    auto V4 = Vertex4_2.data().as_matrix();

    std::vector<real_type> out;

    for (const auto& bubble_in : Bubbles_in) {

        bubble.copy_interpolate(bubble_in);

        auto dual_bubble = bubble+T*gw*gw;
        auto dual_bubble_matrix = dual_bubble.data().as_diagonal_matrix();
        auto FullVertex = Diagrams::BS(dual_bubble_matrix, V4, true);
        complex_type susc = 0.0;
        complex_type bare_susc = 0.0;

        /** Vertex expansion. */
        for (auto w1: gridF.points()) { 
            bare_susc+=bubble(w1);
            for (auto w2: gridF.points()) {
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
real_type getStaticLatticeDMFTSusceptibility(const SolverType& Solver, const typename SolverType::GFType& Bubble_in, const fmatsubara_grid& gridF)
{
    std::vector<typename SolverType::GFType> bubbles;
    bubbles.push_back(Bubble_in);
    return getStaticLatticeDMFTSusceptibility(Solver,bubbles,gridF)[0];
}


template <typename SolverType>
std::vector<std::array<real_type,3>> getStaticLatticeDMFTSkeletonSusceptibility(const SolverType& Solver, const std::vector<typename SolverType::GFType>& Bubbles_in, const fmatsubara_grid& gridF)
{
    typename SolverType::GFType gw(gridF), Sigma(gridF), K0(gridF), K1(gridF);
    gw.copy_interpolate(Solver.gw);
    Sigma.copy_interpolate(Solver.Sigma);
    K0.copy_interpolate(Solver.K0);
    K1.copy_interpolate(Solver.K1);
    auto U = Solver.U;
    auto w_0 = Solver.w_0;
    auto w_1 = Solver.w_1;
    //auto mu = Solver.mu;
    real_type beta = Solver.beta;
    real_type T=1.0/beta;
    typename SolverType::GFType bubble(gridF);

    auto Lambda = 1.0 + gw*(2.0*Sigma - U); 
    auto Lambda2 = (1. + gw*Sigma)*(1.+gw*(Sigma-U));

    std::vector<std::array<real_type,3>> out;

    for (const auto& bubble_in : Bubbles_in) {

        bubble.copy_interpolate(bubble_in);
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
std::array<real_type,3> getStaticLatticeDMFTSkeletonSusceptibility(const SolverType& Solver, const typename SolverType::GFType& Bubble_in, const fmatsubara_grid& gridF)
{
    std::vector<typename SolverType::GFType> bubbles;
    bubbles.push_back(Bubble_in);
    return getStaticLatticeDMFTSkeletonSusceptibility(Solver,bubbles,gridF)[0];
}


} // end of namespace FK
#endif // endif :: #ifndef ___FK_DMFT_HPP___
