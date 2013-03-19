#ifndef __FK_DF_HPP__
#define __FK_DF_HPP__

#include "DF.h"

// Some hints
// GLocalType - g(w)
// GKType - g(w,k...)
// GLocalType::FunctionType g(w) - analytic
// GKType::FunctionType g(w,k...) - analytic
// __repeater<KMesh,D>::get_tuple(_kGrid) - returns a tuple of D kmeshes

namespace FK {

template <size_t D>
template <typename KPoint>
ComplexType DFLadder<D>::getStaticLatticeSusceptibility(const std::array<KPoint, D>& q, const FMatsubaraGrid& gridF)
{

    auto grids = std::tuple_cat(std::make_tuple(gridF),__repeater<KMesh,D>::get_tuple(_kGrid));
    ComplexType susc=0.0;
    INFO_NONEWLINE("Evaluation static susceptibility for q=["); for (int i=0; i<D; ++i) INFO_NONEWLINE(RealType(q[i])<<" "); INFO("]");
    GKType GD0_interp (grids), GD_interp(grids), GLat_interp(grids);
    GD0_interp.copyInterpolate(GD0);
    GD_interp.copyInterpolate(GD);
    GLat_interp.copyInterpolate(GLat);
    auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);
    auto GD_shift = GD_interp.shift(Wq_args_static);
    auto LatticeBubble = Diagrams::getBubble(GLat_interp, Wq_args_static);

    GKType Lwk = this->getGLatDMFT(gridF)/GD0_interp*(-1.0);
    Lwk._f = __fun_traits<typename GKType::FunctionType>::constant(-1.0);
    auto GDL = GD_interp*Lwk;
    __tuple_print<decltype(Wq_args_static)>::print(Wq_args_static);
    auto GDL_bubble = Diagrams::getBubble(GDL, Wq_args_static);
    auto GDL_bubble_vector = GDL_bubble.getData().getAsVector();

    // Prepare static vertex
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> StaticVertex4(std::forward_as_tuple(gridF,gridF)); 
    decltype(StaticVertex4)::PointFunctionType VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return _S.getVertex4(0.0, w1,w2);};
    StaticVertex4.fill(VertexF2);
    auto StaticV4 = StaticVertex4.getData().getAsMatrix();
    auto dual_bubble = Diagrams::getBubble(GD_interp, Wq_args_static);
    auto dual_bubble_matrix = dual_bubble.getData().getAsDiagonalMatrix();
    const auto FullStaticV4 = Diagrams::BS(dual_bubble_matrix, StaticV4, true);

    auto susc1 = (GDL_bubble_vector.transpose()*FullStaticV4*GDL_bubble_vector);
    susc = susc1(0,0);

    susc+=LatticeBubble.sum();
    INFO("Static susceptibility at q=" << q[0] << " = " << susc);

    return susc;
}


template <size_t D>
template <typename KPoint>
ComplexType DFLadder<D>::getStaticLatticeSusceptibility(const std::array<KPoint, D>& q)
{
    return getStaticLatticeSusceptibility(q,_fGrid);
}

} // end of namespace FK
#endif // endif :: #ifndef __FK_DF_H__
