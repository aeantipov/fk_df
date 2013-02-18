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
ComplexType DFLadder<D>::getStaticLatticeSusceptibility(const std::array<KPoint, D>& q)
{
    ComplexType susc=0.0;
    INFO_NONEWLINE("Evaluation static susceptibility for q=["); for (int i=0; i<D; ++i) INFO_NONEWLINE(RealType(q[i])<<" "); INFO("]");
    auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);
    auto GD_shift = GD.shift(Wq_args_static);
    auto LatticeBubble = Diagrams::getBubble(GLat, Wq_args_static);

    GKType Lwk = this->getGLatDMFT()/GD0*(-1.0);
    Lwk._f = __fun_traits<typename GKType::FunctionType>::constant(-1.0);
    auto GDL = GD*Lwk;
    __tuple_print<decltype(Wq_args_static)>::print(Wq_args_static);
    auto GDL_bubble = Diagrams::getBubble(GDL, Wq_args_static);
    auto GDL_bubble_vector = GDL_bubble.getData().getAsVector();

    // Prepare static vertex
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> StaticVertex4(std::forward_as_tuple(_fGrid,_fGrid)); 
    decltype(StaticVertex4)::PointFunctionType VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return _S.getVertex4(0.0, w1,w2);};
    StaticVertex4.fill(VertexF2);
    auto StaticV4 = StaticVertex4.getData().getAsMatrix();
    auto dual_bubble = Diagrams::getBubble(this->GD, Wq_args_static);
    auto dual_bubble_matrix = dual_bubble.getData().getAsDiagonalMatrix();
    const auto FullStaticV4 = Diagrams::BS(dual_bubble_matrix, StaticV4, true);

    auto susc1 = (GDL_bubble_vector.transpose()*FullStaticV4*GDL_bubble_vector);
    susc = susc1(0,0);

    susc+=LatticeBubble.sum();
    INFO("Static susceptibility at q=" << q[0] << " = " << susc);

    return susc;
}

} // end of namespace FK
#endif // endif :: #ifndef __FK_DF_H__
