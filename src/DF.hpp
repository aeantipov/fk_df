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
std::vector<ComplexType> DFLadder<D>::getStaticLatticeSusceptibility(const std::vector<std::array<KPoint, D>>& qpts, const FMatsubaraGrid& gridF)
{

    auto grids = std::tuple_cat(std::make_tuple(gridF),__repeater<KMesh,D>::get_tuple(_kGrid));

    // Prepare static vertex (outdated)
    //GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> FullStaticV4(std::forward_as_tuple(gridF,gridF)); 
    //decltype(StaticVertex4)::PointFunctionType VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return _S.getVertex4(0.0, w1,w2);};
    //StaticVertex4.fill(VertexF2);
    //auto StaticV4 = StaticVertex4.getData().getAsMatrix();

    // Prepate interpolated Green's functions
    GKType GD0_interp (grids), GD_interp(grids), GLat_interp(grids);
    GD0_interp.copyInterpolate(GD0);
    GD_interp.copyInterpolate(GD);
    GLat_interp.copyInterpolate(GLat);
    GLocalType Lambda(gridF);
    Lambda.copyInterpolate(_S.getLambda());
    GKType Lwk = this->getGLatDMFT(gridF)/GD0_interp*(-1.0);
    Lwk._f = __fun_traits<typename GKType::FunctionType>::constant(-1.0);
    auto GDL = GD_interp*Lwk;

    // Prepare output
    size_t nqpts = qpts.size();
    std::vector<ComplexType> out;
    out.reserve(nqpts);

    for (auto q : qpts) {

        ComplexType susc=0.0;
        INFO_NONEWLINE("Evaluation static susceptibility for q=["); for (int i=0; i<D; ++i) INFO_NONEWLINE(RealType(q[i])<<" "); INFO("]");
        auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);
        auto GD_shift = GD_interp.shift(Wq_args_static);
        auto LatticeBubble = Diagrams::getBubble(GLat_interp, Wq_args_static);

        auto GDL_bubble = Diagrams::getBubble(GDL, Wq_args_static);
        auto GDL_bubble_vector = GDL_bubble.getData().getAsVector();

        auto dual_bubble = Diagrams::getBubble(GD_interp, Wq_args_static);
        auto dual_bubble_matrix = dual_bubble.getData().getAsDiagonalMatrix();

        auto mult = _S.beta*_S.U*_S.U*_S.w_0*_S.w_1;
        auto m1 = mult*dual_bubble*Lambda*Lambda;
        ComplexType B=(m1/(1.0+m1)).getData().sum();
        GLocalType B1=m1*Lambda/(1.0+m1);
    
        ComplexType susc1;
        for (auto w1 : gridF.getPoints()) {
            auto F = mult/(1.0+m1(w1))*Lambda(w1)/(1.0-B);
            for (auto w2 : gridF.getPoints()) {
                RealType kronecker = RealType(w1._index == w2._index);
                susc1+=GDL_bubble(w1)*F*(Lambda(w2)-Lambda(w1)*kronecker+B*Lambda(w1)*kronecker-B1(w2))*GDL_bubble(w2);
            }
        }

        //auto FullStaticV4 = Diagrams::BS(dual_bubble_matrix, StaticV4, true);
        //auto susc1 = (GDL_bubble_vector.transpose()*FullStaticV4*GDL_bubble_vector);
        //susc = susc1(0,0);
        
        susc = susc1;

        susc+=LatticeBubble.sum();
        INFO("Static susceptibility at q=" << q[0] << " = " << susc);
        out.push_back(susc);
        };

    return out;
}



template <size_t D>
template <typename KPoint>
ComplexType DFLadder<D>::getStaticLatticeSusceptibility(const std::array<KPoint, D>& q, const FMatsubaraGrid& gridF)
{
    std::vector<std::array<KPoint,D>> qpts;
    qpts.push_back(q);
    return getStaticLatticeSusceptibility(qpts,gridF)[0];
}


template <size_t D>
template <typename KPoint>
ComplexType DFLadder<D>::getStaticLatticeSusceptibility(const std::array<KPoint, D>& q)
{
    return getStaticLatticeSusceptibility(q,_fGrid);
}

} // end of namespace FK
#endif // endif :: #ifndef __FK_DF_H__
