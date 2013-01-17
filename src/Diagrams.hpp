#ifndef __FK_DFDIAGRAMS_HPP__
#define __FK_DFDIAGRAMS_HPP__

#include "Diagrams.h"

namespace FK {

template <typename GKType, typename ... ArgTypes>
inline typename Diagrams::GLocalType Diagrams::getBubble(const GKType &GF, ArgTypes... args)
{
    return getBubble(GF,std::forward_as_tuple(args...));
}

template <typename GKType, typename ... ArgTypes>
inline typename Diagrams::GLocalType Diagrams::getBubble(const GKType &GF, const std::tuple<ArgTypes...>& args)
{
    static_assert(GKType::N==sizeof...(ArgTypes), "Argument size mismatch");
    const auto _fGrid = std::get<0>(GF.getGrids()); 
    GLocalType out(_fGrid);
    GKType GF_shifted(GF.getGrids());
    GF_shifted = GF.shift(args);
    GF_shifted*=GF;

    RealType knorm = 1;
    for (auto val : GF._dims) knorm*=val; 
    knorm/=GF._dims[0];

    for (auto iw: _fGrid.getVals()) { 
        size_t iwn = size_t(iw);
        out[iwn] = GF_shifted[iwn].sum()/RealType(knorm); 
    }

    RealType T = 1.0/(_fGrid._beta);
    return (-T)*out;
}

template <typename ArgType>
inline typename Diagrams::GLocalType Diagrams::getBubble(const GLocalType &GF, ArgType arg)
{
    GLocalType out = GF.shift(arg)*GF;
    RealType T = 1.0/(std::get<0>(GF.getGrids())._beta);
    return (-T)*out;
}



template <typename VertexType>
inline VertexType Diagrams::BS(const VertexType& Chi0, const VertexType &IrrVertex4, bool eval_SC, size_t n_iter, RealType mix)
{
    const auto _fGrid = std::get<0>(IrrVertex4.getGrids());
    VertexType Vertex4_out(IrrVertex4);
    GridObject<RealType,FMatsubaraGrid> EVCheck(_fGrid); 
    //GridObject<RealType,FMatsubaraGrid> EVCheckRe(_fGrid); 
    //GridObject<RealType,FMatsubaraGrid> EVCheckIm(_fGrid); 
    std::function<RealType(FMatsubaraGrid::point)> absEVf = [&](FMatsubaraGrid::point w)->RealType{return std::abs(Chi0(w)*IrrVertex4(w)); };
    //std::function<RealType(FMatsubaraGrid::point)> reEVf = [&](FMatsubaraGrid::point w)->RealType{return std::real(Chi0(w)*IrrVertex4(w)); };
    //std::function<RealType(FMatsubaraGrid::point)> imEVf = [&](FMatsubaraGrid::point w)->RealType{return std::imag(Chi0(w)*IrrVertex4(w)); };
    EVCheck.fill(absEVf);
    //EVCheckRe.fill(reEVf);
    //EVCheckIm.fill(imEVf);
    RealType max_ev = *std::max_element(EVCheck.getData().begin(), EVCheck.getData().end());
    //RealType max_ev_re = *std::max_element(EVCheckRe.getData().begin(), EVCheckRe.getData().end());
    //RealType max_ev_im = *std::max_element(EVCheckIm.getData().begin(), EVCheckIm.getData().end());
    INFO("Maximum EV of Chi0*gamma = " << max_ev); // << "|" << max_ev_re << "|" << max_ev_im);
    if (std::abs(max_ev-1.0) < 1e-6 || eval_SC) {
        VertexType Vertex4_old(IrrVertex4);
        INFO2 ("Caught divergence, evaluating BS equation self_consistently. ");
        RealType diffBS = 1.0;
        for (size_t n=0; n<n_iter && diffBS > 1e-8; ++n) { 
            //INFO("BS iteration " << n << " for iW = " << ComplexType(iW) << ", (qx,qy) = (" << RealType(qx) << "," << RealType(qy) << ").");
            Vertex4_out = IrrVertex4 + IrrVertex4*Chi0*Vertex4_old;
            auto diffBS = IrrVertex4.diff(Vertex4_old);
            INFO("vertex diff = " << diffBS);
            Vertex4_old = Vertex4_out*mix+(1.0-mix)*Vertex4_old;
            }
        }
    else { 
        #ifndef NDEBUG
        DEBUG("Evaluating BS equation using inversion");
        #endif
        Vertex4_out = IrrVertex4/(1.0 - Chi0 * IrrVertex4);
        }
    return Vertex4_out;
}

template <typename VertexType>
inline VertexType Diagrams::getSusc(const VertexType& Chi0, const VertexType &FullVertex4)
{
    return Chi0 + Chi0*FullVertex4*Chi0;
}

} // end of namespace FK 

#endif // endif :: #ifndef __FK_DFDIAGRAMS_HPP__
