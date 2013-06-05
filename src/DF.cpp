#include "DF.h"

namespace FK {

// Implementations to compile

template struct DFLadder<1>;  
template struct DFLadder<2>;  
template struct DFLadder<3>;  
template struct DFLadder<4>;  

/*
            for (auto iW : _bGrid.getPoints()) {
                INFO2("iW = " << iW);
                auto Wq_args = std::tuple_cat(std::make_tuple(iW),q);
                typename GLocalType::PointFunctionType VertexFillf = std::bind(
                    &FKImpuritySolver::getVertex4<BMatsubaraGrid::point, FMatsubaraGrid::point, FMatsubaraGrid::point>, std::cref(_S), iW, std::placeholders::_1, std::placeholders::_1); 
                typename GLocalType::FunctionType Vertexf = 
                    std::bind(&FKImpuritySolver::getVertex4<ComplexType, ComplexType, ComplexType>, std::cref(_S), ComplexType(iW), std::placeholders::_1, std::placeholders::_1); 
                DynVertex4.fill(VertexFillf);
                DynVertex4._f = Vertexf;
                // Bethe-Salpeter vertex
                DualDynBubble = Diagrams::getBubble(GD, Wq_args);
                FullDualDynVertex4 = Diagrams::BS(DualDynBubble, DynVertex4, true, _eval_BS_SC, (_n_BS_iter>0?_n_BS_iter-1:0), _BSmix);
                // Sigma
                auto GD_shift_dyn = GD.shift(Wq_args);
                typename GKType::PointFunctionType SigmaF;
                auto SigmaF2 = [&](wkPointTupleType in)->ComplexType { 
                    auto w = std::get<0>(in);
                    return DynVertex4(w)*DualDynBubble(w)*GD_shift(in)*FullDualDynVertex4(w);
                    };
                SigmaF = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(SigmaF2);
                addSigma.fill(SigmaF);
                addSigma*=T;
                INFO3("Sigma contribution diff = " << addSigma.diff(addSigma*0));
                SigmaD+=addSigma/totalqpts;
                };
            */

}; 
