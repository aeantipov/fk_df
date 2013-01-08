#ifndef __FK_DF_HPP__
#define __FK_DF_HPP__

#include "DF.h"

// Some hints
// GLocalType - g(w)
// GKType - g(w,k...)
// GLocalType::FunctionType g(w) - analytic
// GKType::FunctionType g(w,k...) - analytic
// CubicTraits<D,ksize>::getTuples(_kGrid) - returns a tuple of D kmeshes

namespace FK {

template <class Solver, size_t D, size_t ksize>
DFLadder<Solver,D,ksize>::DFLadder(const Solver &S, const FMatsubaraGrid& fGrid, const BMatsubaraGrid& bGrid, const std::array<typename DFLadder<Solver,D,ksize>::qGridType,D>& qGrids, RealType t):
    CubicDMFTSC<Solver,D,ksize>(S,t),
    _fGrid(fGrid),
    _bGrid(bGrid),
    _qGrids(qGrids),
    GD0(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kGrid))),
    SigmaD(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kGrid)))
{
};

template <class Solver, size_t D, size_t ksize>
template <typename PointType>
std::array<KMesh::point, D> DFLadder<Solver,D,ksize>::shiftKPoint(const std::array<KMesh::point, D> &in, const std::array<PointType, D> &shift) const
{
    std::array<KMesh::point,D> out_k;
    for (size_t t=0; t<D; ++t) { 
        out_k[t]._val = (in[t]._val + RealType(shift[t])); 
        out_k[t]._val-= int(out_k[t]._val/(2.0*PI))*2.0*PI; 
        out_k[t]._index = std::get<1>(_kGrid.find(out_k[t]._val)); 
        };
    return out_k;
};


template <class Solver, size_t D, size_t ksize>
template <typename PointType>
FMatsubaraGrid::point DFLadder<Solver,D,ksize>::shiftMatsubara(const FMatsubaraGrid::point &in, const PointType &shift) const
{
    FMatsubaraGrid::point w2; w2._index=0;
    w2._val=in._val+shift;
    auto find_result = _fGrid.find(w2._val);
    if (!std::get<0>(find_result)) { throw (FMatsubaraGrid::exWrongIndex()); } 
    w2._index = std::get<1>(find_result);
    return w2;
}

template <class Solver, size_t D, size_t ksize>
FMatsubaraGrid::point DFLadder<Solver,D,ksize>::shiftMatsubara(const FMatsubaraGrid::point &in, const BMatsubaraGrid::point &shift) const
{
    int n=_bGrid.getNumber(shift);
    FMatsubaraGrid::point w2; w2._index=0;
    w2._val=in._val+shift._val;
    if (-n<int(in)) w2._index=in._index+n;
    return w2;
} 

template <class Solver, size_t D, size_t ksize>
template <typename ...OrigArgTypes, typename ...ShiftArgTypes> 
std::tuple<OrigArgTypes...> DFLadder<Solver,D,ksize>::shiftArgs(const std::tuple<OrigArgTypes...>& in, const std::tuple<ShiftArgTypes...>& shift)
{
    std::tuple<OrigArgTypes...> out(in);
    return out;
}

/*
template <class Solver, size_t D, size_t ksize>
template <typename ...KP>
GKType DFLadder<Solver,D,ksize>::shiftGF(const GKType &in, BMatsubaraGrid::point W, KP...kpoints) const
{
    GFType GD_shifted(in);
    typename GFType::PointFunctionType ShiftFunction = [&](FMatsubaraGrid::point w1, KP...k)->ComplexType { 
        //DEBUG(w1);
        FMatsubaraGrid::point w2 = this->shiftMatsubara(w1,W);
        std::array<KMesh::point, D> K = {{ k... }};
        std::array<KMesh::point, D> Q = {{ q... }};
        std::array<KMesh::point, D> KpQ = shiftKPoint(K,Q); 
        //DEBUG(kx << ";" << ky << "||" << KpQ[0] << ";" << KpQ[1]);
        return this->GD0(std::tuple_cat(std::forward_as_tuple(w2), KpQ));
        };

    GD_shifted.fill(ShiftFunction);

}
*/

template <class Solver, size_t D, size_t ksize>
template <typename ...KP>
typename DFLadder<Solver,D,ksize>::GLocalType DFLadder<Solver,D,ksize>::getBubble(BMatsubaraGrid::point W, KP... q) const
{
    GLocalType out(this->_fGrid);
    GKType GD_shifted(GD0.shift(W,q...)*GD0); // G(w+W,k+Q)

    for (auto iw: _fGrid.getVals()) { 
        size_t iwn = size_t(iw);
        out[iwn] = GD_shifted[iwn].sum()/RealType(__power<ksize,D>::value);
    }

    RealType T = 1.0/(this->_fGrid._beta);
    out*=(-T);
    return out;
}


template <class Solver, size_t D, size_t ksize>
template <typename ...KP>
ComplexType DFLadder<Solver,D,ksize>::getBubble2(BMatsubaraGrid::point W, FMatsubaraGrid::point w1, KP... q) const
{
    typename EkStorage::PointFunctionType GD0Function = [&](KP...k)->ComplexType { return this->GD0(w1,k...); };
    typename EkStorage::PointFunctionType ShiftFunction = [&](KP...k)->ComplexType { return this->GD0.shift(W,q...)(w1,k...); };
 
    EkStorage G1(_ek.getGrids()), G2(_ek.getGrids());
    G1.fill(GD0Function);
    G2.fill(ShiftFunction);
    G1*=G2;
    return G1.sum()/RealType(__power<ksize,D>::value)/(-this->_fGrid._beta);
}

//
// DFLadder2d
//


template <class Solver, size_t ksize>
typename DFLadder2d<Solver,ksize>::GLocalType DFLadder2d<Solver,ksize>::operator()(bool eval_BS_SC) 
{
    INFO("Using DF Ladder self-consistency in 2 dimensions on a cubic lattice of " << ksize << "^2" <<" atoms.");
    SigmaD = 0.0;
    RealType beta = _fGrid._beta;
    RealType T = 1.0/beta;
    GLocalType gw(_fGrid); 
    gw = _S.gw;
    GLocalType Delta(_fGrid); 
    GLocalType GDLoc(_fGrid); 
    GLocalType GLatLoc(_fGrid); 
    Delta = _S.Delta;
    GLocalType Delta_out(_fGrid); Delta_out=0.0;
    auto wkgrids = std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kGrid));
    GKType GLatDMFT(wkgrids);
    GKType GD(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kGrid)));
    GKType GLat(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kGrid)));

    GLocalType GDsum(_fGrid);

    for (auto iw : _fGrid.getVals()) {
        size_t iwn = size_t(iw);
        GLatDMFT[iwn] = 1.0/(1.0/gw(iw)+Delta(iw)-_ek.getData());
        GD0[iwn] = GLatDMFT[iwn] - gw(iw);
        GDsum[iwn] = GD0[iwn].sum()/RealType(__power<ksize,D>::value);
    };

    auto _f1 = [&](ComplexType w, RealType kx, RealType ky){return (_S.mu - _S.Sigma._f(w)-_ek(kx,ky))/std::abs(w*w)+1.0/w;};
    GLatDMFT._f = _f1;
    auto _f2 = [&](ComplexType w, RealType kx, RealType ky){
            return (-_ek(kx,ky))/std::abs(w*w) + I*imag(Delta._f(w))/std::abs(w*w)-(_ek(kx,ky)*_ek(kx,ky) - 2.0*_ek(kx,ky)*(_S.mu - _S.Sigma._f(w)))/w/std::abs(w*w); };
    GD0._f = _f2;
    GDsum.savetxt("GDsum.dat");

    GD = GD0;
    
/*
    // Prepare Chi0.
    GridObject<ComplexType, BMatsubaraGrid, FMatsubaraGrid, qGridType, qGridType> Chi0full(std::forward_as_tuple(_bGrid, _fGrid, _qGrids[0], _qGrids[1]));
    typename GridObject<ComplexType, BMatsubaraGrid, FMatsubaraGrid, qGridType, qGridType>::PointFunctionType Chi0FillFunction;
    Chi0FillFunction = std::bind(&DFLadder2d::template getBubble2<KMesh::point, KMesh::point>, this, 
        std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    Chi0full.fill(Chi0FillFunction);

    // Prepare reducible vertex.
    GridObject<ComplexType, BMatsubaraGrid, FMatsubaraGrid> Vertex4Full(std::forward_as_tuple(_bGrid,_fGrid));
    typename decltype(Vertex4Full)::PointFunctionType Vertex4FillFunction;
    Vertex4FillFunction = std::bind(&FKImpuritySolver::getVertex4<BMatsubaraGrid::point, FMatsubaraGrid::point>, std::cref(_S), std::placeholders::_1, std::placeholders::_2);
    Vertex4Full.fill(Vertex4FillFunction);
*/

    GLocalType Vertex4(_fGrid);
    GLocalType Chi0(_fGrid);
    INFO("Evaluating Bethe-Salpeter equation");

    for (auto iW : _bGrid.getVals()) {
        INFO("iW = " << iW);
        std::function<ComplexType(FMatsubaraGrid::point)> f1 = std::bind(&FKImpuritySolver::getVertex4<BMatsubaraGrid::point, FMatsubaraGrid::point>, std::cref(_S), iW, std::placeholders::_1);
        Vertex4.fill(f1);
        for (auto qx : _qGrids[0].getVals()) { 
            for (auto qy : _qGrids[1].getVals()) { 
                //DEBUG("qx = " << qx << ", qy = " << qy);
                Chi0 = this->getBubble(iW,qx,qy);
                //DEBUG(Chi0);
                GLocalType IrrVertex4(Vertex4);
                GridObject<RealType,FMatsubaraGrid> EVCheck(_fGrid); 
                std::function<RealType(FMatsubaraGrid::point)> f1 = [&](FMatsubaraGrid::point w)->RealType{return std::abs(Chi0(w)*Vertex4(w)); };
                EVCheck.fill(f1);
                RealType max_ev = *std::max_element(EVCheck.getData().begin(), EVCheck.getData().end());
                INFO("Maximum EV of Chi0*gamma = " << max_ev);
                if (std::abs(max_ev-1.0) < 1e-6 || eval_BS_SC) {
                    GLocalType IrrVertex4_old(Vertex4);
                    INFO ("Caught divergence, evaluating BS equation self_consistently ");
                    RealType diffBS = 1.0;
                    size_t niter = 10;
                    RealType bs_mix = 0.5; 
                    for (size_t n=0; n<niter && diffBS > 1e-8; ++n) { 
                        INFO("BS iteration " << n << " for iW = " << ComplexType(iW) << ", (qx,qy) = (" << RealType(qx) << "," << RealType(qy) << ").");
                        IrrVertex4 = Vertex4 + Vertex4*Chi0*IrrVertex4_old;
                        auto diffBS = IrrVertex4.diff(IrrVertex4_old);
                        INFO("vertex diff = " << diffBS);
                        IrrVertex4_old = IrrVertex4*bs_mix+(1.0-bs_mix)*IrrVertex4_old;
                        }
                    }
                else { 
                    INFO("Evaluating BS equation using inversion");
                    IrrVertex4 = Vertex4/(1.0 - Chi0 * Vertex4);
                    //DEBUG(IrrVertex4);
                     }
                auto GD_shift = GD.shift(iW,qx,qy);
                typename GKType::PointFunctionType SigmaF = [&](FMatsubaraGrid::point w, KMesh::point kx, KMesh::point ky)->ComplexType { 
                    return Vertex4(w)*Chi0(w)*GD_shift(w,kx,ky)*IrrVertex4(w);
                    };
                GKType tmp(this->SigmaD.getGrids());
                tmp.fill(SigmaF);
                //DEBUG("Vertex4 = " << Vertex4);
                //DEBUG("IrrVertex4 = " << IrrVertex4);
                //DEBUG("Chi0 = " << Chi0);
                SigmaD+=tmp*T/2.0/std::pow(_qGrids[0].getSize(),2);
                }
            }
        };

    GD = 1.0/(1.0/GD0 - SigmaD); // Dyson eq;

    for (auto iw : _fGrid.getVals()) {
        size_t iwn = size_t(iw);
        GLat[iwn] = 1.0/(Delta(iw) - _ek.getData()) + 1.0/(Delta(iw) - _ek.getData())/gw(iw)*GD[iwn]/gw(iw)/(Delta(iw) - _ek.getData());

        GDLoc[iwn] = GD[iwn].sum()/RealType(__power<ksize,D>::value);
        GLatLoc[iwn] = GLat[iwn].sum()/RealType(__power<ksize,D>::value);
    }
//    DEBUG(GD0);
    //DEBUG("e_k = " << _ek);
    //DEBUG("e_k.sum = " << _ek.sum());
    //DEBUG("Delta = " << Delta);
    //DEBUG("GLatDMFT = " << GLatDMFT);
    //DEBUG("GLatDMFT.sum(0) = " << GLatDMFT[0].sum()/RealType(__power<ksize,D>::value));
    //DEBUG("GDLoc = " << GDLoc);
//    DEBUG(GD0-SigmaD);

    Delta_out = Delta + 1.0/gw * GDLoc / GLatLoc;
    //DEBUG(Delta_out-Delta);
    //exit(0);
    return Delta_out;
}

    /*
    DEBUG(GLatDMFT(FMatsubara(_fGrid._w_max-2, beta),PI/4.,PI/4.));
    DEBUG(GLatDMFT(FMatsubara(_fGrid._w_max-1, beta),PI/4.,PI/4.));
    DEBUG(GLatDMFT(FMatsubara(_fGrid._w_max,   beta),PI/4.,PI/4.0));
    DEBUG(GLatDMFT(FMatsubara(_fGrid._w_max+1,   beta),PI/4.,PI/4.0));
    DEBUG("-------");
    DEBUG(GD0(FMatsubara(_fGrid._w_max-2, beta),PI/4.,PI/4.));
    DEBUG(GD0(FMatsubara(_fGrid._w_max-1, beta),PI/4.,PI/4.));
    DEBUG(GD0(FMatsubara(_fGrid._w_max,   beta),PI/4.,PI/4.0));
    DEBUG(GD0(FMatsubara(_fGrid._w_max+1,   beta),PI/4.,PI/4.0));
    exit(0);
    */
//    DEBUG(GD0(FMatsubara(_fGrid._w_max-2,_fGrid._beta),0,0));
//    DEBUG(GD0(FMatsubara(_fGrid._w_max-1,_fGrid._beta),0,0));
//    DEBUG(GD0(FMatsubara(_fGrid._w_max,_fGrid._beta),0,0));
//    DEBUG(GD0(FMatsubara(_fGrid._w_max+1,_fGrid._beta),0,0));
    



} // end of namespace FK
#endif // endif :: #ifndef __FK_DF_H__
