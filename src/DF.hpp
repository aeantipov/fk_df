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
template <typename ...KP>
typename DFLadder<Solver,D,ksize>::GLocalType DFLadder<Solver,D,ksize>::getBubble(BMatsubaraGrid::point W, KP... q) const
{
    return this->getBubble(std::forward_as_tuple(W,q...));
}

template <class Solver, size_t D, size_t ksize>
typename DFLadder<Solver,D,ksize>::GLocalType DFLadder<Solver,D,ksize>::getBubble(const WQTupleType &in) const
{
    GLocalType out(this->_fGrid);
    GKType GD_shifted(GD0.shift(in)*GD0); // G(w+W,k+Q)

    for (auto iw: _fGrid.getVals()) { 
        size_t iwn = size_t(iw);
        out[iwn] = GD_shifted[iwn].sum()/RealType(__power<ksize,D>::value);
    }

    RealType T = 1.0/(this->_fGrid._beta);
    out*=(-T);
    return out;
}

/*
template <class Solver, size_t D, size_t ksize>
template <typename ...KP>
ComplexType DFLadder<Solver,D,ksize>::getBubble2(BMatsubaraGrid::point W, KP... q, FMatsubaraGrid::point w1) const
{
    static auto pts = std::make_tuple(W,q...);
    static auto GD_shift = this->GD0.shift(W,q...);
    auto pts2 = std::make_tuple(W,q...);
    if (pts2!=pts) {pts=pts2; GD_shift = this->GD0.shift(W,q...); };
    typename EkStorage::PointFunctionType GD0Function = [&](KP...k)->ComplexType { return GD0(w1,k...); };
    //typename EkStorage::PointFunctionType ShiftFunction = [&](KP...k)->ComplexType { return GD0(w1,k...); };
    typename EkStorage::PointFunctionType ShiftFunction = [&](KP...k)->ComplexType { return GD_shift(w1,k...); };
 
    EkStorage G1(_ek.getGrids()), G2(_ek.getGrids());
    G1.fill(GD0Function);
    G2.fill(ShiftFunction);
    G1*=G2;
    return G1.sum()/RealType(__power<ksize,D>::value)/(-this->_fGrid._beta);
}*/

template <class Solver, size_t D, size_t ksize>
typename DFLadder<Solver,D,ksize>::GLocalType DFLadder<Solver,D,ksize>::operator()(bool eval_BS_SC)
{
    INFO("Using DF Ladder self-consistency in " << D << " dimensions on a cubic lattice of " << ksize << "^" << D <<" atoms.");
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

    typedef typename GKType::ArgTupleType wkTupleType;
    auto glatdmft_f = [&](const wkTupleType &in)->ComplexType{
        ComplexType w = std::get<0>(in);
        auto ktuple = __tuple_tail(in);
        return (_S.mu - _S.Sigma._f(w)-_ek(ktuple))/std::abs(w*w)+1.0/w;
        };
    GLatDMFT._f = __fun_traits<typename GKType::FunctionType>::getFromTupleF(glatdmft_f);
    auto gd_f = [&](const wkTupleType &in)->ComplexType{
            ComplexType w = std::get<0>(in);
            auto ktuple = __tuple_tail(in);
            return (-_ek(ktuple))/std::abs(w*w) + I*imag(Delta._f(w))/std::abs(w*w)-(_ek(ktuple)*_ek(ktuple) - 2.0*_ek(ktuple)*(_S.mu - _S.Sigma._f(w)))/w/std::abs(w*w); };
    GD0._f = __fun_traits<typename GKType::FunctionType>::getFromTupleF(gd_f);

    GDsum.savetxt("GDsum.dat");
    GD = GD0;

    // Put here operations with GD
    GLocalType Vertex4(_fGrid);
    GLocalType Chi0(_fGrid);
    INFO("Starting ladder dual fermion calculations")
    RealType diffGD = 1.0;
    for (size_t nd_iter=0; nd_iter<_n_GD_iter && diffGD > 1e-8; ++nd_iter) { 
        INFO("DF iteration " << nd_iter);
        INFO("Evaluating Bethe-Salpeter equation");

        for (auto iW : _bGrid.getVals()) {
            INFO("iW = " << iW);
            std::function<ComplexType(FMatsubaraGrid::point)> f1 = std::bind(&FKImpuritySolver::getVertex4<BMatsubaraGrid::point, FMatsubaraGrid::point>, std::cref(_S), iW, std::placeholders::_1);
            Vertex4.fill(f1);

            size_t nqpoints = _qGrids[0].getSize();
            size_t totalqpts = int(pow(nqpoints,D));
            std::array<KMesh::point, D> q;
            for (size_t nq=0; nq<totalqpts; ++nq) { // iterate over all qpoints
                size_t offset = 0;
                for (size_t i=0; i<D; ++i) { q[D-1-i]=_qGrids[i][(nq-offset)/(int(pow(nqpoints,i)))%nqpoints]; offset+=(int(pow(nqpoints,i)))*size_t(q[D-1-i]); };
                WQTupleType Wq_args = std::tuple_cat(std::forward_as_tuple(iW),q);
                
                INFO_NONEWLINE(nq << "/" << totalqpts<< ". ");
                //DEBUG("qx = " << qx << ", qy = " << qy);
                //typename ArgBackGenerator<D,KMesh::point,__caller,GLocalType,BMatsubaraGrid::point>::type chi0_caller; // aka __caller<GLocalType,w,q...>
                //chi0_caller._params = Wq_args;
                //chi0_caller._f = std::bind(&DFLadder2d::getBubble, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3); 
                //Chi0 = this->getBubble(iW,qx,qy);
                Chi0 = this->getBubble(Wq_args);
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
                    for (size_t n=0; n<_n_BS_iter && diffBS > 1e-8; ++n) { 
                        //INFO("BS iteration " << n << " for iW = " << ComplexType(iW) << ", (qx,qy) = (" << RealType(qx) << "," << RealType(qy) << ").");
                        IrrVertex4 = Vertex4 + Vertex4*Chi0*IrrVertex4_old;
                        auto diffBS = IrrVertex4.diff(IrrVertex4_old);
                        INFO("vertex diff = " << diffBS);
                        IrrVertex4_old = IrrVertex4*_BSmix+(1.0-_BSmix)*IrrVertex4_old;
                        }
                    }
                else { 
                    INFO("Evaluating BS equation using inversion");
                    IrrVertex4 = Vertex4/(1.0 - Chi0 * Vertex4);
                    //DEBUG(IrrVertex4);
                    }
                auto GD_shift = GD0.shift(Wq_args);
                typename GKType::PointFunctionType SigmaF;
                auto SigmaF2 = [&](wkTupleType in)->ComplexType { 
                    ComplexType w = std::get<0>(in);
                    return Vertex4(w)*Chi0(w)*GD_shift(in)*IrrVertex4(w);
                    };
                SigmaF = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(SigmaF2);
                GKType tmp(this->SigmaD.getGrids());
                tmp.fill(SigmaF);
                //DEBUG("Vertex4 = " << Vertex4);
                //DEBUG("IrrVertex4 = " << IrrVertex4);
                //DEBUG("Chi0 = " << Chi0);
                SigmaD+=tmp*T/2.0/totalqpts;
                } // end of q loop
        }; //end of iW loop

        GD = 1.0/(1.0/GD0 - SigmaD); // Dyson eq;
        diffGD = GD.diff(GD0);
        INFO("DF diff = " << diffGD);
        GD=GD*_GDmix + GD0*(1.0-_GDmix);
    };

    // Finish - prepare all lattice quantities
    for (auto iw : _fGrid.getVals()) {
        size_t iwn = size_t(iw);
        GLat[iwn] = 1.0/(Delta(iw) - _ek.getData()) + 1.0/(Delta(iw) - _ek.getData())/gw(iw)*GD[iwn]/gw(iw)/(Delta(iw) - _ek.getData());
        GDLoc[iwn] = GD[iwn].sum()/RealType(__power<ksize,D>::value);
        GLatLoc[iwn] = GLat[iwn].sum()/RealType(__power<ksize,D>::value);
    }
    Delta_out = Delta + 1.0/gw * GDLoc / GLatLoc;
    return Delta_out;
}

} // end of namespace FK
#endif // endif :: #ifndef __FK_DF_H__
