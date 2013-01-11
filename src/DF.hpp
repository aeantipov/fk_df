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
    GD(GD0.getGrids()),
    SigmaD(GD0.getGrids()), 
    GLat(GD0.getGrids())
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
    GKType GD_shifted(GD.getGrids());
    GD_shifted = GD.shift(in);
    GD_shifted*=GD;
    //GKType GD_shifted(GD.shift(in)*GD); // G(w+W,k+Q)

    for (auto iw: _fGrid.getVals()) { 
        size_t iwn = size_t(iw);
        out[iwn] = GD_shifted[iwn].sum()/RealType(__power<ksize,D>::value);
    }

    RealType T = 1.0/(this->_fGrid._beta);
    return (-T)*out;
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
typename DFLadder<Solver,D,ksize>::GKType DFLadder<Solver,D,ksize>::getGLat(const FMatsubaraGrid &gridF ) const
{
    GKType out(std::tuple_cat(std::forward_as_tuple(gridF), CubicTraits<D,ksize>::getTuples(_kGrid)));
    auto f1 = [&](const typename GKType::PointTupleType &in){return this->GLat(in);};
    auto f2 = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(f1);
    out.fill(f1);
    out._f = GLat._f;
    return out;
}

template <class Solver, size_t D, size_t ksize>
typename DFLadder<Solver,D,ksize>::GLocalType DFLadder<Solver,D,ksize>::operator()()
{
    INFO("Using DF Ladder self-consistency in " << D << " dimensions on a cubic lattice of " << ksize << "^" << D <<" atoms.");
    SigmaD = 0.0;
    RealType beta = _fGrid._beta;
    RealType T = 1.0/beta;
    GLocalType gw(_fGrid); // non-const method. Better copy.
    gw = _S.gw;
    GLocalType Delta(_fGrid); // non-const method. Better copy. 
    GLocalType GDLoc(_fGrid); 
    GLocalType GLatLoc(_fGrid); 
    Delta = _S.Delta;
    GLocalType Delta_out(_fGrid); Delta_out=0.0;
    auto wkgrids = std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kGrid));
    GKType GLatDMFT = CubicDMFTSC<Solver, D, ksize>::getGLat(_fGrid);
    GLocalType GDsum(_fGrid);

    for (auto iw : _fGrid.getVals()) {
        size_t iwn = size_t(iw);
        GD0[iwn] = GLatDMFT[iwn] - gw(iw);
        GDsum[iwn] = GD0[iwn].sum()/RealType(__power<ksize,D>::value);
    };

    typedef typename GKType::ArgTupleType wkTupleType;
    auto gd_f = [&](const wkTupleType &in)->ComplexType{
            ComplexType w = std::get<0>(in);
            auto ktuple = __tuple_tail(in);
            ComplexType e = _ek(ktuple);
            auto tmp = (-e)/std::abs(w*w) -(e*e-_t*_t*2*RealType(D) - 2.0*e*(_S.mu - _S.Sigma._f(w)))/w/std::abs(w*w); 
            return tmp;};
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
            std::function<ComplexType(FMatsubaraGrid::point)> Vertexf = 
                std::bind(&FKImpuritySolver::getVertex4<BMatsubaraGrid::point, FMatsubaraGrid::point>, std::cref(_S), iW, std::placeholders::_1);
            Vertex4.fill(Vertexf);

            size_t nqpoints = _qGrids[0].getSize();
            size_t totalqpts = int(pow(nqpoints,D));
            std::array<KMesh::point, D> q;
            for (size_t nq=0; nq<totalqpts; ++nq) { // iterate over all qpoints
                size_t offset = 0;
                for (size_t i=0; i<D; ++i) { q[D-1-i]=_qGrids[i][(nq-offset)/(int(pow(nqpoints,i)))%nqpoints]; offset+=(int(pow(nqpoints,i)))*size_t(q[D-1-i]); };
                WQTupleType Wq_args = std::tuple_cat(std::forward_as_tuple(iW),q);
                
                INFO_NONEWLINE("\t" << nq << "/" << totalqpts<< ". ");
                //DEBUG("qx = " << qx << ", qy = " << qy);
                Chi0 = this->getBubble(Wq_args);
                GLocalType IrrVertex4(Vertex4);
                GridObject<RealType,FMatsubaraGrid> EVCheck(_fGrid); 
                GridObject<RealType,FMatsubaraGrid> EVCheckRe(_fGrid); 
                GridObject<RealType,FMatsubaraGrid> EVCheckIm(_fGrid); 
                std::function<RealType(FMatsubaraGrid::point)> absEVf = [&](FMatsubaraGrid::point w)->RealType{return std::abs(Chi0(w)*Vertex4(w)); };
                std::function<RealType(FMatsubaraGrid::point)> reEVf = [&](FMatsubaraGrid::point w)->RealType{return std::real(Chi0(w)*Vertex4(w)); };
                std::function<RealType(FMatsubaraGrid::point)> imEVf = [&](FMatsubaraGrid::point w)->RealType{return std::imag(Chi0(w)*Vertex4(w)); };
                EVCheck.fill(absEVf);
                EVCheckRe.fill(reEVf);
                EVCheckIm.fill(imEVf);
                RealType max_ev = *std::max_element(EVCheck.getData().begin(), EVCheck.getData().end());
                RealType max_ev_re = *std::max_element(EVCheckRe.getData().begin(), EVCheckRe.getData().end());
                RealType max_ev_im = *std::max_element(EVCheckIm.getData().begin(), EVCheckIm.getData().end());
                INFO("Maximum EV of Chi0*gamma = " << max_ev << "|" << max_ev_re << "|" << max_ev_im);
                if (std::abs(max_ev-1.0) < 1e-6 || _eval_BS_SC) {
                    GLocalType IrrVertex4_old(Vertex4);
                    INFO2 ("Caught divergence, evaluating BS equation self_consistently. ");
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
                    #ifndef NDEBUG
                    DEBUG("Evaluating BS equation using inversion");
                    #endif
                    IrrVertex4 = Vertex4/(1.0 - Chi0 * Vertex4);
                    //DEBUG(IrrVertex4);
                    }
                auto GD_shift = GD0.shift(Wq_args);
                typename GKType::PointFunctionType SigmaF;
                auto SigmaF2 = [&](wkTupleType in)->ComplexType { 
                    ComplexType w = std::get<0>(in);
                    return 0.5*Vertex4(w)*Chi0(w)*GD_shift(in)*IrrVertex4(w);
                    };
                SigmaF = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(SigmaF2);
                GKType tmp(this->SigmaD.getGrids());
                tmp.fill(SigmaF);
                //DEBUG("Vertex4 = " << Vertex4);
                //DEBUG("IrrVertex4 = " << IrrVertex4);
                SigmaD+=tmp*T/2.0/totalqpts;
                } // end of q loop
        }; //end of iW loop

        GD = 1.0/(1.0/GD0 - SigmaD); // Dyson eq;
        diffGD = GD.diff(GD0);
        INFO("DF diff = " << diffGD);
        GD=GD*_GDmix + GD0*(1.0-_GDmix);
        GD._f = GD0._f; // assume DMFT asymptotics are good 
    };
    INFO("Finished DF iterations");

    // Finish - prepare all lattice quantities
    for (auto iw : _fGrid.getVals()) {
        size_t iwn = size_t(iw);
        GLat[iwn] = 1.0/(Delta(iw) - _ek.getData()) + 1.0/(Delta(iw) - _ek.getData())/gw(iw)*GD[iwn]/gw(iw)/(Delta(iw) - _ek.getData());
        GDLoc[iwn] = GD[iwn].sum()/RealType(__power<ksize,D>::value);
        GLatLoc[iwn] = GLat[iwn].sum()/RealType(__power<ksize,D>::value);
    }
    Delta_out = Delta + 1.0/gw * GDLoc / GLatLoc;
    // Assume DMFT asymptotics
    Delta_out._f = std::bind([&](ComplexType w)->ComplexType{return _t*_t*2*RealType(D)*((_S.mu-_S.w_1*_S.U)/std::abs(w*w) + 1.0/w);}, std::placeholders::_1);
    return Delta_out;
}

} // end of namespace FK
#endif // endif :: #ifndef __FK_DF_H__
