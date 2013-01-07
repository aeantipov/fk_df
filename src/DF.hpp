#ifndef __FK_DF_HPP__
#define __FK_DF_HPP__

#include "DF.h"

// Some hints
// GLocalType - g(w)
// GKType - g(w,k...)
// GLocalType::FunctionType g(w) - analytic
// GKType::FunctionType g(w,k...) - analytic
// CubicTraits<D,ksize>::getTuples(_kgrid) - returns a tuple of D kmeshes

namespace FK {

template <class Solver, size_t D, size_t ksize>
DFLadder<Solver,D,ksize>::DFLadder(const Solver &S, const FMatsubaraGrid& fGrid, const BMatsubaraGrid& bGrid, const std::array<typename DFLadder<Solver,D,ksize>::qGridType,D>& qGrids, RealType t):
    CubicDMFTSC<Solver,D,ksize>(S,t),
    _fGrid(fGrid),
    _bGrid(bGrid),
    _qGrids(qGrids),
    GD0(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kgrid))),
    SigmaD(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kgrid)))
{
};

template <class Solver, size_t D, size_t ksize>
template <typename PointType>
std::array<KMesh::point, D> DFLadder<Solver,D,ksize>::_shift_point(const std::array<KMesh::point, D> &in, const std::array<PointType, D> &shift) const
{
    std::array<KMesh::point,D> out_k;
    for (size_t t=0; t<D; ++t) { 
        out_k[t]._val = (in[t]._val + RealType(shift[t])); 
        out_k[t]._val-= int(out_k[t]._val/(2.0*PI))*2.0*PI; 
        out_k[t]._index = std::get<1>(_kgrid.find(out_k[t]._val)); 
        };
    return out_k;
};



template <class Solver, size_t D, size_t ksize>
template <typename ...KP>
typename DFLadder<Solver,D,ksize>::GLocalType DFLadder<Solver,D,ksize>::getBubble(BMatsubaraGrid::point W, KP... q) const
{

    GLocalType out(this->_fGrid);
    GKType GD_shifted(this->GD0); // G(w+W,k+Q)


    typename GKType::PointFunctionType ShiftFunction = [&](FMatsubaraGrid::point w1, KP...k)->ComplexType { 
        //DEBUG(w1);
        int n=_bGrid.getNumber(W);
        FMatsubaraGrid::point w2; w2._index=0;
        w2._val=w1._val+W._val;
        if (-n<int(w1)) w2._index=w1._index+n;
        std::array<KMesh::point, D> K = {{ k... }};
        std::array<KMesh::point, D> Q = {{ q... }};
        std::array<KMesh::point, D> KpQ = _shift_point(K,Q); 
        //DEBUG(kx << ";" << ky << "||" << KpQ[0] << ";" << KpQ[1]);
        return this->GD0(std::tuple_cat(std::forward_as_tuple(w2), KpQ));
        };

    GD_shifted.fill(ShiftFunction);
    
    //DEBUG(W);
    //DEBUG(GD0);
    //DEBUG(GD_shifted);

    GD_shifted*=GD0;
    
    for (auto iw: _fGrid.getVals()) { 
        size_t iwn = size_t(iw);
        out[iwn] = GD_shifted[iwn].sum()/RealType(__power<ksize,D>::value);
    }

    RealType T = 1.0/(this->_fGrid._beta);
    out*=(-T);
    return out;

}


//
// DFLadder2d
//


template <class Solver, size_t ksize>
typename DFLadder2d<Solver,ksize>::GLocalType DFLadder2d<Solver,ksize>::operator()() 
{
    INFO("Using DF Ladder self-consistency in 2 dimensions on a cubic lattice of " << ksize << "^2" <<" atoms.");
    SigmaD = 0.0;
    RealType beta = _fGrid._beta;
    GLocalType gw(_fGrid); 
    gw = _S.gw;
    GLocalType Delta(_fGrid); 
    GLocalType GDLoc(_fGrid); 
    GLocalType GLatLoc(_fGrid); 
    Delta = _S.Delta;
    GLocalType Delta_out(_fGrid); Delta_out=0.0;
    auto wkgrids = std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kgrid));
    GKType GLatDMFT(wkgrids);
    GKType GD(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kgrid)));
    GKType GLat(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kgrid)));

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
    GDsum.savetxt("GDsum.dat");
//    DEBUG(GD0(FMatsubara(_fGrid._w_max-2,_fGrid._beta),0,0));
//    DEBUG(GD0(FMatsubara(_fGrid._w_max-1,_fGrid._beta),0,0));
//    DEBUG(GD0(FMatsubara(_fGrid._w_max,_fGrid._beta),0,0));
//    DEBUG(GD0(FMatsubara(_fGrid._w_max+1,_fGrid._beta),0,0));
    
    GLocalType Vertex4(_fGrid);
    GLocalType Chi0(_fGrid);
    INFO("Evaluating Bethe-Salpeter equation");
    size_t niter = 10;
    for (auto iW : _bGrid.getVals()) {
        DEBUG("iW = " << iW);
        std::function<ComplexType(FMatsubaraGrid::point)> f1 = std::bind(&FKImpuritySolver::getVertex4<FMatsubaraGrid::point, BMatsubaraGrid::point>, std::cref(_S), std::placeholders::_1, iW);
        Vertex4.fill(f1);
        //DEBUG(Vertex4);
        for (auto qx : _qGrids[0].getVals()) { 
            for (auto qy : _qGrids[1].getVals()) { 
                GLocalType IrrVertex4_old(Vertex4);
                GLocalType IrrVertex4_new(Vertex4);
                //DEBUG(this->getBubble(iW,_qGrid[0],_qGrid[0]));
                Chi0 = this->getBubble(iW,qx,qy);
                //typename GLocalType::PointFunctionType
                GridObject<RealType,FMatsubaraGrid> EVCheck(_fGrid); 
                std::function<RealType(FMatsubaraGrid::point)> f1 = [&](FMatsubaraGrid::point w)->RealType{return std::abs(Chi0(w)*Vertex4(w)); };
                EVCheck.fill(f1);
                RealType max_ev = *std::max_element(EVCheck.getData().begin(), EVCheck.getData().end());
                INFO("Maximum EV of Chi0*gamma = " << max_ev);
                for (size_t n=0; n<niter; ++n) { 
                        INFO("BS iteration " << n << " for iW = " << ComplexType(iW) << ", (qx,qy) = (" << RealType(qx) << "," << RealType(qy) << ").");
                        IrrVertex4_new = Vertex4 + Vertex4*Chi0*IrrVertex4_old;

                        auto diffV = IrrVertex4_new - IrrVertex4_old;
                        auto diff = std::real(_fGrid.integrate(diffV.conj()*diffV));
                        INFO("vertex diff = " << diff);
                        IrrVertex4_old = IrrVertex4_new;
                    }
                auto IrrVertex4_direct = Vertex4/(1.0 - Chi0 * Vertex4);
                DEBUG(IrrVertex4_direct);
                DEBUG(IrrVertex4_new - IrrVertex4_direct);
                //DEBUG(Chi0);
                }
            }
        };

    for (auto iw : _fGrid.getVals()) {
        size_t iwn = size_t(iw);
        GD[iwn] = 1.0/(1.0/GD0[iwn] - SigmaD[iwn]);
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

} // end of namespace FK
#endif // endif :: #ifndef __FK_DF_H__
