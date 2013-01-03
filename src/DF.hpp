#ifndef __FK_DF_HPP__
#define __FK_DF_HPP__

#include "DF.h"

namespace FK {

template <class Solver, size_t D, size_t ksize>
DFLadder<Solver,D,ksize>::DFLadder(const Solver &S, const FMatsubaraGrid& fGrid, const BMatsubaraGrid& bGrid, RealType t):
    CubicDMFTSC<Solver,D,ksize>(S,t),
    _fGrid(fGrid),
    _bGrid(bGrid),
    _qGrid(_kgrid),
    GD0(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kgrid))),
    SigmaD(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kgrid)))
{
    //auto d3 = CubicTraits<2,ksize>::getTuples(_kgrid);
};

template <class Solver, size_t D, size_t ksize>
std::array<KMesh::point, D> DFLadder<Solver,D,ksize>::_shift_point(const std::array<KMesh::point, D> &in, const std::array<KMesh::point, D> &shift)
{
    std::array<KMesh::point,D> out_k;
    for (size_t t=0; t<D; ++t) { 
        out_k[t]._index = (in[t]._index + shift[t]._index)%ksize; 
        out_k[t]._val = out_k[t]._index*2.0*PI/ksize; // %ksize;
        };
    return out_k;
};


template <class Solver, size_t D, size_t ksize>
template <typename ...KP>
typename DFLadder<Solver,D,ksize>::GLocalType DFLadder<Solver,D,ksize>::getBubble(BMatsubaraGrid::point W, KP... Q) const
{
    std::array<KMesh::point,D> arrQ = {{ Q...} };
    return this->getBubble(W,arrQ);
}

template <class Solver, size_t D, size_t ksize>
typename DFLadder<Solver,D,ksize>::GLocalType DFLadder<Solver,D,ksize>::getBubble(BMatsubaraGrid::point W, const std::array<KMesh::point,D>& Q) const
{
    GLocalType out(this->_fGrid);
    
    for (auto w1 : _fGrid.getVals()) {
        int n=_bGrid.getNumber(W);
        FMatsubaraGrid::point w2; w2._index=0;
        w2._val=w1._val+W._val;
        if (-n<int(w1)) w2._index=w1._index+n;
        //auto f2 = [&](FMatsubaraPoint::w1, KPoint
        //out[size_t(w1)]=GD0[size_t(w1)] * GD0(w)
        };
/*
    //typename GKType::FunctionType f1;
    auto shift_point = [&](const std::array<KMesh::point,D>& arrK)->std::array<KMesh::point,D>{ 
        std::array<KMesh::point,D> out_k;
        for (size_t t=0; t<D; ++t) { 
                out_k[t]._index = (arrK[t]._index + Q[t]._index)%ksize; 
                out_k[t]._val = out_k[t]._index*2.0*PI/ksize; // %ksize;
            };
        return out_k;
        };
        FK::additional::__caller<ComplexType,GKType,FMatsubara::point,KMesh::point,KMesh::point> t1(GD0, std::tuple_cat(std::make_tuple(w2),shift_point({{_kgrid[0],_kgrid[0]}})));
    }*/
    return out;
}

template <class Solver, size_t D, size_t ksize>
typename DFLadder<Solver,D,ksize>::GLocalType DFLadder<Solver,D,ksize>::operator()() 
{
    INFO("Using DF Ladder self-consistency in " << D << " dimensions on a cubic lattice of " << ksize << "^" << D <<" atoms.");
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

    DEBUG(Delta(FMatsubara(_fGrid._w_max-1, beta)));
    DEBUG(Delta(FMatsubara(_fGrid._w_max, beta)));

    for (auto iw : _fGrid.getVals()) {
        size_t iwn = size_t(iw);
        //GLatDMFT[iwn] = 1.0/(1.0/_S.gw(iw)+_S.Delta(iw)-_ek.getData());
        GLatDMFT[iwn] = 1.0/(1.0/gw(iw)+Delta(iw)-_ek.getData());
        //SigmaD[size_t(iw)] = -gw(iw)/( 1.0/(Delta(iw) - _ek.getData()) + gw(iw))*gw(iw);
        GD0[iwn] = GLatDMFT[iwn] - gw(iw);
        GDsum[iwn] = GD0[iwn].sum()/RealType(__power<ksize,D>::value);
    };

    if (D==2) {  // Change it later
        auto _f1 = [&](ComplexType w, RealType kx, RealType ky){return (_S.mu - _S.Sigma._f(w)-_ek(kx,ky))/std::abs(w*w)+1.0/w;};
       // auto _f12 = [&](std::tuple<ComplexType> w_t, std::array<RealType,D> kpoints)->ComplexType 
       //     { ComplexType w = std::get<0>(w_t);
       //       return (_S.mu - _S.Sigma._f(w)-_ek(kpoints))/std::abs(w*w)+1.0/w;};
        GLatDMFT._f = _f1;
        auto _f2 = [&](ComplexType w, RealType kx, RealType ky){
            //return (-_ek(kx,ky))/std::abs(w*w) + I*imag(Delta._f(w))/std::abs(w*w)-(_ek(kx,ky)*_ek(kx,ky) - 2.0*_ek(kx,ky)*(_S.mu - _S.Sigma._f(w)))/w/std::abs(w*w); };
            return (-_ek(kx,ky))/std::abs(w*w) + I*imag(Delta._f(w))/std::abs(w*w)-(_ek(kx,ky)*_ek(kx,ky) - 2.0*_ek(kx,ky)*(_S.mu - _S.Sigma._f(w)))/w/std::abs(w*w); };
        GD0._f = _f2;
        };
    DEBUG(GLatDMFT(FMatsubara(_fGrid._w_max-2, beta),PI/4.,PI/4.));
    DEBUG(GLatDMFT(FMatsubara(_fGrid._w_max-1, beta),PI/4.,PI/4.));
    DEBUG(GLatDMFT(FMatsubara(_fGrid._w_max,   beta),PI/4.,PI/4.0));
    DEBUG(GLatDMFT(FMatsubara(_fGrid._w_max+1,   beta),PI/4.,PI/4.0));
    DEBUG("-------");
    DEBUG(GD0(FMatsubara(_fGrid._w_max-2, beta),PI/4.,PI/4.));
    DEBUG(GD0(FMatsubara(_fGrid._w_max-1, beta),PI/4.,PI/4.));
    DEBUG(GD0(FMatsubara(_fGrid._w_max,   beta),PI/4.,PI/4.0));
    DEBUG(GD0(FMatsubara(_fGrid._w_max+1,   beta),PI/4.,PI/4.0));

    //auto fw3 = [this, &Delta](ComplexType w)->ComplexType{return Delta._f(w)/w/w;};
    //GDsum._f = fw3; 
    GDsum.savetxt("GDsum.dat");
    DEBUG(GD0(FMatsubara(_fGrid._w_max-2,_fGrid._beta),0,0));
    DEBUG(GD0(FMatsubara(_fGrid._w_max-1,_fGrid._beta),0,0));
    DEBUG(GD0(FMatsubara(_fGrid._w_max,_fGrid._beta),0,0));
    DEBUG(GD0(FMatsubara(_fGrid._w_max+1,_fGrid._beta),0,0));
    
    exit(0);

    //typename GridPointTypeExtractor<ComplexType, decltype(wkgrids)>::type t1;

    typename GKType::FunctionType t1;
    t1 = [](ComplexType w, RealType k1, RealType k2){return 1.0/w + k1+k2;};
    GLocalType Vertex4(_fGrid);
    GKType Chi0(wkgrids);
    //GKType Chi0Gamma(wkgrids);
    INFO("Evaluating Bethe-Salpeter equation");
    for (auto iW : _bGrid.getVals()) {
        std::function<ComplexType(FMatsubaraGrid::point)> f1 = std::bind(&FKImpuritySolver::getVertex4<FMatsubaraGrid::point, BMatsubaraGrid::point>, std::ref(_S), std::placeholders::_1, iW);
        Vertex4.fill(f1);
        DEBUG(Vertex4);
        DEBUG(this->getBubble(_bGrid[0],_qGrid[0],_qGrid[0]));
        exit(0);
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
