#ifndef __FK_DF_HPP__
#define __FK_DF_HPP__

#include "DF.h"

namespace FK {

template <class Solver, size_t D, size_t ksize>
DFLadder<Solver,D,ksize>::DFLadder(const Solver &S, const FMatsubaraGrid& fGrid, const BMatsubaraGrid& bGrid, RealType t):
    CubicDMFTSC<Solver,D,ksize>(S,t),
    _fGrid(fGrid),
    _bGrid(bGrid),
    GD0(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kgrid))),
    SigmaD(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kgrid)))
{
    //auto d3 = CubicTraits<2,ksize>::getTuples(_kgrid);
};

template <class Solver, size_t D, size_t ksize>
typename DFLadder<Solver,D,ksize>::GLocalType DFLadder<Solver,D,ksize>::operator()() 
{
    SigmaD = 0.0;
    GLocalType gw(_fGrid); 
    gw = _S.gw;
    GLocalType Delta(_fGrid); 
    GLocalType GDLoc(_fGrid); 
    GLocalType GLatLoc(_fGrid); 
    Delta = _S.Delta;
    GLocalType Delta_out(_fGrid); Delta_out=0.0;
    GKType GLatDMFT(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kgrid)));
    GKType GD(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kgrid)));
    GKType GLat(std::tuple_cat(std::make_tuple(_fGrid),CubicTraits<D,ksize>::getTuples(_kgrid)));
    for (auto iw : _fGrid.getVals()) {
        size_t iwn = size_t(iw);
        DEBUG(iw);
        //GLatDMFT[iwn] = 1.0/(1.0/_S.gw(iw)+_S.Delta(iw)-_ek.getData());
        GLatDMFT[iwn] = 1.0/(1.0/gw(iw)+Delta(iw)-_ek.getData());
        //SigmaD[size_t(iw)] = -gw(iw)/( 1.0/(Delta(iw) - _ek.getData()) + gw(iw))*gw(iw);
        GD0[iwn] = GLatDMFT[iwn] - gw(iw);

        GD[iwn] = GD0[iwn];
        GLat[iwn] = GLatDMFT[iwn];

        GDLoc[iwn] = GD[iwn].sum()/RealType(__power<ksize,D>::value);
        GLatLoc[iwn] = GLat[iwn].sum()/RealType(__power<ksize,D>::value);
    }
//    DEBUG(GD0);
//    DEBUG(GDLoc);
//    DEBUG(GD0-SigmaD);

    Delta_out = Delta + 1.0/gw * GDLoc / GLatLoc;
    DEBUG(Delta_out);
    //exit(0);
    return Delta_out;
}

} // end of namespace FK
#endif // endif :: #ifndef __FK_DF_H__
