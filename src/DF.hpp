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
typename DFLadder<Solver,D,ksize>::GLocalType DFLadder<Solver,D,ksize>::operator()() 
{
    INFO("Using DF Ladder self-consistency in " << D << " dimensions on a cubic lattice of " << ksize << "^" << D <<" atoms.");
    SigmaD = 0.0;
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
    for (auto iw : _fGrid.getVals()) {
        size_t iwn = size_t(iw);
        //GLatDMFT[iwn] = 1.0/(1.0/_S.gw(iw)+_S.Delta(iw)-_ek.getData());
        GLatDMFT[iwn] = 1.0/(1.0/gw(iw)+Delta(iw)-_ek.getData());
        //SigmaD[size_t(iw)] = -gw(iw)/( 1.0/(Delta(iw) - _ek.getData()) + gw(iw))*gw(iw);
        GD0[iwn] = GLatDMFT[iwn] - gw(iw);
    };
    
    typename GridPointTypeExtractor<ComplexType, decltype(wkgrids)>::type t1;
    t1 = [](ComplexType w, RealType k1, RealType k2){return 1.0/w + k1+k2;};
    GLocalType Vertex4(_fGrid);
    GKType Chi0(wkgrids);
    Chi0 = t1;
    DEBUG(Chi0);
    exit(0);
    GKType Chi0Gamma(wkgrids);
    INFO("Evaluating Bethe-Salpeter equation");
    for (auto iW : _bGrid.getVals()) {
        std::function<ComplexType(FMatsubaraGrid::point)> f1 = std::bind(&FKImpuritySolver::getVertex4<FMatsubaraGrid::point, BMatsubaraGrid::point>, std::ref(_S), std::placeholders::_1, iW);
        Vertex4.fill(f1);
        DEBUG(Vertex4);
        exit(0);
        for (auto q : _qGrid.getVals()) {
            //GLocalType Vertex = 
            //DEBUG("iW = " << iW << ", q = " << q);
            
        }
    }

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
