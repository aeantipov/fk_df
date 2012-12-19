#ifndef __FK_DF_H__
#define __FK_DF_H__

#include "FKCommon.h"
#include "SelfConsistency.h"

namespace FK {

template <class Solver, size_t D, size_t ksize=32>
struct DFLadder : CubicDMFTSC<Solver,D,ksize> {
    using typename CubicDMFTSC<Solver,D,ksize>::EkStorage;
    typedef typename Solver::GFType GLocalType;
    typedef typename ArgBackGenerator<D,KMesh,GridObject,ComplexType,FMatsubaraGrid>::type GKType;
    using CubicDMFTSC<Solver,D,ksize>::_S;
    using CubicDMFTSC<Solver,D,ksize>::_kgrid;
    using CubicDMFTSC<Solver,D,ksize>::_ek;
    const FMatsubaraGrid _fGrid;
    const BMatsubaraGrid _bGrid;
    GKType GD0;
    GKType SigmaD;
private:
    GKType& _gen_gdual();
public:
    DFLadder(const Solver &S, const FMatsubaraGrid& fGrid, const BMatsubaraGrid& bGrid, RealType t);
    GLocalType operator ()();

};

} // end of namespace FK

#include "DF.hpp"

#endif // endif :: #ifndef __FK_DF_H__
