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
    using CubicDMFTSC<Solver,D,ksize>::_t;
    using CubicDMFTSC<Solver,D,ksize>::_kgrid;
    using CubicDMFTSC<Solver,D,ksize>::_ek;
    const FMatsubaraGrid _fGrid;
    const BMatsubaraGrid _bGrid;
    const KMesh _qGrid;
    GKType GD0;
    GKType SigmaD;
protected:
    static std::array<KMesh::point, D> _shift_point(const std::array<KMesh::point, D> &in, const std::array<KMesh::point, D> &shift);
    template <typename wtype, typename ... kpoints> std::function<ComplexType(wtype,kpoints...)> _ektoiw(){ return [&](wtype w, kpoints ... k){return _ek(k...)/w;};};

public:
    template <typename ...KP> GLocalType getBubble(BMatsubaraGrid::point W, KP...kpoints) const;
    GLocalType getBubble(BMatsubaraGrid::point W, const std::array<KMesh::point,D>& q) const;
    DFLadder(const Solver &S, const FMatsubaraGrid& fGrid, const BMatsubaraGrid& bGrid, RealType t);
    GLocalType operator ()();

};

} // end of namespace FK

#include "DF.hpp"

#endif // endif :: #ifndef __FK_DF_H__
