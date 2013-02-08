#pragma once
#ifndef __FK_DF_H__
#define __FK_DF_H__

#include "FKCommon.h"
#include "SelfConsistency.h"
#include "Diagrams.h"

namespace FK {

struct DFBase
{
    RealType _GDmix = 0.5;
    size_t _n_GD_iter = 100;

    bool _eval_BS_SC = false;
    RealType _BSmix = 1.0;
    size_t _n_BS_iter = 100;

    bool _EvaluateStaticDiagrams = true;
    bool _EvaluateDynamicDiagrams = true;
};

template <size_t D>
struct DFLadder : CubicDMFTSC<D>, DFBase {

    using CubicDMFTSC<D>::NDim;
    using typename CubicDMFTSC<D>::EkStorage;
    typedef typename FKImpuritySolver::GFType GLocalType;
    typedef typename ArgBackGenerator<D,KMesh,GridObject,ComplexType,FMatsubaraGrid>::type GKType;
    typedef typename GKType::ArgTupleType wkArgTupleType;
    typedef typename GKType::PointTupleType wkPointTupleType;
    typedef decltype(std::tuple_cat(std::make_tuple(BMatsubaraGrid::point()), std::array<KMesh::point, D>())) WQTupleType; 
    typedef typename ArgBackGenerator<D,KMeshPatch,GridObject,ComplexType,BMatsubaraGrid>::type SuscType;
    using CubicDMFTSC<D>::_S;
    using CubicDMFTSC<D>::_t;
    using CubicDMFTSC<D>::_kGrid;
    using CubicDMFTSC<D>::_ek;
    const FMatsubaraGrid _fGrid;
    const BMatsubaraGrid _bGrid;
    GKType GD0;
    GKType GD;
    GKType SigmaD;
    GKType GLat;
    GLocalType GLatLoc; 
private:
    void _initialize();
public:
    DFLadder(const FKImpuritySolver &S, const FMatsubaraGrid& fGrid, KMesh kGrid, const BMatsubaraGrid& bGrid, RealType t);
    template <typename ...KP> GLocalType getBubble(const typename DFLadder<D>::GKType& GF, BMatsubaraGrid::point W, KP...kpoints) const;
    GLocalType getBubble(const GKType& GF, const WQTupleType& in) const;
    GKType getGLatDMFT() const { return CubicDMFTSC<D>::getGLat(_fGrid); };
    GKType getGLat() const { return GLat; };
    GKType getGLat(const FMatsubaraGrid &gridF ) const;
    //template <typename ...KP> ComplexType getBubble2(BMatsubaraGrid::point W, KP...kpoints, FMatsubaraGrid::point w1) const;
    GLocalType operator()();
    std::tuple<SuscType> calculateLatticeData(const BMatsubaraGrid& gridB);
    std::tuple<SuscType> calculateLatticeData(const BMatsubaraGrid& gridB, const std::array<KMeshPatch, D>& kpoints);
};

} // end of namespace FK

#include "DF.hpp"

#endif // endif :: #ifndef __FK_DF_H__
