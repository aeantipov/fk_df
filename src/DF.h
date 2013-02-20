#pragma once
#ifndef __FK_DF_H__
#define __FK_DF_H__

#include "SelfConsistency.h"
#include "Diagrams.h"

namespace FK {

using namespace GFTools;

struct DFBase
{
    typedef typename FKImpuritySolver::GFType GLocalType;
    RealType _GDmix = 0.5;
    RealType _SC_cutoff = 1e-8;
    size_t _n_GD_iter = 100;

    bool _eval_BS_SC = false;
    RealType _BSmix = 1.0;
    size_t _n_BS_iter = 1000;

    bool _EvaluateStaticDiagrams = true;
    bool _EvaluateDynamicDiagrams = true;

    virtual GLocalType getGLoc() = 0;
};

template <size_t D>
struct DFLadder : CubicDMFTSC<D>, DFBase {

    using CubicDMFTSC<D>::NDim;
    using typename CubicDMFTSC<D>::EkStorage;
    using typename DFBase::GLocalType;
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
    GKType getGLatDMFT() const ;
    GKType getGLat() const ;
    GKType getGLat(const FMatsubaraGrid &gridF ) const;
    //template <typename ...KP> ComplexType getBubble2(BMatsubaraGrid::point W, KP...kpoints, FMatsubaraGrid::point w1) const;
    GLocalType operator()();
    std::tuple<SuscType> calculateLatticeData(const BMatsubaraGrid& gridB);
    std::tuple<SuscType> calculateLatticeData(const BMatsubaraGrid& gridB, const std::array<KMeshPatch, D>& kpoints);
    template <typename KPoint> ComplexType getStaticLatticeSusceptibility(const std::array<KPoint, D>& q);
    GLocalType getGLoc();
};

} // end of namespace FK

#include "DF.hpp"

#endif // endif :: #ifndef __FK_DF_H__
