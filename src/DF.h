#pragma once
#ifndef __FK_DF_H__
#define __FK_DF_H__

#include "DMFT.h"
#include "Diagrams.h"

namespace FK {

using namespace GFTools;

struct DFBase
{
    typedef typename FKImpuritySolver::GFType GLocalType;
    RealType _GDmix = 1.0;
    RealType _SC_cutoff = 1e-7;
    size_t _n_GD_iter = 200;

    bool _eval_BS_SC = false;
    RealType _BSmix = 0.1;
    size_t _n_BS_iter = 10;

    bool _EvaluateStaticDiagrams = true;
    bool _EvaluateDynamicDiagrams = false;

    bool _update_mixing = true;

    virtual GLocalType getGLoc() = 0;
};

template <typename LatticeT, size_t D>
struct DFLadder : LatticeDMFTSC<LatticeT>, DFBase {

    using LatticeDMFTSC<LatticeT>::_D;
    using typename LatticeDMFTSC<LatticeT>::lattice_traits;
    using LatticeDMFTSC<LatticeT>::lattice;
    using typename LatticeDMFTSC<LatticeT>::EkStorage;
    using typename DFBase::GLocalType;
    typedef typename ArgBackGenerator<D,KMesh,GridObject,ComplexType,FMatsubaraGrid>::type GKType;
    typedef typename ArgBackGenerator<D,KMesh,GridObject,ComplexType,FMatsubaraGrid,FMatsubaraGrid>::type FullVertexType;
    typedef typename GKType::ArgTupleType wkArgTupleType;
    typedef typename GKType::PointTupleType wkPointTupleType;
    typedef decltype(std::tuple_cat(std::make_tuple(BMatsubaraGrid::point()), std::array<KMesh::point, D>())) WQTupleType; 
    typedef typename ArgBackGenerator<D,KMeshPatch,GridObject,ComplexType,BMatsubaraGrid>::type SuscType;
    using LatticeDMFTSC<LatticeT>::_S;
    using LatticeDMFTSC<LatticeT>::_kGrid;
    using LatticeDMFTSC<LatticeT>::_ek;
    const FMatsubaraGrid _fGrid;
    GKType GD0;
    GKType GD;
    GKType SigmaD;
    GKType GLat;
    GLocalType GLatLoc; 
private:
    void _initialize();
public:
    template <typename ...LatticeParams> 
        DFLadder(const FKImpuritySolver &S, const FMatsubaraGrid& fGrid, KMesh kGrid, LatticeParams ... lattice_p);
    template <typename ...KP> GLocalType getBubble(const typename DFLadder<LatticeT,D>::GKType& GF, BMatsubaraGrid::point W, KP...kpoints) const;
    GLocalType getBubble(const GKType& GF, const WQTupleType& in) const;
    GKType getGLatDMFT(const FMatsubaraGrid& gridF) const ;
    GKType getGLat() const ;
    GKType getGLat(const FMatsubaraGrid &gridF ) const;
    //template <typename ...KP> ComplexType getBubble2(BMatsubaraGrid::point W, KP...kpoints, FMatsubaraGrid::point w1) const;
    GLocalType operator()();
    std::tuple<SuscType> calculateLatticeData(const BMatsubaraGrid& gridB);
    std::tuple<SuscType> calculateLatticeData(const BMatsubaraGrid& gridB, const std::array<KMeshPatch, D>& kpoints);
    template <typename KPoint> std::vector<ComplexType> getStaticLatticeSusceptibility(const std::vector<std::array<KPoint, D>>& q, const FMatsubaraGrid& fGrid);
    template <typename KPoint> ComplexType getStaticLatticeSusceptibility(const std::array<KPoint, D>& q, const FMatsubaraGrid& fGrid);
    template <typename KPoint> ComplexType getStaticLatticeSusceptibility(const std::array<KPoint, D>& q);
    GLocalType getGLoc();
    struct exRuntimeError : public std::runtime_error { exRuntimeError(const std::string &s):std::runtime_error(s){};};
};

template <size_t D>
using DFLadderCubic = DFLadder<CubicTraits<D>,D>;

} // end of namespace FK

#include "DF.hpp"

#endif // endif :: #ifndef __FK_DF_H__
