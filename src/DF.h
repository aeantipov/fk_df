#pragma once
#ifndef __FK_DF_H__
#define __FK_DF_H__

#include "DMFT.h"
#include "Diagrams.h"

namespace FK {

using namespace gftools;

struct DFBase
{
    typedef typename FKImpuritySolver::GFType GLocalType;
    real_type _GDmix = 1.0;
    real_type _SC_cutoff = 1e-7;
    size_t _n_GD_iter = 200;

    bool _eval_BS_SC = false;
    real_type _BSmix = 0.1;
    size_t _n_BS_iter = 10;
    bool _bs_evaluate_only_order_n = false;

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
    typedef typename tools::ArgBackGenerator<D,kmesh,grid_object,complex_type,fmatsubara_grid>::type GKType;
    typedef typename tools::ArgBackGenerator<D,kmesh,grid_object,complex_type,fmatsubara_grid,fmatsubara_grid>::type FullVertexType;
    typedef typename GKType::arg_tuple wkarg_tuple;
    typedef typename GKType::point_tuple wkpoint_tuple;
    typedef typename Diagrams::WQTupleType<D> WQTupleType;
    typedef typename tools::ArgBackGenerator<D,kmesh_patch,grid_object,complex_type,bmatsubara_grid>::type SuscType;
    using LatticeDMFTSC<LatticeT>::_S;
    using LatticeDMFTSC<LatticeT>::_kGrid;
    using LatticeDMFTSC<LatticeT>::_ek;
    const fmatsubara_grid _fGrid;
    GKType GD0;
    GKType GD;
    GKType SigmaD;
    GKType GLat;
    GLocalType GLatLoc; 
private:
    void _initialize();
public:
    template <typename ...LatticeParams> 
        DFLadder(const FKImpuritySolver &S, const fmatsubara_grid& fGrid, kmesh kGrid, LatticeParams ... lattice_p);
    template <typename ...KP> GLocalType getBubble(const typename DFLadder<LatticeT,D>::GKType& GF, bmatsubara_grid::point W, KP...kpoints) const;
    GLocalType getBubble(const GKType& GF, const WQTupleType& in) const;
    GKType getGLatDMFT(const fmatsubara_grid& gridF) const ;
    GKType getGLat() const ;
    GKType getGLat(const fmatsubara_grid &gridF ) const;
    //template <typename ...KP> complex_type getBubble2(bmatsubara_grid::point W, KP...kpoints, fmatsubara_grid::point w1) const;
    GLocalType operator()();
    std::tuple<SuscType> calculateLatticeData(const bmatsubara_grid& gridB);
    std::tuple<SuscType> calculateLatticeData(const bmatsubara_grid& gridB, const std::array<kmesh_patch, D>& kpoints);
    template <typename KPoint> std::vector<complex_type> getStaticLatticeSusceptibility(const std::vector<std::array<KPoint, D>>& q, const fmatsubara_grid& fGrid);
    template <typename KPoint> complex_type getStaticLatticeSusceptibility(const std::array<KPoint, D>& q, const fmatsubara_grid& fGrid);
    template <typename KPoint> complex_type getStaticLatticeSusceptibility(const std::array<KPoint, D>& q);
    GLocalType getGLoc();
    struct exRuntimeError : public std::runtime_error { exRuntimeError(const std::string &s):std::runtime_error(s){};};
};

template <size_t D>
using DFLadderCubic = DFLadder<CubicTraits<D>,D>;

} // end of namespace FK

#include "DF.hpp"

#endif // endif :: #ifndef __FK_DF_H__
