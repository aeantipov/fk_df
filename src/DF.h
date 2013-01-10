#ifndef __FK_DF_H__
#define __FK_DF_H__

#include "FKCommon.h"
#include "SelfConsistency.h"

namespace FK {

template <class Solver, size_t D, size_t ksize=32>
struct DFLadder : CubicDMFTSC<Solver,D,ksize> {

    using typename CubicDMFTSC<Solver,D,ksize>::EkStorage;
    typedef typename Solver::GFType GLocalType;
    typedef KMeshPatch qGridType;
    typedef typename ArgBackGenerator<D,KMesh,GridObject,ComplexType,FMatsubaraGrid>::type GKType;
    typedef decltype(std::tuple_cat(std::make_tuple(BMatsubaraGrid::point()), std::array<KMesh::point, D>())) WQTupleType; 
    using CubicDMFTSC<Solver,D,ksize>::_S;
    using CubicDMFTSC<Solver,D,ksize>::_t;
    using CubicDMFTSC<Solver,D,ksize>::_kGrid;
    using CubicDMFTSC<Solver,D,ksize>::_ek;
    const FMatsubaraGrid _fGrid;
    const BMatsubaraGrid _bGrid;
    const std::array<qGridType,D> _qGrids;
    GKType GD0;
    GKType GD;
    GKType SigmaD;
    GKType GLat;
    RealType _GDmix = 1.0;
    size_t _n_GD_iter = 100;
    RealType _BSmix = 1.0;
    size_t _n_BS_iter = 100;
public:
    DFLadder(const Solver &S, const FMatsubaraGrid& fGrid, const BMatsubaraGrid& bGrid, const std::array<qGridType,D>& qGrids, RealType t);
    template <typename ...KP> GLocalType getBubble(BMatsubaraGrid::point W, KP...kpoints) const;
    GKType getGLatDMFT() const { return CubicDMFTSC<Solver,D,ksize>::getGLat(_fGrid); };
    GKType getGLat() const { return GLat; };
    GLocalType getBubble(const WQTupleType& in) const;
    //template <typename ...KP> ComplexType getBubble2(BMatsubaraGrid::point W, KP...kpoints, FMatsubaraGrid::point w1) const;
    GLocalType operator()(bool eval_BS_SC = false);
};

} // end of namespace FK

#include "DF.hpp"

#endif // endif :: #ifndef __FK_DF_H__
