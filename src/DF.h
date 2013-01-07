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
    typedef KMeshPatch qGridType;
    using CubicDMFTSC<Solver,D,ksize>::_S;
    using CubicDMFTSC<Solver,D,ksize>::_t;
    using CubicDMFTSC<Solver,D,ksize>::_kgrid;
    using CubicDMFTSC<Solver,D,ksize>::_ek;
    const FMatsubaraGrid _fGrid;
    const BMatsubaraGrid _bGrid;
    const std::array<qGridType,D> _qGrids;
    GKType GD0;
    GKType SigmaD;
protected:
    template <typename PointType>
        std::array<KMesh::point, D> _shift_point(const std::array<KMesh::point, D> &in, const std::array<PointType, D> &shift) const;
    //template <typename funType> funType _ektoiw();
    //template <typename wtype, typename ... kpoints> std::function<ComplexType(ComplexType,kpoints...)> _ektoiw<std::function<ComplexType(ComplexType,kpoints...)>>(){ return [&](ComplexType w, kpoints ... k){return _ek(k...)/w;};};
    DFLadder(const Solver &S, const FMatsubaraGrid& fGrid, const BMatsubaraGrid& bGrid, const std::array<qGridType,D>& qGrids, RealType t);

public:
    template <typename ...KP> GLocalType getBubble(BMatsubaraGrid::point W, KP...kpoints) const;
    GLocalType operator ()();
};

template <class Solver, size_t ksize>
struct DFLadder2d : DFLadder<Solver, 2, ksize>
{
    static const int D=2;
    using typename DFLadder<Solver,2,ksize>::EkStorage;
    using typename DFLadder<Solver,2,ksize>::GLocalType;
    using typename DFLadder<Solver,2,ksize>::GKType;
    using typename DFLadder<Solver,2,ksize>::qGridType;
    using DFLadder<Solver,2,ksize>::_S;
    using DFLadder<Solver,2,ksize>::_t;
    using DFLadder<Solver,2,ksize>::_kgrid;
    using DFLadder<Solver,2,ksize>::_ek;
    using DFLadder<Solver,2,ksize>::_fGrid;
    using DFLadder<Solver,2,ksize>::_bGrid;
    using DFLadder<Solver,2,ksize>::_qGrids;
    using DFLadder<Solver,2,ksize>::_shift_point;
    
    using DFLadder<Solver,2,ksize>::GD0;
    using DFLadder<Solver,2,ksize>::SigmaD;

    DFLadder2d(const Solver &S, const FMatsubaraGrid& fGrid, const BMatsubaraGrid& bGrid, const std::array<qGridType,2>& qGrids, RealType t):
        DFLadder<Solver,2,ksize>(S,fGrid,bGrid,qGrids,t){};
    //GLocalType getBubble(BMatsubaraGrid::point W, KMesh::point q1, KMesh::point q2) const;
    GLocalType operator ()();
    
};

} // end of namespace FK

#include "DF.hpp"

#endif // endif :: #ifndef __FK_DF_H__
