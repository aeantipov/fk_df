#ifndef __FK_DF_H__
#define __FK_DF_H__

#include "FKCommon.h"
#include "SelfConsistency.h"
#include "DFDiagrams.h"

namespace FK {

struct DFBase
{
    RealType _GDmix = 1.0;
    size_t _n_GD_iter = 100;

    bool _eval_BS_SC = false;
    RealType _BSmix = 1.0;
    size_t _n_BS_iter = 100;

    virtual void calculateLatticeData(const BMatsubaraGrid& gridB)=0;
};

template <class Solver, size_t D, size_t ksize=32>
struct DFLadder : CubicDMFTSC<Solver,D,ksize>, DFBase {

    using CubicDMFTSC<Solver,D,ksize>::NDim;
    using typename CubicDMFTSC<Solver,D,ksize>::EkStorage;
    typedef typename Solver::GFType GLocalType;
    typedef typename ArgBackGenerator<D,KMesh,GridObject,ComplexType,FMatsubaraGrid>::type GKType;
    typedef typename GKType::ArgTupleType wkTupleType;
    typedef decltype(std::tuple_cat(std::make_tuple(BMatsubaraGrid::point()), std::array<KMesh::point, D>())) WQTupleType; 
    typedef typename ArgBackGenerator<D,KMesh,GridObject,ComplexType,BMatsubaraGrid>::type SuscType;
    using CubicDMFTSC<Solver,D,ksize>::_S;
    using CubicDMFTSC<Solver,D,ksize>::_t;
    using CubicDMFTSC<Solver,D,ksize>::_kGrid;
    using CubicDMFTSC<Solver,D,ksize>::_ek;
    const FMatsubaraGrid _fGrid;
    const BMatsubaraGrid _bGrid;
    const DFDiagrams<D> Diagrams;
    GKType GD0;
    GKType GD;
    GKType SigmaD;
    GKType GLat;
    GLocalType GLatLoc; 
    SuscType LatticeSusc;
private:
    void _initialize();
public:
    DFLadder(const Solver &S, const FMatsubaraGrid& fGrid, const BMatsubaraGrid& bGrid, RealType t);
    template <typename ...KP> GLocalType getBubble(const typename DFLadder<Solver,D,ksize>::GKType& GF, BMatsubaraGrid::point W, KP...kpoints) const;
    GLocalType getBubble(const GKType& GF, const WQTupleType& in) const;
    GKType getGLatDMFT() const { return CubicDMFTSC<Solver,D,ksize>::getGLat(_fGrid); };
    GKType getGLat() const { return GLat; };
    GKType getGLat(const FMatsubaraGrid &gridF ) const;
    //template <typename ...KP> ComplexType getBubble2(BMatsubaraGrid::point W, KP...kpoints, FMatsubaraGrid::point w1) const;
    GLocalType operator()();
    void calculateLatticeData(const BMatsubaraGrid& gridB);
};

} // end of namespace FK

#include "DF.hpp"

#endif // endif :: #ifndef __FK_DF_H__
