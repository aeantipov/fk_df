#ifndef ___FK_DMFT_H___
#define ___FK_DMFT_H___

#include <KMesh.hpp>
#include <RealGrid.hpp>
#include <GridObject.hpp>
#include "Solver.h"
#include "LatticeTraits.hpp"

namespace FK {

using namespace GFTools;

struct DMFTBase 
{
    typedef typename FKImpuritySolver::GFType GFType;
    const FKImpuritySolver &_S;
    DMFTBase(const FKImpuritySolver &S):_S(S){};
    virtual GFType operator()() = 0;
    template <typename MPoint> GFType getBubblePI(MPoint in) const;
    template <typename MPoint> GFType getBubble0(MPoint in) const;
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> getBubblePI() const;
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> getBubble0() const;
};

struct BetheSC : public DMFTBase
{
    typedef typename FKImpuritySolver::GFType GFType;
    const RealType _t;
    BetheSC(const FKImpuritySolver &S, RealType t);
    GFType operator()();
};

template <typename LatticeT, size_t D> struct LatticeDMFTSC : public DMFTBase
{
    using DMFTBase::_S;
    typedef typename FKImpuritySolver::GFType GFType;
    typedef LatticeT lattice_traits;
    //typedef Eigen::Array<ComplexType,__power<ksize,D>::value,1,Eigen::ColMajor> EkStorage;
    typedef typename ArgBackGenerator<D,KMesh,GridObject,ComplexType>::type EkStorage;
    typedef typename ArgBackGenerator<D,KMesh,GridObject,ComplexType,FMatsubaraGrid>::type GKType; // G(w,kx...)
    //EkStorage _ek_vals;
protected:
    typedef typename ArgFunGenerator<D,GFType,RealType>::type ArgFunType;
public:
    static const size_t NDim = D;
    const RealType _t;
    const KMesh _kGrid;
    mutable EkStorage _ek;
    GFType _gloc;

    LatticeDMFTSC(const FKImpuritySolver &S, KMesh kGrid, RealType t);
    template <typename ...ArgTypes> RealType dispersion(ArgTypes... kpoints) const;
    template <typename ...ArgTypes> RealType dispersion(const std::tuple<ArgTypes...>& kpoints) const;
    template <typename ...ArgTypes> GFType glat(ArgTypes... kpoints) const;
    template <typename MPoint, typename ...ArgTypes> ComplexType glat_analytic(MPoint w, ArgTypes... kpoints) const;
    template <typename MPoint, typename ...ArgTypes> ComplexType glat_analytic(std::tuple<MPoint,ArgTypes...> in) const;
    GKType getGLat(const FMatsubaraGrid& in) const;
    GFType operator()();
    
    template <typename MPoint> GFType getBubble0(MPoint in) const;
    template <typename MPoint> GFType getBubblePI(MPoint in) const;
    template <typename MPoint, typename KPoint> GFType getBubble(MPoint in, std::array<KPoint,D> q) const;
};

template <size_t D>
using CubicDMFTSC = LatticeDMFTSC<CubicTraits<D>,D>;


struct CubicInfDMFTSC : public DMFTBase
{
    typedef typename FKImpuritySolver::GFType GFType;
    typedef GridObject<ComplexType, RealGrid> ComplW;
    using DMFTBase::_S;
public:
    const RealType _t;
    const RealGrid _realgrid;
    ComplW _nominator;

    CubicInfDMFTSC(const FKImpuritySolver &S, RealType t, const RealGrid& realgrid);
    GFType operator()();
    using DMFTBase::getBubblePI;
    template <typename MPoint> GFType getBubble0(MPoint in) const;
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> getBubble0() const;
};

/** Extra tools. */
template <typename SolverType>
std::vector<RealType> getStaticLatticeDMFTSusceptibility(const SolverType& Solver, const std::vector<GFWrap>& Bubbles_in, const FMatsubaraGrid& gridF);
template <typename SolverType>
RealType getStaticLatticeDMFTSusceptibility(const SolverType& Solver, const GFWrap& Bubble_in, const FMatsubaraGrid& gridF);
template <typename SolverType>
std::vector<std::array<RealType,3>> getStaticLatticeDMFTSkeletonSusceptibility(const SolverType& Solver, const std::vector<GFWrap>& Bubbles_in, const FMatsubaraGrid& gridF);
template <typename SolverType>
std::array<RealType,3> getStaticLatticeDMFTSkeletonSusceptibility(const SolverType& Solver, const GFWrap& Bubble_in, const FMatsubaraGrid& gridF);

} // end of namespace FK
#endif // endif :: ifndef ___FK_DMFT_H___

#include "DMFT.hpp"

