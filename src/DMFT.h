#ifndef ___FK_DMFT_H___
#define ___FK_DMFT_H___

#include <gftools/kmesh.hpp>
#include <gftools/real_grid.hpp>
#include "Solver.h"
#include "LatticeTraits.hpp"

namespace FK {

using namespace gftools;

struct DMFTBase 
{
    typedef typename FKImpuritySolver::GFType GFType;
    const FKImpuritySolver &_S;
    DMFTBase(const FKImpuritySolver &S):_S(S){};
    virtual GFType operator()() = 0;
    template <typename MPoint> GFType getBubblePI(MPoint in) const;
    template <typename MPoint> GFType getBubble0(MPoint in) const;
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> getBubblePI() const;
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> getBubble0() const;
};

struct BetheSC : public DMFTBase
{
    typedef typename FKImpuritySolver::GFType GFType;
    const real_type _t;
    BetheSC(const FKImpuritySolver &S, real_type t);
    GFType operator()();
};

template <typename LatticeT> struct LatticeDMFTSC : public DMFTBase
{
    using DMFTBase::_S;
    typedef typename FKImpuritySolver::GFType GFType;
    typedef LatticeT lattice_traits;
    static const size_t _D = lattice_traits::_D;
    //typedef Eigen::Array<complex_type,__power<ksize,_D>::value,1,Eigen::ColMajor> EkStorage;
    typedef typename tools::ArgBackGenerator<_D,kmesh,grid_object,complex_type>::type EkStorage;
    typedef typename tools::ArgBackGenerator<_D,kmesh,grid_object,complex_type,fmatsubara_grid>::type GKType; // G(w,kx...)
    //EkStorage _ek_vals;
protected:
    typedef typename tools::ArgFunGenerator<_D,GFType,real_type>::type ArgFunType;
public:
    lattice_traits lattice;
    const kmesh _kGrid;
    mutable EkStorage _ek;
    GFType _gloc;

    template <typename ...LatticeParams> 
        LatticeDMFTSC(const FKImpuritySolver &S, kmesh kGrid, LatticeParams ... ls);
    template <typename ...ArgTypes> real_type dispersion(ArgTypes... kpoints) const;
    template <typename ...ArgTypes> real_type dispersion(const std::tuple<ArgTypes...>& kpoints) const;
    template <typename ...ArgTypes> GFType glat(ArgTypes... kpoints) const;
    template <typename MPoint, typename ...ArgTypes> complex_type glat_analytic(MPoint w, ArgTypes... kpoints) const;
    template <typename MPoint, typename ...ArgTypes> complex_type glat_analytic(std::tuple<MPoint,ArgTypes...> in) const;
    GKType getGLat(const fmatsubara_grid& in) const;
    GFType operator()();
    
    template <typename MPoint> GFType getBubble0(MPoint in) const;
    template <typename MPoint> GFType getBubblePI(MPoint in) const;
    template <typename MPoint, typename KPoint> GFType getBubble(MPoint in, std::array<KPoint,_D> q) const;
};

template <size_t D>
using CubicDMFTSC = LatticeDMFTSC<CubicTraits<D>>;
typedef LatticeDMFTSC<TriangularTraits> TriangularDMFT;


struct CubicInfDMFTSC : public DMFTBase
{
    typedef typename FKImpuritySolver::GFType GFType;
    typedef grid_object<complex_type, real_grid> ComplW;
    using DMFTBase::_S;
public:
    const real_type _t;
    const real_grid _realgrid;
    ComplW _nominator;

    CubicInfDMFTSC(const FKImpuritySolver &S, real_type t, const real_grid& realgrid);
    GFType operator()();
    using DMFTBase::getBubblePI;
    template <typename MPoint> GFType getBubble0(MPoint in) const;
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> getBubble0() const;
};

/** Extra tools. */
template <typename SolverType>
std::vector<real_type> getStaticLatticeDMFTSusceptibility(const SolverType& Solver, const std::vector<typename SolverType::GFType>& Bubbles_in, const fmatsubara_grid& gridF);
template <typename SolverType>
real_type getStaticLatticeDMFTSusceptibility(const SolverType& Solver, const typename SolverType::GFType& Bubble_in, const fmatsubara_grid& gridF);
template <typename SolverType>
std::vector<std::array<real_type,3>> getStaticLatticeDMFTSkeletonSusceptibility(const SolverType& Solver, const std::vector<typename SolverType::GFType>& Bubbles_in, const fmatsubara_grid& gridF);
template <typename SolverType>
std::array<real_type,3> getStaticLatticeDMFTSkeletonSusceptibility(const SolverType& Solver, const typename SolverType::GFType& Bubble_in, const fmatsubara_grid& gridF);

} // end of namespace FK
#endif // endif :: ifndef ___FK_DMFT_H___

#include "DMFT.hpp"

