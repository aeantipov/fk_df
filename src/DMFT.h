#ifndef ___FK_DMFT_H___
#define ___FK_DMFT_H___

#include <KMesh.hpp>
#include <RealGrid.hpp>
#include <GridObject.hpp>
#include "Solver.h"

namespace FK {

using namespace GFTools;

struct DMFT 
{
    typedef typename FKImpuritySolver::GFType GFType;
    const FKImpuritySolver &_S;
    DMFT(const FKImpuritySolver &S):_S(S){};
    virtual GFType operator()() = 0;
    template <typename MPoint> GFType getBubblePI(MPoint in) const;
    template <typename MPoint> GFType getBubble0(MPoint in) const;
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> getBubblePI() const;
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> getBubble0() const;
};

struct BetheSC : public DMFT
{
    typedef typename FKImpuritySolver::GFType GFType;
    const RealType _t;
    BetheSC(const FKImpuritySolver &S, RealType t);
    GFType operator()();
};

template <size_t D> struct CubicDMFTSC : public DMFT
{
    using DMFT::_S;
    typedef typename FKImpuritySolver::GFType GFType;
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

    CubicDMFTSC(const FKImpuritySolver &S, RealType t, KMesh kGrid);
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

/** A typedef for a point in the Brillouin zone. */
template <size_t D>
using BZPoint = typename std::array<KMesh::point,D>;

/** Stream the BZPoint. */
template <size_t D>
std::ostream& operator<<(std::ostream& out, BZPoint<D> in)
{for (size_t i=0; i<D; ++i) out << RealType(in[i]) << " "; return out; } 
    
template <size_t M> 
struct CubicTraits{ 
    /** Returns an analytic std::function of the dispersion. */
    template <typename FunctionType, typename ...ArgTypes> static FunctionType get_dispersion(RealType t);
    /** Finds the equivalent point, which is used is calculations. */
    static std::array<RealType,M> findSymmetricBZPoint(const std::array<RealType,M>& in);
    static BZPoint<M> findSymmetricBZPoint(const BZPoint<M>& in, const KMesh& kGrid);
    /** Returns a vector of D-dimensional arrays of points on the KMesh, which covers the Brillouin zone */
    static std::vector<BZPoint<M>> getAllBZPoints(const KMesh& in);
    /** Returns a vector of a pair of a D-dimensional arrays of points on the KMesh and the amount of points that can be obtained from a symmetry operation in the lattice. */
    static std::map<BZPoint<M>, std::vector<BZPoint<M>>> getUniqueBZPoints(const KMesh& kGrid);
};

template <>
struct CubicTraits<0>{ 
    /** Actual dispersion relation. */
    template <typename ArgType1, typename ...ArgTypes> static RealType ek(RealType t, ArgType1 kpoint1, ArgTypes... kpoints); 
    template <typename ArgType1> static RealType ek(RealType t, ArgType1 kpoint1); 
    /** Returns an analytic std::function of the dispersion. */
    template <typename FunctionType, typename ...ArgTypes> static FunctionType get_dispersion(RealType t) { return [t](ArgTypes ... in){return ek(t,in...);};};
};


struct CubicInfDMFTSC : public DMFT
{
    typedef typename FKImpuritySolver::GFType GFType;
    typedef GridObject<ComplexType, RealGrid> ComplW;
    using DMFT::_S;
public:
    const RealType _t;
    const RealGrid _realgrid;
    ComplW _nominator;

    CubicInfDMFTSC(const FKImpuritySolver &S, RealType t, const RealGrid& realgrid);
    GFType operator()();
    using DMFT::getBubblePI;
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

