#ifndef ___FK_SELFCONSISTENCY_H___
#define ___FK_SELFCONSISTENCY_H___

#include "FKCommon.h"
#include "GridObject.h"
#include "Solver.h"

namespace FK {

struct SelfConsistency 
{
    typedef typename FKImpuritySolver::GFType GFType;
    const FKImpuritySolver &_S;
    SelfConsistency(const FKImpuritySolver &S):_S(S){};
    virtual GFType operator()() = 0;
    template <typename MPoint> GFType getBubblePI(MPoint in) const;
    template <typename MPoint> GFType getBubble0(MPoint in) const;
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> getBubblePI() const;
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> getBubble0() const;
};

struct BetheSC : public SelfConsistency
{
    typedef typename FKImpuritySolver::GFType GFType;
    const RealType _t;
    BetheSC(const FKImpuritySolver &S, RealType t);
    GFType operator()();
};

template <size_t D> struct CubicDMFTSC : public SelfConsistency
{
    using SelfConsistency::_S;
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
};
    
template <size_t M> 
struct CubicTraits{ 
    //const KMesh grid = KMesh(ksize);
    //typedef Eigen::Matrix<RealType,__power<ksize,M>::value,1,Eigen::ColMajor> EkStorage;
    //typedef typename ArgBackGenerator<M,KMesh,std::tuple>::type KMeshTupleType;
    typedef typename std::array<KMesh,M> KMeshTupleType;
    //template <class ... GridTypes> static void gen_tuples(std::tuple<GridTypes...> grids){gen_tuples<GridTypes...,KMesh>(std::tuple_cat(grids,KMesh(ksize)));}
    static KMeshTupleType getTuples(const KMesh &grid){std::array<KMesh,M> out; out.fill(grid); return out;};
    template <class IteratorType, typename ...ArgTypes> static void fill(IteratorType in, RealType t, const KMesh& grid, ArgTypes... other_pos); 
    template <class ContainerType, typename ...ArgTypes> static void fillContainer(ContainerType &in, RealType t, const KMesh& grid, ArgTypes... other_pos); 
    template <class ObjType, typename ...ArgTypes> static void setAnalytics(ObjType &in, RealType t);
    template <typename FunctionType, typename ...ArgTypes> static FunctionType get_dispersion(RealType t)
        { return CubicTraits<M-1>::template get_dispersion<FunctionType,ArgTypes...,RealType>(t); }
    //static RealType ek(RealType t, ArgTypes... kpoints); 
};

template <>
struct CubicTraits<0>{ 
    template <class IteratorType, typename ...ArgTypes> static void fill(IteratorType in, RealType t, const KMesh& grid, ArgTypes... other_pos); 
    template <class DataType, typename ...ArgTypes> static void fillContainer(DataType &in, RealType t, const KMesh& grid, ArgTypes... other_pos); 
    template <typename ...GridTypes> static std::tuple<GridTypes...> getKMeshTuple(GridTypes ... extra_grids)
        {return std::make_tuple(extra_grids...);};
    //template <typename ...ArgTypes> static RealType ek(RealType t, ArgTypes... kpoints) { return 1.0; } ; 
    template <typename ArgType1, typename ...ArgTypes> static RealType ek(RealType t, ArgType1 kpoint1, ArgTypes... kpoints); 
    template <typename ArgType1> static RealType ek(RealType t, ArgType1 kpoint1); 
    template <typename FunctionType, typename ...ArgTypes> static FunctionType get_dispersion(RealType t) { return [t](ArgTypes ... in){return ek(t,in...);};};
};


struct CubicInfDMFTSC : public SelfConsistency
{
    typedef typename FKImpuritySolver::GFType GFType;
    typedef GridObject<ComplexType, RealGrid> ComplW;
    using SelfConsistency::_S;
public:
    const RealType _t;
    const RealGrid _realgrid;
    ComplW _nominator;

    CubicInfDMFTSC(const FKImpuritySolver &S, RealType t, const RealGrid& realgrid);
    GFType operator()();
    using SelfConsistency::getBubblePI;
    template <typename MPoint> GFType getBubble0(MPoint in) const;
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> getBubble0() const;
};

} // end of namespace FK
#endif // endif :: ifndef ___FK_SELFCONSISTENCY_H___

#include "SelfConsistency.hpp"

