#ifndef ___FK_SELFCONSISTENCY_H___
#define ___FK_SELFCONSISTENCY_H___

#include "FKCommon.h"
#include "GridObject.h"
#include "Solver.h"

namespace FK {

template <class Solver>
struct SelfConsistency 
{
    typedef typename Solver::GFType GFType;
    const Solver &_S;
    SelfConsistency(const Solver &S):_S(S){};
    virtual GFType operator()() = 0;
    template <typename MPoint> GFType getLatticeDMFTVertex4(MPoint in) const;
    template <typename MPoint> GFType getBubblePI(MPoint in) const;
    template <typename MPoint> GFType getBubble0(MPoint in) const;
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> getBubblePI() const;
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> getBubble0() const;

    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> getStaticLatticeDMFTVertex4() const;
};

template <class Solver>
struct BetheSC : public SelfConsistency<Solver>
{
    typedef typename Solver::GFType GFType;
    const RealType _t;
    BetheSC(const Solver &S, RealType t);
    GFType operator()();
};

template <class Solver, size_t D, size_t ksize=32> struct CubicDMFTSC : public SelfConsistency<Solver>
{
    using SelfConsistency<Solver>::_S;
    typedef typename Solver::GFType GFType;
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

    CubicDMFTSC(const Solver &S, RealType t);
    template <typename ...ArgTypes> RealType dispersion(ArgTypes... kpoints) const;
    template <typename ...ArgTypes> RealType dispersion(const std::tuple<ArgTypes...>& kpoints) const;
    template <typename ...ArgTypes> GFType glat(ArgTypes... kpoints) const;
    template <typename MPoint, typename ...ArgTypes> ComplexType glat_val(MPoint w, ArgTypes... kpoints) const;
    GKType getGLat(const FMatsubaraGrid& in) const;
    GFType operator()();
    
    /*
    template <typename MPoint, typename ...QPoints>
    RealType getLatticeSusceptibility(MPoint w, QPoints ...qpts) const;
    template <typename MPoint, typename ...QPoints>
    RealType getLatticeSusceptibility(const std::tuple<MPoint,QPoints...> & qpts) const;*/
    
    template <typename MPoint> GFType getBubble0(MPoint in) const;
    template <typename MPoint> GFType getBubblePI(MPoint in) const;
};
    
template <size_t M, size_t ksize> 
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
        { return CubicTraits<M-1, ksize>::template get_dispersion<FunctionType,ArgTypes...,RealType>(t); }
    //static RealType ek(RealType t, ArgTypes... kpoints); 
};

template <size_t ksize>
struct CubicTraits<0, ksize>{ 
    template <class IteratorType, typename ...ArgTypes> static void fill(IteratorType in, RealType t, const KMesh& grid, ArgTypes... other_pos); 
    template <class DataType, typename ...ArgTypes> static void fillContainer(DataType &in, RealType t, const KMesh& grid, ArgTypes... other_pos); 
    template <typename ...GridTypes> static std::tuple<GridTypes...> getKMeshTuple(GridTypes ... extra_grids)
        {return std::make_tuple(extra_grids...);};
    //template <typename ...ArgTypes> static RealType ek(RealType t, ArgTypes... kpoints) { return 1.0; } ; 
    template <typename ArgType1, typename ...ArgTypes> static RealType ek(RealType t, ArgType1 kpoint1, ArgTypes... kpoints); 
    template <typename ArgType1> static RealType ek(RealType t, ArgType1 kpoint1); 
    template <typename FunctionType, typename ...ArgTypes> static FunctionType get_dispersion(RealType t) { return [t](ArgTypes ... in){return ek(t,in...);};};
};


template <class Solver>
struct CubicInfDMFTSC : public SelfConsistency<Solver>
{
    typedef typename Solver::GFType GFType;
    typedef GridObject<ComplexType, RealGrid> ComplW;
    using SelfConsistency<Solver>::_S;
public:
    const RealType _t;
    const RealGrid _realgrid;
    ComplW _nominator;

    CubicInfDMFTSC(const Solver &S, RealType t, const RealGrid& realgrid);
    GFType operator()();
        //template <typename MPoint, typename ...QPoints>
    //RealType getLatticeSusceptibility(MPoint w, bool zeroOrPi) const;
    using SelfConsistency<Solver>::getBubblePI;
    template <typename MPoint> GFType getBubble0(MPoint in) const;
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> getBubble0() const;
};

} // end of namespace FK
#endif // endif :: ifndef ___FK_SELFCONSISTENCY_H___

#include "SelfConsistency.hpp"

