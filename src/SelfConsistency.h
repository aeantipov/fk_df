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
    GFType operator()() const;
};

template <class Solver>
struct BetheSC : public SelfConsistency<Solver>
{
    typedef typename Solver::GFType GFType;
    const RealType _t;
    BetheSC(const Solver &S, RealType t);
    GFType operator()() const;
};

template <class Solver, size_t D> struct CubicDMFTSC : public SelfConsistency<Solver>
{
    typedef typename Solver::GFType GFType;
    static const size_t _ksize = 32;
    typedef Eigen::Array<ComplexType,__power<_ksize,D>::value,1,Eigen::ColMajor> EkStorage;
    EkStorage _ek_vals;
private:
    typedef typename ArgFunGenerator<D,GFType,RealType>::type ArgFunType;
public:
    const RealType _t;
    const size_t _npoints;
    const KMesh _kgrid;
    GFType _gloc;

    CubicDMFTSC(const Solver &S, RealType t, size_t npoints);
    template <typename ...ArgTypes> RealType dispersion(const ArgTypes&... kpoints);
    template <typename ...ArgTypes> GFType glat(ArgTypes... kpoints) const;
    GFType operator()();
};
    
template <size_t M, size_t ksize> 
struct __fill_ek{ 
    //typedef Eigen::Matrix<RealType,__power<ksize,M>::value,1,Eigen::ColMajor> EkStorage;
    template <class IteratorType, typename ...ArgTypes> static void fill(IteratorType in, RealType t, const KMesh& grid, ArgTypes... other_pos); 
    //static RealType ek(RealType t, ArgTypes... kpoints); 
};

template <size_t ksize>
struct __fill_ek<0, ksize>{ 
    template <class IteratorType, typename ...ArgTypes> static void fill(IteratorType in, RealType t, const KMesh& grid, ArgTypes... other_pos); 
    //template <typename ...ArgTypes> static RealType ek(RealType t, ArgTypes... kpoints) { return 1.0; } ; 
    template <typename ArgType1, typename ...ArgTypes> static RealType ek(RealType t, ArgType1 kpoint1, ArgTypes... kpoints); 
    template <typename ArgType1> static RealType ek(RealType t, ArgType1 kpoint1); 
};


template <class Solver>
struct CubicInfDMFTSC : public SelfConsistency<Solver>
{
    typedef typename Solver::GFType GFType;
    typedef GridObject<ComplexType, RealGrid> ComplW;
public:
    const RealType _t;
    const RealGrid _realgrid;
    ComplW _nominator;

    CubicInfDMFTSC(const Solver &S, RealType t, const RealGrid& realgrid);
    GFType operator()() const;
};

} // end of namespace FK
#endif // endif :: ifndef ___FK_SELFCONSISTENCY_H___

#include "SelfConsistency.hpp"

