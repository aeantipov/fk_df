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
    typedef Eigen::Matrix<RealType,__power<_ksize,D>::value,1,Eigen::ColMajor> EkStorage;
    EkStorage _ek_vals;
private:
    template <size_t M, typename ...ArgTypes> struct fill_ek{ static void fill(EkStorage &in, size_t pos, ArgTypes... other_pos); };
    template <typename ArgType1, typename ...ArgTypes> static RealType ek(RealType t, ArgType1 kpoint1, ArgTypes... kpoints); 
    template <typename ArgType1> static RealType ek(RealType t, ArgType1 kpoint1); 
    typedef typename ArgFunGenerator<D,GFType,RealType>::type ArgFunType;
public:
    const RealType _t;
    const size_t _npoints;
    const KMesh _kgrid;

    CubicDMFTSC(const Solver &S, RealType t, size_t npoints);
    template <typename ...ArgTypes> RealType dispersion(const ArgTypes&... kpoints);
    template <typename ...ArgTypes> GFType glat(const GFType& gw, const GFType& Delta, ArgTypes... kpoints) const;
    template <typename ...ArgTypes> ComplexType glat_val(const ComplexType &gw, const ComplexType &Delta, ArgTypes... kpoints) const;
    GFType operator()() const;
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

