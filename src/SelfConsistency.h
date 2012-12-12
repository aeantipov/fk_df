#ifndef ___FK_SELFCONSISTENCY_H___
#define ___FK_SELFCONSISTENCY_H___

#include "FKCommon.h"
#include "GridObject.h"
#include "Solver.h"

namespace FK {

struct SelfConsistency 
{
    typedef GridObject<ComplexType,FMatsubaraGrid> GFType;
protected:
    SelfConsistency(){};
};

struct BetheSC : public SelfConsistency
{
    const RealType _t;
    BetheSC(RealType t);
    GFType operator()(const GFType& gw) const; 
    GFType operator()(const GFType& gw, const GFType &Delta) const;
};

template <size_t D> struct CubicDMFTSC : public SelfConsistency
{
private:
    template <typename ArgType1, typename ...ArgTypes> static RealType ek(RealType t, ArgType1 kpoint1, ArgTypes... kpoints); 
    template <typename ArgType1> static RealType ek(RealType t, ArgType1 kpoint1); 
    typedef typename ArgFunGenerator<D,GFType,RealType>::type ArgFunType;
public:
    const RealType _t;
    const size_t _npoints;
    const KMesh _kgrid;

    CubicDMFTSC(RealType t, size_t npoints);
    template <typename ...ArgTypes> RealType dispersion(const ArgTypes&... kpoints);
    template <typename ...ArgTypes> GFType glat(const GFType& gw, const GFType& Delta, ArgTypes... kpoints) const;
    template <typename ...ArgTypes> ComplexType glat_val(const ComplexType &gw, const ComplexType &Delta, ArgTypes... kpoints) const;
    GFType operator()(const GFType& gw, const GFType &Delta) const;
};

struct CubicInfDMFTSC : public SelfConsistency
{
    typedef GridObject<RealType, RealGrid> RealW;
    typedef GridObject<ComplexType, RealGrid> ComplW;
public:
    const RealType _t;
    const RealGrid _realgrid;
    RealW _nominator;
    ComplW _denominator;

    CubicInfDMFTSC(RealType t, const RealGrid realgrid);
    GFType operator()(const GFType& gw, const GFType &Delta) const;
};

} // end of namespace FK
#endif // endif :: ifndef ___FK_SELFCONSISTENCY_H___

#include "SelfConsistency.hpp"

