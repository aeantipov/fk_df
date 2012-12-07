#ifndef ___FK_SELFCONSISTENCY_H___
#define ___FK_SELFCONSISTENCY_H___

#include "FKCommon.h"
#include "GridObject.h"
#include "Solver.h"

namespace FK {

class SelfConsistency 
{
protected:
    SelfConsistency(){};
public:
    typedef GridObject<ComplexType,FMatsubaraGrid> GFType;
    virtual GFType operator()(const GFType &gw)=0; 
};

struct BetheSC : public SelfConsistency
{
    const RealType _t;
    BetheSC(RealType t);
    GFType operator()(const GFType &gw);
};

template <size_t D> struct CubicDMFTSC : public SelfConsistency
{
private:
    template <typename ArgType1, typename ...ArgTypes> RealType ek(const ArgType1& kpoint1, const ArgTypes&... kpoints); 
    template <typename ArgType1> RealType ek(const ArgType1& kpoint1); 
public:
    const RealType _t;
    const size_t _npoints;
    const KMesh _kgrid;
    //const FMatsubaraGrid _wgrid;

    CubicDMFTSC(RealType t, size_t npoints);
    template <typename ...ArgTypes> RealType dispersion(const ArgTypes&... kpoints);
    template <typename ...ArgTypes> ComplexType glat(const GFType&gw, const GFType&Delta, const ArgTypes&... kpoints);
    GFType operator()(const GFType &gw);
};

} // end of namespace FK
#endif // endif :: ifndef ___FK_SELFCONSISTENCY_H___

#include "SelfConsistency.hpp"

