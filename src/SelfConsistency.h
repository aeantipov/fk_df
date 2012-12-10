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
public:
   // GFType operator()(const GFType& gw)=0; 
};

struct BetheSC : public SelfConsistency
{
    const RealType _t;
    BetheSC(RealType t);
    GFType operator()(const GFType& gw); 
};

template <size_t D> struct CubicDMFTSC : public SelfConsistency
{
private:
    template <typename ArgType1, typename ...ArgTypes> RealType ek(const ArgType1& kpoint1, const ArgTypes&... kpoints); 
    template <typename ArgType1> RealType ek(const ArgType1& kpoint1); 
    //template <> auto getSC<1>(const GFType& gw, const GFType& Delta) const -> GFType;
    public:
    const RealType _t;
    const size_t _npoints;
    const KMesh _kgrid;
    //const FMatsubaraGrid _wgrid;
    //GFType &_Delta_old;
    //template <typename ArgType1, typename ...ArgTypes> auto getSC(const GFType& gw, const GFType& Delta) -> decltype(getSC<ArgTypes...>(gw,Delta));

    CubicDMFTSC(RealType t, size_t npoints);
    template <typename ...ArgTypes> RealType dispersion(const ArgTypes&... kpoints);
    template <typename ...ArgTypes> GFType glat(const GFType& gw, const GFType& Delta, const ArgTypes&... kpoints);
    GFType operator()(const GFType& gw, const GFType &Delta);
};

template <size_t N, class Obj,typename ...OtherArgs> struct RecursiveGridIntegratorType;

template <class Obj,typename ...OtherArgs> struct RecursiveGridIntegratorType<1,Obj,OtherArgs...>
{
    typedef decltype(std::declval<Obj>()(1.0,std::declval<OtherArgs>()...)) type;
};

template<size_t ...> struct seq { };
template<size_t N, size_t...S> struct gens : gens<N-1, N-1, S...> {};
template<size_t...S> struct gens<0, S...>{ typedef seq<S...> type; };

template <size_t N, class Obj,typename ...OtherArgs> struct RecursiveGridIntegratorType
{
    template <typename KArg1=RealType, typename ...KArgs> struct mytype { 
        typedef std::result_of<Obj(KArg1, KArgs..., OtherArgs...)> type;
    };
    static auto t1 = [&](RealType k, OtherArgs... args){return std::declval<Obj>()(k,args...);};
    typedef decltype(t1) mytype1;
};

template <size_t N, class Obj,typename ...OtherArgs> struct RecursiveGridIntegrator {
    inline auto integrate(const KMesh& kgrid, Obj& in, const OtherArgs&... Args)
        -> RecursiveGridIntegratorType<N,Obj,OtherArgs...> 
           
    { 
        //auto f1 = [&in,&Args...](RealType k){return in(k,Args...);};
        //return kgrid.integrate(f1);
    } 
};

template <class Obj, typename ...OtherArgs> struct RecursiveGridIntegrator<1,Obj,OtherArgs...> {
    typedef decltype(std::declval<Obj>()(0.0,std::declval<OtherArgs>()...)) type;
    inline type integrate(const KMesh& kgrid, const Obj& in, const OtherArgs&... Args){ 
        DEBUG("Using this method");
        auto f1 = [&in,&Args...](RealType k){return in(k,Args...);};
        return kgrid.integrate(f1);
    }
};
} // end of namespace FK
#endif // endif :: ifndef ___FK_SELFCONSISTENCY_H___

#include "SelfConsistency.hpp"

