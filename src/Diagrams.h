#ifndef __FK_DFDIAGRAMS_H__
#define __FK_DFDIAGRAMS_H__

#include "FKCommon.h"
#include "GridObject.h"
#include "GFWrap.h"

namespace FK { 

struct Diagrams
{
    typedef GFWrap GLocalType;
    //template <size_t D> 
    //using GKType = typename ArgBackGenerator<D,KMesh,GridObject,ComplexType,FMatsubaraGrid>::type;
    template <size_t D> 
    using wkTupleType = decltype(std::tuple_cat(std::make_tuple(FMatsubaraGrid::point()), std::array<KMesh::point, D>()));
    template <size_t D> 
    using WQTupleType = decltype(std::tuple_cat(std::make_tuple(BMatsubaraGrid::point()), std::array<KMesh::point, D>())); 
    //const FMatsubaraGrid& _fGrid;
    //const KMesh& _kGrid;
    //const size_t _ksize;
    //const size_t _knorm;

    //Diagrams(const FMatsubaraGrid& fGrid, const KMesh& kGrid);

    template <typename GKType, typename ... ArgTypes>
    static GLocalType getBubble(const GKType &GF, ArgTypes... args);
    template <typename GKType, typename ... ArgTypes>
    static GLocalType getBubble(const GKType &GF, const std::tuple<ArgTypes...>& args);
    template <typename ArgType>
    static GLocalType getBubble(const GLocalType& GF, ArgType arg);

    template <typename VertexType> 
    static VertexType BS(const VertexType &Chi0, const VertexType &IrrVertex4, bool eval_SC = false, size_t n_iter = 100, RealType mix = 1.0);

    template <typename VertexType> 
    static VertexType getSusc(const VertexType &Chi0, const VertexType &FullVertex4);

};

} // end of namespace FK 

#include "Diagrams.hpp"

#endif // endif :: #ifndef __FK_DFDIAGRAMS_H__