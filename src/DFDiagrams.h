#ifndef __FK_DFDIAGRAMS_H__
#define __FK_DFDIAGRAMS_H__

#include "FKCommon.h"
#include "GridObject.h"
#include "GFWrap.h"

namespace FK { 

template <size_t D>
struct DFDiagrams
{
    typedef GFWrap GLocalType;
    typedef typename ArgBackGenerator<D,KMesh,GridObject,ComplexType,FMatsubaraGrid>::type GKType;
    typedef decltype(std::tuple_cat(std::make_tuple(FMatsubaraGrid::point()), std::array<KMesh::point, D>())) wkTupleType; 
    typedef decltype(std::tuple_cat(std::make_tuple(BMatsubaraGrid::point()), std::array<KMesh::point, D>())) WQTupleType; 
    const FMatsubaraGrid& _fGrid;
    const KMesh& _kGrid;
    const size_t _ksize;
    const size_t _knorm;

    DFDiagrams(const FMatsubaraGrid& fGrid, const KMesh& kGrid);

    template <typename ...KP> GLocalType getBubble(const typename DFDiagrams<D>::GKType& GF, BMatsubaraGrid::point W, KP...kpoints) const;
    GLocalType getBubble(const GKType& GF, const WQTupleType& args) const;

    template <typename VertexType> 
    VertexType BS(const VertexType &Chi0, const VertexType &Vertex4, bool eval_SC = false, size_t n_iter = 100, RealType mix = 1.0) const;

    template <typename VertexType> 
    static VertexType getSusc(const VertexType &Chi0, const VertexType &FullVertex4);

};

} // end of namespace FK 

#include "DFDiagrams.hpp"

#endif // endif :: #ifndef __FK_DFDIAGRAMS_H__
