#pragma once
#ifndef __FK_DFDIAGRAMS_H__
#define __FK_DFDIAGRAMS_H__

#include <gftools/kmesh.hpp>
#include <gftools/grid_object.hpp>

namespace FK { 

using namespace gftools;

struct Diagrams
{
	//template <typename V> using MatrixType = typename container<V,2>::MatrixType;
    typedef grid_object<complex_type, fmatsubara_grid> GLocalType;
    //template <size_t D> 
    //using GKType = typename tools::ArgBackGenerator<D,kmesh,grid_object,complex_type,fmatsubara_grid>::type;
    template <size_t D> 
    using wkTupleType = decltype(std::tuple_cat(std::make_tuple(std::declval<typename fmatsubara_grid::point>()),
    							 std::declval<typename tuple_tools::repeater<typename kmesh::point, D>::tuple_type>()));
    template <size_t D> 
    using WQTupleType = decltype(std::tuple_cat(std::make_tuple(std::declval<typename bmatsubara_grid::point>()),
    							 std::declval<typename tuple_tools::repeater<typename kmesh::point, D>::tuple_type>()));
    //const fmatsubara_grid& _fGrid;
    //const kmesh& _kGrid;
    //const size_t _ksize;
    //const size_t _knorm;

    //Diagrams(const fmatsubara_grid& fGrid, const kmesh& kGrid);

    template <typename GKType>
    static GLocalType getBubble(const GKType &GF, const GKType &GF_shift);
    template <typename GKType, typename ... ArgTypes>
    static GLocalType getBubble(const GKType &GF, ArgTypes... args);
    template <typename GKType, typename ... ArgTypes>
    static GLocalType getBubble(const GKType &GF, const std::tuple<ArgTypes...>& args);
    template <typename ArgType>
    static GLocalType getBubble(const GLocalType& GF, ArgType arg);

    template <typename G>
    using IsGridObject1 = typename std::enable_if<(G::N==1) && sizeof(typename G::value_type), G>::type;

    template <typename G>
    static IsGridObject1<G> BS(const G &Chi0, const G &IrrVertex4,
            bool forward, bool eval_SC = false, size_t n_iter = 100, real_type mix = 1.0);

    template <typename G>
    using IsGridObject = typename std::enable_if<(G::N>1) && sizeof(typename G::value_type), G>::type;

    template <typename G>
    static IsGridObject<G> BS(
            const G &Chi0,
            const G &IrrVertex4,
            bool forward, bool eval_SC = false, size_t n_iter = 100, real_type mix = 1.0);

    template <typename MatrixType>
    using IsMatrix = typename std::enable_if<std::is_convertible<decltype(std::declval<MatrixType>().rows()), int>::value, MatrixType>::type;

    template <typename MatrixType>
    static IsMatrix<MatrixType> BS(
        const MatrixType &Chi0,
        const MatrixType &IrrVertex4,
        bool forward, bool eval_SC = false, size_t n_iter = 100, real_type mix= 1.0);
//    template <typename VertexType> 
//    static VertexType BS(const VertexType &Chi0, const VertexType &IrrVertex4, bool eval_SC = false, size_t n_iter = 100, real_type mix = 1.0);

    template <typename VertexType> 
    static VertexType getSusc(const VertexType &Chi0, const VertexType &FullVertex4);

};

} // end of namespace FK 

#include "Diagrams.hpp"

#endif // endif :: #ifndef __FK_DFDIAGRAMS_H__
