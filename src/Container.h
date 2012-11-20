#ifndef ___FK_CONTAINER_H___
#define ___FK_CONTAINER_H___

#include "FKCommon.h"
#include "Logger.h"
#include "Grid.h"

namespace FK { 

/** This class is a N-dimensional container to store the data of
 * the ValueType type. */
template <size_t N, typename ValueType>
class Container {
static_assert(N>1, "N=1 has a specialized variant");
protected:
    /** The data is stored as a recursive vector<vector ... > > structure. */
    std::vector<Container<N-1, ValueType> >  _vals;
public:
    typedef typename std::vector<Container<N-1, ValueType>>::iterator iterator;
    typedef typename std::vector<Container<N-1, ValueType>>::const_iterator const_iterator;
    /* Constructor from the tuple of Grids. */
    template <class ...T> Container ( const std::tuple<T...> &in);
    /* Constructor from the std::array of size_t. */
    template <size_t M> Container ( const std::array<size_t, M> &in);
    /* Returns the value at point i. */
    auto operator[](size_t i)->decltype(_vals[0]);
    /* Begin iterator. */
    typename Container<N,ValueType>::iterator begin();
    /* End iterator. */
    typename Container<N,ValueType>::iterator end();
    /** Make the object streamable. */
    template <size_t M, typename ValType> friend std::ostream& operator<<(std::ostream& lhs, const Container<M,ValType> &in);
};

/* A specialization for the 1d Container. */
template <typename ValueType>
class Container<1,ValueType> {
protected:
    std::vector<ValueType>  _vals;
public:
    typedef typename std::vector<ValueType>::iterator iterator;
    typedef typename std::vector<ValueType>::const_iterator const_iterator;
    /* Constructor from the tuple of Grids. The last grid in tuple is used. */
    template <class ...T> Container ( const std::tuple<T...> &in);
    /* Constructor from the std::array of size_t. The last value is used. */
    template <size_t M> Container ( const std::array<size_t, M> &in);
    /* Returns the value at index i. */
    ValueType& operator[](size_t i);
    /* Begin iterator. */
    typename Container<1,ValueType>::iterator begin();
    /* End iterator. */
    typename Container<1,ValueType>::iterator end();
    /** Make the object streamable. */
    template <typename ValType> friend std::ostream& operator<<(std::ostream& lhs, const Container<1,ValType> &in);
};

}; // end of namespace FK
#endif

#include "Container.hpp"
