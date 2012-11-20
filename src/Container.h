#ifndef ___FK_CONTAINER_H___
#define ___FK_CONTAINER_H___

#include "FKCommon.h"
#include "Logger.h"
#include "Grid.h"

namespace FK { 

template <size_t N, typename ValueType>
class Container {
static_assert(N>1, "N=1 has a specialized variant");
protected:
    std::vector<Container<N-1, ValueType> >  _vals;
public:
    typedef typename std::vector<Container<N-1, ValueType>>::iterator iterator;
    typedef typename std::vector<Container<N-1, ValueType>>::const_iterator const_iterator;
    template <class ...T> Container ( const std::tuple<T...> &in); 
    template <size_t M> Container ( const std::array<size_t, M> &in); 
    auto operator[](size_t i)->decltype(_vals[0]);
    typename Container<N,ValueType>::iterator begin();
    typename Container<N,ValueType>::iterator end();
    template <size_t M, typename ValType> friend std::ostream& operator<<(std::ostream& lhs, const Container<M,ValType> &in);
};

template <typename ValueType>
class Container<1,ValueType> {
protected:
    std::vector<ValueType>  _vals;
public:
    typedef typename std::vector<ValueType>::iterator iterator;
    typedef typename std::vector<ValueType>::const_iterator const_iterator;
    template <class ...T> Container ( const std::tuple<T...> &in); 
    template <size_t M> Container ( const std::array<size_t, M> &in); 
    ValueType& operator[](size_t i);
    typename Container<1,ValueType>::iterator begin();
    typename Container<1,ValueType>::iterator end();
    template <typename ValType> friend std::ostream& operator<<(std::ostream& lhs, const Container<1,ValType> &in);
};

}; // end of namespace FK
#endif

#include "Container.hpp"
