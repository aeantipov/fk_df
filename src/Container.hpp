#ifndef ___FK_CONTAINER_HPP___
#define ___FK_CONTAINER_HPP___

#include "Container.h"

namespace FK {

//
// Container, N=1
//

template <typename ValueType>
template <class ...T> 
inline Container<1,ValueType>::Container ( const std::tuple<T...> &in) : 
    _vals(std::vector<ValueType>( std::get<std::tuple_size<std::tuple<T...> >::value-1>(in).size()))
{
    DEBUG("Constructing from grid size.");
};

template <typename ValueType>
template <size_t M>
inline Container<1,ValueType>::Container ( const std::array<size_t, M> &in):_vals(std::vector<ValueType>(std::get<M-1>(in)))
{
    DEBUG("Constructing from size_t.");
}

template <typename ValueType>
inline ValueType& Container<1,ValueType>::operator[] (size_t i)
{
    return _vals[i];
}

template <typename ValueType> 
std::ostream& operator<<(std::ostream& lhs, const Container<1,ValueType> &in)
{
    lhs << "[";
    std::ostream_iterator<ValueType> out_it (lhs,", ");
    std::copy(index_begin<ValueType>(in._vals),index_end<ValueType>(in._vals) , out_it);
    lhs << "]";
    return lhs;
}

template <typename ValueType> 
typename Container<1,ValueType>::iterator Container<1,ValueType>::begin()
{
    return _vals.begin();
}

template <typename ValueType> 
typename Container<1,ValueType>::iterator Container<1,ValueType>::end()
{
    return _vals.end();
}

//
// Container, N!=1
//

template <size_t N, typename ValueType>
template <class ...T> //, typename = typename std::enable_if<std::is_function<std::tuple_element<0, std::tuple<T...> >.getGrid()>::value> 
inline Container<N,ValueType>::Container ( const std::tuple<T...> &in) : 
    _vals(std::vector<Container<N-1, ValueType>>( 
        (std::get<std::tuple_size<std::tuple<T...> >::value-N>(in), Container<N-1, ValueType>(in)).size()
        )
    )
{
};

template <size_t N, typename ValueType>
template <size_t M>
inline Container<N,ValueType>::Container ( const std::array<size_t, M> &in):
    _vals(std::vector<Container<N-1,ValueType> >(std::get<M-1>(in), Container<N-1,ValueType>(in)))
{
}



template <size_t N, typename ValueType>
inline auto Container<N,ValueType>::operator[](size_t i)->decltype(_vals[0])
{
    return _vals[i];
}

template <size_t N, typename ValueType> 
std::ostream& operator<<(std::ostream& lhs, const Container<N,ValueType> &in)
{
    lhs << "[";
    std::ostream_iterator<Container<N-1,ValueType> > out_it (lhs,", ");
    std::copy(in._vals.begin(),in._vals.end(), out_it);
    lhs << "]";
    return lhs;
}

template <size_t N, typename ValueType> 
typename Container<N,ValueType>::iterator Container<N,ValueType>::begin()
{
    return _vals.begin();
}

template <size_t N, typename ValueType> 
typename Container<N,ValueType>::iterator Container<N,ValueType>::end()
{
    return _vals.end();
}

} // end of namespace FK
#endif
