#ifndef ___FK_GRIDOBJECT_HPP___
#define ___FK_GRIDOBJECT_HPP___

#include "GridObject.h"

namespace FK {
//
// GridObject::ContainerExtractor
//

template< typename ValueType, typename ...GridType> 
template <size_t Nc, typename ArgType1, typename ...ArgTypes>
inline ValueType& GridObject<ValueType,GridType...>::ContainerExtractor<Nc,ArgType1,ArgTypes...>::get(
    Container<Nc, ValueType> &data, const std::tuple<GridType...> &grids, 
    const ArgType1& arg1, const ArgTypes&... args) 
{
    const auto & grid=std::get<N-Nc>(grids);
    auto &tmp = grid.getValue(data, arg1);
    return ContainerExtractor<Nc-1,ArgTypes...>::get(tmp,grids,args...);
};

template< typename ValueType, typename ...GridType>
template <typename ArgType1> 
inline ValueType& GridObject<ValueType,GridType...>::ContainerExtractor<1,ArgType1>::get(
    Container<1, ValueType> &data, const std::tuple<GridType...> &grids, const ArgType1& arg1)
{
    const auto & grid=std::get<N-1>(grids);
    auto &tmp = grid.getValue(data, arg1);
    return tmp;
}

//
// GridObject
//

template <typename ValueType, typename ...GridType> 
GridObject<ValueType,GridType...>::GridObject( const std::tuple<GridType...> &in):
    _grids(in)
{
    GetGridSizes<N>::TupleSizeToArray(_grids,_dims);   
    _data.reset(new Container<sizeof...(GridType),ValueType>(_dims));
}

template <typename ValueType, typename ...GridType> 
template <typename ...ArgTypes> 
inline ValueType& GridObject<ValueType,GridType...>::operator()(const ArgTypes&... in)
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridType), "GridObject call number of input parameters mismatch."); 
    return ContainerExtractor<sizeof...(GridType), ArgTypes...>::get(*_data,_grids,in...);
}

template <typename ValueType, typename ...GridType> 
std::ostream& operator<<(std::ostream& lhs, const GridObject<ValueType,GridType...> &in)
{
    lhs << *(in._data);
    return lhs;
}


template <typename ValueType, typename ...GridType> 
inline void GridObject<ValueType,GridType...>::fill(const std::function<ValueType(GridType...)> & in)
{
    _f = in;
}

} // end of namespace FK
#endif
