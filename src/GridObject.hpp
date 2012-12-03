#ifndef ___FK_GRIDOBJECT_HPP___
#define ___FK_GRIDOBJECT_HPP___

#include "GridObject.h"

namespace FK {
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
    //return _get_value_from_container(*_data, in);
}


/*
template <typename ValueType, typename ...GridType> 
template <typename ...ArgTypes> 
inline ValueType& GridObject<ValueType,GridType...>::operator()(const std::tuple<ArgTypes...>& in)
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridType), "GridObject call number of input parameters mismatch."); 
    return _get_value_from_container(*_data, in);
}

template <typename ValueType, typename ...GridType> 
template <size_t Nc, typename ...ArgTypes> 
inline ValueType& GridObject<ValueType, GridType...>::_get_value_from_container(
    Container<Nc, ValueType> &data, 
    const std::tuple<ArgTypes...>& args)
{
    static_assert(Nc>1,"");
    auto & grid=std::get<N-Nc>(_grids);
    auto & x = std::get<N-Nc>(args);
    //return _get_value_from_container(grid.getValue(data, x), args); 
    auto &tmp = grid.getValue(data, x);
    DEBUG("!" << Nc);
    _get_value_from_container(tmp, args);
    //if (Nc>1) return _get_value_from_container(tmp, args); 
    //else return tmp;
}

*/
template <typename ValueType, typename ...GridType> 
std::ostream& operator<<(std::ostream& lhs, const GridObject<ValueType,GridType...> &in)
{
    lhs << *(in._data);
    return lhs;
}

} // end of namespace FK
#endif
