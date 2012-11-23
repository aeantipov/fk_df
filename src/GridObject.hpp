#ifndef ___FK_GRIDOBJECT_HPP___
#define ___FK_GRIDOBJECT_HPP___

#include "GridObject.h"

namespace FK {
//
// GridObject
//


template <typename ValueType, class ...GridType> 
GridObject<ValueType,GridType...>::GridObject( const std::tuple<GridType...> &in):
    _grids(in)
{
    GetGridSizes<N>::TupleSizeToArray(_grids,_dims);   
    _data.reset(new Container<sizeof...(GridType),ValueType>(_dims));
}

template <typename ValueType, class ...GridType> 
std::ostream& operator<<(std::ostream& lhs, const GridObject<ValueType,GridType...> &in)
{
}

} // end of namespace FK
#endif
