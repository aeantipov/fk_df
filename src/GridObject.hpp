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
//    _data(Container<sizeof...(GridType),ValueType>(in)) 
{
    //_dims = std::array<unsigned int, 3>({1,2,3});
}

template <typename ValueType, class ...GridType> 
std::ostream& operator<<(std::ostream& lhs, const GridObject<ValueType,GridType...> &in)
{
}

} // end of namespace FK
#endif
