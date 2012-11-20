#ifndef ___FK_GRID_OBJECT_H___
#define ___FK_GRID_OBJECT_H___

#include "FKCommon.h"
#include "Logger.h"
#include "Grid.h"
#include "Container.h"

namespace FK {

/** A GridObject is a heterogeneous container class, that stores data, 
 * defined on multiple different grids, which can be of different types.
 */
template< typename ValueType, typename ...GridType> 
class GridObject 
{
protected:
    static const int _N = sizeof...(GridType);
    static void _fill_dims(const std::tuple<GridType...>& in, std::array<size_t,_N> &out); 
    std::tuple<GridType...> _grids;
    std::unique_ptr<Container<_N, ValueType> > _data;
public:
    /** Constructs a grid object out of a tuple containing various grids. */
    GridObject( const std::tuple<GridType...> &in);

    /** Returns element number i, which corresponds to (*_grid)[i]. */
    auto operator[](unsigned int i)->decltype((*_data)[0]);
    /** Returns the associated grid. */
    //auto getGrid()->decltype(std::get<0>(_grids)) const; 
//    ValueType operator()(decltype((*_grid)[0]) in ) const;
    template <typename ValType, class ...GridType2> friend std::ostream& operator<<(std::ostream& lhs, const GridObject<ValType,GridType2...> &in);
};

} // end of namespace FK
#endif // endif::ifndef ___FK_GRID_OBJECT_H___

#include "GridObject.hpp"

//
// N=1
//

/*
template<typename ValueType, class GridType>
inline
GridObject<ValueType,GridType>::GridObject()
{
//    DEBUG("Default, N=1");
}

template<typename ValueType, class GridType>
inline
GridObject<ValueType,GridType>::GridObject(const GridType &grid):
    _grid(std::make_shared<GridType>(grid)),
    _vals(std::vector<ValueType>  (grid.getSize()))
{
//    DEBUG("Initializing from grid, N=1");
}

template<typename ValueType, class GridType>
template <class ...T> 
inline
GridObject<ValueType,GridType>::GridObject(const std::tuple<T...> &in):GridObject(std::get<std::tuple_size<std::tuple<T...> >::value-1>(in))
{
}

template<typename ValueType, class GridType>
inline
const GridType & GridObject<ValueType,GridType>::getGrid() const
{
    return (*_grid);
} 

template< typename ValueType, class GridType>
inline
void GridObject<ValueType,GridType>::set(std::function<ValueType (decltype((*_grid)[0]))> f)
{
    assert(_vals.size() == (*_grid).getSize());
    for (unsigned int i=0; i<_vals.size(); ++i) (_vals)[i]=f((*_grid)[i]);
}

//
// N!=1
//

template<typename ValueType, typename GridTypeN, typename ...GridType> 
inline
GridObject<ValueType,GridTypeN, GridType...>::GridObject()
{
//    DEBUG("Default, N=" << sizeof...(GridType)+1);
}

template<typename ValueType, typename GridTypeN, typename ...GridType> 
template <class ...T>
inline
GridObject<ValueType,GridTypeN, GridType...>::GridObject( const std::tuple<T...> &in):
    _grid(std::make_shared<GridTypeN>(std::get<0>(in))),
    //_vals(std::vector<GridObject<ValueType,GridType...> > (std::get<0>(in).getSize())) 
    _vals(std::vector<GridObject<ValueType,GridType...> >(std::get<0>(in).getSize(), std::move(GridObject<ValueType,GridType...>(in)))) 
{

    static_assert(std::tuple_size<std::tuple<T...> >::value >= sizeof...(GridType)+1,"Too small elements in the tuple");
//    for (int i=0; i<_vals.size(); ++i) {
//            (_vals)[i]=std::move(GridObject<ValueType,GridType...>(in));
//        }
//    std::fill(index_begin<GridObject<ValueType,GridType...>>(_vals),index_end<GridObject<ValueType,GridType...>>(_vals),std::move(GridObject<ValueType,GridType...>(in)));
//    DEBUG("Grid constructor N = " << sizeof...(GridType)+1);
//
          //std::vector<GridObject<ValueType,GridType...> >(
           // std::vector< GridObject<ValueType,GridType...> >::Constant(std::get<0>(in).getSize(),std::move(GridObject<ValueType,GridType...>(in)) ) 
}

template< typename ValueType, typename GridTypeN, typename ...GridType>
auto GridObject<ValueType,GridTypeN, GridType...>::operator[](unsigned int i)->decltype((_vals)[0])
{
    return (_vals)[i];
}

template< typename ValueType, typename GridTypeN, typename ...GridType>

const GridTypeN & GridObject<ValueType,GridTypeN, GridType...>::getGrid() const
{
    return (*_grid);
}


template <typename ValType, class GridType1, class ...GridType2> 
std::ostream& operator<<(std::ostream& lhs, const GridObject<ValType,GridType1, GridType2...> &in)
{
    if (sizeof...(GridType2)>0) {
        std::ostream_iterator<GridObject<ValType,GridType2...> > out_it (lhs,", ");
        std::copy(in._vals.begin(), in._vals.end(), out_it);
        }
    else {
        std::ostream_iterator<ValType> out_it (lhs,", ");
        std::copy(in._vals.begin(), in._vals.end(), out_it);
         }
    return lhs;
}

*/

