#ifndef ___FK_GRID_OBJECT_H___
#define ___FK_GRID_OBJECT_H___

#include "FKCommon.h"
#include "Logger.h"
#include "Grid.h"

namespace FK { 

template< typename ValueType, typename GridTypeN, typename ...GridType> 
class GridObject : public GridObject<ValueType,GridType...>
{
static_assert(sizeof...(GridType)>0, "N=1 has a specialized variant");
protected:
    std::shared_ptr<GridTypeN> _grid;
    std::shared_ptr<VectorType<GridObject<ValueType,GridType...> > > _vals;
public:
    /** Default constructor. Required by Eigen::Vector. Does nothing. */
    GridObject(); 
    /** Constructs a grid object out of a tuple containing various grids. */
    template <class ...T> GridObject( const std::tuple<T...> &in);
    /** Returns element number i, which corresponds to (*_grid)[i]. */
    auto operator[](unsigned int i)->decltype((*_vals)[0]);
    const GridTypeN & getGrid() const; 
//    void set(std::function<ValueType (decltype((*_grid)[0]))> f);
//    ValueType operator()(decltype((*_grid)[0]) in ) const;
//    template <unsigned int M, typename ValType, class GridType2> friend std::ostream& operator<<(std::ostream& lhs, const GridObject<M,ValType,GridType2> &in);
};

template <typename ValueType, class GridType>
class GridObject<ValueType, GridType> 
{
protected:
    std::shared_ptr<GridType> _grid;
    std::shared_ptr<VectorType<ValueType> > _vals;
public:
    /** Default constructor. Required by Eigen::Vector. Does nothing. Couldn't find how to make it protected. 
     * A link - http://accu.org/index.php/journals/296. */
    GridObject(); 
    /** Generate an object from a grid. */
    GridObject(const GridType& grid);
    template <class ...T> GridObject(const std::tuple<T...> &in);
    const GridType & getGrid() const; 
}; 

//
// N=1
//

template<typename ValueType, class GridType>
GridObject<ValueType,GridType>::GridObject()
{
    DEBUG("Default, N=1");
}

template<typename ValueType, class GridType>
GridObject<ValueType,GridType>::GridObject(const GridType &grid):
    _grid(std::make_shared<GridType>(grid)),
    _vals(std::make_shared<VectorType<ValueType> > (grid.getSize()))
{
    DEBUG("Initializing from grid, N=1");
}

template<typename ValueType, class GridType>
template <class ...T> 
GridObject<ValueType,GridType>::GridObject(const std::tuple<T...> &in):GridObject(std::get<std::tuple_size<std::tuple<T...> >::value-1>(in))
{
}

template<typename ValueType, class GridType>
const GridType & GridObject<ValueType,GridType>::getGrid() const
{
    return (*_grid);
} 
//
// N!=1
//

template<typename ValueType, typename GridTypeN, typename ...GridType> 
GridObject<ValueType,GridTypeN, GridType...>::GridObject()
{
    DEBUG("Default, N=" << sizeof...(GridType)+1);
}

template<typename ValueType, typename GridTypeN, typename ...GridType> 
template <class ...T>
GridObject<ValueType,GridTypeN, GridType...>::GridObject( const std::tuple<T...> &in):
    _grid(std::make_shared<GridTypeN>(std::get<0>(in))),
    _vals(std::make_shared<VectorType<GridObject<ValueType,GridType...> > >( 
          VectorType<GridObject<ValueType,GridType...> >(
            VectorType< GridObject<ValueType,GridType...> >::Constant(std::get<0>(in).getSize(),GridObject<ValueType,GridType...>(in) ) 
            )
          )
         )
{
    static_assert(std::tuple_size<std::tuple<T...> >::value >= sizeof...(GridType)+1,"Too small elements in the tuple");
    DEBUG("Grid constructor N = " << sizeof...(GridType)+1);
}

template< typename ValueType, typename GridTypeN, typename ...GridType>
auto GridObject<ValueType,GridTypeN, GridType...>::operator[](unsigned int i)->decltype((*_vals)[0])
{
    return (*_vals)[i];
}

template< typename ValueType, typename GridTypeN, typename ...GridType>

const GridTypeN & GridObject<ValueType,GridTypeN, GridType...>::getGrid() const
{
    return (*_grid);
}

/*
template< typename ValueType, class GridType>
GridObject<ValueType,GridType>::GridObject(std::shared_ptr<GridType> grid,  std::function<ValueType (decltype((*grid)[0]))> f):GridObject<ValueType,GridType>(grid)
{
    this->set(f);
}

template< typename ValueType, class GridType>
void GridObject<ValueType,GridType>::set(std::function<ValueType (decltype((*_grid)[0]))> f)
{
    assert(_vals->size() == (*_grid).getSize());
    for (unsigned int i=0; i<_vals->size(); ++i) (*_vals)[i]=f((*_grid)[i]);
}

template <unsigned int M, typename ValType, class GridType2> 
std::ostream& operator<<(std::ostream& lhs, const GridObject<M,ValType,GridType2> &in)
{
    lhs << *(in._vals); 
    return lhs;
}
*/
}
#endif // endif::ifndef ___FK_GRID_OBJECT_H___

