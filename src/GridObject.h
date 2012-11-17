#ifndef ___FK_GRID_OBJECT_H___
#define ___FK_GRID_OBJECT_H___

#include "FKCommon.h"
#include "Logger.h"
#include "Grid.h"

namespace FK { 

template<unsigned int N, typename ValueType, typename GridTypeN, typename ...GridType> 
class Grid1dObject : public Grid1dObject<N-1,ValueType,GridType...>
{
static_assert(N>1, "N=1 has a specialized variant");
protected:
    std::shared_ptr<GridTypeN> _grid;
    std::shared_ptr<VectorType<Grid1dObject<N-1,ValueType,GridType...> > > _vals;
    Grid1dObject();
public:
//    Grid1dObject(std::shared_ptr<GridType> grid,  std::function<ValueType (decltype((*grid)[0]))> f);
//    void set(std::function<ValueType (decltype((*_grid)[0]))> f);
    Grid1dObject( std::tuple<GridTypeN, GridType...> in);
//    ValueType operator()(decltype((*_grid)[0]) in ) const;
//    ValueType& operator[](unsigned int i);
//    template <unsigned int M, typename ValType, class GridType2> friend std::ostream& operator<<(std::ostream& lhs, const Grid1dObject<M,ValType,GridType2> &in);
};

template <typename ValueType, class GridType>
class Grid1dObject<1,ValueType, GridType> 
{
protected:
    std::shared_ptr<GridType> _grid;
    std::shared_ptr<VectorType<ValueType> > _vals;
    Grid1dObject();
public:
    Grid1dObject(GridType grid);
}; 

//
// N=1
//

template<typename ValueType, class GridType>
Grid1dObject<1,ValueType,GridType>::Grid1dObject()
{
    DEBUG("N=1");
    //_grid.reset(new GridType);
    _vals.reset(new VectorType<ValueType>);
}

template<typename ValueType, class GridType>
Grid1dObject<1,ValueType,GridType>::Grid1dObject(GridType grid):_grid(std::make_shared<GridType>(grid))
{
    DEBUG("Initializing from grid, N=1");
    _vals.reset(new VectorType<ValueType>);
}

//
// N!=1
//

template<unsigned int N,typename ValueType, typename GridTypeN, typename ...GridType> 
Grid1dObject<N,ValueType,GridTypeN, GridType...>::Grid1dObject()
{
    std::cout << "Default: N = " << sizeof...(GridType) << std::endl;
    //_grid.reset(new GridType);
    //_vals.reset(new VectorType<ValueType>);
}


template<unsigned int N,typename ValueType, typename GridTypeN, typename ...GridType> 
Grid1dObject<N,ValueType,GridTypeN, GridType...>::Grid1dObject( std::tuple<GridTypeN, GridType...> in)
{
    DEBUG("Grid constructor N = " << sizeof...(GridType)+1);
    int tuple_size = sizeof...(GridType)+1;
    //_vals.reset(new VectorType<Grid1dObject<ValueType,GridType...> >::constant(std::get<0>(in).getSize(),  ));
}

/*
template<unsigned int N, typename ValueType, class GridType>
Grid1dObject<N,ValueType,GridType>::Grid1dObject(std::shared_ptr<GridType> grid):_grid(grid)
{
    DEBUG("N=" << N);
    _vals.reset(new VectorType<ValueType>(_grid->getSize()));
}

template<unsigned int N, typename ValueType, class GridType>
ValueType& Grid1dObject<N,ValueType,GridType>::operator[](unsigned int i)
{
    return (*_vals)[i];
}

template<unsigned int N, typename ValueType, class GridType>
Grid1dObject<N,ValueType,GridType>::Grid1dObject(std::shared_ptr<GridType> grid,  std::function<ValueType (decltype((*grid)[0]))> f):Grid1dObject<N,ValueType,GridType>(grid)
{
    this->set(f);
}

template<unsigned int N, typename ValueType, class GridType>
void Grid1dObject<N,ValueType,GridType>::set(std::function<ValueType (decltype((*_grid)[0]))> f)
{
    assert(_vals->size() == (*_grid).getSize());
    for (unsigned int i=0; i<_vals->size(); ++i) (*_vals)[i]=f((*_grid)[i]);
}

template <unsigned int M, typename ValType, class GridType2> 
std::ostream& operator<<(std::ostream& lhs, const Grid1dObject<M,ValType,GridType2> &in)
{
    lhs << *(in._vals); 
    return lhs;
}
*/
}
#endif // endif::ifndef ___FK_GRID_OBJECT_H___

