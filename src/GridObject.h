#ifndef ___FK_GRID_OBJECT_H___
#define ___FK_GRID_OBJECT_H___

#include "FKCommon.h"
#include "Logger.h"
#include "Grid.h"

namespace FK { 

template<typename ValueType, class GridType>
class Grid1dObject {
protected:
    GridType _grid;
    VectorType<ValueType> _vals;
public:
    Grid1dObject(const GridType& grid);
    Grid1dObject(const GridType& grid,  std::function<ValueType (decltype(grid[0]))> f);
    void set(std::function<ValueType (decltype(_grid[0]))> f);
    ValueType operator()(decltype(_grid[0]) in ) const;
    ValueType operator[](unsigned int i) const;
    template <typename ValType, class GridType2> friend std::ostream& operator<<(std::ostream& lhs, const Grid1dObject<ValType,GridType2> &in);
};

template<typename ValueType, class GridType>
Grid1dObject<ValueType,GridType>::Grid1dObject(const GridType& grid,  std::function<ValueType (decltype(grid[0]))> f):_grid(grid)
{
    _vals.resize(grid.getSize());
}

template<typename ValueType, class GridType>
void Grid1dObject<ValueType,GridType>::set(std::function<ValueType (decltype(_grid[0]))> f)
{
    assert(_vals.size() == _grid.getSize());
    for (unsigned int i=0; i<_vals.size(); ++i) _vals[i]=f(_grid[i]);
}

template <typename ValType, class GridType2> 
std::ostream& operator<<(std::ostream& lhs, const Grid1dObject<ValType,GridType2> &in)
{
    lhs << in._vals; 
    return lhs;
}

}
#endif // endif::ifndef ___FK_GRID_OBJECT_H___

