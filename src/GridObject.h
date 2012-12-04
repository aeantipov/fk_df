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
template< typename ValueType, typename ...GridTypes> 
class GridObject 
{
protected:
    static const size_t N = sizeof...(GridTypes);
    std::tuple<GridTypes...> _grids;
    std::array<size_t, N> _dims;
    template <size_t Nc, typename ArgType1, typename ...ArgTypes> struct ContainerExtractor {
        static ValueType& get(Container<Nc, ValueType> &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1, const ArgTypes&... args);
        static void set(Container<Nc, ValueType> &data, const std::tuple<GridTypes...> &grids, const std::function<ValueType(ArgType1, ArgTypes...)> &f);
             };
    template <typename ArgType1> struct ContainerExtractor<1,ArgType1> {
        static ValueType& get(Container<1, ValueType> &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1); 
        static void set(Container<1, ValueType> &data, const std::tuple<GridTypes...> &grids, const std::function<ValueType(ArgType1)> &f);
        };

    std::shared_ptr<Container<N, ValueType>> _data;
public:

    /** Constructs a grid object out of a tuple containing various grids. */
    GridObject( const std::tuple<GridTypes...> &grids);
    /** Constructor of grids and data. */
    GridObject( const std::tuple<GridTypes...> &grids, const Container<sizeof...(GridTypes), ValueType>& data):
        _grids(grids), 
        _data(std::make_shared<Container<sizeof...(GridTypes), ValueType>> (data)){};
    /** Copy constructor. */
    GridObject( const GridObject<ValueType, GridTypes...>& rhs):_grids(rhs._grids), _data(rhs._data){}; 
    /** Move constructor. */
    GridObject( GridObject<ValueType, GridTypes...>&& rhs);

    /** Returns element number i, which corresponds to (*_grid)[i]. */
    //auto operator[](unsigned int i)->decltype((*_data)[0]);
    Container<sizeof...(GridTypes), ValueType>& getData(){return *_data;};
    template <typename ...ArgTypes> void fill(const std::function<ValueType(ArgTypes...)> &);

    template <int M> ValueType operator[](const std::array<size_t,M>& in) const;
    template <typename ...ArgTypes> ValueType& operator()(const ArgTypes&... in);
    template <typename ...ArgTypes> GridObject& operator= (const std::function<ValueType(ArgTypes...)> &);
    GridObject& operator= (const GridObject & rhs);
    GridObject& operator*= (const GridObject & rhs);
    GridObject operator* (const GridObject & rhs) const;
    GridObject& operator+= (const GridObject & rhs);
    GridObject operator+ (const GridObject & rhs) const;
    GridObject& operator-= (const GridObject & rhs);
    GridObject operator- (const GridObject & rhs) const;

    template <typename ValType, class ...GridTypes2> friend std::ostream& operator<<(std::ostream& lhs, const GridObject<ValType,GridTypes2...> &in);
};

} // end of namespace FK
#endif // endif::ifndef ___FK_GRID_OBJECT_H___

#include "GridObject.hpp"

