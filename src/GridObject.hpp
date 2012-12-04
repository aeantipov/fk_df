#ifndef ___FK_GRIDOBJECT_HPP___
#define ___FK_GRIDOBJECT_HPP___

#include "GridObject.h"

namespace FK {
//
// GridObject::ContainerExtractor
//

template< typename ValueType, typename ...GridTypes> 
template <size_t Nc, typename ArgType1, typename ...ArgTypes>
inline ValueType& GridObject<ValueType,GridTypes...>::ContainerExtractor<Nc,ArgType1,ArgTypes...>::get(
    Container<Nc, ValueType> &data, const std::tuple<GridTypes...> &grids, 
    const ArgType1& arg1, const ArgTypes&... args) 
{
    const auto & grid=std::get<N-Nc>(grids);
    auto &tmp = grid.getValue(data, arg1);
    return ContainerExtractor<Nc-1,ArgTypes...>::get(tmp,grids,args...);
};

template< typename ValueType, typename ...GridTypes>
template <typename ArgType1> 
inline ValueType& GridObject<ValueType,GridTypes...>::ContainerExtractor<1,ArgType1>::get(
    Container<1, ValueType> &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1)
{
    const auto & grid=std::get<N-1>(grids);
    auto &tmp = grid.getValue(data, arg1);
    return tmp;
}

template< typename ValueType, typename ...GridTypes> 
template <size_t Nc, typename ArgType1, typename ...ArgTypes>
inline void GridObject<ValueType,GridTypes...>::ContainerExtractor<Nc,ArgType1,ArgTypes...>::set(
    Container<Nc, ValueType> &data, 
    const std::tuple<GridTypes...> &grids, 
    const std::function<ValueType(ArgType1, ArgTypes...)> &f)
{
    const auto & grid=std::get<N-Nc>(grids);
    const size_t grid_size = grid.getSize();
    const auto& grid_vals = grid.getVals();
    static_assert(std::is_convertible<decltype(grid_vals[0]), ArgType1>::value, "!");
    for (size_t i=0; i<grid_size; ++i) { 
        const auto& cur_val = grid_vals[i];
        const auto f1 = [&f,&cur_val](const ArgTypes&... Args){return f(cur_val,Args...);};
        ContainerExtractor<Nc-1, ArgTypes...>::set(data[i],grids,f1);
    }
}
 
template< typename ValueType, typename ...GridTypes> 
template <typename ArgType1> 
inline void GridObject<ValueType,GridTypes...>::ContainerExtractor<1,ArgType1>::set(
    Container<1, ValueType> &data, 
    const std::tuple<GridTypes...> &grids, 
    const std::function<ValueType(ArgType1)> &f)
{
    const auto & grid=std::get<N-1>(grids);
    const size_t grid_size = grid.getSize();
    const auto& grid_vals = grid.getVals();
    static_assert(std::is_convertible<decltype(grid_vals[0]), ArgType1>::value, "!");
    for (size_t i=0; i<grid_size; ++i) { 
        data[i]=f(grid_vals[i]);
        };
}
//
// GridObject
//

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>::GridObject( const std::tuple<GridTypes...> &in):
    _grids(in)
{
    GetGridSizes<N>::TupleSizeToArray(_grids,_dims);   
    _data.reset(new Container<sizeof...(GridTypes),ValueType>(_dims));
}


template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>::GridObject( GridObject<ValueType,GridTypes...> && rhs):_grids(rhs._grids)
{
    _data.swap(rhs._data);
}

template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline ValueType& GridObject<ValueType,GridTypes...>::operator()(const ArgTypes&... in)
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "GridObject call, number of input parameters mismatch."); 
    return ContainerExtractor<sizeof...(GridTypes), ArgTypes...>::get(*_data,_grids,in...);
}

template <typename ValueType, typename ...GridTypes> 
std::ostream& operator<<(std::ostream& lhs, const GridObject<ValueType,GridTypes...> &in)
{
    lhs << *(in._data);
    return lhs;
}


template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline void GridObject<ValueType,GridTypes...>::fill(const std::function<ValueType(ArgTypes...)> & in)
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "GridObject fill, number of input parameters mismatch."); 
    ContainerExtractor<sizeof...(GridTypes), ArgTypes...>::set(*_data,_grids,in);
}

//
// Operators
//

template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator= (
    const std::function<ValueType(ArgTypes...)> & in)
{
    this->fill(in);
    return *this;
}
template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    *_data=*(rhs._data);
    return *this;
}


template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator+= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    *_data+=*(rhs._data);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator*= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    *_data*=*(rhs._data);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator-= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    *_data-=*(rhs._data);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator+ (
    const GridObject<ValueType,GridTypes...>& rhs) const
{
    GridObject out(*this);
    out+=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator* (
    const GridObject<ValueType,GridTypes...>& rhs) const
{
    GridObject out(*this);
    out*=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator- (
    const GridObject<ValueType,GridTypes...>& rhs) const
{
    GridObject out(*this);
    out-=rhs;
    return out;
}

} // end of namespace FK
#endif
