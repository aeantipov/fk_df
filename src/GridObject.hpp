#ifndef ___FK_GRIDOBJECT_HPP___
#define ___FK_GRIDOBJECT_HPP___

#include "GridObject.h"
#include <iomanip>

namespace FK {
//
// Specification of integration with GridObject
//


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
auto GridObject<ValueType,GridTypes...>::operator[](size_t i)->decltype((*_data)[i])
{
    return (*_data)[i];
}

/*
template <typename ValueType, typename ...GridTypes> 
template <int M> 
ValueType& GridObject<ValueType,GridTypes...>::operator[](const std::array<size_t,M>& in)
{
}*/

template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline ValueType& GridObject<ValueType,GridTypes...>::get(const ArgTypes&... in)
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "GridObject call, number of input parameters mismatch."); 
    return ContainerExtractor<sizeof...(GridTypes), ArgTypes...>::get(*_data,_grids,in...);
}

template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline ValueType GridObject<ValueType,GridTypes...>::operator()(const ArgTypes&... in) const
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
template <size_t M>
inline auto GridObject<ValueType,GridTypes...>::getGrid() const -> const decltype(std::get<M>(_grids))
{
    return std::get<M>(_grids);
}



template <typename ValueType, typename ...GridTypes> 
inline auto GridObject<ValueType,GridTypes...>::getGrid() const -> const decltype(std::get<0>(_grids))
{
    return std::get<0>(_grids);
}


template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline void GridObject<ValueType,GridTypes...>::fill(const std::function<ValueType(ArgTypes...)> & in)
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "GridObject fill, number of input parameters mismatch."); 
    ContainerExtractor<sizeof...(GridTypes), ArgTypes...>::set(*_data,_grids,in);
}

template <typename ValueType, typename ...GridTypes> 
template <template <typename, class> class Filler, typename ...ArgTypes> 
inline void GridObject<ValueType,GridTypes...>::fill(const Filler<ValueType, ArgTypes...> &in)
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "GridObject fill, number of input parameters mismatch."); 
    ContainerExtractor<sizeof...(GridTypes), ArgTypes...>::set(*_data,_grids,in);
}


template <typename ValueType, typename ...GridTypes> 
template <typename U, typename std::enable_if<std::is_same<U, ComplexType>::value, int>::type>
inline GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::conj()
{
    GridObject<ValueType,GridTypes...> out(*this);
    *(out._data) = out._data->conj();
    return out;
}

template <typename ValueType, typename ...GridTypes> 
inline ValueType GridObject<ValueType,GridTypes...>::sum()
{
    return _data->sum();
}

template <typename ValueType, typename ...GridTypes> 
inline void GridObject<ValueType,GridTypes...>::savetxt(const std::string& fname) const
{
    std::ofstream out;
    out.open(fname.c_str());
    for (auto x : std::get<0>(_grids).getVals())
        {
            out << std::scientific << __num_format<decltype(x)>(x) << "    " << __num_format<ValueType>((*this)(x)) << std::endl;
        }
    out.close();
}

template <typename ValueType, typename ...GridTypes> 
inline void GridObject<ValueType,GridTypes...>::loadtxt(const std::string& fname)
{
    std::ifstream in;
    static const RealType read_tol = 1e-3;
    in.open(fname.c_str());
    if (!in.good()) { throw exIOProblem(); }

    for (auto x : std::get<0>(_grids).getVals())
        {
            __num_format<decltype(x)> tmp(x);
            in >> tmp;
            if (std::abs(tmp._v._val-ValueType(x))>read_tol) { ERROR("loadtxt - grid mismatch"); throw exIOProblem(); };
             __num_format<ValueType> tmp2(this->get(x));
            in >> tmp2;
            this->get(x) = tmp2._v;
        }
    in.close();
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
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator= (
    const ValueType& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    *_data=rhs;
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
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator+= (
    const ValueType & rhs)
{
    *_data+=rhs;
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
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator*= (
    const ValueType & rhs)
{
    *_data*=rhs;
    return *this;
}


template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator/= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    *_data/=*(rhs._data);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator/= (
    const ValueType & rhs)
{
    *_data/=rhs;
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
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator-= (
    const ValueType & rhs)
{
    *_data-=rhs;
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
inline GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator+ (
    const ValueType & rhs) const
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
inline GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator* (
    const ValueType & rhs) const
{
    GridObject out(*this);
    out*=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator/ (
    const GridObject<ValueType,GridTypes...>& rhs) const
{
    GridObject out(*this);
    out/=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator/ (
    const ValueType & rhs) const
{
    GridObject out(*this);
    out/=rhs;
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

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator- (
    const ValueType & rhs) const
{
    GridObject out(*this);
    out-=rhs;
    return out;
}


} // end of namespace FK
#endif
