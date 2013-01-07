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
    _grids(in),
    _f(__fun_traits<FunctionType>::constant(0.0))
{
    GetGridSizes<N>::TupleSizeToArray(_grids,_dims);   
    _data.reset(new Container<sizeof...(GridTypes),ValueType>(_dims));
}

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>::GridObject( const GridObject<ValueType, GridTypes...>& rhs):
    _grids(rhs._grids), 
    _data(new Container<sizeof...(GridTypes), ValueType>(*(rhs._data))),
    _f(rhs._f)
{
}; 

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>::GridObject( GridObject<ValueType,GridTypes...> && rhs):
    _grids(rhs._grids),
    _f(rhs._f)
{
    _data.swap(rhs._data);
}

template <typename ValueType, typename ...GridTypes> 
const std::tuple<GridTypes...> GridObject<ValueType,GridTypes...>::getGrids() const 
{ 
    return _grids; 
};

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
    try { return ContainerExtractor<sizeof...(GridTypes), ArgTypes...>::get(*_data,_grids,in...); }
    //catch (std::exception &e) { 
    catch (...) { 
        #ifndef NDEBUG
        DEBUG("Using analytical expression");
        #endif
        return _f(in...);
     };

}

template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline ValueType GridObject<ValueType,GridTypes...>::operator()(const std::tuple<ArgTypes...>& in) const
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "GridObject call, number of input parameters mismatch."); 
    std::function<ValueType(ArgTypes...)> f1 = [&](ArgTypes... in1)->ValueType{return this->operator()(in1...); };// ContainerExtractor<sizeof...(GridTypes), ArgTypes...>::get(*_data,_grids,in...);};
    __caller<ValueType,ArgTypes...> t = {in,f1};

    return t.call();
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

template <class F> struct PointFToF;
template <typename ValueType, typename ... PointTypes > 
struct PointFToF<std::function<ValueType(PointTypes...)>> {
    template <typename ...ArgTypes> static std::function<ValueType(ArgTypes...)> convert(const std::function<ValueType(PointTypes...)>& in);
};
/*
template <size_t N, template <size_t, typename ...> class, typename ... > struct __genContainerExtractor;
template <size_t N, template <size_t, typename ...> class T, template <typename ...> class B, typename GridType1, typename ... GridTypes, typename ...ArgTypes> 
struct __genContainerExtractor<N,T,B<GridType1, GridTypes...>, ArgTypes... > : __genContainerExtractor<N-1,T, B<GridTypes...>, typename GridType1::point, ArgTypes...> {} ;
template <template <size_t, typename ...> class T, template <typename ...> class B, typename ... ArgTypes> 
struct __genContainerExtractor<0,T,B<>, ArgTypes...> { typedef T<sizeof...(ArgTypes),ArgTypes...> type; };

template <typename ValueType, typename ...GridTypes> 
inline void GridObject<ValueType,GridTypes...>::fill(const typename GridObject<ValueType,GridTypes...>::PointFunctionType& in)
{

    __genContainerExtractor<sizeof...(GridTypes), ContainerExtractor, std::tuple<GridTypes...>>::type::set(*_data,_grids,in);

    //ContainerExtractor<sizeof...(GridTypes), ArgTypes...>::set(*_data,_grids,in);
    //_f = in;
}
*/


template <size_t N, template <size_t, typename ...> class, typename ... > struct __genContainerExtractorForFunctionType;
template <size_t N, template <size_t, typename ...> class T, template <typename ...> class B, typename GridType1, typename ... GridTypes, typename ...ArgTypes> 
struct __genContainerExtractorForFunctionType<N,T,B<GridType1, GridTypes...>, ArgTypes... > :
     __genContainerExtractorForFunctionType<N-1,T, B<GridTypes...>, decltype(GridType1::point::_val), ArgTypes...> {} ;
template <template <size_t, typename ...> class T, template <typename ...> class B, typename ... ArgTypes> 
struct __genContainerExtractorForFunctionType<0,T,B<>, ArgTypes...> { typedef T<sizeof...(ArgTypes),ArgTypes...> type; };

template <typename ValueType, typename ...GridTypes> 
inline void GridObject<ValueType,GridTypes...>::fill(const typename GridObject<ValueType,GridTypes...>::FunctionType& in)
{
    __genContainerExtractorForFunctionType<sizeof...(GridTypes), ContainerExtractor, std::tuple<GridTypes...>>::type::set(*_data,_grids,in);
    _f = in;
}



template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline void GridObject<ValueType,GridTypes...>::fill_tuple(const std::function<ValueType(std::tuple<ArgTypes...>)> & in)
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "GridObject fill, number of input parameters mismatch."); 
    std::function<ValueType(ArgTypes...)> f1 = [&](ArgTypes... in1)->ValueType{return in(std::make_tuple(in1...)); };
    this->fill(f1);
}


/*
template <typename ValueType, typename ...GridTypes> 
template <template <typename, class> class Filler, typename ...ArgTypes> 
inline void GridObject<ValueType,GridTypes...>::fill(const Filler<ValueType, ArgTypes...> &in)
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "GridObject fill, number of input parameters mismatch."); 
    ContainerExtractor<sizeof...(GridTypes), ArgTypes...>::set(*_data,_grids,in);
}
*/

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
    INFO("Saving " << typeid(*this).name() << " to " << fname);
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
    static const RealType read_tol = 1e-6;
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

/*
template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator= (
    const std::function<ValueType(ArgTypes...)> & in)
{
    this->fill(in);
    return *this;
}
*/

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    *_data=*(rhs._data);
    _f = rhs._f;
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator= (
    const ValueType& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    *_data=rhs;
    _f = __fun_traits<FunctionType>::constant(rhs);
    return *this;
}



template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator+= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    *_data+=*(rhs._data);
    _f=__fun_traits<FunctionType>::add(_f, rhs._f);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator+= (
    const ValueType & rhs)
{
    *_data+=rhs;
    _f=__fun_traits<FunctionType>::add(_f, __fun_traits<FunctionType>::constant(rhs));
    return *this;
}


template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator*= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    *_data*=*(rhs._data);
    _f=__fun_traits<FunctionType>::multiply(_f, rhs._f);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator*= (
    const ValueType & rhs)
{
    *_data*=rhs;
    _f=__fun_traits<FunctionType>::multiply(_f, __fun_traits<FunctionType>::constant(rhs));
    return *this;
}


template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator/= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    *_data/=*(rhs._data);
    _f=__fun_traits<FunctionType>::divide(_f, rhs._f);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator/= (
    const ValueType & rhs)
{
    *_data/=rhs;
    _f=__fun_traits<FunctionType>::divide(_f, __fun_traits<FunctionType>::constant(rhs));
    return *this;
}


template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator-= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    *_data-=*(rhs._data);
    _f=__fun_traits<FunctionType>::subtract(_f, rhs._f);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator-= (
    const ValueType & rhs)
{
    *_data-=rhs;
    _f=__fun_traits<FunctionType>::subtract(_f, __fun_traits<FunctionType>::constant(rhs));
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
