#ifndef ___FK_FK_H___
#define ___FK_FK_H___

#include<iostream>
#include<complex>
#include<string>
#include<vector>
#include<list>
#include<map>
#include<array>
#include<iomanip>

#include<memory>
#include<utility>
#include<functional>

#include<Eigen/Core>
#include<Eigen/StdVector>

#include "EigenIterator.h"

#define REALTYPE_DOUBLE

namespace FK {

typedef double RealType;
typedef std::complex<RealType> ComplexType;

template <typename T>
using VectorType = Eigen::Matrix<T,Eigen::Dynamic, 1>;

template <typename T>
using MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;

/** A short name for imaginary unit. */
static const ComplexType I = ComplexType(0.0,1.0);    // 'static' to prevent linking problems

const RealType PI = std::atan(1.0)*4;
/** Dense complex matrix. */
typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> ComplexMatrixType;
/** Dense real matrix. */
typedef Eigen::Matrix<RealType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> RealMatrixType;

/** Dense complex vector. */
typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,1,Eigen::AutoAlign> ComplexVectorType;
/** Dense real vector. */
typedef Eigen::Matrix<RealType,Eigen::Dynamic,1,Eigen::AutoAlign> RealVectorType;
/** Dense vector of integers. */
typedef Eigen::Matrix<int,Eigen::Dynamic,1,Eigen::AutoAlign> IntVectorType;

template <bool Fermion> inline ComplexType Matsubara(int n, RealType beta){return PI*I/beta*ComplexType(2*n+Fermion);};

inline ComplexType FMatsubara(int n, RealType beta){return Matsubara<1>(n,beta);};
inline ComplexType BMatsubara(int n, RealType beta){return Matsubara<0>(n,beta);};

template <typename T> 
struct __num_format {
    const int _prec = 12;
    T _v;
    __num_format(T v):_v(v){};
    operator T(){return _v;};
    friend std::ostream& operator<<(std::ostream& lhs, const __num_format<T> &in){lhs << std::setprecision(in._prec) << in._v; return lhs;};
    friend std::istream& operator>>(std::istream& lhs, __num_format<T> &out){lhs >> out._v; return lhs;};
};

inline std::ostream& operator<<(std::ostream& lhs, const __num_format<ComplexType> &in){lhs << std::setprecision(in._prec) << real(in._v) << " " << imag(in._v); return lhs;};
inline std::istream& operator>>(std::istream& lhs, __num_format<ComplexType> &out){RealType re,im; lhs >> re; lhs >> im; out._v = re+I*im; return lhs;};

/** A tool to generate a type for an object T of D Args. */
template <size_t N, typename ArgType1, template <typename ...> class T, typename ...ArgTypes>  
struct ArgBackGenerator:ArgBackGenerator<N-1,ArgType1,T,ArgTypes...,ArgType1 > {};

template <typename ArgType1, template <typename ...> class T, typename ...ArgTypes> 
struct ArgBackGenerator<1, ArgType1, T, ArgTypes...> { typedef T<ArgTypes..., ArgType1> type; };

/** A tool to generate a type for a function of N Args. */
template <size_t N, typename ValueType, typename ArgType1, typename ...ArgTypes>  
struct ArgFunGenerator : ArgFunGenerator<N-1,ValueType,ArgType1, ArgTypes...,ArgType1 >{static_assert(N>1,"");};

template <typename ValueType, typename ArgType1, typename ...ArgTypes> 
struct ArgFunGenerator<1, ValueType, ArgType1, ArgTypes...> { 
    typedef std::function<ValueType(ArgTypes..., ArgType1)> type; };

/** A tool to calc an integer power function of an int. */
template<int base, unsigned exponent >
struct __power {
    enum { value = base * __power<base, exponent - 1>::value };
};
template< int base >
struct __power<base,0> {
    enum { value = 1 };
};

/** Function traits. */
template <typename FunctionType> struct __fun_traits;
template <typename ValType, typename ... ArgTypes> 
struct __fun_traits<std::function<ValType(ArgTypes...)> >
{
    static std::function<ValType(ArgTypes...)> constant(const ValType &c) 
    { return [c](ArgTypes...in){return c;};}
    static std::function<ValType(ArgTypes...)> add(std::function<ValType(ArgTypes...)> f1, std::function<ValType(ArgTypes...)> f2)
    { return [f1,f2](ArgTypes... in){return f1(in...)+f2(in...);}; }
    static std::function<ValType(ArgTypes...)> multiply(std::function<ValType(ArgTypes...)> f1, std::function<ValType(ArgTypes...)> f2)
    { return [f1,f2](ArgTypes... in){return f1(in...)*f2(in...);}; }
    static std::function<ValType(ArgTypes...)> subtract(std::function<ValType(ArgTypes...)> f1, std::function<ValType(ArgTypes...)> f2)
    { return [f1,f2](ArgTypes... in){return f1(in...)-f2(in...);}; }
    static std::function<ValType(ArgTypes...)> divide(std::function<ValType(ArgTypes...)> f1, std::function<ValType(ArgTypes...)> f2)
    { return [f1,f2](ArgTypes... in){return f1(in...)*f2(in...);}; }
    static std::function<ValType(ArgTypes...)> getFromTupleF(const std::function<ValType(std::tuple<ArgTypes...>)>& f1)
    { return [f1](ArgTypes...in){return f1(in...);};}
};


/** A tool to wrap a call a class method from a tuple. */
template<int ...> struct __seq {};
template<int N, int ...S> struct __gens : __gens<N-1, N-1, S...> {};
template<int ...S> struct __gens<0, S...>{ typedef __seq<S...> type; };

template <typename ReturnType, typename ...Args> struct __caller { 
    std::tuple<Args...> _params;
    std::function<ReturnType(Args...)> _f;
    template<int ...S> ReturnType _callf(__seq<S...>) { return _f(std::get<S>(_params)...); };
    ReturnType call(){ return _callf(typename __gens<sizeof...(Args)>::type()); };
};

} // end namespace FK

#endif // endif::ifndef ___FK_FK_H___
