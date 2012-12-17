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

inline ComplexType FMatsubara(int n, RealType beta){return PI*I/beta*ComplexType(2*n+1);};
inline ComplexType BMatsubara(int n, RealType beta){return PI*I/beta*ComplexType(2*n);};

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

/* A tool to generate an array of grid sizes from a given tuple of grids. */
template <size_t N>
struct GetGridSizes {
    template <typename... GridType, size_t M>
    static inline void TupleSizeToArray( const std::tuple<GridType...>& in, std::array<size_t, M> &out ) {
        static_assert(N>1,"!");
        std::get<N-1>(out) = std::get<N-1>(in).getSize();
        GetGridSizes<N-1>::TupleSizeToArray( in, out );
    }
};

template <>
template <typename... GridType, size_t M>
inline void GetGridSizes<1>::TupleSizeToArray( const std::tuple<GridType...>& in, std::array<size_t, M> &out ) {
    std::get<0>(out) = std::get<0>(in).getSize();
}

template<typename T> 
struct function_traits;  

template<typename R, typename ...Args> 
struct function_traits<std::function<R(Args...)>>
{
    static const size_t nargs = sizeof...(Args);

    typedef R result_type;

    template <size_t i>
    struct arg
    {
        typedef typename std::tuple_element<i, std::tuple<Args...>>::type type;
    };
};

/** A tool to generate a type for an object T of D Args. */
    template <size_t N, typename ...> struct ArgGenerator;
    template <size_t N, typename ArgType1, template <typename ...> class T, typename ...ArgTypes>  
        struct ArgGenerator<N,ArgType1,T<ArgTypes...>,ArgTypes...>:
            ArgGenerator<N-1,ArgType1,T<ArgTypes...>,ArgTypes...,ArgType1 >{static_assert(N>1,"");};

    template <typename ArgType1, template <typename ...> class T, typename ...ArgTypes> 
        struct ArgGenerator<1, ArgType1, T<ArgTypes...>, ArgTypes...> 
        { typedef T<ArgType1, ArgTypes...> type; };

/** A tool to generate a type for a function of N Args. */
    //template <size_t N, typename ...> struct ArgGenerator;
    template <size_t N, typename ValueType, typename ArgType1, typename ...ArgTypes>  
        struct ArgFunGenerator : ArgFunGenerator<N-1,ValueType,ArgType1, ArgTypes...,ArgType1 >{static_assert(N>1,"");};

    template <typename ValueType, typename ArgType1, typename ...ArgTypes> struct ArgFunGenerator<1, ValueType, ArgType1, ArgTypes...> 
        { typedef std::function<ValueType(ArgTypes..., ArgType1)> type; };

/** A tool to calc an integer power function of an int. */
template<int base, unsigned exponent >
struct __power {
    enum { value = base * __power<base, exponent - 1>::value };
};
template< int base >
struct __power<base,0> {
    enum { value = 1 };
};
} // end namespace FK

#endif // endif::ifndef ___FK_FK_H___
