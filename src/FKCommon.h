#ifndef ___FK_FK_H___
#define ___FK_FK_H___

#include<iostream>
#include<complex>
#include<string>
#include<vector>
#include<list>
#include<map>
#include<array>

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

} // end namespace FK

#endif // endif::ifndef ___FK_FK_H___
