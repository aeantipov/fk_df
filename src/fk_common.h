#ifndef ___FK_FK_H___
#define ___FK_FK_H___

#include<iostream>
#include<complex>
#include<string>
#include<vector>
#include<list>
#include<map>

#include<memory>
#include<utility>

#include<Eigen/Core>

#define REALTYPE_DOUBLE

namespace FK {

typedef double RealType;
typedef std::complex<RealType> ComplexType;

/** A short name for imaginary unit. */
static const ComplexType I = ComplexType(0.0,1.0);    // 'static' to prevent linking problems

/** Dense complex matrix. */
typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> MatrixType;
/** Dense real matrix. */
typedef Eigen::Matrix<RealType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> RealMatrixType;

/** Dense complex vector. */
typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,1,Eigen::AutoAlign> ComplexVectorType;
/** Dense real vector. */
typedef Eigen::Matrix<RealType,Eigen::Dynamic,1,Eigen::AutoAlign> RealVectorType;
/** Dense vector of integers. */
typedef Eigen::Matrix<int,Eigen::Dynamic,1,Eigen::AutoAlign> IntVectorType;


} // end namespace FK

#endif // endif::ifndef ___FK_FK_H___
