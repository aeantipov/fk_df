#ifndef __fk_DF_h_213al5
#define __fk_DF_h_213al5

#include <Container.hpp>
#include <fftw3.h>

namespace FK {

template <typename ValueType, size_t D, typename BC, typename std::enable_if<D==2, bool>::type=0> 
Container<ValueType,D> run_fft (const ContainerBase<ValueType,D,BC> &in)
{
    MatrixType<ComplexType> kdata = in.getAsMatrix();
    MatrixType<ComplexType> out(kdata);
    fftw_plan p;
    p = fftw_plan_dft_2d(kdata.rows(), kdata.cols(), reinterpret_cast<fftw_complex*>(kdata.data()), reinterpret_cast<fftw_complex*>(out.data()), FFTW_BACKWARD, FFTW_ESTIMATE); 
    fftw_execute(p);
    out/=(kdata.rows()*kdata.cols());
    return Container<ValueType,2>(out);
}

template <typename ValueType, size_t D, typename BC, typename std::enable_if<D==3, bool>::type=0> 
Container<ValueType,D> run_fft (const ContainerBase<ValueType,D,BC> &in)
{
    Container<ValueType,D> out(in);
    fftw_plan p;
    std::array<size_t,3> shape = {{ out._data.shape()[0], out._data.shape()[1], out._data.shape()[2] }};
    p = fftw_plan_dft_3d(shape[0], shape[1], shape[2], 
                         reinterpret_cast<fftw_complex*>( in._data.origin()), 
                         reinterpret_cast<fftw_complex*>(out._data.origin()), 
                         FFTW_BACKWARD, FFTW_ESTIMATE); 
    RealType norm=1.0*shape[0]*shape[1]*shape[2];
    fftw_execute(p);
    out/=norm;
    return out;
}


template <typename ValueType, size_t D, typename BC, typename std::enable_if<D>=4 || D==1, bool>::type=0> 
Container<ValueType,D> run_fft (const ContainerBase<ValueType,D,BC> &in)
{
    ERROR("No FFT defined for D="<<D);
    return in; 
}

} // end of namespace FK
#endif
