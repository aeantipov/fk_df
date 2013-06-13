#ifndef __fk_DF_h_213al5
#define __fk_DF_h_213al5

#include <Container.hpp>
#include <fftw3.h>

namespace FK {

template <size_t D, typename BC, typename std::enable_if<D==1, bool>::type=0> 
Container<ComplexType,D> run_fft (const ContainerBase<ComplexType,D,BC> &in, int direction)
{
    Container<ComplexType,D> out(in);
    fftw_plan p;
    std::array<size_t,1> shape = {{ out._data.shape()[0] }}; 
    p = fftw_plan_dft_1d(shape[0], 
                         reinterpret_cast<fftw_complex*>( in._data.origin()), 
                         reinterpret_cast<fftw_complex*>(out._data.origin()), 
                         direction, FFTW_ESTIMATE); 
    RealType norm=1.0*shape[0];
    fftw_execute(p);
    if (direction == FFTW_BACKWARD) out/=norm;
    return out;
}


template <size_t D, typename BC, typename std::enable_if<D==2, bool>::type=0> 
Container<ComplexType,D> run_fft (const ContainerBase<ComplexType,D,BC> &in, int direction)
{
    Container<ComplexType,D> out(in);
    fftw_plan p;
    std::array<size_t,2> shape = {{ out._data.shape()[0], out._data.shape()[1] }};
    p = fftw_plan_dft_2d(shape[0], shape[1],  
                         reinterpret_cast<fftw_complex*>( in._data.origin()), 
                         reinterpret_cast<fftw_complex*>(out._data.origin()), 
                         direction, FFTW_ESTIMATE); 
    RealType norm=1.0*shape[0]*shape[1];
    fftw_execute(p);
    if (direction == FFTW_BACKWARD) out/=norm;
    return out;
}

template <size_t D, typename BC, typename std::enable_if<D==3, bool>::type=0> 
Container<ComplexType,D> run_fft (const ContainerBase<ComplexType,D,BC> &in, int direction)
{
    Container<ComplexType,D> out(in);
    fftw_plan p;
    std::array<size_t,3> shape = {{ out._data.shape()[0], out._data.shape()[1], out._data.shape()[2] }};
    p = fftw_plan_dft_3d(shape[0], shape[1], shape[2], 
                         reinterpret_cast<fftw_complex*>( in._data.origin()), 
                         reinterpret_cast<fftw_complex*>(out._data.origin()), 
                         direction, FFTW_ESTIMATE); 
    RealType norm=1.0*shape[0]*shape[1]*shape[2];
    fftw_execute(p);
    if (direction == FFTW_BACKWARD) out/=norm;
    return out;
}

template <size_t D, typename BC, typename std::enable_if<D==4, bool>::type=0> 
Container<ComplexType,D> run_fft (const ContainerBase<ComplexType,D,BC> &in, int direction)
{
    Container<ComplexType,D> out(in);
    fftw_plan p;
    const std::array<int,4> shape = {{ out._data.shape()[0], out._data.shape()[1], out._data.shape()[2], out._data.shape()[3] }};
    p = fftw_plan_dft(4, shape.data(),
                         reinterpret_cast<fftw_complex*>( in._data.origin()), 
                         reinterpret_cast<fftw_complex*>(out._data.origin()), 
                         direction, FFTW_ESTIMATE); 
    RealType norm=1.0*shape[0]*shape[1]*shape[2]*shape[3];
    fftw_execute(p);
    if (direction == FFTW_BACKWARD) out/=norm;
    return out;
}



template <size_t D, typename BC, typename std::enable_if<D>=5, bool>::type=0> 
Container<ComplexType,D> run_fft (const ContainerBase<ComplexType,D,BC> &in, int direction)
{
    ERROR("No FFT defined for D="<<D);
    return in; 
}

} // end of namespace FK
#endif
