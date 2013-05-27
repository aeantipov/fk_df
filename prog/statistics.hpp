#include <fftw3.h>

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}

void sighandler(int signal)
{
    static size_t count = 0;
    count++;
    INFO("Caught INTERRUPT, signal " << signal <<" " << count << " times. ")
    INTERRUPT = true;
    if (count >= 3) { INFO("Force exiting"); exit(signal); }
}

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

template <typename ValueType, size_t D, typename BC, typename std::enable_if<D!=2, bool>::type=0> 
Container<ValueType,D> run_fft (const ContainerBase<ValueType,D,BC> &in)
{
    ERROR("No FFT defined for D="<<D);
    return in; 
}


