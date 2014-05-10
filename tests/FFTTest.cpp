#include <numeric>

//#include <kmesh.hpp>
#include <gftools/enum_grid.hpp>
//#include <grid_object.hpp>
#include "DMFT.h"
#include "FFT.hpp"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;
typedef grid_object<complex_type,fmatsubara_grid> GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, real_type tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}


int main()
{
    std::srand(std::time(0));
    size_t ksize=rand()%10+12;
    INFO("Ksize = " << ksize << std::endl);
    kmesh k1(ksize);
    enum_grid rgrid(0,k1.size());

    typedef grid_object<complex_type, kmesh, kmesh> EK_complex;
    typedef grid_object<real_type, kmesh, kmesh> EK_real;

    EK_complex disp2d_complex(std::forward_as_tuple(k1,k1));
    EK_real disp2d_real(std::forward_as_tuple(k1,k1));
    disp2d_complex.tail_ = CubicTraits<2>(1.0).get_dispersion(); disp2d_complex.fill(disp2d_complex.tail_);
    disp2d_real.fill(CubicTraits<2>(1.0).get_dispersion());
    EK_complex disp2d_complex_backup(disp2d_complex);
    
    DEBUG(disp2d_real);

    grid_object<complex_type, enum_grid, enum_grid> t_matrix_complex(std::forward_as_tuple(rgrid,rgrid));
    auto t_matrix_complex2(t_matrix_complex);
    grid_object<int, enum_grid, enum_grid> t_matrix_int(std::forward_as_tuple(rgrid,rgrid));
    grid_object<int, enum_grid, enum_grid> t_matrix_int2(std::forward_as_tuple(rgrid,rgrid));

    t_matrix_complex.data() = run_fft(disp2d_complex.data(), FFTW_BACKWARD);
    t_matrix_complex2.data() = run_fft(disp2d_complex.data(), FFTW_FORWARD)/ksize/ksize;
    disp2d_complex.data() = run_fft(t_matrix_complex.data(), FFTW_FORWARD);

    if (!is_equal(disp2d_complex.diff(disp2d_complex_backup),0)) return EXIT_FAILURE;
    INFO("PASSED FFT1");

    t_matrix_int.fill(grid_object<int, enum_grid, enum_grid> ::function_type([&](int a1, int a2){return std::round(std::real(t_matrix_complex(a1,a2)));}));
    t_matrix_int2.fill(grid_object<int, enum_grid, enum_grid>::function_type([&](int a1, int a2){return std::round(std::real(t_matrix_complex2(a1,a2)));}));
    DEBUG(t_matrix_int);
    if ((t_matrix_int - t_matrix_int2).sum()!=0) return EXIT_FAILURE;
    INFO ("PASSED FFT2");

    if (t_matrix_int.sum()!=-4) return EXIT_FAILURE;
    INFO ("PASSED SUM TEST");

    if (t_matrix_int[0][1]!=-1) return EXIT_FAILURE;
    if (t_matrix_int[0][ksize-1]!=-1) return EXIT_FAILURE;
    if (t_matrix_int[1][0]!=-1) return EXIT_FAILURE;
    if (t_matrix_int[ksize-1][0]!=-1) return EXIT_FAILURE;
    INFO ("PASSED ELEMENTS TEST");


    grid_object<complex_type, kmesh> disp1d (std::forward_as_tuple(k1));
    disp1d.tail_ = CubicTraits<1>(1.0).get_dispersion(); disp1d.fill(disp1d.tail_);
    auto disp1d_backup(disp1d);
    disp1d.data() = run_fft(run_fft(disp1d.data(), FFTW_BACKWARD), FFTW_FORWARD);
    if (!is_equal(disp1d_backup.diff(disp1d),0)) return EXIT_FAILURE;
    INFO("PASSED FFT1d")

    disp2d_complex.data() = run_fft(run_fft(disp2d_complex.data(), FFTW_BACKWARD), FFTW_FORWARD);
    if (!is_equal(disp2d_complex.diff(disp2d_complex_backup),0)) return EXIT_FAILURE;
    INFO("PASSED FFT2d")

    grid_object<complex_type, kmesh, kmesh, kmesh> disp3d (std::forward_as_tuple(k1,k1,k1));
    disp3d.tail_ = CubicTraits<3>(1.0).get_dispersion(); disp3d.fill(disp3d.tail_);
    auto disp3d_backup(disp3d);
    if (!is_equal(disp3d_backup.diff(disp3d),0)) return EXIT_FAILURE;
    INFO("PASSED FFT3d")

    grid_object<complex_type, kmesh, kmesh, kmesh, kmesh> disp4d (std::forward_as_tuple(k1,k1,k1,k1));
    disp4d.tail_ = CubicTraits<4>(1.0).get_dispersion(); disp4d.fill(disp4d.tail_);
    auto disp4d_backup(disp4d);
    if (!is_equal(disp4d_backup.diff(disp4d),0)) return EXIT_FAILURE;
    INFO("PASSED FFT4d")

    return EXIT_SUCCESS;
}
