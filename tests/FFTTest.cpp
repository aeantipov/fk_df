#include <numeric>

#include <KMesh.hpp>
#include <EnumerateGrid.hpp>
#include <GridObject.hpp>
#include "DMFT.h"
#include "FFT.hpp"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;
typedef GFWrap GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}


int main()
{
    std::srand(std::time(0));
    size_t ksize=rand()%10+12;
    INFO("Ksize = " << ksize << std::endl);
    KMesh k1(ksize);
    EnumerateGrid rgrid(0,k1.getSize());

    typedef GridObject<ComplexType, KMesh, KMesh> EK_complex;
    typedef GridObject<RealType, KMesh, KMesh> EK_real;

    EK_complex disp2d_complex(std::forward_as_tuple(k1,k1));
    EK_real disp2d_real(std::forward_as_tuple(k1,k1));
    disp2d_complex.fill(CubicTraits<2>::template get_dispersion<typename decltype(disp2d_complex)::FunctionType> (1.0));
    disp2d_real.fill(CubicTraits<2>::template get_dispersion<typename decltype(disp2d_real)::FunctionType> (1.0));
    EK_complex disp2d_complex_backup(disp2d_complex);
    
    DEBUG(disp2d_real);

    GridObject<ComplexType, EnumerateGrid, EnumerateGrid> t_matrix_complex(std::forward_as_tuple(rgrid,rgrid));
    auto t_matrix_complex2(t_matrix_complex);
    GridObject<int, EnumerateGrid, EnumerateGrid> t_matrix_int(std::forward_as_tuple(rgrid,rgrid));
    GridObject<int, EnumerateGrid, EnumerateGrid> t_matrix_int2(std::forward_as_tuple(rgrid,rgrid));

    t_matrix_complex.getData() = run_fft(disp2d_complex.getData(), FFTW_BACKWARD);
    t_matrix_complex2.getData() = run_fft(disp2d_complex.getData(), FFTW_FORWARD)/ksize/ksize;
    disp2d_complex.getData() = run_fft(t_matrix_complex.getData(), FFTW_FORWARD);

    if (!is_equal(disp2d_complex.diff(disp2d_complex_backup),0)) return EXIT_FAILURE;
    INFO("PASSED FFT1");

    t_matrix_int.fill(typename decltype(t_matrix_int)::FunctionType([&](int a1, int a2){return std::round(std::real(t_matrix_complex(a1,a2)));}));
    t_matrix_int2.fill(typename decltype(t_matrix_int)::FunctionType([&](int a1, int a2){return std::round(std::real(t_matrix_complex2(a1,a2)));}));
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


    GridObject<ComplexType, KMesh> disp1d (std::forward_as_tuple(k1));
    disp1d.fill(CubicTraits<1>::template get_dispersion<typename decltype(disp1d)::FunctionType> (1.0));
    auto disp1d_backup(disp1d);
    disp1d.getData() = run_fft(run_fft(disp1d.getData(), FFTW_BACKWARD), FFTW_FORWARD);
    if (!is_equal(disp1d_backup.diff(disp1d),0)) return EXIT_FAILURE;
    INFO("PASSED FFT1d")

    disp2d_complex.getData() = run_fft(run_fft(disp2d_complex.getData(), FFTW_BACKWARD), FFTW_FORWARD);
    if (!is_equal(disp2d_complex.diff(disp2d_complex_backup),0)) return EXIT_FAILURE;
    INFO("PASSED FFT2d")

    GridObject<ComplexType, KMesh, KMesh, KMesh> disp3d (std::forward_as_tuple(k1,k1,k1));
    disp3d.fill(CubicTraits<3>::template get_dispersion<typename decltype(disp3d)::FunctionType> (1.0));
    auto disp3d_backup(disp3d);
    if (!is_equal(disp3d_backup.diff(disp3d),0)) return EXIT_FAILURE;
    INFO("PASSED FFT3d")

    GridObject<ComplexType, KMesh, KMesh, KMesh, KMesh> disp4d (std::forward_as_tuple(k1,k1,k1,k1));
    disp4d.fill(CubicTraits<4>::template get_dispersion<typename decltype(disp4d)::FunctionType> (1.0));
    auto disp4d_backup(disp4d);
    if (!is_equal(disp4d_backup.diff(disp4d),0)) return EXIT_FAILURE;
    INFO("PASSED FFT4d")

    return EXIT_SUCCESS;
}
