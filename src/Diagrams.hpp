#ifndef __FK_DFDIAGRAMS_HPP__
#define __FK_DFDIAGRAMS_HPP__

#include "Diagrams.h"
#include <Eigen/LU>
#include <Eigen/Dense>
#include "FFT.hpp"

namespace FK {

template <typename GKType>
inline GKType Diagrams::getStaticBubbles(const GKType &GF)
{
    GKType out(GF.grids());
    const auto& fgrid = std::get<0>(GF.grids());
    int knorm = GF[0].size();
    for (fmatsubara_grid::point iw1 : fgrid.points())  {
        auto g1 = run_fft(GF[iw1], FFTW_FORWARD);
        out[iw1] = run_fft(g1*g1.conj(), FFTW_BACKWARD)/knorm;
        };
    return out / (-fgrid.beta());
} 

template <typename GKType, typename ... ArgTypes>
inline typename Diagrams::GLocalType Diagrams::getBubble(const GKType &GF, const std::tuple<ArgTypes...>& args)
{
    static_assert(GKType::N==sizeof...(ArgTypes), "Argument size mismatch");
    GKType GF_shifted = GF.shift(args);
    return getBubble(GF,GF_shifted);
}

template <typename GKType, typename ... ArgTypes>
inline typename Diagrams::GLocalType Diagrams::getBubble(const GKType &GF, ArgTypes... args)
{
    return getBubble(GF,std::forward_as_tuple(args...));
}

template <typename GKType>
inline typename Diagrams::GLocalType Diagrams::getBubble(const GKType &GF, const GKType &GF_shift)
{
    const auto _fGrid = std::get<0>(GF.grids()); 
    GLocalType out(_fGrid);
    auto GF_shifted=GF*GF_shift;

    real_type knorm = 1;
    for (auto val : GF.data().shape()) knorm*=val;
    knorm/=GF.data().shape()[0];

    for (auto iw: _fGrid.points()) { 
        size_t iwn = size_t(iw);
        out[iwn] = GF_shifted[iwn].sum()/real_type(knorm); 
    }

    real_type T = 1.0/(_fGrid.beta());
    return (-T)*out;
}

template <typename ArgType>
inline typename Diagrams::GLocalType Diagrams::getBubble(const GLocalType &GF, ArgType arg)
{
    GLocalType out = GF.shift(arg)*GF;
    real_type T = 1.0/(std::get<0>(GF.grids()).beta());
    return (-T)*out;
}


template <typename MatrixType>
inline typename Diagrams::IsMatrix<MatrixType>
Diagrams::BS(const MatrixType &Chi0, const MatrixType &IrrVertex4, bool forward, bool eval_SC, size_t n_iter, real_type mix)
{
    INFO_NONEWLINE("\tRunning" << ((!forward)?" inverse ":" ") << "matrix BS equation...");
    size_t size = IrrVertex4.rows(); 

    MatrixType V4Chi;
    if (forward)
        V4Chi = MatrixType::Identity(size,size) - IrrVertex4*Chi0;
    else
        V4Chi = MatrixType::Identity(size,size) + Chi0*IrrVertex4;
    //auto ((IrrVertex4*Chi0).eigenvalues())
    auto D1 = V4Chi.determinant();
    if (std::imag(D1)>1e-7) { ERROR("Determinant : " << D1); throw (std::logic_error("Complex determinant in BS. Exiting.")); };
    if (std::real(D1)<1e-2) INFO3("Determinant : " << D1);

    if (!eval_SC && std::real(D1)>std::numeric_limits<real_type>::epsilon()) {
        try {
            if (forward) {
                V4Chi = V4Chi.colPivHouseholderQr().solve(IrrVertex4);
                //V4Chi = V4Chi.inverse()*IrrVertex4; 
                }
            else
                V4Chi=IrrVertex4*V4Chi.inverse();
            INFO("done.");
            return V4Chi;
        }
        catch (std::exception &e) {
            ERROR("Couldn't invert the vertex");
        }
    }; // From here solver by iterations
    auto V4 = IrrVertex4;
    auto V4_old = IrrVertex4;
    V4Chi=IrrVertex4*Chi0;
    INFO_NONEWLINE("\tEvaluating BS self-consistently. Making " << n_iter << " iterations.");
    real_type diffBS = 1.0;
    for (size_t n=0; n<n_iter && diffBS > 1e-8; ++n) { 
        INFO_NONEWLINE("\t\t" << n+1 << "/" << n_iter<< ". ")
        if (forward)
            V4 = (IrrVertex4 + V4Chi*V4_old)*mix + (1.0-mix)*V4_old;
        else 
            V4 = (IrrVertex4 - V4_old*V4Chi)*mix + (1.0-mix)*V4_old;
        diffBS = (V4-V4_old).norm();
        INFO("vertex diff = " << diffBS);
        V4_old = V4;
        }
    return V4;
}

template <typename G>
inline Diagrams::IsGridObject1<G> Diagrams::BS(const G &Chi0, const G &IrrVertex4,
            bool forward, bool eval_SC, size_t n_iter, real_type mix)
{
	typedef typename std::tuple_element<0, typename G::grid_tuple>::type GridType;
	typedef typename G::value_type ValueType;
    INFO_NONEWLINE("\tRunning " << ((!forward)?"inverse":"") << " BS equation...");
    const auto _fGrid = std::get<0>(IrrVertex4.grids());
    grid_object<ValueType,GridType> Vertex4_out(IrrVertex4);
    grid_object<real_type,GridType> EVCheck(_fGrid); 
    std::function<real_type(typename GridType::point)> absEVf = [&](typename GridType::point w)->real_type{return std::abs(Chi0(w)*IrrVertex4(w)); };
    EVCheck.fill(absEVf);
    real_type max_ev = *std::max_element(EVCheck.data().begin(), EVCheck.data().end());
    INFO_NONEWLINE("\tMaximum EV of Chi0*gamma = " << max_ev << " ."); // << "|" << max_ev_re << "|" << max_ev_im);
    if (std::abs(max_ev-1.0) < 1e-6 || eval_SC) {
        grid_object<ValueType,GridType> Vertex4_old(IrrVertex4);
        INFO2 ("Caught divergence, evaluating BS equation self_consistently. ");
        real_type diffBS = 1.0;
        for (size_t n=0; n<n_iter && diffBS > 1e-8; ++n) { 
            //INFO("BS iteration " << n << " for iW = " << complex_type(iW) << ", (qx,qy) = (" << real_type(qx) << "," << real_type(qy) << ").");
            Vertex4_out = IrrVertex4 + IrrVertex4*Chi0*Vertex4_old;
            diffBS = Vertex4_out.diff(Vertex4_old);
            INFO("vertex diff = " << diffBS);
            Vertex4_old = Vertex4_out*mix+(1.0-mix)*Vertex4_old;
            }
        }
    else { 
        #ifndef NDEBUG
        DEBUG("Evaluating BS equation using inversion");
        #endif
        Vertex4_out = IrrVertex4/(1.0 - (2*forward-1)*Chi0 * IrrVertex4);
        }
    INFO("done");
    return Vertex4_out;
}

template <typename G>
inline Diagrams::IsGridObject<G> Diagrams::BS(const G &Chi0, const G &IrrVertex4,
            bool forward, bool eval_SC, size_t n_iter, real_type mix)
{
    INFO_NONEWLINE("\tRunning" << ((!forward)?" inverse ":" ") << "BS equation...");
    G Vertex4_out(IrrVertex4);
    if (!eval_SC) {
        try {
            Vertex4_out = IrrVertex4/(1.0 - (2*forward-1)*Chi0 * IrrVertex4);
            INFO("done.");
            return Vertex4_out;
            }
         catch (std::exception &e) {
                ERROR("Couldn't invert the vertex");
            };
        };
    G Vertex4_old(IrrVertex4);
    real_type diffBS = 1.0;
    INFO_NONEWLINE("\tEvaluating BS self-consistently. ");
    for (size_t n=0; n<n_iter && diffBS > 1e-8; ++n) { 
        Vertex4_out = IrrVertex4 + IrrVertex4*Chi0*Vertex4_old*(2*forward-1);
        diffBS = Vertex4_out.diff(Vertex4_old);
        INFO("vertex diff = " << diffBS);
        Vertex4_old = Vertex4_out*mix+(1.0-mix)*Vertex4_old;
        };
    INFO("done");
    return Vertex4_out;
}




template <typename VertexType>
inline VertexType Diagrams::getSusc(const VertexType& Chi0, const VertexType &FullVertex4)
{
    return Chi0 + Chi0*FullVertex4*Chi0;
}

} // end of namespace FK 

#endif // endif :: #ifndef __FK_DFDIAGRAMS_HPP__
