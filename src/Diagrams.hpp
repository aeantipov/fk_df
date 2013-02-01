#ifndef __FK_DFDIAGRAMS_HPP__
#define __FK_DFDIAGRAMS_HPP__

#include "Diagrams.h"
#include <Eigen/LU>
#include <Eigen/Dense>

namespace FK {

template <typename GKType, typename ... ArgTypes>
inline typename Diagrams::GLocalType Diagrams::getBubble(const GKType &GF, ArgTypes... args)
{
    return getBubble(GF,std::forward_as_tuple(args...));
}

template <typename GKType, typename ... ArgTypes>
inline typename Diagrams::GLocalType Diagrams::getBubble(const GKType &GF, const std::tuple<ArgTypes...>& args)
{
    static_assert(GKType::N==sizeof...(ArgTypes), "Argument size mismatch");
    const auto _fGrid = std::get<0>(GF.getGrids()); 
    GLocalType out(_fGrid);
    static GKType GF_shifted(GF.getGrids());
    GF_shifted = GF.shift(args);
    GF_shifted*=GF;

    RealType knorm = 1;
    for (auto val : GF._dims) knorm*=val; 
    knorm/=GF._dims[0];

    for (auto iw: _fGrid.getPoints()) { 
        size_t iwn = size_t(iw);
        out[iwn] = GF_shifted[iwn].sum()/RealType(knorm); 
    }

    RealType T = 1.0/(_fGrid._beta);
    return (-T)*out;
}

template <typename ArgType>
inline typename Diagrams::GLocalType Diagrams::getBubble(const GLocalType &GF, ArgType arg)
{
    GLocalType out = GF.shift(arg)*GF;
    RealType T = 1.0/(std::get<0>(GF.getGrids())._beta);
    return (-T)*out;
}

/*
template <typename VertexType>
inline VertexType Diagrams::BS(const VertexType &Chi0, const VertexType &IrrVertex4, bool eval_SC, size_t n_iter, RealType mix)
{
    VertexType one(IrrVertex4), Vertex4(IrrVertex4);
    one = 1.0;
    try {
        return IrrVertex4/(one - Chi0 * IrrVertex4); 
    }
    catch (std::exception &e) {
        ERROR("Couldn't invert the vertex");
        exit(1);
    }
}*/

template <typename ValueType>
inline MatrixType<ValueType> Diagrams::BS(const MatrixType<ValueType> &Chi0, const MatrixType<ValueType> &IrrVertex4, bool forward, bool eval_SC, size_t n_iter, RealType mix)
{
    INFO_NONEWLINE("\tRunning" << ((!forward)?" inverse ":" ") << "BS equation...");
    size_t size = IrrVertex4.rows(); 

    MatrixType<ValueType> V4Chi;
    if (forward)
        V4Chi = MatrixType<ValueType>::Identity(size,size) - IrrVertex4*Chi0;
    else
        V4Chi = MatrixType<ValueType>::Identity(size,size) + Chi0*IrrVertex4;
    auto D1 = V4Chi.determinant();

    if (!eval_SC && std::abs(D1)>std::numeric_limits<RealType>::epsilon()) {
        try {
            #ifdef NDEBUG
            #undef NDEBUG
            #define NDEBUG1
            #endif
            //Eigen::ColPivHouseholderQR<MatrixType<ValueType>> Solver(MatrixType<ValueType>::Identity(size,size) - (2*forward-1)*IrrVertex4*Chi0);
            //Eigen::ColPivHouseholderQR<MatrixType<ValueType>> Solver(MatrixType<ValueType>::Identity(size,size) - (IrrVertex4*Chi0));
            //auto V4 = Solver.solve(IrrVertex4); 
            if (forward) {
                V4Chi = V4Chi.colPivHouseholderQr().solve(IrrVertex4);
                }
            else
                V4Chi=IrrVertex4*V4Chi.inverse();
            INFO("done.");
            #ifdef NDEBUG1
            #define NDEBUG
            #undef NDEBUG1
            #endif
            return V4Chi;
        }
        catch (std::exception &e) {
            ERROR("Couldn't invert the vertex");
        }
    }; // From here solver by iterations
    V4Chi=IrrVertex4*Chi0;
    auto V4 = IrrVertex4;
    auto V4_old = IrrVertex4;
    INFO_NONEWLINE("\tEvaluating BS self-consistently. ");
    RealType diffBS = 1.0;
    for (size_t n=0; n<n_iter && diffBS > 1e-8; ++n) { 
        INFO_NONEWLINE("\t\t" << n+1 << "/" << n_iter<< ". ")
        if (forward)
            V4 = IrrVertex4 + V4Chi*V4_old;
        else 
            V4 = IrrVertex4 - V4_old*V4Chi;
        diffBS = (V4-V4_old).norm();
        INFO("vertex diff = " << diffBS);
        V4_old = V4*mix+(1.0-mix)*V4_old;
        }
    return V4;
}


template <typename ValueType, typename GridType>
inline GridObject<ValueType,GridType> Diagrams::BS (const GridObject<ValueType,GridType>& Chi0, const GridObject<ValueType,GridType> &IrrVertex4, bool forward, bool eval_SC, size_t n_iter, RealType mix)
{
    //INFO_NONEWLINE("\tRunning " << ((!forward)?"inverse":"") << " BS equation...");
    const auto _fGrid = std::get<0>(IrrVertex4.getGrids());
    GridObject<ValueType,GridType> Vertex4_out(IrrVertex4);
    GridObject<RealType,GridType> EVCheck(_fGrid); 
    std::function<RealType(typename GridType::point)> absEVf = [&](typename GridType::point w)->RealType{return std::abs(Chi0(w)*IrrVertex4(w)); };
    EVCheck.fill(absEVf);
    RealType max_ev = *std::max_element(EVCheck.getData().begin(), EVCheck.getData().end());
    INFO_NONEWLINE("\tMaximum EV of Chi0*gamma = " << max_ev << " ."); // << "|" << max_ev_re << "|" << max_ev_im);
    if (std::abs(max_ev-1.0) < 1e-6 || eval_SC) {
        GridObject<ValueType,GridType> Vertex4_old(IrrVertex4);
        INFO2 ("Caught divergence, evaluating BS equation self_consistently. ");
        RealType diffBS = 1.0;
        for (size_t n=0; n<n_iter && diffBS > 1e-8; ++n) { 
            //INFO("BS iteration " << n << " for iW = " << ComplexType(iW) << ", (qx,qy) = (" << RealType(qx) << "," << RealType(qy) << ").");
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
    //INFO("done");
    return Vertex4_out;
}

template <typename VertexType>
inline VertexType Diagrams::getSusc(const VertexType& Chi0, const VertexType &FullVertex4)
{
    return Chi0 + Chi0*FullVertex4*Chi0;
}

} // end of namespace FK 

#endif // endif :: #ifndef __FK_DFDIAGRAMS_HPP__
