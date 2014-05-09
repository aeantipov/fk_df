#include "Solver.h"

namespace FK { 

FKImpuritySolver::FKImpuritySolver(real_type U, real_type mu, real_type e_d, GFType& Delta):
    _v_mult(0.0),
    U(U), mu(mu), e_d(e_d),
    w_grid(Delta.grid()),
    half_grid(0,std::max(w_grid.w_max_*2,int(w_grid._beta)*10),w_grid._beta),
    beta(w_grid._beta), 
    Delta(Delta), 
    gw(GFType(w_grid)), K0(GFType(w_grid)), K1(GFType(w_grid)), Lambda(GFType(w_grid)), Sigma(GFType(w_grid))
{
};

void FKImpuritySolver::run(bool calc_weight)
{
    INFO("Running FK Solver, beta = " << beta << ", U = " << U << ", mu = " << mu << ", e_d = " << e_d << ". Using " << std::max(w_grid.w_max_,int(beta)*10) << " freqs.");
    GFType iw(w_grid);
    std::function<complex_type(complex_type)> iwtail = [](complex_type w){return w;};
    iw.fill(iwtail);
    K0 = 1.0/(iw+mu-Delta);
    K1 = 1.0/(iw+mu-Delta-U);
    K0.tail() = [=](complex_type w)->complex_type{return 1.0/(w+mu-Delta.tail_eval(w));};
    K1.tail() = [=](complex_type w)->complex_type{return 1.0/(w+mu-Delta.tail_eval(w)-U);};
   
    if (calc_weight) {
        real_type alpha = beta*(mu-e_d-U/2);
        std::function<real_type(fmatsubara_grid::point)> tempF = [this](fmatsubara_grid::point w){
            return 2.*log(std::abs(K0(w))) - 2.*log(std::abs(K1(w))); 
            };
        alpha+= half_grid.integrate(tempF)*beta;
        alpha=exp(alpha);
        w_1 = alpha / (1.0 + alpha);
        w_0 = 1.0-w_1;
        };
    
    _v_mult = beta*U*U*w_0*w_1;
    INFO("w_0 = " << w_0 << "; w_1 = " << w_1 );

    gw = K0*w_0 + K1*w_1;
    Lambda = K0*K1/gw/gw;
    Sigma = U*U*w_1*w_0/(1.0/K0-w_0*U) + w_1*U; 

    Sigma.tail() = std::bind([=](complex_type w){return w_1*U + w_0*w_1*U*U/w + w_0*w_1*U*U*(mu-w_0*U)/std::abs(w*w);}, std::placeholders::_1);
    //gw.tail() = std::bind([=](complex_type w){return (mu-w_1*U-real(Delta.tail()(w)))/std::abs(w*w)+1.0/w;}, std::placeholders::_1);
    gw.tail() = std::bind([=](complex_type w){return (mu-w_1*U)/std::abs(w*w)+1.0/w;}, std::placeholders::_1);
    Lambda.tail() = std::bind([=](complex_type w){return 1.0 - U*(w_1*w_1 - w_0*w_0)/w;}, std::placeholders::_1);
    //DEBUG("Sigma = " << Sigma);
}

grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> FKImpuritySolver::getBubble() const
{
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> out(std::forward_as_tuple(w_grid, w_grid));
    real_type T = 1.0/w_grid._beta;
    typedef decltype(out) tmp1;
    tmp1::point_function_type f = [&](fmatsubara_grid::point w1, fmatsubara_grid::point w2) {
        return -T*(gw(w1)*gw(w2));
    };
    out.fill(f);
    out.tail() = [&](complex_type w1, complex_type w2){return -T*gw.tail_eval(w1)*gw.tail_eval(w2);};
    return out;
}

typename FKImpuritySolver::GFType FKImpuritySolver::getLambda() const
{
    return Lambda;
}

} // end of namespace FK

