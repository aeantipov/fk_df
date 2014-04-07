#include "Solver.h"

namespace FK { 

FKImpuritySolver::FKImpuritySolver(RealType U, RealType mu, RealType e_d, GFType& Delta): 
    _v_mult(0.0),
    U(U), mu(mu), e_d(e_d),
    w_grid(Delta.getGrid()), 
    half_grid(0,std::max(w_grid._w_max*2,int(beta)*10),w_grid._beta),
    beta(w_grid._beta), 
    Delta(Delta), 
    gw(GFType(w_grid)), K0(GFType(w_grid)), K1(GFType(w_grid)), Lambda(GFType(w_grid)), Sigma(GFType(w_grid))
{
};

void FKImpuritySolver::run(bool calc_weight)
{
    INFO("Running FK Solver, beta = " << beta << ", U = " << U << ", mu = " << mu << ", e_d = " << e_d << ". Using " << std::max(w_grid._w_max,int(beta)*10) << " freqs.");
    GFType iw(w_grid);
    std::function<ComplexType(ComplexType)> iw_f = [](ComplexType w){return w;};
    iw.fill(iw_f);
    K0 = 1.0/(iw+mu);
    K1 = 1.0/(iw+mu-U);
    K0._f = [=](ComplexType w)->ComplexType{return 1.0/(w+mu);};
    K1._f = [=](ComplexType w)->ComplexType{return 1.0/(w+mu-U);};
   
    if (calc_weight) {
        RealType alpha = beta*(mu-e_d-U/2);
        std::function<RealType(FMatsubaraGrid::point)> tempF = [this](FMatsubaraGrid::point w){
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

    Sigma._f = std::bind([=](ComplexType w){return w_1*U + w_0*w_1*U*U/w + w_0*w_1*U*U*(mu-w_0*U)/std::abs(w*w);}, std::placeholders::_1);
    //gw._f = std::bind([=](ComplexType w){return (mu-w_1*U-real(Delta._f(w)))/std::abs(w*w)+1.0/w;}, std::placeholders::_1);
    gw._f = std::bind([=](ComplexType w){return (mu-w_1*U)/std::abs(w*w)+1.0/w;}, std::placeholders::_1);
    Lambda._f = std::bind([=](ComplexType w){return 1.0 - U*(w_1*w_1 - w_0*w_0)/w;}, std::placeholders::_1);
    //DEBUG("Sigma = " << Sigma);
}

GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> FKImpuritySolver::getBubble() const
{
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> out(std::forward_as_tuple(w_grid, w_grid));
    RealType T = 1.0/w_grid._beta;
    typedef decltype(out) tmp1;
    tmp1::PointFunctionType f = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2) { 
        return -T*(gw(w1)*gw(w2));
    };
    out.fill(f);
    out._f = [&](ComplexType w1, ComplexType w2){return -T*gw._f(w1)*gw._f(w2);};
    return out;
}

typename FKImpuritySolver::GFType FKImpuritySolver::getLambda() const
{
    return Lambda;
}

} // end of namespace FK

