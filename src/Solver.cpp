#include "Solver.h"

namespace FK { 

FKImpuritySolver::FKImpuritySolver(RealType U, RealType mu, RealType e_d, GFType& Delta): 
    _v_mult(0.0),
    U(U), mu(mu), e_d(e_d),
    w_grid(Delta.getGrid()), 
    beta(w_grid._beta), 
    Delta(Delta), 
    gw(GFType(w_grid)), K0(GFType(w_grid)), K1(GFType(w_grid)), Sigma(GFType(w_grid))
{
};

void FKImpuritySolver::run()
{
    INFO("Running FK Solver, beta = " << beta << ", U = " << U << ", mu = " << mu << ", e_d = " << e_d);
    std::function<ComplexType(ComplexType)> K0f, K1f;
    K0f = [this](ComplexType w){return 1.0/(w+mu-Delta(w));};
    K1f = [this](ComplexType w){return 1.0/(w+mu-Delta(w)-U);};
    auto K0finv = [this](ComplexType w){return w+mu-Delta(w);};
    K0 = K0f;
    K1 = K1f;
    auto tempF = [this](ComplexType w){return K1(w)/K0(w)*(1.0+(U+e_d-2*mu)/w);};
    auto Z1toZ0 = w_grid.prod(tempF);
    w_0 = 1.0/(1.0+Z1toZ0);
    w_1 = 1.0-w_0;

    /*auto Zf1 = [this](ComplexType w, RealType mu1){return 1.0+(mu1 - Delta(w))/w;};
    std::function<ComplexType(RealType)> Zfprod = [this,Zf1](RealType mu1){return w_grid.prod(std::bind(Zf1, std::placeholders::_1, mu1));};
    ComplexType Z0 = Zfprod(mu);
    ComplexType Z1 = Zfprod(mu-U)*std::exp(beta*(mu-e_d-U/2));
    ComplexType Z=Z0+Z1;
    DEBUG("Z0 = " << Z0 << ", Z1 = " << Z1);
    w_0 = Z0/Z;
    w_1 = Z1/Z;*/
    
    _v_mult = beta*U*U*w_0*w_1;
    INFO("w_0 = " << w_0 << "; w_1 = " << w_1 );
    gw = K0*w_0 + K1*w_1;
    DEBUG("gw = " << gw);
    std::function<ComplexType(ComplexType)> Sigmaf = [this,K0finv](ComplexType w){return U*U*w_1*w_0/(K0finv(w)+w_1*U) + w_1*U;};
    Sigma = Sigmaf;
    DEBUG("Sigma = " << Sigma);
    }

} // end of namespace FK

