#include "Solver.h"

namespace FK { 

FKImpuritySolver::FKImpuritySolver(RealType U, RealType mu, RealType e_d, GFType& Delta): 
    _v_mult(0.0),
    U(U), mu(mu), e_d(e_d),
    w_grid(Delta.getGrid()), 
    half_grid(0,w_grid._w_max,w_grid._beta),
    beta(w_grid._beta), 
    Delta(Delta), 
    gw(GFType(w_grid)), K0(GFType(w_grid)), K1(GFType(w_grid)), Sigma(GFType(w_grid))
{
};

void FKImpuritySolver::run(bool calc_weight)
{
    INFO("Running FK Solver, beta = " << beta << ", U = " << U << ", mu = " << mu << ", e_d = " << e_d);
    GFType iw(w_grid);
    std::function<ComplexType(ComplexType)> iw_f = [](ComplexType w){return w;};
    iw.fill(iw_f);
    K0 = 1.0/(iw+mu-Delta);
    K1 = 1.0/(iw+mu-Delta-U);
    //std::function<ComplexType(FMatsubaraGrid::point)> K0f, K1f;
    //K0f = [this](FMatsubaraGrid::point w){return 1.0/(ComplexType(w)+mu-Delta(w));};
    //K1f = [this](FMatsubaraGrid::point w){return 1.0/(ComplexType(w)+mu-Delta(w)-U);};
    //auto tailf = [this](FMatsubaraGrid::point w)->ComplexType{return 1.0/(ComplexType(w));};
    //K0.fill(K0f);
    //K1 = K1f;
    //K0._f = tailf;
    //K1._f = tailf; 
    if (calc_weight) {
        //auto tempF = [this](FMatsubaraGrid::point w){return K1(w)/K0(w)*(1.0+(U+e_d-2*mu)/ComplexType(w));};
        std::function<RealType(FMatsubaraGrid::point)> tempF = [this](ComplexType w){return std::pow(std::abs(K1(w))/std::abs(K0(w)),2)*(1.0+std::pow(U+e_d-2*mu,2)/std::pow(imag(ComplexType(w)),2));};
        //auto tempF = [this](FMatsubaraGrid::point w){return K1(w)/K0(w)*(1.0+(U+e_d-2*mu)/ComplexType(w));};
        RealType Z1toZ0 = half_grid.prod(tempF);
        w_0 = 1.0/(1.0+Z1toZ0);
        w_1 = 1.0-w_0;
        };

    
    /*auto Zf1 = [this](FMatsubaraGrid::point w, RealType mu1){return 1.0+(mu1 - Delta(w))/ComplexType(w);};
    std::function<ComplexType(RealType)> Zfprod = [this,Zf1](RealType mu1){return w_grid.prod(std::bind(Zf1, std::placeholders::_1, mu1));};
    ComplexType Z0 = Zfprod(mu);
    ComplexType Z1 = Zfprod(mu-U)*std::exp(beta*(mu-e_d-U/2));
    ComplexType Z=Z0+Z1;
    DEBUG("Z0 = " << Z0 << ", Z1 = " << Z1);
    w_0 = Z0/Z;
    w_1 = Z1/Z;
    */
    
    _v_mult = beta*U*U*w_0*w_1;
    INFO("w_0 = " << w_0 << "; w_1 = " << w_1 );
    gw = K0*w_0 + K1*w_1;
    Sigma = U*U*w_1*w_0/(1.0/K0-w_0*U) + w_1*U; 
    Sigma._f = std::bind([this](ComplexType w){return w_1*U + w_0*w_1*U*U/w;}, std::placeholders::_1);
    //gw._f = std::bind([=](ComplexType w){return (mu-w_1*U-real(Delta._f(w)))/std::abs(w*w)+1.0/w;}, std::placeholders::_1);
    gw._f = std::bind([this](ComplexType w){return (mu-w_1*U)/std::abs(w*w)+1.0/w;}, std::placeholders::_1);
    //DEBUG("Sigma = " << Sigma);
}

template <> 
ComplexType FKImpuritySolver::getVertex4<FMatsubaraGrid::point, FMatsubaraGrid::point> (FMatsubaraGrid::point w1, FMatsubaraGrid::point w2) const 
{
    return _v_mult*K0(w1)*K0(w2)*K1(w1)*K1(w2);
}
 
template <> 
ComplexType FKImpuritySolver::getVertex4<BMatsubaraGrid::point, FMatsubaraGrid::point> (BMatsubaraGrid::point wB, FMatsubaraGrid::point wF) const 
{
    int bindex=BMatsubaraGrid(0,0,w_grid._beta).getNumber(wB);
    FMatsubaraGrid::point w2;
    w2._val=wF._val+wB._val;
    if (-bindex<int(wF)) {
        w2._index=wF._index+bindex;
        return -this->getVertex4<FMatsubaraGrid::point, FMatsubaraGrid::point>(wF,w2);
        }
    else return -_v_mult*K0(wF)*K0(ComplexType(wF)+ComplexType(wB))*K1(wF)*K1(ComplexType(wF)+ComplexType(wB));
}


} // end of namespace FK

