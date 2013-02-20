#include <numeric>

#include "Solver.h"
#include "SelfConsistency.h"

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
    //Log.setDebugging(true);

    INFO("Hi! Doing Falicov-Kimball. ");
    RealType U = 2.0;
    RealType mu = U/2;
    RealType e_d = 0.0;
    RealType beta = 1./0.14752; //0.07759;
    RealType T = 1.0/beta;
    RealType t = 1.0; 
    size_t maxit = 1000;
    RealType mix = 0.5;
    
    size_t n_freq = 256;
    size_t n_b_freq = 15;
    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    BMatsubaraGrid gridB(-n_b_freq, n_b_freq+1, beta);
    GF iwn(gridF);
    iwn.fill([](ComplexType w){return w;});
    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*2.0*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    RealType diff=1.0;
    CubicInfDMFTSC SC(Solver, t, RealGrid(-6.0*t,6.0*t,1024));

    for (int i=0; i<maxit && diff>1e-8; ++i) {
        INFO("Iteration " << i);
        if (i<3) Solver.run(true);
        else Solver.run(false);
        auto Delta_new = SC();
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);
        Delta = Delta_new*mix + (1.0-mix)*Solver.Delta;
        Solver.Delta = Delta; 
        }

    bool success = false;

    auto gw = Solver.gw;
    auto Sigma = Solver.Sigma;
    auto K0 = Solver.K0;
    auto K1 = Solver.K1; 
    auto w_0 = Solver.w_0;
    auto w_1 = Solver.w_1;

    INFO("Some relations check");

    auto test00 = K0 - K1;
    auto test01 = (-1.0)*K0*K1*U;

    //DEBUG(test00.diff(test01));
    if (!is_equal(test00.diff(test01),0)) return EXIT_FAILURE;

    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> Vertex4_out(std::forward_as_tuple(gridF,gridF)); 
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> Vertex4_2(std::forward_as_tuple(gridF,gridF)); 
    decltype(Vertex4_out)::PointFunctionType VertexF = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return Solver.getFVertex4(w1,w2);};
    decltype(Vertex4_out)::PointFunctionType VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return Solver.getVertex4(0.0, w1,w2);};
    Vertex4_out.fill(VertexF);
    Vertex4_2.fill(VertexF2);
    typename GF::PointFunctionType VertexDiagF = [&](FMatsubaraGrid::point w1){return Solver.getFVertex4(w1,w1);};
    GFWrap Vertex4_diag(gridF);
    Vertex4_diag.fill(VertexDiagF);
    

    auto test02 = (1.0+gw*(2*Sigma-U))*gw*gw/(1+gw*Sigma)/(1.+gw*(Sigma-U));
    auto test03 = gw*gw + T*gw*gw*Vertex4_diag*gw*gw;
    //DEBUG(test02.diff(test03));
    if (!is_equal(test02.diff(test03),0)) return EXIT_FAILURE;
    auto local_bubble = -T*gw*gw;
    ComplexType s1=0.0;
    ComplexType s2=0.0;
    ComplexType s3=0.0;

    for (auto w1: gridF.getPoints()) { 
        s1+=local_bubble(w1);
        s3+=local_bubble(w1);
        s1-=local_bubble(w1)*Vertex4_out(w1,w1)*local_bubble(w1); 
        s2-=T*(w_0*K0(w1)*K0(w1)+w_1*K1(w1)*K1(w1));
        for (auto w2: gridF.getPoints()) {
            s1+=local_bubble(w1)*Vertex4_out(w1,w2)*local_bubble(w2); 
            s2+=T*(w_0*K0(w1)*K0(w2)+w_1*K1(w1)*K1(w2));
            s3+=local_bubble(w1)*Vertex4_2(w1,w2)*local_bubble(w2);
            }
        };

    INFO("Local susceptibility (from vertex): " << s1);
    INFO("Local susceptibility (exact): " << s2);
    INFO("Local susceptibility (from vertex 2 ): " << s3);
    if (!(is_equal(s1,s2)) || !(is_equal(s1,s3))) return EXIT_FAILURE;
    

    INFO("T = " << T << " | U =" << U);
    INFO("Values from skeleton expansion");
    auto Bubbleq0 = SC.getBubble0(0.0);
    auto BubbleqPI = SC.getBubblePI(0.0); 
    auto Lambda = 1.0 + gw*(2.0*Sigma - U); 
    auto Lambda2 = (1. + gw*Sigma)*(1.+gw*(Sigma-U));

    if(!is_equal((gw*gw/Lambda2).diff(K0*K1),0)) return EXIT_FAILURE;

    auto St = Lambda/Lambda2*gw*gw;
    //DEBUG(St.diff(test03));
    if (!is_equal(St.diff(test03),0)) return EXIT_FAILURE;
    std::vector<std::string> names = {"local", "pi", "zero"};
    std::vector<GF> bubbles = { -T*gw*gw, BubbleqPI, Bubbleq0 };
    std::vector<GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid>> long_bubbles = {Solver.getBubble(), SC.getBubblePI(), SC.getBubble0()}; 
    std::map<std::string,ComplexType> susc_vals;
     
    for (size_t i=0; i<bubbles.size(); ++i) { 
        auto bubble = bubbles[i];
        auto dual_bubble = bubbles[i]+T*gw*gw;
        auto long_bubble = long_bubbles[i];
        auto long_dual_bubble = long_bubble - Solver.getBubble();
        auto nu = gw*(-1.0/gw/gw - T/bubble);
        auto d1 = Lambda*gw*nu + Lambda2;
        auto d1_2 = (-1.0)*Lambda2 * T * gw*gw/bubble * (1.0 + Vertex4_diag * dual_bubble);
    //    DEBUG(d1);
    //    DEBUG(d1_2);
        if (!is_equal(d1.diff(d1_2),0)) return EXIT_FAILURE;
        auto ugamma_down = 1.0 - (w_0*w_1*U*U*gw*gw*gw*nu/(Lambda2 * (d1))).sum();
        auto sigma_l = ((Vertex4_diag * dual_bubble)/(1+Vertex4_diag * dual_bubble)).sum();
        auto ugamma_down2 = 1.0 - sigma_l ; 
        if (!is_equal(ugamma_down,ugamma_down2)) return EXIT_FAILURE;
        auto ugamma = (w_0*w_1*U*U*gw*gw/(d1)).sum() / (ugamma_down); 
        //DEBUG("U*gamma = " << ugamma); 
        auto chi_v = -T*((Lambda - ugamma)*gw*gw/d1).sum();
        auto chi_part1 = 1.0 * ( (1+T*Vertex4_diag*gw*gw)*bubble / (1+Vertex4_diag*dual_bubble)).sum();
        auto chi_part2 = T*(ugamma*gw*gw/d1).sum(); 
        ComplexType chi_part2_2 = 0.0;
        ComplexType sigma_l_2 = 0.0;
        for (auto wn : gridF.getPoints()) {
            for (auto wm : gridF.getPoints()) {
                chi_part2_2+=Vertex4_out(wn,wm)*bubble(wn)*bubble(wm)/(1.0 + Vertex4_out(wm,wm)*dual_bubble(wm))/(1.0 + Vertex4_out(wn,wn)*dual_bubble(wn));
                sigma_l_2+=(Vertex4_out(wn,wm)*dual_bubble(wm))/(1.0 + Vertex4_out(wn,wm)*dual_bubble(wm));
                }
            }; 


        if (!is_equal(chi_part2_2/(1.0-sigma_l),chi_part2)) return EXIT_FAILURE;;

        susc_vals[names[i]] = chi_v;
        INFO2(names[i] << " susceptibility = " << chi_v);
        DEBUG(chi_part1 << "+" <<  chi_part2 << " = " << chi_part1 + chi_part2);
    };

    if (!is_equal(susc_vals["local"],s1)) return EXIT_FAILURE; 
    
    INFO("BS static susc");
    std::map<std::string,ComplexType> susc_vals2;

    auto V4 = Vertex4_2.getData().getAsMatrix();

    //auto FullVertexqPI = StaticV4.getData().getAsMatrix();
    //auto FullVertexq0 = StaticV4.getData().getAsMatrix();

    

    for (size_t i=0; i<bubbles.size(); ++i) { 
        auto long_bubble = long_bubbles[i];
        auto dual_long_bubble = long_bubble - Solver.getBubble();
        auto FullVertex2 = Vertex4_out/(1.0 - Vertex4_out*dual_long_bubble);
        DEBUG(FullVertex2.diff(Vertex4_out));

        auto bubble = bubbles[i];
        auto dual_bubble = bubbles[i]+T*gw*gw;
        auto dual_bubble_matrix = dual_bubble.getData().getAsDiagonalMatrix();
        
        //auto FullVertex = V4 / (1.0 - V4*dual_bubble_matrix);
        auto FullVertex = Diagrams::BS(dual_bubble_matrix, V4, true);

        ComplexType susc = 0.0;
        ComplexType susc_cor = 0.0;
        //DEBUG(long_bubble.sum());
        //DEBUG(bubble.sum());

        for (auto w1: gridF.getPoints()) { 
            susc+=bubble(w1);
            auto v1 = bubble(w1)*FullVertex(size_t(w1),size_t(w1))*bubble(w1);
            auto v2 = long_bubble(w1,w1)*FullVertex2(w1,w1)*long_bubble(w1,w1);
            //DEBUG(v1 << " <---> " << v2);
            //susc-=v1+v2; 
            //susc_cor-=v2;
            for (auto w2: gridF.getPoints()) {
                //susc_cor+=long_bubble(w1,w2);
                susc+=bubble(w1)*FullVertex(size_t(w1),size_t(w2))*bubble(w2); 
                //susc_cor+=long_bubble(w1,w2)*FullVertex2(w1,w2)*long_bubble(w1,w2);
                }
            };
        susc_vals[names[i]] = susc;
        INFO2("Static susc " << names[i] <<" (bs) = " << susc);
        //DEBUG(susc_cor);
        }
    INFO("T = " << T << " | U =" << U);
    INFO("Static bare q=0 susc = " << Bubbleq0.sum());
    INFO("Static bare q=PI susc = " << BubbleqPI.sum());
    INFO("Static q=0 susc = " << susc_vals["zero"]);
    INFO("Static q=PI susc = " << susc_vals["pi"]);
    success = true;
    return EXIT_SUCCESS;
}
