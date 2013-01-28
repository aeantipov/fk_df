#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"
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
    Log.setDebugging(true);

    INFO("Hi! Doing Falicov-Kimball. ");
    RealType U = 2.0;
    RealType mu = U/2;
    RealType e_d = 0.0;
    RealType beta = 1./0.1262009;
    RealType T = 1.0/beta;
    RealType t = 1.0; 
    size_t maxit = 1000;
    RealType mix = 0.5;
    constexpr size_t KPOINTS=16;
    constexpr size_t D=2;
    
    size_t n_freq = 512;
    FMatsubaraGrid gridF(-n_freq, n_freq, beta);

    GF iwn(gridF);
    iwn.fill([](ComplexType w){return w;});
    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1;
    f1 = [t](ComplexType w) -> ComplexType {return t*2.0*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    RealType diff=1.0;
    CubicDMFTSC<FKImpuritySolver, D> SC(Solver, t, KMesh(KPOINTS));

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

    auto gw = Solver.gw;
    auto Sigma = Solver.Sigma;
    auto K0 = Solver.K0;
    auto K1 = Solver.K1; 
    auto w_0 = Solver.w_0;
    auto w_1 = Solver.w_1;

    INFO("Calculating additional statistics.");
    INFO("Static susceptibility");

    size_t nkpoints = 1;
    size_t npaths = 2;
    std::vector<std::vector<std::array<RealType, D>>> paths(npaths); // Here generate the path in BZ
    
    std::function<RealType(size_t)> zero_pi_fill = [nkpoints](int x)->RealType{RealType y = 6*(RealType(x+1)/(nkpoints)); return PI*(1.0 - (pow(0.1,y) - pow(0.1,6)));}; 
    auto kloggrid = RealGrid(0,nkpoints,zero_pi_fill);

    for (size_t i=0; i<nkpoints; ++i) { 
        std::array<RealType, D> path1_pts; path1_pts.fill(kloggrid[i]);
        std::array<RealType, D> path2_pts; path2_pts.fill(kloggrid[i]); path1_pts[0]=PI;
        paths[0].push_back(path1_pts);
        paths[1].push_back(path2_pts);
        };

    
    auto Lambda = 1.0 + gw*(2.0*Sigma - U); 
    auto Lambda2 = (1. + gw*Sigma)*(1.+gw*(Sigma-U));

    for (size_t p=0; p<npaths; ++p) {
        Container<1, ComplexType> susc_cc_vals(nkpoints), susc_cf_vals(nkpoints), susc_ff_vals(nkpoints);
        for (size_t i=0; i<nkpoints; ++i) { 

            INFO(i << "/" << nkpoints);
            std::array<RealType, D> path_pts = paths[p][i];
            
            auto glat = SC.glat(path_pts[0], path_pts[1]);
            //auto bubble = Diagrams::getBubble(glat,0.0);
            auto bubble = SC.getBubblePI(0.0); 
            auto bubble2 = static_cast<SelfConsistency<FKImpuritySolver>&>(SC).getBubblePI(0);
            DEBUG(bubble.diff(bubble2));
            exit(0);

            //DEBUG(glat(_[0]));
            //DEBUG(SC.getGLat()(-);

            //DEBUG(bubble);
            //DEBUG(BubbleqPI);
            //exit(0);

            // Skeleton expansion.
            auto nu = gw*(-1.0/gw/gw - T/bubble);
            auto d1 = Lambda*gw*nu + Lambda2;
            auto ugamma_down = 1.0 - (w_0*w_1*U*U*gw*gw*gw*nu/(Lambda2 * (d1))).sum();
            auto ugamma = (w_0*w_1*U*U*gw*gw/(d1)).sum() / (ugamma_down); 

            auto chi_cf = ugamma/U;
            auto chi_ff = w_0*w_1/T/ugamma_down;
            auto chi_cc = -T*((Lambda - ugamma)*gw*gw/d1).sum();
            //susc_vals2[names[i]] = chi_cc;
            
            
            INFO2("Static cc susc " << path_pts[0] << "," << path_pts[1] <<" (bs) = " << chi_cc);
            INFO2("Static cf susc " << path_pts[0] << "," << path_pts[1] <<" (bs) = " << chi_cf);
            INFO2("Static ff susc " << path_pts[0] << "," << path_pts[1] <<" (bs) = " << chi_ff);
            };
        };

    /*__num_format<ComplexType>(bare_susc_vals["zero"]).savetxt("StaticChi0q0.dat");
    __num_format<ComplexType>(bare_susc_vals["pi"]).savetxt("StaticChi0qPI.dat");
    __num_format<ComplexType>(bare_susc_vals["local"]).savetxt("StaticChi0Local.dat");
    __num_format<ComplexType>(susc_vals["zero"]).savetxt("StaticChiq0.dat");
    __num_format<ComplexType>(susc_vals["pi"]).savetxt("StaticChiqPI.dat");
    __num_format<ComplexType>(susc_vals["local"]).savetxt("StaticChiLocal.dat");
*/
    }
