#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "container.h"
#include "grid_object.h"
#include "Solver.h"
#include "DMFT.h"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;
typedef typename Solver::GFType GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, real_type tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}


int main()
{
    Log.setDebugging(true);

    INFO("Hi! Doing Falicov-Kimball. ");
    real_type U = 2.0;
    real_type mu = U/2;
    real_type e_d = 0.0;
    real_type T = 0.1325348;
    real_type beta = 1./T;
    real_type t = 1.0; 
    size_t maxit = 1000;
    real_type mix = 0.5;
    constexpr size_t KPOINTS=32;
    constexpr size_t D=2;
    
    size_t n_freq = 1024;
    fmatsubara_grid gridF(-n_freq, n_freq, beta);

    GF iwn(gridF);
    iwn.fill([](complex_type w){return w;});
    GF Delta(gridF);
    std::function<complex_type(complex_type)> f1;
    f1 = [t](complex_type w) -> complex_type {return t*2.0*t/w;};
    Delta.fill(f1);
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    real_type diff=1.0;
    auto kGridRegular = kmesh(KPOINTS);
    CubicDMFTSC<D> SC(Solver, t, kGridRegular);

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

    size_t nkpoints = 128;
    size_t npaths = 2;
    std::vector<std::vector<std::array<real_type, D>>> paths(npaths); // Here generate the path in BZ
    std::vector<std::string> path_names = {"Path1", "Path2"}; // name of paths
    
    std::function<real_type(size_t)> zero_pi_fill = [nkpoints](int x)->real_type{real_type y = 6*(real_type(x)/(nkpoints-1)); return PI*(1.0 - (pow(0.1,y) - pow(0.1,6)));}; 
    auto kloggrid = real_grid(0,nkpoints,zero_pi_fill);

    for (size_t i=0; i<nkpoints; ++i) { 
        std::array<real_type, D> path1_pts; path1_pts.fill(kloggrid[i]);
        std::array<real_type, D> path2_pts; path2_pts.fill(kloggrid[i]); path1_pts[0]=PI;
        paths[0].push_back(path1_pts);
        paths[1].push_back(path2_pts);
        };

    
    auto Lambda = 1.0 + gw*(2.0*Sigma - U); 
    auto Lambda2 = (1. + gw*Sigma)*(1.+gw*(Sigma-U));

    auto glatglat = SC.getGLat(gridF);
    auto bubble2 = static_cast<DMFT&>(SC).getBubblePI(0);

    for (size_t p=0; p<npaths; ++p) {
    INFO("Path " << p);
        grid_object<complex_type, real_grid> susc_cc_vals(kloggrid), susc_cf_vals(kloggrid), susc_ff_vals(kloggrid);
        for (size_t i=0; i<nkpoints; ++i) { 

            INFO_NONEWLINE("\t" << i << "/" << nkpoints << " : k= ");
            std::array<real_type, D> path_pts = paths[p][i];
            __tuple_print<std::tuple<real_type,real_type>>::print(path_pts);
            
            typedef decltype(SC)::GKType::arg_tuple argwktuple;
            auto f1 = [&](argwktuple in)->complex_type
                {argwktuple shift = std::tuple_cat(std::make_tuple(0.0), path_pts);
                 argwktuple out = glatglat._shiftArgs(in,shift);
                  return -T*SC.glat_analytic(in)*SC.glat_analytic(out); 
                };
            decltype(glatglat)::function_type f2 = tools::fun_traits<decltype(glatglat)::function_type>::getFromTupleF(f1);
            glatglat.fill(f2);
            GF bubble(gridF);
            typename GF::point_function_type f3 = [&](fmatsubara_grid::point w)->complex_type{return glatglat[size_t(w)].sum()/pow(KPOINTS,D);};
            bubble.fill(f3);
            
            
            INFO("Diff with q=pi bare cc susc" << bubble.diff(bubble2));

            // Skeleton expansion.
            auto nu = gw*(-1.0/gw/gw - T/bubble);
            auto d1 = Lambda*gw*nu + Lambda2;
            auto ugamma_down = 1.0 - (w_0*w_1*U*U*gw*gw*gw*nu/(Lambda2 * (d1))).sum();
            auto ugamma = (w_0*w_1*U*U*gw*gw/(d1)).sum() / (ugamma_down); 

            auto chi_cf = ugamma/U;
            auto chi_ff = w_0*w_1/T/ugamma_down;
            auto chi_cc = -T*((Lambda - ugamma)*gw*gw/d1).sum();
            //susc_vals2[names[i]] = chi_cc;
            
            
            INFO2("Static cc susc (bs) = " << chi_cc);
            INFO2("Static cf susc (bs) = " << chi_cf);
            INFO2("Static ff susc (bs) = " << chi_ff);
            susc_cc_vals.get(path_pts[1]) = chi_cc;
            susc_cf_vals.get(path_pts[1]) = chi_cf;
            susc_ff_vals.get(path_pts[1]) = chi_ff;
            };
        susc_cc_vals.savetxt("ChargeCC"+path_names[p]+".dat");
        susc_cf_vals.savetxt("ChargeCF"+path_names[p]+".dat");
        susc_ff_vals.savetxt("ChargeFF"+path_names[p]+".dat");
        };

    /*num_io<complex_type>(bare_susc_vals["zero"]).savetxt("StaticChi0q0.dat");
    num_io<complex_type>(bare_susc_vals["pi"]).savetxt("StaticChi0qPI.dat");
    num_io<complex_type>(bare_susc_vals["local"]).savetxt("StaticChi0Local.dat");
    num_io<complex_type>(susc_vals["zero"]).savetxt("StaticChiq0.dat");
    num_io<complex_type>(susc_vals["pi"]).savetxt("StaticChiqPI.dat");
    num_io<complex_type>(susc_vals["local"]).savetxt("StaticChiLocal.dat");
*/
    }
