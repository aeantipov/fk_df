#include <numeric>

#include <MatsubaraGrid.hpp>
#include <KMesh.hpp>
#include <GridObject.hpp>
#include "Solver.h"
#include "DMFT.h"
#include "FFT.hpp"

#ifdef LATTICE_bethe
    typedef FK::BetheSC sc_type;
    static constexpr size_t D=0;
#elif LATTICE_cubicinfd
    typedef FK::CubicInfDMFTSC sc_type;
    static constexpr size_t D=0;
#elif LATTICE_cubic1d
    typedef FK::CubicDMFTSC<1> sc_type;
    static constexpr size_t D=1;
    #define _calc_extra_stats
#elif LATTICE_cubic2d
    typedef FK::CubicDMFTSC<2> sc_type;
    static constexpr size_t D=2;
    #define _calc_extra_stats
#elif LATTICE_cubic3d
    typedef FK::CubicDMFTSC<3> sc_type;
    static constexpr size_t D=3;
    #define _calc_extra_stats
#elif LATTICE_cubic4d
    typedef FK::CubicDMFTSC<4> sc_type;
    static constexpr size_t D=4;
    #define _calc_extra_stats
#elif LATTICE_triangular
    typedef FK::TriangularDMFT sc_type;
    static constexpr size_t D=2;
    #define _calc_extra_stats
#endif

#include "FKOptionParserDMFT.h"

#include <iostream>
#include <ctime>
#include <array>
#include <unordered_map>
#include <csignal>
#include <bitset>

#include <fftw3.h>

using namespace GFTools;
using namespace FK;

int __get_n_lines(const std::string& fname);
int __get_min_number(const std::string& fname, RealType beta);

RealType beta;
size_t extraops;
bool INTERRUPT = false;
 
typedef GridObject<ComplexType,FMatsubaraGrid> GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}

void sighandler(int signal)
{
    static size_t count = 0;
    count++;
    INFO("Caught INTERRUPT, signal " << signal <<" " << count << " times. ")
    INTERRUPT = true;
    if (count >= 3) { INFO("Force exiting"); exit(signal); }
}


#ifdef _calc_extra_stats
void getExtraDMFTData(const sc_type& SC);
#endif

int main(int argc, char *argv[])
{
  // Catch CTRL-C
  std::signal(SIGABRT, &sighandler);
  std::signal(SIGTERM, &sighandler);
  std::signal(SIGINT , &sighandler);

  FKOptionParserDMFT opt;
   try {
        opt.parse(&argv[1], argc-1); // Skip argv[0].
        INFO("Hi! Doing Falicov-Kimball. ");
        std::cout << "FK. Parameters " << std::endl;
        std::cout << "beta                 : " << opt.beta << std::endl;
        std::cout << "U                    : " << opt.U    << std::endl;
        std::cout << "t                    : " << opt.t    << std::endl;
        std::cout << "tp                   : " << opt.tp   << std::endl;
        std::cout << "mu                   : " << opt.mu   << std::endl;
        std::cout << "e_d                  : " << opt.e_d << std::endl;
        std::cout << "Number Of Matsubaras : " << opt.n_freq << std::endl;
        std::cout << "Max number of iterations : " << opt.n_iter << std::endl;
        std::cout << "Convergence cutoff   : " << opt.cutoff << std::endl;
    } catch (const optparse::unrecognized_option& e) {
        std::cout << "unrecognized option: " << e.what() << std::endl;
        return 1;
    } catch (const optparse::invalid_value& e) {
        std::cout << "invalid value: " << e.what() << std::endl;
        return 1;
    }
    
    RealType U = opt.U;
    RealType mu = opt.mu;
    RealType e_d = opt.e_d;
    beta = opt.beta;
    RealType t = opt.t; 
    RealType tp = opt.tp; 
    size_t n_freq = opt.n_freq;
    size_t maxit = opt.n_iter;
    RealType mix = opt.mix;
    extraops = opt.extra_ops;
    size_t kpoints = opt.kpts;

    KMesh kgrid(kpoints);
    
    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1, f2;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};

    try { 
        std::string fname = "Delta_full.dat";
        GF Delta2(std::make_tuple(FMatsubaraGrid(__get_min_number(fname,beta), __get_min_number(fname,beta)+__get_n_lines(fname), beta)));
        Delta2.loadtxt(fname);
        Delta.copyInterpolate(Delta2);
        Delta._f = __fun_traits<decltype(Delta._f)>::constant(0.0);
        } 

    catch (std::exception &e) { Delta.fill(f1); };

    Delta.savetxt("Delta_0.dat");
    
    FKImpuritySolver Solver(U,mu,e_d,Delta);

    #ifdef LATTICE_bethe
    sc_type SC = BetheSC(Solver,t);
    #elif LATTICE_cubicinfd
    sc_type SC = CubicInfDMFTSC(Solver,t,RealGrid(-6.0*t,6.0*t,1024));
    #elif LATTICE_cubic1d
    sc_type SC = CubicDMFTSC<1>(Solver, kgrid, t);
    #elif LATTICE_cubic2d
    sc_type SC = CubicDMFTSC<2>(Solver, kgrid, t);
    #elif LATTICE_cubic3d
    sc_type SC = CubicDMFTSC<3>(Solver, kgrid, t);
    #elif LATTICE_cubic4d
    sc_type SC = CubicDMFTSC<4>(Solver, kgrid, t);
    #elif LATTICE_triangular
    sc_type SC = sc_type(Solver, kgrid, t, tp);
    #endif 

    bool update_weights = opt.update_weights;
    Solver.w_0 = opt.w_0;
    Solver.w_1 = opt.w_1;

    RealType diff=1.0;
    for (int i=0; i<maxit && diff>opt.cutoff &&!INTERRUPT; ++i) {
        INFO("Iteration " << i <<". Mixing = " << mix);
        update_weights = update_weights && diff/mix>1e-3;
        Solver.run(update_weights);
        Delta = SC();
        auto Delta_new = Delta*mix+(1.0-mix)*Solver.Delta;
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);
        Solver.Delta = Delta_new;
        }
   
    FMatsubaraGrid gridF_half(0, std::max(n_freq*3,size_t(100)), beta);
    GF Delta_half(gridF_half); Delta_half.copyInterpolate(Delta);
    GF gw_half(gridF_half); gw_half.copyInterpolate(Solver.gw);
    GF sigma_half(gridF_half); sigma_half.copyInterpolate(Solver.Sigma);
    sigma_half.savetxt("Sigma.dat");
    gw_half.savetxt("Gw.dat");
    Delta_half.savetxt("Delta.dat");
    int __msize = std::max(n_freq*5,size_t(1024));
    FMatsubaraGrid gridF_large(-__msize, __msize, beta);
    GF Delta_large(gridF_large); Delta_large.copyInterpolate(Delta);
    Delta_large.savetxt("Delta_full.dat");

    #ifdef _calc_extra_stats
    getExtraDMFTData(SC);
    #endif
}

#ifdef _calc_extra_stats
void getExtraDMFTData(const sc_type& SC)
{
    INFO("\nCalculating extra statistics");
    std::bitset<10> flags(extraops);

    const auto &Solver = SC._S;
    RealType T=1.0/Solver.beta;
    RealType U = Solver.U;

    if (flags[0]) {
        INFO("\nCalculating static cc, cf, ff susceptibilities at q=0, q=pi and r=0");
        size_t n_freq = std::max(int(beta*2), 512);
        FMatsubaraGrid gridF(-n_freq, n_freq, beta);
        GF gw_interp(gridF);
        gw_interp.copyInterpolate(Solver.gw);
        auto Bubbleq0 = SC.getBubble0(0.0);
        auto BubbleqPI = SC.getBubblePI(0.0); 
        std::vector<std::string> names = {"local", "pi", "zero"};
        std::vector<GF> bubbles = { -T*gw_interp*gw_interp, BubbleqPI, Bubbleq0 };
        
        auto skeleton_vals = getStaticLatticeDMFTSkeletonSusceptibility(Solver,bubbles,gridF); 
        //auto bs_vals = getStaticLatticeDMFTSusceptibility(Solver,bubbles,gridF);

        for (size_t i=0; i<bubbles.size(); ++i) { 

            /** Skeleton expansion. */
            auto chi_cc = skeleton_vals[i][0];
            auto chi_cf = skeleton_vals[i][1];
            auto chi_ff = skeleton_vals[i][2];
            
            /** Vertex expansion. */
            //auto susc = bs_vals[i];

            //INFO2("Static cc susc " << names[i] <<" (bs) = " << susc);
            INFO2("Static cc susc " << names[i] <<" (exact) = " << chi_cc);
            INFO2("Static cf susc " << names[i] <<" (exact) = " << chi_cf);
            INFO2("Static ff susc " << names[i] <<" (exact) = " << chi_ff);

            //__num_format<RealType>(susc).savetxt("StaticChiCC_" + names[i] + ".dat");
            __num_format<RealType>(chi_cc).savetxt("StaticChiCC_" + names[i] + ".dat");
            __num_format<RealType>(chi_cc).savetxt("StaticChiCC_" + names[i] + "_skeleton.dat");
            __num_format<RealType>(chi_cf).savetxt("StaticChiCF_" + names[i] + "_skeleton.dat");
            __num_format<RealType>(chi_ff).savetxt("StaticChiFF_" + names[i] + "_skeleton.dat");
            //__num_format<ComplexType>(bubbles[i].).savetxt("StaticChi0CC_" + names[i] + ".dat");
        };
    };

    if (flags[1]) {
        INFO ("Calculating static susceptibility along (pi,0)->(pi,2*pi) direction");
        size_t n_freq = 256; 
        FMatsubaraGrid local_grid(-n_freq,n_freq,beta);
        GridObject<ComplexType, FMatsubaraGrid> Lambda(local_grid);
        Lambda.copyInterpolate(Solver.getLambda());
        auto glat = SC.getGLat(local_grid); 
        GridObject<RealType,KMesh> out_b(SC._kGrid), out_chi(SC._kGrid), out_dual_bubble(SC._kGrid), out_lattice_bubble(SC._kGrid);
        std::array<KMesh::point, D> q;
        q.fill(SC._kGrid.findClosest(PI));
        for (auto q_pt : SC._kGrid.getPoints()) {
            INFO2("q = " << RealType(q_pt));
            std::get<D-1>(q) = q_pt;
            auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);

            auto lattice_bubble = SC.getBubble(0.0,q);
            //RealType susc_val = getStaticLatticeDMFTSkeletonSusceptibility(Solver, lattice_bubble, FMatsubaraGrid(-n_freq,n_freq,beta))[0];
            RealType susc_val = getStaticLatticeDMFTSkeletonSusceptibility(Solver, lattice_bubble, Solver.w_grid)[0];
            out_chi[size_t(q_pt)]= susc_val;
            
            auto dual_bubble = lattice_bubble + T*Solver.gw*Solver.gw; 
            auto Bw1 = beta*Solver.w_0*Solver.w_1*Solver.U*Solver.U*Solver.getLambda()*Solver.getLambda()*dual_bubble;
            auto Bw = Bw1/(1.0+Bw1);
            ComplexType B = Bw.sum();
            out_b[size_t(q_pt)]= std::real(B);
            
            RealType dual_bubble_val = std::real(dual_bubble.sum());
            out_dual_bubble[size_t(q_pt)]= dual_bubble_val;

            RealType lattice_bubble_val =  std::real(Diagrams::getBubble(glat, Wq_args_static).sum());
            out_lattice_bubble[size_t(q_pt)] = lattice_bubble_val;
            };
        out_chi.savetxt("Susc_dir.dat");
        out_b.savetxt("B_dir.dat");
        out_dual_bubble.savetxt("dual_bubble_dir.dat");
        out_lattice_bubble.savetxt("lattice_bubble_dir.dat");
        }
    
    if (flags[2]) {
        INFO2("Calculating B(q=pi)");
        auto BubbleqPI = SC.getBubblePI(0.0); 
        auto dual_bubble_pi = BubbleqPI + T*Solver.gw*Solver.gw;
        auto Bw1 = beta*Solver.w_0*Solver.w_1*U*U*Solver.getLambda()*Solver.getLambda()*dual_bubble_pi;
        auto Bw = Bw1/(1.0+Bw1);
        dual_bubble_pi.savetxt("DualBubbleCC_pi.dat");
        Bw1.savetxt("BwNominator_pi.dat");
        Bw.savetxt("Bw_pi.dat");
        ComplexType B = Bw.sum();
        INFO("B(pi) = " << B);
        __num_format<ComplexType>(B).savetxt("B_pi.dat");
    }

    if (flags[3]) {
        INFO2("Saving G(w,k=pi)");
        GF glat_pi(Solver.w_grid);
        auto glat = SC.getGLat(Solver.w_grid); 
        std::array<RealType, D> q;
        q.fill(PI);
        typename GF::PointFunctionType f = [&](FMatsubaraGrid::point w){return glat(std::tuple_cat(std::make_tuple(w),q));};
        glat_pi.fill(f);
        glat_pi.savetxt("glat_pi.dat");
    }

    if (flags[4]) {
        INFO2("Saving Green's functions - doing FFT");
        auto glat_k = SC.getGLat(Solver.w_grid); 
        typedef typename ArgBackGenerator<D,RealGrid,GridObject,ComplexType,FMatsubaraGrid>::type glat_r_type;
        auto grid_r = RealGrid(0,SC._kGrid.getSize(),SC._kGrid.getSize(),false);
        glat_r_type glat_r(std::tuple_cat(std::forward_as_tuple(Solver.w_grid),__repeater<RealGrid,D>::get_tuple(grid_r)));
        for (auto w : Solver.w_grid.getPoints()) { 
            glat_r[size_t(w)] = run_fft(glat_k[size_t(w)],FFTW_BACKWARD);
            }

        size_t distance = 4;
        GF glat_rp(Solver.w_grid);
        for (size_t i=0; i<distance; ++i) {
            std::array<typename decltype(grid_r)::point, D> r_p; r_p.fill(grid_r[0]); r_p[D-1]=grid_r[i]; // Makes (i,0...) point in real space
            typename GF::PointFunctionType f = [&](FMatsubaraGrid::point w){return glat_r(std::tuple_cat(std::make_tuple(w),r_p));};
            glat_rp.fill(f);
            std::stringstream fname_stream;
            fname_stream << "glat_r"; 
            for (auto rr : r_p) fname_stream << size_t(rr);
            fname_stream << ".dat";
            std::string fname; fname_stream >> fname;
            glat_rp.savetxt(fname);
            };
        };

    if (flags[5]) {
        INFO("\nCalculating static cc susceptibility(q)");
        auto bzpoints_map = CubicTraits<D>::getUniqueBZPoints(SC._kGrid); 
        std::vector<BZPoint<D>> bzpoints;
        std::vector<GF> bubbles;
        INFO2("Preparing bare bubbles");
        for (auto map_it = bzpoints_map.begin(); map_it!=bzpoints_map.end(); map_it++){
            INFO3(std::distance(bzpoints_map.begin(),map_it)+1<<" / "<< bzpoints_map.size());
            bzpoints.push_back(map_it->first);
            bubbles.push_back(SC.getBubble(0.0,map_it->first));
            }
        INFO2("done.");
        size_t n_freq = std::max(int(beta*2), 1024);
        auto stat_susc_bz = getStaticLatticeDMFTSkeletonSusceptibility(Solver, bubbles, FMatsubaraGrid(-n_freq,n_freq,beta));
        std::map<BZPoint<D>, RealType> susc_map;
        for (size_t nq=0; nq<stat_susc_bz.size(); ++nq) susc_map[bzpoints[nq]] = std::real(stat_susc_bz[nq][0]);
        auto all_bz_points = CubicTraits<D>::getAllBZPoints(SC._kGrid);
        size_t nqpts = all_bz_points.size();
        size_t dimsize = SC._kGrid.getSize();
        std::ofstream out;
        out.open("StaticChiDMFTCC.dat");
        for (size_t nq=0; nq<nqpts; ++nq) {
            BZPoint<D> current_point = all_bz_points[nq];
            BZPoint<D> sym_point = CubicTraits<D>::findSymmetricBZPoint(current_point,SC._kGrid);
            out << all_bz_points[nq] << susc_map[sym_point] << std::endl;
            if ((nq+1)%dimsize==0) out << std::endl;
        }
        //__num_format<ComplexType>(stat_susc_pi).savetxt("StaticChiDFCC_pi.dat");
        out.close();
    };

    if (flags[6]) {
        INFO2("Saving G(w,k)");
        auto glat = SC.getGLat(Solver.w_grid); 
        GF glat_k(Solver.w_grid);
        auto all_bz_points = CubicTraits<D>::getAllBZPoints(SC._kGrid);
        size_t nqpts = all_bz_points.size();
        for (size_t nq=0; nq<nqpts; ++nq) {
            BZPoint<D> q = all_bz_points[nq];
            typename GF::PointFunctionType f = [&](FMatsubaraGrid::point w){return glat(std::tuple_cat(std::make_tuple(w),q));};
            glat_k.fill(f);
            std::stringstream s1;
            s1 << "glat_k" << RealType(std::get<0>(q)) << "_" << RealType(std::get<1>(q)) << ".dat";
            std::string s; s1>>s;
            glat_k.savetxt(s);
        }
    };



/*    if (extraops>=2) { 
    INFO("Dynamic susceptibility");
    size_t n_b_freq = std::max(std::min(Solver.w_grid._w_max/2,int(2*beta)),10);
    BMatsubaraGrid gridB(-n_b_freq, n_b_freq+1, beta);
    GridObject<RealType,BMatsubaraGrid> chi0_q0_vals(gridB), chi0_qPI_vals(gridB);
    GridObject<RealType,BMatsubaraGrid> chi_q0_vals(gridB), chi_qPI_vals(gridB), chi_q0_dmft_vals(gridB), chi_qPI_dmft_vals(gridB);

    for (auto iW : gridB.getPoints()) {
        if (INTERRUPT) exit(0);
        INFO("iW = " << iW);
        size_t iWn = size_t(iW);
        GF Vertex4(gridF);
        Vertex4.fill(typename GF::PointFunctionType([&](FMatsubaraGrid::point w){return Solver.getBVertex4(iW,w);}));
        auto gw_bubble = Diagrams::getBubble(gw, iW);
        auto Chi0q0 = SC.getBubble0(iW);
        auto Chiq0 = Diagrams::getSusc<GF>(Chi0q0, Diagrams::BS(Chi0q0 - gw_bubble, Vertex4 , true));
        auto Chi0qPI = SC.getBubblePI(iW);
        auto ChiqPI = Diagrams::getSusc<GF>(Chi0qPI, Diagrams::BS(Chi0qPI - gw_bubble, Vertex4, true)); 
        chi_q0_vals[iWn] = std::real(Chiq0.sum());
        chi0_q0_vals[iWn] = std::real(Chi0q0.sum());
        chi0_qPI_vals[iWn] = std::real(Chi0qPI.sum());
        chi_qPI_vals[iWn] = std::real(ChiqPI.sum());
        auto chiq0_dmft = -T/ComplexType(iW)*(gw-gw.shift(iW)).sum();
        chi_q0_dmft_vals[size_t(iW)] = std::real(chiq0_dmft); 

        if (is_equal(ComplexType(iW),0.0)) { 
            INFO("Static val = " << chi_q0_vals[iWn])
            chi_q0_dmft_vals[iWn] = chi_q0_vals[iWn];
            };
        };

    INFO("Chi0[q=0]     = " << chi0_q0_vals);
    INFO("Chi0[q=PI]     = " << chi0_qPI_vals);

    INFO("Chi[q=0]     = " << chi_q0_vals);
    INFO("ChiDMFT[q=0] = " << chi_q0_dmft_vals);
    INFO("Full Chi, q=0 diff = " << chi_q0_vals.diff(chi_q0_dmft_vals));

    chi0_q0_vals.savetxt("DynamicChi0q0.dat");
    chi0_qPI_vals.savetxt("DynamicChi0qPI.dat");
    chi_q0_vals.savetxt("DynamicChiq0.dat");
    chi_q0_dmft_vals.savetxt("DynamicChiq0_DMFT.dat");
    chi_qPI_vals.savetxt("DynamicChiqPI.dat");
    };
*/
}

#endif

int __get_n_lines(const std::string& fname)
{   
    std::ifstream f(fname.c_str()); 
    std::string line; 
    int nlines = 0; 
    while (std::getline(f, line)) ++nlines; 
    f.close(); 
    return nlines;
} 

int __get_min_number(const std::string& fname, RealType beta)
{
    std::ifstream f(fname.c_str());
    __num_format<FMatsubaraGrid::point> tmp(FMatsubaraGrid::point(0,0));
    f >> tmp;
    f.close();
    ComplexType w_min = tmp._v;
    int w_min_index = FMatsubaraIndex(w_min,beta);
    return w_min_index;
}


