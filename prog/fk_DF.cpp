#include <numeric>

#include <gftools/enum_grid.hpp>
#include "Solver.h"
#include "DF.h"
#include "FFT.hpp"

#ifdef LATTICE_cubic1d
    typedef FK::CubicDMFTSC<1> dmft_sc_type;
    typedef FK::DFLadderCubic<1> df_sc_type;
    static constexpr size_t D=1;
    #define _calc_extra_stats
#elif LATTICE_cubic2d
    typedef FK::CubicDMFTSC<2> dmft_sc_type;
    typedef FK::DFLadderCubic<2> df_sc_type;
    static constexpr size_t D=2;
    #define _calc_extra_stats
#elif LATTICE_cubic3d
    typedef FK::CubicDMFTSC<3> dmft_sc_type;
    typedef FK::DFLadderCubic<3> df_sc_type;
    static constexpr size_t D=3;
    #define _calc_extra_stats
#elif LATTICE_cubic4d
    typedef FK::CubicDMFTSC<4> dmft_sc_type;
    typedef FK::DFLadderCubic<4> df_sc_type;
    static constexpr size_t D=4;
    #define _calc_extra_stats
#elif LATTICE_triangular
    typedef FK::TriangularDMFT dmft_sc_type;
    typedef FK::DFLadder<FK::TriangularTraits,2> df_sc_type;
    static constexpr size_t D=2;
    #define _calc_extra_stats
#endif

#include "FKOptionParserDF.h"

#include <iostream>
#include <ctime>
#include <array>
#include <unordered_map>
#include <csignal>
#include <bitset>

using namespace gftools;
using namespace FK;

int __get_n_lines(const std::string& fname);
int __get_min_number(const std::string& fname, real_type beta);

//extern template grid_object<complex_type,fmatsubara_grid>;

real_type beta;
size_t extraops=0;
bool INTERRUPT = false;
std::string dual_sigma_file;

typedef grid_object<complex_type,fmatsubara_grid> GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, real_type tolerance = 1e-7)
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
template <class SCType> void getExtraData(SCType& SC, const fmatsubara_grid& gridF);
#endif


int main(int argc, char *argv[])
{
  // Catch CTRL-C
  std::signal(SIGABRT, &sighandler);
  std::signal(SIGTERM, &sighandler);
  std::signal(SIGINT , &sighandler);

  FKOptionParserDF opt;
   try {
        opt.parse(&argv[1], argc-1); // Skip argv[0].
        INFO("Hi! Doing Falicov-Kimball. ");
        std::cout << "FK. Parameters " << std::endl;
        std::cout << "beta                 : " << opt.beta << std::endl;
        std::cout << "T                    : " << 1.0/opt.beta << std::endl;
        std::cout << "U                    : " << opt.U    << std::endl;
        std::cout << "t                    : " << opt.t    << std::endl;
        std::cout << "tp                   : " << opt.tp   << std::endl;
        std::cout << "mu                   : " << opt.mu   << std::endl;
        std::cout << "e_d                  : " << opt.e_d << std::endl;
        std::cout << "Number Of Matsubaras : " << opt.n_freq << std::endl;
        std::cout << "Max number of DMFT iterations : " << opt.NDMFTRuns << std::endl;
        std::cout << "Max number of DF   iterations : " << opt.NDFRuns << std::endl;
        INFO("Extra options flag: " << opt.extraops);
    } catch (const optparse::unrecognized_option& e) {
        std::cout << "unrecognized option: " << e.what() << std::endl;
        return 1;
    } catch (const optparse::invalid_value& e) {
        std::cout << "invalid value: " << e.what() << std::endl;
        return 1;
    }
    
    real_type U = opt.U;
    real_type mu = opt.mu;
    real_type e_d = opt.e_d;
    beta = opt.beta;
    real_type t = opt.t; 
    real_type tp = opt.tp; 
    size_t n_freq = opt.n_freq;
    size_t ksize = opt.kpts;
    //size_t n_dual_freq = opt.n_dual_freq;
    real_type mix = opt.mix;
    extraops = opt.extraops;
    size_t NDMFTRuns = opt.NDMFTRuns;
    size_t NDFRuns = opt.NDFRuns;
    real_type DFCutoff = opt.DFCutoff;
    bool update_mixing = opt.update_mixing;
    bool read_dual_sigma = opt.read_dual_sigma;
    dual_sigma_file = opt.dual_sigma_file;

    kmesh kGrid(ksize);

    fmatsubara_grid gridF(-n_freq, n_freq, beta);
    fmatsubara_grid gridF_half(0, 2*n_freq, beta);
    int __msize = std::max(n_freq*5,size_t(1024));
    fmatsubara_grid gridF_large(-__msize, __msize, beta);

    GF Delta(gridF);
    std::function<complex_type(complex_type)> f1, f2;
    f1 = [t](complex_type w) -> complex_type {return t*t/w;};

    try { 
        std::string fname = "Delta_full.dat";
        GF Delta2(std::make_tuple(fmatsubara_grid(__get_min_number(fname,beta), __get_min_number(fname,beta)+__get_n_lines(fname), beta)));
        Delta2.loadtxt(fname);
        Delta.copy_interpolate(Delta2);
        Delta.tail_ = tools::fun_traits<typename GF::function_type>::constant(0.0);
        } 
    catch (std::exception &e) { Delta.fill(f1); };
    Delta.savetxt("Delta_0.dat");
    
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    
    kmesh_patch qGrid(kGrid);

    #ifdef LATTICE_triangular
    dmft_sc_type SC_DMFT(Solver, kGrid, t, tp);
    df_sc_type SC_DF(Solver, gridF, kGrid, t, tp);
    std::cout << "==> Triangular lattice" << std::endl;
    #else
    dmft_sc_type SC_DMFT(Solver, kGrid, t);
    df_sc_type SC_DF(Solver, gridF, kGrid, t);
    std::cout << "==> Hypercubic lattice" << std::endl;
    #endif 

    SC_DF._n_GD_iter = opt.DFNumberOfSelfConsistentIterations;
    SC_DF._GDmix = opt.DFSCMixing;
    SC_DF._SC_cutoff = opt.DFSCCutoff;
    SC_DF._eval_BS_SC = opt.DFEvaluateBSSelfConsistent;
    SC_DF._bs_evaluate_only_order_n = opt.DFEvaluateBSOnlyN;
    SC_DF._n_BS_iter = opt.DFNumberOfBSIterations;
    SC_DF._BSmix = opt.DFBSMixing;
    SC_DF._EvaluateStaticDiagrams = opt.DFEvaluateStaticDiagrams;
    SC_DF._EvaluateDynamicDiagrams = opt.DFEvaluateDynamicDiagrams;

    bool update_weights = opt.update_weights;
    bool read_weights = opt.read_weights;
    Solver.w_0 = opt.w_0;
    Solver.w_1 = opt.w_1;

    if (read_weights) {
        try {
        		double s=0.0;
                num_io<real_type> wf(s);
                wf.loadtxt("w_0.dat");
                Solver.w_0 = wf;
                wf.loadtxt("w_1.dat");
                Solver.w_1 = wf;
                if (std::abs(Solver.w_0 + Solver.w_1 - 1.0)<std::numeric_limits<real_type>::epsilon()) throw(1);
            }
        catch (...)
            {
                INFO("Couldn't load weights from file.");
            }
        };

    // read dual self-energy
    if (read_dual_sigma) {
        try { 
            INFO("Loading dual self-energy from "<<dual_sigma_file);
            SC_DF.SigmaD.loadtxt(dual_sigma_file);
            } 
        catch (std::exception &e) { 
            INFO("Starting from zero dual self-energy");
            SC_DF.SigmaD = 0.0; 
            };
    };

    real_type diff=1.0, diff_min = 1.0;
    size_t diff_min_count = 0;
    std::ofstream diff_stream("diff.dat",std::ios::out);
    diff_stream.close();
    bool calc_DMFT = (NDMFTRuns>0);

    size_t i_dmft = 0; 
    size_t i_df = 0;

    for (; i_dmft<=NDMFTRuns-calc_DMFT && i_df<=NDFRuns && diff>1e-8+(1-calc_DMFT)*(DFCutoff-1e-8) &&!INTERRUPT; (calc_DMFT)?i_dmft++:i_df++) {
        INFO("Iteration " << i_dmft+i_df <<". Mixing = " << mix);

        update_weights = update_weights && diff/mix>1e-3 && calc_DMFT;
        Solver.run(update_weights);

        if (calc_DMFT) {  
            Delta = SC_DMFT();
            }
        else { 
            Delta = SC_DF();
            std::stringstream tmp;
            tmp << "GwDF" << i_df << ".dat";
            Solver.gw.savetxt(tmp.str());
            tmp.str( std::string() ); tmp.clear();
            tmp << "DeltaDF" << i_df << ".dat";
            Delta.savetxt(tmp.str());
            tmp.str( std::string() ); tmp.clear();
            tmp << "glocDF" << i_df << ".dat";
            SC_DF.getGLoc().savetxt(tmp.str());
             }
        auto Delta_new = Delta*mix+(1.0-mix)*Solver.Delta;
        diff = Delta_new.diff(Solver.Delta);
        INFO("diff = " << diff);

        if (diff<diff_min) { diff_min = diff; diff_min_count = 0; }
            else diff_min_count++;
        if (diff_min_count > 5 ) {
            ERROR("\n\tCaught loop cycle. Reducing main loop mixing to " << mix/2. << " .\n");
            mix=std::max(mix/2., 0.01);
            diff_min_count = 0;
            diff_min = diff;
            };
        if (!calc_DMFT) {
            diff_stream.open("diff.dat",std::ios::app);
            diff_stream << diff << "  " << mix << std::endl;
            diff_stream.close();
        };
        
        Solver.Delta = Delta_new;
        if (diff<=1e-8 && calc_DMFT) { 
            GF Delta_large(gridF_large); Delta_large.copy_interpolate(Solver.Delta);
            Delta_large.savetxt("DeltaDMFT.dat");
            Solver.Sigma.savetxt("SigmaDMFT.dat"); 
            Solver.gw.savetxt("GwDMFT.dat");
            num_io<real_type>(Solver.w_0).savetxt("w_0DMFT.dat");
            num_io<real_type>(Solver.w_1).savetxt("w_1DMFT.dat");
            diff = 1.0; calc_DMFT = false; }; // now continue with DF 
        }
   
    GF Delta_half(gridF_half); Delta_half.copy_interpolate(Delta);
    GF gw_half(gridF_half); gw_half.copy_interpolate(Solver.gw);
    GF sigma_half(gridF_half); sigma_half.copy_interpolate(Solver.Sigma);
    gw_half.savetxt("Gw.dat");
    //Delta_large.savetxt("Delta_full.dat");
    num_io<real_type>(Solver.w_0).savetxt("w_0.dat");
    num_io<real_type>(Solver.w_1).savetxt("w_1.dat");

    fmatsubara_grid stat_grid(-128, 128, beta);
    #ifdef _calc_extra_stats
    switch (D){
        case 1: getExtraData(SC_DF,fmatsubara_grid(-256, 256, beta)); break;
        case 2: getExtraData(SC_DF,fmatsubara_grid(-512, 512, beta)); break;
        case 3: getExtraData(SC_DF,fmatsubara_grid(-128, 128, beta)); break;
        case 4: getExtraData(SC_DF,gridF); break; 
    }
    #endif
}



#ifdef _calc_extra_stats
template <class SCType> void getExtraData(SCType& SC, const fmatsubara_grid& gridF)
{
    constexpr size_t D = SCType::_D;
    INFO("Calculating additional statistics.");
    const auto &Solver = SC._S;
    real_type beta = Solver.beta;
    real_type T=1.0/beta;
    SC.GLatLoc.savetxt("gloc.dat");
    const auto &kgrid = SC._kGrid;

    std::bitset<10> flags(extraops);

    typename SCType::GKType dual_bubbles = Diagrams::getStaticBubbles(SC.GD); 
    typename SCType::GKType dual_bubbles0 = Diagrams::getStaticBubbles(SC.GD0); 
    typename SCType::GKType lattice_bubbles = Diagrams::getStaticBubbles(SC.getGLat()); 
    typedef typename SCType::GLocalType GLocalType;

    if (flags[0]) {
        std::array<real_type, D> q;
        q.fill(PI);
        size_t n_freq;
        complex_type stat_susc_pi;
        switch(D)
        {
            case 1: 
                n_freq = std::max(int(beta*2), 1024); 
                stat_susc_pi = SC.getStaticLatticeSusceptibility(q, fmatsubara_grid(-n_freq,n_freq,beta));
                break;
            case 2: 
                stat_susc_pi = SC.getStaticLatticeSusceptibility(q, Solver.w_grid);
                //n_freq = std::max(int(beta*2), 512); 
                //stat_susc_pi = SC.getStaticLatticeSusceptibility(q, fmatsubara_grid(-n_freq,n_freq,beta));
                break;
            case 3: 
                stat_susc_pi = SC.getStaticLatticeSusceptibility(q, Solver.w_grid);
                break;
            case 4: 
                stat_susc_pi = SC.getStaticLatticeSusceptibility(q, Solver.w_grid);
                break;
        };
        num_io<complex_type>(stat_susc_pi).savetxt("StaticChiDFCC_pi.dat");
    };

    if (flags[1]) {
        INFO ("\nCalculating static susceptibility along (pi,0)->(pi,2*pi) direction");
        size_t n_freq = 256; 
        fmatsubara_grid local_grid(-n_freq,n_freq,beta);
        grid_object<complex_type, fmatsubara_grid> Lambda(local_grid);
        Lambda.copy_interpolate(Solver.getLambda());
        auto glat = SC.getGLat(); 
        grid_object<real_type,kmesh> out_b(SC._kGrid), out_chi(SC._kGrid), out_dual_bubble(SC._kGrid), out_lattice_bubble(SC._kGrid), out_dual_bubble0(SC._kGrid);
        std::array<typename kmesh::point, D> q = tuple_tools::repeater<typename kmesh::point, D>::get_array(SC._kGrid.find_nearest(PI));
        q.fill(SC._kGrid.find_nearest(PI));

        GLocalType dual_bubble(dual_bubbles.grid());
        GLocalType dual_bubble0(dual_bubbles.grid());
        GLocalType lattice_bubble(dual_bubbles.grid());

        for (auto q_pt : SC._kGrid.points()) {
            INFO2("q = " << real_type(q_pt));
            std::get<D-1>(q) = q_pt;
            auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);

            complex_type susc_val = SC.getStaticLatticeSusceptibility(q,Solver.w_grid);
            out_chi[size_t(q_pt)]= std::real(susc_val);
            
            dual_bubble.fill([&](typename fmatsubara_grid::point w){return dual_bubbles(std::tuple_cat(std::make_tuple(w), q)); });
            dual_bubble0.fill([&](typename fmatsubara_grid::point w){return dual_bubbles0(std::tuple_cat(std::make_tuple(w), q)); });
            auto Bw1 = beta*Solver.w_0*Solver.w_1*Solver.U*Solver.U*Solver.getLambda()*Solver.getLambda()*dual_bubble;
            auto Bw = Bw1/(1.0+Bw1);
            complex_type B = Bw.sum();
            out_b[size_t(q_pt)]= std::real(B);
            
            real_type dual_bubble_val = std::real(dual_bubble.sum());
            out_dual_bubble[size_t(q_pt)]= dual_bubble_val;
            out_dual_bubble0[size_t(q_pt)]= std::real(dual_bubble0.sum());

            lattice_bubble.fill([&](typename fmatsubara_grid::point w){return lattice_bubbles(std::tuple_cat(std::make_tuple(w), q)); });
            real_type lattice_bubble_val =  std::real(lattice_bubble.sum());
            out_lattice_bubble[size_t(q_pt)] = lattice_bubble_val;
            };
        out_chi.savetxt("Susc_dir.dat");
        out_b.savetxt("B_dir.dat");
        out_dual_bubble.savetxt("dual_bubble_dir.dat");
        out_dual_bubble0.savetxt("dual_bubble0_dir.dat");
        out_lattice_bubble.savetxt("lattice_bubble_dir.dat");
        }
        
    
    if (flags[2]) {
        INFO2("Calculating B(q=pi)");
        std::array<real_type, D> q;
        q.fill(PI);
        auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);
        GLocalType dual_bubble_pi (SC.GD.grid());
        dual_bubble_pi.fill([&](typename fmatsubara_grid::point w){return dual_bubbles(std::tuple_cat(std::make_tuple(w), q)); });
        GLocalType dual_bubble0_pi (SC.GD0.grid());
        dual_bubble0_pi.fill([&](typename fmatsubara_grid::point w){return dual_bubbles0(std::tuple_cat(std::make_tuple(w), q)); });
        auto Bw1 = beta*Solver.w_0*Solver.w_1*Solver.U*Solver.U*Solver.getLambda()*Solver.getLambda()*dual_bubble_pi;
        auto Bw = Bw1/(1.0+Bw1);
        complex_type B = Bw.sum();
        dual_bubble_pi.savetxt("dual_bubble_pi_w.dat");
        dual_bubble0_pi.savetxt("dual_bubble0_pi_w.dat");
        Bw1.savetxt("BwNominator_pi.dat");
        Bw.savetxt("Bw_pi.dat");
        INFO("B(pi) = " << B);
        num_io<complex_type>(B).savetxt("B_pi.dat");
    }

    if (flags[3]) {
        INFO2("Saving G(w,k=pi)");
        GF glat_pi(Solver.w_grid);
        auto glat = SC.getGLat(); 
        std::array<real_type, D> q;
        q.fill(PI);
        typename GF::point_function_type f = [&](fmatsubara_grid::point w){return glat(std::tuple_cat(std::make_tuple(w),q));};
        glat_pi.fill(f);
        glat_pi.savetxt("glat_pi.dat");
    }

    if (flags[4]) {
        INFO2("Saving Green's functions - doing FFT");
        auto glat_k = SC.getGLat(); 
        typedef typename tools::ArgBackGenerator<D,enum_grid,grid_object,complex_type,fmatsubara_grid>::type glat_r_type;
        auto grid_r = enum_grid(0,SC._kGrid.size(),false);
        glat_r_type glat_r(std::tuple_cat(std::forward_as_tuple(SC._fGrid),tuple_tools::repeater<enum_grid,D>::get_tuple(grid_r)));
        for (auto w : SC._fGrid.points()) { 
            glat_r[size_t(w)] = run_fft(glat_k[size_t(w)],FFTW_BACKWARD);
            }

        size_t distance = std::min(5,int(kgrid.size()));
        GF glat_rp(Solver.w_grid);
        for (size_t i=0; i<distance; ++i) {
            std::array<typename enum_grid::point, D> r_p = tuple_tools::repeater<typename enum_grid::point, D>::get_array(grid_r[0]);
            r_p[D-1]=grid_r[i]; // Makes (i,0...) point in real space
            typename GF::point_function_type f = [&](fmatsubara_grid::point w){return glat_r(std::tuple_cat(std::make_tuple(w),r_p));};
            glat_rp.fill(f);
            std::stringstream fname_stream;
            fname_stream << "glat_r" << i << ".dat";
            //for (auto rr : r_p) fname_stream << size_t(rr);
            std::string fname; fname_stream >> fname;
            glat_rp.savetxt(fname);
            };
    }



    if (flags[5]) {
        auto glat = SC.getGLat(); 
        auto bzpoints_map = CubicTraits<D>::getUniqueBZPoints(SC._kGrid); 
        std::vector<BZPoint<D>> bzpoints;
        for (auto map_it = bzpoints_map.begin(); map_it!=bzpoints_map.end(); map_it++){
            bzpoints.push_back(map_it->first);
            }
        size_t n_freq = Solver.gw.grid().size(); //std::max(int(beta*2), 256);
        auto stat_susc_bz = SC.getStaticLatticeSusceptibility(bzpoints, fmatsubara_grid(-n_freq,n_freq,beta));
        std::map<BZPoint<D>, real_type> susc_map, lattice_bubble_map, dual_bubble_map;
        INFO_NONEWLINE("\tGetting bubbles [" << stat_susc_bz.size() <<"] : ");

        GLocalType dual_bubble(dual_bubbles.grid());
        GLocalType lattice_bubble(lattice_bubbles.grid());

        for (size_t nq=0; nq<stat_susc_bz.size(); ++nq) { 
            INFO_NONEWLINE(nq << " "); 
            dual_bubble.fill([&](typename fmatsubara_grid::point w){return dual_bubbles(std::tuple_cat(std::make_tuple(w), bzpoints[nq])); });
            lattice_bubble.fill([&](typename fmatsubara_grid::point w){return lattice_bubbles(std::tuple_cat(std::make_tuple(w), bzpoints[nq])); });
            susc_map[bzpoints[nq]] = std::real(stat_susc_bz[nq]);
            lattice_bubble_map[bzpoints[nq]] = std::real(lattice_bubble.sum());
            dual_bubble_map[bzpoints[nq]] = std::real(dual_bubble.sum()); 
        };
        INFO("");

        auto all_bz_points = CubicTraits<D>::getAllBZPoints(SC._kGrid);
        size_t nqpts = all_bz_points.size();

        typedef typename tools::ArgBackGenerator<D,kmesh        ,grid_object,complex_type>::type susc_k_type;
        typedef typename tools::ArgBackGenerator<D,enum_grid,grid_object,complex_type>::type susc_r_type;

        auto kmeshes = tuple_tools::repeater<kmesh,D>::get_tuple(SC._kGrid);
        auto rgrid = enum_grid(0,SC._kGrid.size(),false);
        auto rmeshes = tuple_tools::repeater<enum_grid,D>::get_tuple(rgrid);

        susc_k_type susc_full(kmeshes), susc_lattice_bubbles(kmeshes), gd_bubbles(kmeshes);
        susc_r_type susc_full_r(rmeshes), susc_lattice_bubbles_r(rmeshes), gd_bubbles_r(rmeshes);

        for (size_t nq=0; nq<nqpts; ++nq) {
            BZPoint<D> current_point = all_bz_points[nq];
            typename susc_k_type::point_tuple current_point_tuple = std::tuple_cat(current_point);
            BZPoint<D> sym_point = CubicTraits<D>::findSymmetricBZPoint(current_point,SC._kGrid);
            susc_full.get(current_point_tuple) = susc_map[sym_point];
            susc_lattice_bubbles.get(current_point_tuple) = lattice_bubble_map[sym_point];
            gd_bubbles.get(current_point_tuple) = dual_bubble_map[sym_point];
            };

        susc_full.savetxt("StaticChiDFCC.dat");
        susc_lattice_bubbles.savetxt("StaticLatticeBubbleDFCC.dat");
        gd_bubbles.savetxt("StaticDualBubbleCC.dat");

        INFO2("Doing FFT of susceptibilities");

        susc_full_r.data() = run_fft(susc_full.data(),FFTW_BACKWARD);
        susc_lattice_bubbles_r.data() = run_fft(susc_lattice_bubbles.data(),FFTW_BACKWARD);
        gd_bubbles_r.data() = run_fft(gd_bubbles.data(),FFTW_BACKWARD);

        susc_full_r.savetxt("StaticChiDFCC_r.dat");
        susc_lattice_bubbles_r.savetxt("StaticLatticeBubbleDFCC_r.dat");
        gd_bubbles_r.savetxt("StaticDualBubbleCC_r.dat");

        INFO2("Saving susceptibilities in different directions");

        grid_object<real_type,enum_grid> cc_x (rgrid), cc_y(rgrid), cc_xy(rgrid); 
        for (auto p1:rgrid.points()) {
            auto zero_pt = rgrid.find_nearest(0);
            std::array<typename enum_grid::point, D> xpt = tuple_tools::repeater<typename enum_grid::point, D>::get_array(zero_pt);
            std::array<typename enum_grid::point, D> ypt = xpt;
            std::array<typename enum_grid::point, D> xypt = xpt;
            xpt.fill(zero_pt); ypt.fill(zero_pt); xypt.fill(zero_pt);
            xpt[D-1]=p1; 
            xypt[D-1]=p1;
            ypt[std::max(int(D),2)-2]=p1;
            xypt[std::max(int(D),2)-2]=p1;

            cc_x.get(p1) = std::real(susc_full_r(typename susc_r_type::point_tuple(std::tuple_cat(xpt))));
            cc_y.get(p1) = std::real(susc_full_r(typename susc_r_type::point_tuple(std::tuple_cat(ypt))));
            cc_xy.get(p1) = std::real(susc_full_r(typename susc_r_type::point_tuple(std::tuple_cat(xypt))));
            }

        cc_x.savetxt("StaticChiDFCC_dir_x.dat");
        cc_y.savetxt("StaticChiDFCC_dir_y.dat");
        cc_xy.savetxt("StaticChiDFCC_dir_xy.dat");
        };


    if (flags[6]) {
        INFO2("Saving dual self-energy");
        SC.SigmaD.savetxt(dual_sigma_file);
        };

    if (flags[7]){
        INFO2("G(k)");
        auto glat = SC.getGLat(); 
        auto glatdmft = SC.getGLatDMFT(SC._fGrid);

        auto w = Solver.w_grid.find_nearest(I*PI/beta);
        auto rgrid = enum_grid(0,SC._kGrid.size(),false);

        typedef typename tools::ArgBackGenerator<D,kmesh        ,grid_object,complex_type>::type gd_k_type;
        typedef typename tools::ArgBackGenerator<D,enum_grid,grid_object,complex_type>::type gd_r_type;

        auto kmeshes = tuple_tools::repeater<kmesh,D>::get_tuple(kgrid);
        auto rmeshes = tuple_tools::repeater<enum_grid,D>::get_tuple(rgrid);

        gd_k_type glat_k(kmeshes), gd_k(kmeshes), gd0_k(kmeshes), glat_dmft_k(kmeshes), sigmad_k(kmeshes);
        gd_r_type glat_r(rmeshes), gd_r(rmeshes), gd0_r(rmeshes), glat_dmft_r(rmeshes), sigmad_r(rmeshes);

        glat_k.data() = glat[w.index_];
        gd_k.data() = SC.GD[w.index_];
        gd0_k.data() = SC.GD0[w.index_];
        glat_dmft_k.data() = glatdmft[w.index_];
        sigmad_k.data() = SC.SigmaD[w.index_];

        glat_k.savetxt("glat_k_w0.dat");
        gd_k.savetxt("gd_k_w0.dat");
        gd0_k.savetxt("gd0_k_w0.dat");
        glat_dmft_k.savetxt("glat_dmft_k_w0.dat");
        sigmad_k.savetxt("sigmad_k_w0.dat");

        glat_r.data() = run_fft(glat_k.data(), FFTW_BACKWARD);
        gd_r.data() = run_fft(gd_k.data(), FFTW_BACKWARD);
        gd0_r.data() = run_fft(gd0_k.data(), FFTW_BACKWARD);
        glat_dmft_r.data() = run_fft(glat_dmft_k.data(), FFTW_BACKWARD);
        sigmad_r.data() = run_fft(sigmad_k.data(), FFTW_BACKWARD);

        glat_r.savetxt("glat_r_w0.dat");
        gd_r.savetxt("gd_r_w0.dat");
        gd0_r.savetxt("gd0_r_w0.dat");
        glat_dmft_r.savetxt("glat_dmft_r_w0.dat");
        sigmad_r.savetxt("sigmad_r_w0.dat");
        };

    if (flags[8]){
        INFO2("Saving G(w,k)");
        auto glat = SC.getGLat(); 
        glat.savetxt("glat_k.dat");
        SC.GD0.savetxt("gd0_k.dat");
        SC.GD.savetxt("gd_k.dat");
        SC.getGLatDMFT(SC._fGrid).savetxt("glat_dmft_k.dat");
        };

    if (flags[9]){
        INFO2("Saving bare 4-point vertex (gamma4)");
        typedef grid_object<std::complex<double>,fmatsubara_grid,fmatsubara_grid> VertexType;
        VertexType StaticVertex4(std::forward_as_tuple(SC._fGrid,SC._fGrid)); 
        grid_object<std::complex<double>,fmatsubara_grid,fmatsubara_grid>::point_function_type VertexF2 = [&](fmatsubara_grid::point w1, fmatsubara_grid::point w2){return Solver.getVertex4(0.0, w1,w2);};
        StaticVertex4.fill(VertexF2);
        StaticVertex4.savetxt("gamma4.dat");
        }

 }
#endif // endif :: #ifdef _calc_extra_stats

int __get_n_lines(const std::string& fname)
{   
    std::ifstream f(fname.c_str()); 
    std::string line; 
    int nlines = 0; 
    while (std::getline(f, line)) ++nlines; 
    f.close(); 
    return nlines;
} 

int __get_min_number(const std::string& fname, real_type beta)
{
    std::ifstream f(fname.c_str());
    complex_type x;
    num_io<complex_type> tmp(x);
    f >> tmp;
    f.close();
    complex_type w_min = tmp.value_;
    int w_min_index = FMatsubaraIndex(w_min,beta);
    return w_min_index;
}


