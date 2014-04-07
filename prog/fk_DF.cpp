#include <numeric>

#include <MatsubaraGrid.hpp>
#include <KMesh.hpp>
#include <GridObject.hpp>
#include <EnumerateGrid.hpp>
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

using namespace GFTools;
using namespace FK;

int __get_n_lines(const std::string& fname);
int __get_min_number(const std::string& fname, RealType beta);

//extern template GridObject<ComplexType,FMatsubaraGrid>;

RealType beta;
size_t extraops=0;
bool INTERRUPT = false;
std::string dual_sigma_file;

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
template <class SCType> void getExtraData(SCType& SC, const FMatsubaraGrid& gridF);
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
    
    RealType U = opt.U;
    RealType mu = opt.mu;
    RealType e_d = opt.e_d;
    beta = opt.beta;
    RealType t = opt.t; 
    RealType tp = opt.tp; 
    size_t n_freq = opt.n_freq;
    size_t ksize = opt.kpts;
    //size_t n_dual_freq = opt.n_dual_freq;
    RealType mix = opt.mix;
    extraops = opt.extraops;
    size_t NDMFTRuns = opt.NDMFTRuns;
    size_t NDFRuns = opt.NDFRuns;
    RealType DFCutoff = opt.DFCutoff;
    bool update_mixing = opt.update_mixing;
    bool read_dual_sigma = opt.read_dual_sigma;
    dual_sigma_file = opt.dual_sigma_file;

    KMesh kGrid(ksize);

    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    FMatsubaraGrid gridF_half(0, 2*n_freq, beta);
    int __msize = std::max(n_freq*5,size_t(1024));
    FMatsubaraGrid gridF_large(-__msize, __msize, beta);

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
    
    KMeshPatch qGrid(kGrid);

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
                __num_format<RealType> wf(0.0); 
                wf.loadtxt("w_0.dat");
                Solver.w_0 = wf;
                wf.loadtxt("w_1.dat");
                Solver.w_1 = wf;
                if (std::abs(Solver.w_0 + Solver.w_1 - 1.0)<std::numeric_limits<RealType>::epsilon()) throw(1);
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

    RealType diff=1.0, diff_min = 1.0;
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
            GF Delta_large(gridF_large); Delta_large.copyInterpolate(Solver.Delta);
            Delta_large.savetxt("DeltaDMFT.dat");
            Solver.Sigma.savetxt("SigmaDMFT.dat"); 
            Solver.gw.savetxt("GwDMFT.dat");
            __num_format<RealType>(Solver.w_0).savetxt("w_0DMFT.dat");
            __num_format<RealType>(Solver.w_1).savetxt("w_1DMFT.dat");
            diff = 1.0; calc_DMFT = false; }; // now continue with DF 
        }
   
    GF Delta_half(gridF_half); Delta_half.copyInterpolate(Delta);
    GF gw_half(gridF_half); gw_half.copyInterpolate(Solver.gw);
    GF sigma_half(gridF_half); sigma_half.copyInterpolate(Solver.Sigma);
    gw_half.savetxt("Gw.dat");
    //Delta_large.savetxt("Delta_full.dat");
    __num_format<RealType>(Solver.w_0).savetxt("w_0.dat");
    __num_format<RealType>(Solver.w_1).savetxt("w_1.dat");

    FMatsubaraGrid stat_grid(-128, 128, beta);
    #ifdef _calc_extra_stats
    switch (D){
        case 1: getExtraData(SC_DF,FMatsubaraGrid(-256, 256, beta)); break;
        case 2: getExtraData(SC_DF,FMatsubaraGrid(-512, 512, beta)); break;
        case 3: getExtraData(SC_DF,FMatsubaraGrid(-128, 128, beta)); break;
        case 4: getExtraData(SC_DF,gridF); break; 
    }
    #endif
}



#ifdef _calc_extra_stats
template <class SCType> void getExtraData(SCType& SC, const FMatsubaraGrid& gridF)
{
    constexpr size_t D = SCType::_D;
    INFO("Calculating additional statistics.");
    const auto &Solver = SC._S;
    RealType beta = Solver.beta;
    RealType T=1.0/beta;
    SC.GLatLoc.savetxt("gloc.dat");

    std::bitset<10> flags(extraops);

    if (flags[0]) {
        std::array<RealType, D> q;
        q.fill(PI);
        size_t n_freq;
        ComplexType stat_susc_pi;
        switch(D)
        {
            case 1: 
                n_freq = std::max(int(beta*2), 1024); 
                stat_susc_pi = SC.getStaticLatticeSusceptibility(q, FMatsubaraGrid(-n_freq,n_freq,beta));
                break;
            case 2: 
                stat_susc_pi = SC.getStaticLatticeSusceptibility(q, Solver.w_grid);
                //n_freq = std::max(int(beta*2), 512); 
                //stat_susc_pi = SC.getStaticLatticeSusceptibility(q, FMatsubaraGrid(-n_freq,n_freq,beta));
                break;
            case 3: 
                stat_susc_pi = SC.getStaticLatticeSusceptibility(q, Solver.w_grid);
                break;
            case 4: 
                stat_susc_pi = SC.getStaticLatticeSusceptibility(q, Solver.w_grid);
                break;
        };
        __num_format<ComplexType>(stat_susc_pi).savetxt("StaticChiDFCC_pi.dat");
    };

    if (flags[1]) {
        INFO ("\nCalculating static susceptibility along (pi,0)->(pi,2*pi) direction");
        size_t n_freq = 256; 
        FMatsubaraGrid local_grid(-n_freq,n_freq,beta);
        GridObject<ComplexType, FMatsubaraGrid> Lambda(local_grid);
        Lambda.copyInterpolate(Solver.getLambda());
        auto glat = SC.getGLat(); 
        GridObject<RealType,KMesh> out_b(SC._kGrid), out_chi(SC._kGrid), out_dual_bubble(SC._kGrid), out_lattice_bubble(SC._kGrid), out_dual_bubble0(SC._kGrid);
        std::array<KMesh::point, D> q;
        q.fill(SC._kGrid.findClosest(PI));
        for (auto q_pt : SC._kGrid.getPoints()) {
            INFO2("q = " << RealType(q_pt));
            std::get<D-1>(q) = q_pt;
            auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);

            ComplexType susc_val = SC.getStaticLatticeSusceptibility(q,Solver.w_grid);
            out_chi[size_t(q_pt)]= std::real(susc_val);
            
            auto dual_bubble = Diagrams::getBubble(SC.GD, Wq_args_static);
            auto dual_bubble0 = Diagrams::getBubble(SC.GD0, Wq_args_static);
            auto Bw1 = beta*Solver.w_0*Solver.w_1*Solver.U*Solver.U*Solver.getLambda()*Solver.getLambda()*dual_bubble;
            auto Bw = Bw1/(1.0+Bw1);
            ComplexType B = Bw.sum();
            out_b[size_t(q_pt)]= std::real(B);
            
            RealType dual_bubble_val = std::real(dual_bubble.sum());
            out_dual_bubble[size_t(q_pt)]= dual_bubble_val;
            out_dual_bubble0[size_t(q_pt)]= std::real(dual_bubble0.sum());

            RealType lattice_bubble_val =  std::real(Diagrams::getBubble(glat, Wq_args_static).sum());
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
        std::array<RealType, D> q;
        q.fill(PI);
        auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);
        auto dual_bubble_pi = Diagrams::getBubble(SC.GD, Wq_args_static);
        auto dual_bubble0_pi = Diagrams::getBubble(SC.GD0, Wq_args_static);
        auto Bw1 = beta*Solver.w_0*Solver.w_1*Solver.U*Solver.U*Solver.getLambda()*Solver.getLambda()*dual_bubble_pi;
        auto Bw = Bw1/(1.0+Bw1);
        ComplexType B = Bw.sum();
        dual_bubble_pi.savetxt("dual_bubble_pi_w.dat");
        dual_bubble0_pi.savetxt("dual_bubble0_pi_w.dat");
        Bw1.savetxt("BwNominator_pi.dat");
        Bw.savetxt("Bw_pi.dat");
        INFO("B(pi) = " << B);
        __num_format<ComplexType>(B).savetxt("B_pi.dat");
    }

    if (flags[3]) {
        INFO2("Saving G(w,k=pi)");
        GF glat_pi(Solver.w_grid);
        auto glat = SC.getGLat(); 
        std::array<RealType, D> q;
        q.fill(PI);
        typename GF::PointFunctionType f = [&](FMatsubaraGrid::point w){return glat(std::tuple_cat(std::make_tuple(w),q));};
        glat_pi.fill(f);
        glat_pi.savetxt("glat_pi.dat");
    }

    if (flags[4]) {
        INFO2("Saving Green's functions - doing FFT");
        auto glat_k = SC.getGLat(); 
        typedef typename ArgBackGenerator<D,EnumerateGrid,GridObject,ComplexType,FMatsubaraGrid>::type glat_r_type;
        auto grid_r = EnumerateGrid(0,SC._kGrid.getSize(),false);
        glat_r_type glat_r(std::tuple_cat(std::forward_as_tuple(SC._fGrid),__repeater<EnumerateGrid,D>::get_tuple(grid_r)));
        for (auto w : SC._fGrid.getPoints()) { 
            glat_r[size_t(w)] = run_fft(glat_k[size_t(w)],FFTW_BACKWARD);
            }

        size_t distance = 5;
        GF glat_rp(Solver.w_grid);
        for (size_t i=0; i<distance; ++i) {
            std::array<typename EnumerateGrid::point, D> r_p; r_p.fill(grid_r[0]); r_p[D-1]=grid_r[i]; // Makes (i,0...) point in real space
            typename GF::PointFunctionType f = [&](FMatsubaraGrid::point w){return glat_r(std::tuple_cat(std::make_tuple(w),r_p));};
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
        size_t n_freq = Solver.gw.getGrid().getSize(); //std::max(int(beta*2), 256);
        auto stat_susc_bz = SC.getStaticLatticeSusceptibility(bzpoints, FMatsubaraGrid(-n_freq,n_freq,beta));
        std::map<BZPoint<D>, RealType> susc_map, lattice_bubble_map, dual_bubble_map;
        INFO_NONEWLINE("\tGetting bubbles [" << stat_susc_bz.size() <<"] : ");
        for (size_t nq=0; nq<stat_susc_bz.size(); ++nq) { 
            INFO_NONEWLINE(nq << " "); 
            susc_map[bzpoints[nq]] = std::real(stat_susc_bz[nq]);
            lattice_bubble_map[bzpoints[nq]] = std::real(std::real(Diagrams::getBubble(glat, std::tuple_cat(std::make_tuple(0.0), bzpoints[nq])).sum()));
            dual_bubble_map[bzpoints[nq]] = std::real(Diagrams::getBubble(SC.GD, std::tuple_cat(std::make_tuple(0.0), bzpoints[nq])).sum());
        };
        INFO("");

        auto all_bz_points = CubicTraits<D>::getAllBZPoints(SC._kGrid);
        size_t nqpts = all_bz_points.size();

        typedef typename ArgBackGenerator<D,KMesh        ,GridObject,ComplexType>::type susc_k_type;
        typedef typename ArgBackGenerator<D,EnumerateGrid,GridObject,ComplexType>::type susc_r_type;

        auto kmeshes = __repeater<KMesh,D>::get_tuple(SC._kGrid);
        auto rgrid = EnumerateGrid(0,SC._kGrid.getSize(),false);
        auto rmeshes = __repeater<EnumerateGrid,D>::get_tuple(rgrid);

        susc_k_type susc_full(kmeshes), susc_lattice_bubbles(kmeshes), dual_bubbles(kmeshes);
        susc_r_type susc_full_r(rmeshes), susc_lattice_bubbles_r(rmeshes), dual_bubbles_r(rmeshes);

        for (size_t nq=0; nq<nqpts; ++nq) {
            BZPoint<D> current_point = all_bz_points[nq];
            typename susc_k_type::PointTupleType current_point_tuple = std::tuple_cat(current_point);
            BZPoint<D> sym_point = CubicTraits<D>::findSymmetricBZPoint(current_point,SC._kGrid);
            susc_full.get(current_point_tuple) = susc_map[sym_point];
            susc_lattice_bubbles.get(current_point_tuple) = lattice_bubble_map[sym_point];
            dual_bubbles.get(current_point_tuple) = dual_bubble_map[sym_point];
            };

        susc_full.savetxt("StaticChiDFCC.dat");
        susc_lattice_bubbles.savetxt("StaticLatticeBubbleDFCC.dat");
        dual_bubbles.savetxt("StaticDualBubbleCC.dat");

        INFO2("Doing FFT of susceptibilities");

        susc_full_r.getData() = run_fft(susc_full.getData(),FFTW_BACKWARD);
        susc_lattice_bubbles_r.getData() = run_fft(susc_lattice_bubbles.getData(),FFTW_BACKWARD);
        dual_bubbles_r.getData() = run_fft(dual_bubbles.getData(),FFTW_BACKWARD);

        susc_full_r.savetxt("StaticChiDFCC_r.dat");
        susc_lattice_bubbles_r.savetxt("StaticLatticeBubbleDFCC_r.dat");
        dual_bubbles_r.savetxt("StaticDualBubbleCC_r.dat");

        INFO2("Saving susceptibilities in different directions");

        GridObject<RealType,EnumerateGrid> cc_x (rgrid), cc_y(rgrid), cc_xy(rgrid); 
        for (auto p1:rgrid.getPoints()) {
            auto zero_pt = rgrid.findClosest(0);
            std::array<typename EnumerateGrid::point, D> xpt,ypt,xypt; 
            xpt.fill(zero_pt); ypt.fill(zero_pt); xypt.fill(zero_pt);
            xpt[D-1]=p1; 
            xypt[D-1]=p1;
            ypt[std::max(int(D),2)-2]=p1;
            xypt[std::max(int(D),2)-2]=p1;

            cc_x.get(p1) = std::real(susc_full_r(typename susc_r_type::PointTupleType(std::tuple_cat(xpt))));
            cc_y.get(p1) = std::real(susc_full_r(typename susc_r_type::PointTupleType(std::tuple_cat(ypt))));
            cc_xy.get(p1) = std::real(susc_full_r(typename susc_r_type::PointTupleType(std::tuple_cat(xypt))));
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

        auto w = Solver.w_grid.findClosest(I*PI/beta);
        const auto &kgrid = SC._kGrid;
        auto rgrid = EnumerateGrid(0,SC._kGrid.getSize(),false);

        typedef typename ArgBackGenerator<D,KMesh        ,GridObject,ComplexType>::type gd_k_type;
        typedef typename ArgBackGenerator<D,EnumerateGrid,GridObject,ComplexType>::type gd_r_type;

        auto kmeshes = __repeater<KMesh,D>::get_tuple(kgrid);
        auto rmeshes = __repeater<EnumerateGrid,D>::get_tuple(rgrid);

        gd_k_type glat_k(kmeshes), gd_k(kmeshes), gd0_k(kmeshes), glat_dmft_k(kmeshes), sigmad_k(kmeshes);
        gd_r_type glat_r(rmeshes), gd_r(rmeshes), gd0_r(rmeshes), glat_dmft_r(rmeshes), sigmad_r(rmeshes);

        glat_k.getData() = glat[w.index_];
        gd_k.getData() = SC.GD[w.index_];
        gd0_k.getData() = SC.GD0[w.index_];
        glat_dmft_k.getData() = glatdmft[w.index_];
        sigmad_k.getData() = SC.SigmaD[w.index_];

        glat_k.savetxt("glat_k_w0.dat");
        gd_k.savetxt("gd_k_w0.dat");
        gd0_k.savetxt("gd0_k_w0.dat");
        glat_dmft_k.savetxt("glat_dmft_k_w0.dat");
        sigmad_k.savetxt("sigmad_k_w0.dat");

        glat_r.getData() = run_fft(glat_k.getData(), FFTW_BACKWARD);
        gd_r.getData() = run_fft(gd_k.getData(), FFTW_BACKWARD);
        gd0_r.getData() = run_fft(gd0_k.getData(), FFTW_BACKWARD);
        glat_dmft_r.getData() = run_fft(glat_dmft_k.getData(), FFTW_BACKWARD);
        sigmad_r.getData() = run_fft(sigmad_k.getData(), FFTW_BACKWARD);

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


