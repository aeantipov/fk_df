#include <numeric>

#include <MatsubaraGrid.hpp>
#include <KMesh.hpp>
#include <GridObject.hpp>
#include <EnumerateGrid.hpp>
#include "GFWrap.h"
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
#endif

#include "FKOptionParserDF.h"

#include <iostream>
#include <ctime>
#include <array>
#include <unordered_map>
#include <csignal>


using namespace GFTools;
using namespace FK;

//extern template GridObject<ComplexType,FMatsubaraGrid>;

RealType beta;
size_t extraops=0;
bool INTERRUPT = false;

typedef GFWrap GF;

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
    size_t n_freq = opt.n_freq;
    size_t ksize = opt.kpts;
    //size_t n_dual_freq = opt.n_dual_freq;
    RealType mix = opt.mix;
    extraops = opt.extraops;
    size_t NDMFTRuns = opt.NDMFTRuns;
    size_t NDFRuns = opt.NDFRuns;
    RealType DFCutoff = opt.DFCutoff;
    bool update_mixing = opt.update_mixing;

    KMesh kGrid(ksize);

    FMatsubaraGrid gridF(-n_freq, n_freq, beta);
    FMatsubaraGrid gridF_half(0, 2*n_freq, beta);
    int __msize = std::max(n_freq*5,size_t(1024));
    FMatsubaraGrid gridF_large(-__msize, __msize, beta);

    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1, f2;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};

    try { GF Delta2("Delta_full.dat", beta); Delta = Delta2;} 
    catch (std::exception &e) { Delta.fill(f1); };
    Delta.savetxt("Delta_0.dat");
    
    FKImpuritySolver Solver(U,mu,e_d,Delta);
    
    KMeshPatch qGrid(kGrid);

    dmft_sc_type SC_DMFT(Solver, t, kGrid);
    df_sc_type SC_DF(Solver, gridF, kGrid, t);

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
            GF Delta_large(gridF_large); Delta_large = Solver.Delta;
            Delta_large.savetxt("DeltaDMFT.dat");
            Solver.Sigma.savetxt("SigmaDMFT.dat"); 
            Solver.gw.savetxt("GwDMFT.dat");
            __num_format<RealType>(Solver.w_0).savetxt("w_0DMFT.dat");
            __num_format<RealType>(Solver.w_1).savetxt("w_1DMFT.dat");
            diff = 1.0; calc_DMFT = false; }; // now continue with DF 
        }
   
    GF Delta_half(gridF_half); Delta_half = Delta;
    GF gw_half(gridF_half); gw_half = Solver.gw;
    gw_half.savetxt("Gw.dat");
    //Delta_half.savetxt("Delta.dat");
    GF Delta_large(gridF_large); Delta_large = Delta;
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
    constexpr size_t D = SCType::NDim;
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
                n_freq = std::max(int(beta*2), 512); 
                stat_susc_pi = SC.getStaticLatticeSusceptibility(q, FMatsubaraGrid(-n_freq,n_freq,beta));
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
        GridObject<RealType,KMesh> out_b(SC._kGrid), out_chi(SC._kGrid), out_dual_bubble(SC._kGrid), out_lattice_bubble(SC._kGrid);
        std::array<KMesh::point, D> q;
        q.fill(SC._kGrid.findClosest(PI));
        for (auto q_pt : SC._kGrid.getPoints()) {
            INFO2("q = " << RealType(q_pt));
            std::get<D-1>(q) = q_pt;
            auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);

            ComplexType susc_val = SC.getStaticLatticeSusceptibility(q,Solver.w_grid);
            out_chi[size_t(q_pt)]= std::real(susc_val);
            
            auto dual_bubble = Diagrams::getBubble(SC.GD, Wq_args_static);
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
        std::array<RealType, D> q;
        q.fill(PI);
        auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);
        auto dual_bubble_pi = Diagrams::getBubble(SC.GD, Wq_args_static);
        auto Bw1 = beta*Solver.w_0*Solver.w_1*Solver.U*Solver.U*Solver.getLambda()*Solver.getLambda()*dual_bubble_pi;
        auto Bw = Bw1/(1.0+Bw1);
        ComplexType B = Bw.sum();
        dual_bubble_pi.savetxt("DualBubbleCC_pi.dat");
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
            std::array<typename decltype(grid_r)::point, D> r_p; r_p.fill(grid_r[0]); r_p[D-1]=grid_r[i]; // Makes (i,0...) point in real space
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
        size_t n_freq = std::max(int(beta*2), 256);
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
            typename susc_k_type::PointTupleType current_point_tuple = current_point;
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

            cc_x.get(p1) = std::real(susc_full_r(typename susc_r_type::PointTupleType(xpt)));
            cc_y.get(p1) = std::real(susc_full_r(typename susc_r_type::PointTupleType(ypt)));
            cc_xy.get(p1) = std::real(susc_full_r(typename susc_r_type::PointTupleType(xypt)));
            }

        cc_x.savetxt("StaticChiDFCC_dir_x.dat");
        cc_y.savetxt("StaticChiDFCC_dir_y.dat");
        cc_xy.savetxt("StaticChiDFCC_dir_xy.dat");
        };


    if (flags[6]) {
        INFO2("Saving dual self-energy");
        SC.SigmaD.savetxt("SigmaDwk.dat");
        };

    if (flags[7]){
        INFO2("Saving G(w,k)");
        auto glat = SC.getGLat(); 
        GF glat_k(Solver.w_grid);
        glat_k.savetxt("glat_k.dat");
        };

/*
    if (dynamic_flag) {

        if (dynamic_flag >= 2) {
            size_t n_b_freq = gridF._w_max/2; // std::max(gridF._w_max/2,int(beta));
             BMatsubaraGrid gridB(-n_b_freq, n_b_freq+1, beta);

             auto glat = SC.getGLat();
             auto gloc = SC.GLatLoc;
             size_t ksize = SC._kGrid.getSize();
             
             GF iw_gf(gridF); 
             iw_gf.fill([](ComplexType w){return w;});
             GridObject<ComplexType,BMatsubaraGrid> chi_q0(gridB), chi_qPI(gridB);


             KMeshPatch grid0PI(SC._kGrid, {{0, SC._kGrid.getSize()/2}} );
             std::array<KMeshPatch, D> qgrids = __repeater<KMeshPatch,D>::get_array(grid0PI); 
             auto data = SC.calculateLatticeData(gridB, qgrids); // Heavy operation

             auto LatticeSusc = std::get<0>(data);

             std::array<KMesh::point,SCType::NDim> q_0, q_PI;
             q_0.fill(SC._kGrid[0]);
             q_PI.fill(SC._kGrid[ksize/2]);

             for (auto iW : gridB.getPoints()) {
                 auto args_0 = std::tuple_cat(std::forward_as_tuple(iW),q_0);
                 auto args_pi = std::tuple_cat(std::forward_as_tuple(iW),q_PI);
                 chi_q0[size_t(iW)] = LatticeSusc(args_0);
                 chi_qPI[size_t(iW)] = LatticeSusc(args_pi);
                 };

             chi_q0.savetxt("Chiq0.dat");
             chi_qPI.savetxt("ChiqPI.dat");
             auto chi_q0_0 = -T*chi_q0.sum();
             auto chi_qPI_0 = -T*chi_qPI.sum();
             INFO("Chi0(q=0) sum  = " << chi_q0_0);
             INFO("Chi0(q=pi) sum = " << chi_qPI_0);
         };
        */
}
#endif

