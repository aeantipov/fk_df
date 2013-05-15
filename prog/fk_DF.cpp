#include <numeric>

#include <MatsubaraGrid.hpp>
#include <KMesh.hpp>
#include <GridObject.hpp>
#include "GFWrap.h"
#include "Solver.h"
#include "DF.h"

#include "FKOptionParserDF.h"

#include <iostream>
#include <ctime>
#include <array>
#include <unordered_map>
#include <csignal>

#include <fftw3.h>

using namespace GFTools;
using namespace FK;

//extern template GridObject<ComplexType,FMatsubaraGrid>;

RealType beta;
size_t D = 0;
size_t extraops=0;

typedef GFWrap GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}

bool INTERRUPT = false;
void sighandler(int signal)
{
    static size_t count = 0;
    count++;
    INFO("Caught INTERRUPT, signal " << signal <<", " << count << "/3 times. ")
    INTERRUPT = true;
    if (count >= 3) { INFO("Force exiting"); exit(signal); }
}

template <typename ValueType, size_t D, typename std::enable_if<D==2, bool>::type=0> 
Container<D,ValueType> run_fft (const Container<D,ValueType> &in)
{
    MatrixType<ComplexType> kdata = in.getAsMatrix();
    MatrixType<ComplexType> out(kdata);
    fftw_plan p;
    p = fftw_plan_dft_2d(kdata.rows(), kdata.cols(), reinterpret_cast<fftw_complex*>(kdata.data()), reinterpret_cast<fftw_complex*>(out.data()), FFTW_BACKWARD, FFTW_ESTIMATE); 
    fftw_execute(p);
    out/=(kdata.rows()*kdata.cols());
    return Container<2,ValueType>(out);
}

template <typename ValueType, size_t D, typename std::enable_if<D!=2, bool>::type=0> 
Container<D,ValueType> run_fft (const Container<D,ValueType> &in)
{
    ERROR("No FFT defined for D="<<D);
    return in; 
}


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
        auto bzpoints_map = CubicTraits<D>::getUniqueBZPoints(SC._kGrid); 
        std::vector<BZPoint<D>> bzpoints;
        for (auto map_it = bzpoints_map.begin(); map_it!=bzpoints_map.end(); map_it++){
            bzpoints.push_back(map_it->first);
            }
        size_t n_freq = std::max(int(beta*2), 256);
        auto stat_susc_bz = SC.getStaticLatticeSusceptibility(bzpoints, FMatsubaraGrid(-n_freq,n_freq,beta));
        std::map<BZPoint<D>, RealType> susc_map;
        for (size_t nq=0; nq<stat_susc_bz.size(); ++nq) susc_map[bzpoints[nq]] = std::real(stat_susc_bz[nq]);
        auto all_bz_points = CubicTraits<D>::getAllBZPoints(SC._kGrid);
        size_t nqpts = all_bz_points.size();
        size_t dimsize = SC._kGrid.getSize();
        std::ofstream out;
        out.open("StaticChiDFCC.dat");
        for (size_t nq=0; nq<nqpts; ++nq) {
            BZPoint<D> current_point = all_bz_points[nq];
            BZPoint<D> sym_point = CubicTraits<D>::findSymmetricBZPoint(current_point,SC._kGrid);
            out << all_bz_points[nq] << susc_map[sym_point] << std::endl;
            if ((nq+1)%dimsize==0) out << std::endl;
        }
        out.close();
    };

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

    if (flags[4] && D==2) {
        INFO2("Saving Green's functions - doing FFT");
        auto glat_k = SC.getGLat(); 
        typedef typename ArgBackGenerator<D,RealGrid,GridObject,ComplexType,FMatsubaraGrid>::type glat_r_type;
        auto grid_r = RealGrid(0,SC._kGrid.getSize(),SC._kGrid.getSize(),false);
        glat_r_type glat_r(std::tuple_cat(std::forward_as_tuple(SC._fGrid),__repeater<RealGrid,D>::get_tuple(grid_r)));
        for (auto w : SC._fGrid.getPoints()) { 
            glat_r[size_t(w)] = run_fft(glat_k[size_t(w)]);
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
    }
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
        std::cout << "Selfconsistency      : " << opt.sc_type << std::endl;
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
    size_t n_dual_freq = opt.n_dual_freq;
    RealType mix = opt.mix;
    auto sc_switch = opt.sc_index;
    extraops = opt.extraops;
    size_t NDMFTRuns = opt.NDMFTRuns;
    size_t NDFRuns = opt.NDFRuns;
    RealType DFCutoff = opt.DFCutoff;

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
    
    std::unique_ptr<SelfConsistency> SC_DF_ptr, SC_DMFT_ptr;
    std::unique_ptr<DFBase> DF_ptr;
    typedef FKOptionParserDF::SC enumSC;
    KMeshPatch qGrid(kGrid);
    switch (sc_switch) {
        case enumSC::DFCubic1d: 
            SC_DMFT_ptr.reset(new CubicDMFTSC<1>(Solver, t, kGrid));
            SC_DF_ptr.reset(new DFLadder<1>(Solver, gridF, kGrid, t)); 
            DF_ptr.reset(static_cast<DFLadder<1>*> (SC_DF_ptr.get()));
            D=1; break;
        case enumSC::DFCubic2d: 
            SC_DMFT_ptr.reset(new CubicDMFTSC<2>(Solver, t, kGrid));
            SC_DF_ptr.reset(new DFLadder<2>(Solver, gridF, kGrid, t)); 
            DF_ptr.reset(static_cast<DFLadder<2>*> (SC_DF_ptr.get()));
            D=2; break;
        case enumSC::DFCubic3d: 
            SC_DMFT_ptr.reset(new CubicDMFTSC<3>(Solver, t, kGrid));
            SC_DF_ptr.reset(new DFLadder<3>(Solver, gridF, kGrid, t)); 
            DF_ptr.reset(static_cast<DFLadder<3>*> (SC_DF_ptr.get()));
            D=3; break;
        case enumSC::DFCubic4d: 
            SC_DMFT_ptr.reset(new CubicDMFTSC<4>(Solver, t, kGrid));
            SC_DF_ptr.reset(new DFLadder<4>(Solver, gridF, kGrid, t)); 
            DF_ptr.reset(static_cast<DFLadder<4>*> (SC_DF_ptr.get()));
            D=4; break;
        default:
            ERROR("Couldn't find the desired SC type. Exiting.");
            exit(0);
    };
    DF_ptr->_n_GD_iter = opt.DFNumberOfSelfConsistentIterations;
    DF_ptr->_GDmix = opt.DFSCMixing;
    DF_ptr->_SC_cutoff = opt.DFSCCutoff;
    DF_ptr->_eval_BS_SC = opt.DFEvaluateBSSelfConsistent;
    DF_ptr->_n_BS_iter = opt.DFNumberOfBSIterations;
    DF_ptr->_BSmix = opt.DFBSMixing;
    DF_ptr->_EvaluateStaticDiagrams = opt.DFEvaluateStaticDiagrams;
    DF_ptr->_EvaluateDynamicDiagrams = opt.DFEvaluateDynamicDiagrams;

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

    auto &SC_DMFT = *SC_DMFT_ptr;
    auto &SC_DF   = *SC_DF_ptr;
  
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
            DF_ptr->getGLoc().savetxt(tmp.str());
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
    GF sigma_half(gridF_half); sigma_half = Solver.Sigma;
    sigma_half.savetxt("Sigma.dat");
    gw_half.savetxt("Gw.dat");
    Delta_half.savetxt("Delta.dat");
    GF Delta_large(gridF_large); Delta_large = Delta;
    Delta_large.savetxt("Delta_full.dat");
    __num_format<RealType>(Solver.w_0).savetxt("w_0.dat");
    __num_format<RealType>(Solver.w_1).savetxt("w_1.dat");

    if (extraops>0) {
        switch (sc_switch) {
            case enumSC::DFCubic1d: 
                getExtraData(*(static_cast<DFLadder<1>*> (SC_DF_ptr.get())), gridF); 
                break;
            case enumSC::DFCubic2d: 
                getExtraData(*(static_cast<DFLadder<2>*> (SC_DF_ptr.get())), gridF); 
                break;
            case enumSC::DFCubic3d: 
                getExtraData(*(static_cast<DFLadder<3>*> (SC_DF_ptr.get())), gridF); 
                break;
            case enumSC::DFCubic4d: 
                getExtraData(*(static_cast<DFLadder<4>*> (SC_DF_ptr.get())), gridF); 
                break;
            default: break;
            }; 
        };


    SC_DMFT_ptr.release(); 
    SC_DF_ptr.release();
    DF_ptr.release();
}

