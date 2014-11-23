#ifndef __FK_DF_HPP__
#define __FK_DF_HPP__

#include "DF.h"
#include "FFT.hpp"

namespace FK {

// Some hints
// GLocalType - g(w)
// GKType - g(w,k...)
// GLocalType::function_type g(w) - analytic
// GKType::function_type g(w,k...) - analytic
// tuple_tools::repeater<kmesh,D>::get_tuple(_kGrid) - returns a tuple of D kmeshes

template <typename LatticeT, size_t D>
template <typename ...LatticeParams> 
DFLadder<LatticeT,D>::DFLadder(const FKImpuritySolver &S, const fmatsubara_grid& fGrid, kmesh kGrid, LatticeParams ... lattice_p):
    LatticeDMFTSC<LatticeT>(S,kGrid,lattice_p...),
    _fGrid(fGrid),
    GD0(std::tuple_cat(std::make_tuple(_fGrid),tuple_tools::repeater<kmesh,D>::get_tuple(_kGrid))),
    GD(GD0.grids()),
    SigmaD(GD0.grids()), 
    GLat(GD0.grids()),
    GLatLoc(_fGrid)
{
    _initialize();
    SigmaD = 0.0;
};


template <typename LatticeT, size_t D>
void DFLadder<LatticeT,D>::_initialize()
{
    GLat = LatticeDMFTSC<LatticeT>::getGLat(_fGrid);
    for (auto iw : _fGrid.points()) {
        size_t iwn = size_t(iw);
        GD0[iwn] = GLat[iwn] - _S.gw(iw);
    };

    std::function<complex_type(wkarg_tuple)> gd_f = [&](const wkarg_tuple &in)->complex_type{
            complex_type w = std::get<0>(in);
            auto ktuple = tuple_tools::tuple_tail(in);
            complex_type e = _ek(ktuple);
            return (-e)/std::abs(w*w) -(e*e-lattice.disp_square_sum() - 2.0*e*(_S.mu - _S.Sigma.tail_eval(w)))/w/std::abs(w*w);
            //return (-e)/std::abs(w*w) -(e*e-_t*_t*2*real_type(D) - 2.0*e*(_S.mu - _S.w_1*_S.U))/w/std::abs(w*w); 
            };
    GD0.tail_ = tools::extract_tuple_f(gd_f);
    GD=1.0/(1.0/GD0 - SigmaD);
    GD.tail_ = GD0.tail_;
};

template <typename LatticeT, size_t D>
inline typename DFLadder<LatticeT,D>::GKType DFLadder<LatticeT,D>::getGLat(const fmatsubara_grid &gridF ) const
{
    GKType out(std::tuple_cat(std::forward_as_tuple(gridF), tuple_tools::repeater<kmesh,D>::get_tuple(_kGrid)));
    auto f1 = [&](const typename GKType::point_tuple &in){return this->GLat(in);};
    auto f2 = tools::fun_traits<typename GKType::point_function_type>::getFromTupleF(f1);
    out.fill(f2);
    out.tail_ = GLat.tail_;
    return out;
}

template <typename LatticeT, size_t D>
inline typename DFLadder<LatticeT,D>::GKType DFLadder<LatticeT,D>::getGLatDMFT(const fmatsubara_grid& gridF) const 
{ 
    return LatticeDMFTSC<LatticeT>::getGLat(gridF); 
};

template <typename LatticeT, size_t D>
inline typename DFLadder<LatticeT,D>::GKType DFLadder<LatticeT,D>::getGLat() const 
{ 
    return GLat; 
};

template <typename LatticeT, size_t D>
typename DFLadder<LatticeT,D>::GLocalType DFLadder<LatticeT,D>::operator()()
{
    INFO("Using DF Ladder self-consistency in " << D << " dimensions on a lattice of " << _kGrid.size() << "^" << D <<" atoms.");
    real_type beta = _fGrid.beta();
    real_type T = 1.0/beta;
    GLocalType gw(_fGrid); // non-const method. Better copy.
    gw = _S.gw;
    GLocalType Delta(_fGrid); // non-const method. Better copy. 
    GLocalType GDLoc(_fGrid); 
    Delta = _S.Delta;
    GLocalType Delta_out(_fGrid); Delta_out=0.0;
    auto wkgrids = std::tuple_cat(std::make_tuple(_fGrid),tuple_tools::repeater<kmesh,D>::get_tuple(_kGrid));
    //typedef typename tools::ArgBackGenerator<D,kmesh_patch,grid_object,complex_type,fmatsubara_grid>::type SigmaPatchType;
    //auto reduced_grids = std::tuple_cat(fmatsubara_grid(0,_fGrid._max,bea), tuple_tools::repeater<kmesh_patch,D>::get_tuple(
    _initialize();
    GKType GD_initial(GD); // copy dressed GD if non-zero Sigma was provided in the beginning

    // Put here operations with GD
    //DEBUG("GD0 = " << GD0);
    //DEBUG("GD  = " << GD);
    //DEBUG("SigmaD = " << SigmaD);
    SigmaD = 0.0;

    // Generate a list of unique q-points
    //size_t ksize = _kGrid.size();
    //const auto all_q_pts = lattice_traits::getAllBZPoints(_kGrid); 
    const auto unique_q_pts = lattice_traits::getUniqueBZPoints(_kGrid);
    size_t totalqpts = size_t(pow(_kGrid.size(),D)); 
    real_type knorm = real_type(totalqpts);

    // Prepare static vertex
    auto mult = _S.beta*_S.U*_S.U*_S.w_0*_S.w_1;
    GLocalType Lambda(_fGrid);
    Lambda.copy_interpolate(_S.getLambda());
    typedef grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> VertexType;
    VertexType StaticVertex4(std::forward_as_tuple(_fGrid,_fGrid)); 
    grid_object<complex_type,fmatsubara_grid,fmatsubara_grid>::point_function_type VertexF2 = [&](fmatsubara_grid::point w1, fmatsubara_grid::point w2){return _S.getVertex4(0.0, w1,w2);};
    StaticVertex4.fill(VertexF2);
    auto StaticV4 = StaticVertex4.data().as_matrix();

    GKType FullStaticVertex(wkgrids);
    
    // Prepare dynamic vertex
   
    std::unique_ptr<FullVertexType> full_dyn_vertex; 
    if (_EvaluateDynamicDiagrams) { 
        INFO2("Allocating memory for dynamic vertex");
        full_dyn_vertex.reset(new FullVertexType(std::tuple_cat(std::make_tuple(_fGrid,_fGrid),tuple_tools::repeater<kmesh,D>::get_tuple(_kGrid))));
    }
    // Miscelanneous - stream for checking convergence
    std::ofstream diffDF_stream("diffDF.dat",std::ios::out);
    diffDF_stream.close();

    INFO("Starting ladder dual fermion calculations")
    GLocalType GDsum(_fGrid);
    for (auto iw : _fGrid.points()) { GDsum[iw.index_] = std::abs(GD[iw.index_].sum())/knorm; }; 
    INFO("Beginning with GD sum = " << std::abs(GDsum.sum())/real_type(_fGrid.size()));

    real_type diffGD = 1.0, diffGD_min = 1.0; 
    size_t diffGD_min_count = 0; // A counter to estimate, how many iterations have passed after the minimal achieved diff
    for (size_t nd_iter=0; nd_iter<_n_GD_iter && diffGD > _SC_cutoff; ++nd_iter) { 
        INFO("DF iteration " << nd_iter << ". Evaluating BS equation.");

        size_t nq = 1;

        GKType bubbles = Diagrams::getStaticBubbles(this->GD); 
        // iterate over all unique kpoints (patch of BZ)
        for (auto pts_it = unique_q_pts.begin(); pts_it != unique_q_pts.end(); pts_it++) { 
            std::array<kmesh::point, D> q = pts_it->first; // point
            real_type q_weight = real_type(pts_it->second.size()); // it's weight
            auto other_pts = pts_it -> second; // other points, equivalent by symmetry
            INFO_NONEWLINE(nq << "/" << unique_q_pts.size() << ": [");
            for (size_t i=0; i<D; ++i) INFO_NONEWLINE(real_type(q[i]) << " "); INFO_NONEWLINE("]. Weight : " << q_weight << ". ");

            auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);
            auto dual_bubble = Diagrams::getBubble(this->GD, Wq_args_static);
            GLocalType db2 (_fGrid); 
            db2.fill([&](typename fmatsubara_grid::point w){return bubbles(std::tuple_cat(std::make_tuple(w), q)); });
            std::cout << db2 << std::endl << dual_bubble << std::endl;
            exit(0);
            INFO("");

            if (_EvaluateStaticDiagrams) {
                INFO_NONEWLINE("\tStatic contribution...");
                #ifdef bs_matrix
                auto dual_bubble_matrix = dual_bubble.data().as_diagonal_matrix();
                auto FullStaticV4_m = Diagrams::BS(dual_bubble_matrix, StaticV4, true, _eval_BS_SC, _n_BS_iter, _BSmix);
                auto FullStaticV4 = FullStaticV4_m.diagonal();
                for (fmatsubara_grid::point iw1 : _fGrid.points())  {
                    auto f_val = FullStaticV4(iw1.index_);//, iw1.index_);
                    for (auto q_pt : other_pts) { 
                        FullStaticVertex.get(std::tuple_cat(std::make_tuple(iw1),q_pt)) = f_val;
                        };
                    };
                #else
                auto m1 = mult*dual_bubble*Lambda*Lambda;
                complex_type B_=(m1/(1.0+m1)).data().sum();
                if (std::imag(B_)>1e-5) throw (exRuntimeError("B is imaginary."));
                real_type B = std::real(B_);
                INFO("\t\tB = "<<B);
                GLocalType B1=m1*Lambda/(1.0+m1);
                GLocalType FullVertex11(_fGrid);
                if (B < 1.0 && (!_eval_BS_SC)) {
                    FullVertex11 = mult*Lambda/(1.0+m1)*(Lambda*B - B1)/(1.0-B); // Diagonal part of vertex 
                    }
                else {
                    //size_t n_iter = 10;
                    auto dual_bubble_matrix = dual_bubble.data().as_diagonal_matrix();
                    ERROR("DF iteration" << nd_iter << ". Evaluating BS equation using iterations.");
                    decltype(StaticV4) FullStaticV4 = Diagrams::BS(dual_bubble_matrix, StaticV4, true, true, _n_BS_iter, _BSmix).diagonal();
                    std::copy(FullStaticV4.data(), FullStaticV4.data()+FullStaticV4.size(), FullVertex11.data().data());
                    }
                for (fmatsubara_grid::point iw1 : _fGrid.points())  {
                    auto f_val = FullVertex11(iw1);
                    for (auto q_pt : other_pts) { 
                        FullStaticVertex.get(std::tuple_cat(std::make_tuple(iw1),q_pt)) = f_val;
                        };
                    };
                #endif
                };

            
            if (_EvaluateDynamicDiagrams) { 
                INFO_NONEWLINE("\tDynamic contribution...");
                decltype(StaticVertex4) DualBubbleDynamic(StaticVertex4.grids());
                VertexType::point_function_type dbfill = [&](fmatsubara_grid::point w1, fmatsubara_grid::point w2){
                    throw(std::logic_error("no dynamic diagrams"));
                    return 1.0;
                    //return -T*(GD[size_t(w1)]*GD_shift[size_t(w2)]).sum()/real_type(totalqpts);
                    };
                decltype(StaticVertex4) DynamicFullVertex4(StaticVertex4.grids());
                DualBubbleDynamic.fill(dbfill);
                DynamicFullVertex4 = Diagrams::BS(DualBubbleDynamic, StaticVertex4*(-1.0), true, _eval_BS_SC,(_n_BS_iter>0?_n_BS_iter-1:0),_BSmix);
                DynamicFullVertex4*=(-1.0)*StaticVertex4*DualBubbleDynamic;
                for (fmatsubara_grid::point iw1 : _fGrid.points())  {
                    for (fmatsubara_grid::point iw2 : _fGrid.points())  {
                        auto f_val = DynamicFullVertex4(iw1,iw2);
                        for (auto q_pt : other_pts) { 
                            full_dyn_vertex->get(std::tuple_cat(std::make_tuple(iw1,iw2),q_pt)) = f_val;
                            };
                        };
                    };
                };
            nq++;
            }; // end of q loop

        // Finished obtaining vertex
        INFO("Updating dual self-energy with the static contribution...");
        if (_EvaluateStaticDiagrams) {
            INFO2("Static diagrams ... Running FFT");
            for (fmatsubara_grid::point iw1 : _fGrid.points())  {
                auto v4r = run_fft(FullStaticVertex[iw1.index_], FFTW_FORWARD)/knorm;
                auto gdr = run_fft(GD[iw1.index_], FFTW_BACKWARD);
                SigmaD[iw1.index_]+= (1.0*T)*run_fft(v4r*gdr, FFTW_FORWARD); 
                };
        };
        if (_EvaluateDynamicDiagrams) {
            INFO2("Dynamic diagrams ... Running FFT");
            for (fmatsubara_grid::point iw1 : _fGrid.points())  {
                for (fmatsubara_grid::point iw2 : _fGrid.points())  {
                    auto v4r = run_fft((*full_dyn_vertex)[iw1.index_][iw2.index_], FFTW_FORWARD)/knorm;
                    auto gdr = run_fft(GD[iw2.index_], FFTW_BACKWARD);
                    SigmaD[iw1.index_]+= (1.0*T)*run_fft(v4r*gdr, FFTW_FORWARD); 
                    };
                };
            };
        INFO("Total sigma diff = " << SigmaD.diff(SigmaD*0));

        auto GD_new = _GDmix/(1.0/GD0 - SigmaD) + GD*(1.0-_GDmix); // Dyson eq;
        diffGD = GD_new.diff(GD);
        if (diffGD<diffGD_min-_SC_cutoff/10.) { diffGD_min = diffGD; diffGD_min_count = 0; }
        else diffGD_min_count++;
        INFO2("DF diff = " << diffGD);
        if (diffGD_min_count > 12 && std::abs(_GDmix-0.05)>1e-3 && _update_mixing) {
            ERROR("\n\tCaught loop cycle. Reducing DF mixing to " << _GDmix/2 << " .\n");
            _GDmix=std::max(_GDmix/1.5, 0.05);
            GD_new = GD_initial;
            SigmaD = 0.0;
            diffGD_min = diffGD;
            diffGD_min_count = 0;
            }

        diffDF_stream.open("diffDF.dat",std::ios::app);
        diffDF_stream << diffGD << "  " << _GDmix << std::endl;
        diffDF_stream.close();
        GD=GD_new;
        GD.tail_ = GD0.tail_; // assume DMFT asymptotics are good
        SigmaD = 0.0;

       for (auto iw : _fGrid.points()) { GDsum[iw.index_] = std::abs(GD[iw.index_].sum())/knorm; }; 
       INFO2("GD sum = " << std::abs(GDsum.sum())/real_type(_fGrid.size()));
    };
    INFO("Finished DF iterations");
        
    SigmaD = 1.0/GD0 - 1.0/GD;
    for (auto iw : _fGrid.points()) {
        size_t iwn = size_t(iw);
        GLat[iwn] = 1.0/(Delta(iw) - _ek.data()) + 1.0/(Delta(iw) - _ek.data())/gw(iw)*GD[iwn]/gw(iw)/(Delta(iw) - _ek.data());
        GDLoc[iwn] = GD[iwn].sum()/knorm; 
        GLatLoc[iwn] = GLat[iwn].sum()/knorm; 
        };

    full_dyn_vertex.release(); 
    // Finish - prepare all lattice quantities
    Delta_out = Delta + 1.0/gw * GDLoc / GLatLoc;
    // Assume DMFT asymptotics
    Delta_out.tail_ = std::bind([&](complex_type w)->complex_type{return lattice.disp_square_sum()*((_S.mu-_S.w_1*_S.U)/std::abs(w*w) + 1.0/w);}, std::placeholders::_1);
    //DEBUG("GD0 = " << GD0);
    //DEBUG("GD  = " << GD);
    //DEBUG("SigmaD = " << SigmaD);

    return Delta_out;
}

template <typename LatticeT, size_t D>
inline typename DFLadder<LatticeT,D>::GLocalType DFLadder<LatticeT,D>::getGLoc()
{
    return GLatLoc; 
}

template <typename LatticeT, size_t D>
std::tuple<typename DFLadder<LatticeT,D>::SuscType> DFLadder<LatticeT,D>::calculateLatticeData(const bmatsubara_grid& gridB)
{
    kmesh_patch fullgrid(_kGrid);
    std::array<kmesh_patch, D> grids = tuple_tools::repeater<kmesh_patch,D>::get_array(fullgrid); 
    auto t1 = tuple_tools::repeater<kmesh_patch,D>::get_tuple(fullgrid); 
    return calculateLatticeData(gridB, grids);
}

template <typename LatticeT, size_t D>
std::tuple<typename DFLadder<LatticeT,D>::SuscType> DFLadder<LatticeT,D>::calculateLatticeData(const bmatsubara_grid& gridB, const std::array<kmesh_patch, D>& qgrids)
{
    SuscType LatticeSusc(std::tuple_cat(std::forward_as_tuple(gridB),qgrids));
    //real_type T = 1.0/_fGrid.beta();
    GKType Lwk = this->getGLatDMFT(_fGrid)/GD0*(-1.0);
    Lwk.tail_ = tools::fun_traits<typename GKType::function_type>::constant(-1.0);
    auto GDL = GD*Lwk;
    GLocalType DynVertex4(_fGrid), Chi0(_fGrid), FullDualDynVertex4(_fGrid), ChiLat(_fGrid), GDL_chi0(_fGrid);
        
    size_t totalqpts = 1;
    for (const auto& qmesh : qgrids) { 
            totalqpts*=qmesh.size();
        }

    bool calc_static = true;
    bmatsubara_grid::point zero_point = gridB.find_nearest(0.0);


    INFO2("Calculating lattice susceptibility");
    for (auto iW : gridB.points()) {
        INFO2("iW = " << iW);
        typename GLocalType::point_function_type VertexFillf = 
            std::bind(&FKImpuritySolver::getBVertex4<bmatsubara_grid::point, fmatsubara_grid::point>, std::cref(_S), iW, std::placeholders::_1); 
        typename GLocalType::function_type Vertexf = 
            std::bind(&FKImpuritySolver::getBVertex4<complex_type, complex_type>, std::cref(_S), complex_type(iW), std::placeholders::_1); 
        DynVertex4.fill(VertexFillf);
        DynVertex4.tail_ = Vertexf;
 
        std::array<kmesh::point, D> q;
        for (size_t nq=0; nq<totalqpts; ++nq) { // iterate over all kpoints
            size_t offset = 0;
            for (size_t i=0; i<D; ++i) { 
                q[D-1-i]=qgrids[i][(nq-offset)/(int(pow(qgrids[i].size(),i)))%qgrids[i].size()]; 
                offset+=(int(pow(qgrids[i].size(),i)))*size_t(q[D-1-i]); 
                };
            WQTupleType Wq_args = std::tuple_cat(std::forward_as_tuple(iW),q);
            
            INFO_NONEWLINE("\t\t" << nq << "/" << totalqpts<< ". "); 
            size_t count=0; for (auto qval : q) { INFO_NONEWLINE("q_" << count++ << "="<<real_type(qval) << " "); }; 
            // Bethe-Salpeter vertex
            Chi0 = Diagrams::getBubble(GD, Wq_args);
            FullDualDynVertex4 = Diagrams::BS(Chi0, DynVertex4, true, _eval_BS_SC, _n_BS_iter, _BSmix);
            
            // Lattice susceptibility
            GLocalType ChiLat0 = Diagrams::getBubble(GLat, Wq_args);
            GDL_chi0 = Diagrams::getBubble(GDL, Wq_args);
            ChiLat = GDL_chi0*FullDualDynVertex4*GDL_chi0 + ChiLat0;
            auto ChiVal = ChiLat.sum();
            LatticeSusc.get(Wq_args) = ChiVal;
            //if (calc_static) LatticeSusc.get(std::tuple_cat(std::forward_as_tuple(zero_point), q))-=ChiVal*T;
        } // end of q loop
    }; //end of iW loop
    return std::forward_as_tuple(LatticeSusc);
}


template <typename LatticeT, size_t D>
template <typename KPoint>
std::vector<complex_type> DFLadder<LatticeT,D>::getStaticLatticeSusceptibility(const std::vector<std::array<KPoint, D>>& qpts, const fmatsubara_grid& gridF)
{

    auto grids = std::tuple_cat(std::make_tuple(gridF),tuple_tools::repeater<kmesh,D>::get_tuple(_kGrid));

    // Prepare static vertex (outdated)
    //grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> FullStaticV4(std::forward_as_tuple(gridF,gridF)); 
    //decltype(StaticVertex4)::point_function_type VertexF2 = [&](fmatsubara_grid::point w1, fmatsubara_grid::point w2){return _S.getVertex4(0.0, w1,w2);};
    //StaticVertex4.fill(VertexF2);
    //auto StaticV4 = StaticVertex4.data().as_matrix();

    // Prepate interpolated Green's functions
    GKType GD0_interp (grids), GD_interp(grids), GLat_interp(grids);
    if (gridF.w_max_ != _fGrid.w_max_ || gridF.w_min_ != _fGrid.w_min_) {
        GD0_interp.copy_interpolate(GD0);
        GD_interp.copy_interpolate(GD);
        GLat_interp.copy_interpolate(GLat);
        }
    else { 
        GD0_interp = GD0;
        GD_interp = GD;
        GLat_interp=GLat;
        };

    auto mult = _S.beta*_S.U*_S.U*_S.w_0*_S.w_1;
    GLocalType Lambda(gridF);
    Lambda.copy_interpolate(_S.getLambda());
    #ifdef bs_matrix
        grid_object<complex_type,fmatsubara_grid,fmatsubara_grid> StaticVertex4(std::forward_as_tuple(gridF,gridF)); 
        decltype(StaticVertex4)::point_function_type VertexF2 = [&](fmatsubara_grid::point w1, fmatsubara_grid::point w2){return _S.getVertex4(0.0, w1,w2);};
    #endif
    // put here custom vertex functions
    #ifdef bs_matrix
        StaticVertex4.fill(VertexF2);
        auto StaticV4 = StaticVertex4.data().as_matrix();
    #endif
    
    GKType Lwk = this->getGLatDMFT(gridF)/GD0_interp*(-1.0);
    Lwk.tail_ = tools::fun_traits<typename GKType::function_type>::constant(-1.0);
    auto GDL = GD_interp*Lwk;

    // Prepare output
    size_t nqpts = qpts.size();
    std::vector<complex_type> out;
    out.reserve(nqpts);

    for (auto q : qpts) {

        complex_type susc=0.0;
        INFO_NONEWLINE("Evaluation of static susceptibility for q=["); for (int i=0; i<D; ++i) INFO_NONEWLINE(real_type(q[i])<<" "); INFO("]");
        auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);
        auto LatticeBubble = Diagrams::getBubble(GLat_interp, Wq_args_static);

        auto GDL_bubble = Diagrams::getBubble(GDL, Wq_args_static);

        auto dual_bubble = Diagrams::getBubble(GD_interp, Wq_args_static);

        #ifdef bs_matrix
        auto GDL_bubble_vector = GDL_bubble.data().as_vector();
        auto dual_bubble_matrix = dual_bubble.data().as_diagonal_matrix();
        auto size = StaticV4.rows();
        auto V4Chi = MatrixType<complex_type>::Identity(size,size) - StaticV4*dual_bubble_matrix;
        auto D1 = V4Chi.determinant();
        if (std::imag(D1)<1e-7 && std::real(D1)>0) { 
            auto FullStaticV4 = Diagrams::BS(dual_bubble_matrix, StaticV4, true, false);
            susc = (GDL_bubble_vector.transpose()*FullStaticV4*GDL_bubble_vector)(0,0);
            }
        else susc = -1.;
        #else
        auto m1 = mult*dual_bubble*Lambda*Lambda;
        complex_type B=(m1/(1.0+m1)).data().sum();
        GLocalType B1=m1*Lambda/(1.0+m1);

        for (auto w1 : gridF.points()) {
            auto F = mult/(1.0+m1(w1))*Lambda(w1)/(1.0-B);
            for (auto w2 : gridF.points()) {
                real_type kronecker = real_type(w1.index_ == w2.index_);
                susc+=GDL_bubble(w1)*F*(Lambda(w2)-Lambda(w1)*kronecker+B*Lambda(w1)*kronecker-B1(w2))*GDL_bubble(w2);
            }
        }
        #endif

        susc+=LatticeBubble.sum();
        INFO("Static susceptibility at q=" << q[0] << " = " << susc);
        out.push_back(susc);
        };

    return out;
}



template <typename LatticeT, size_t D>
template <typename KPoint>
complex_type DFLadder<LatticeT,D>::getStaticLatticeSusceptibility(const std::array<KPoint, D>& q, const fmatsubara_grid& gridF)
{
    std::vector<std::array<KPoint,D>> qpts;
    qpts.push_back(q);
    return getStaticLatticeSusceptibility(qpts,gridF)[0];
}


template <typename LatticeT, size_t D>
template <typename KPoint>
complex_type DFLadder<LatticeT,D>::getStaticLatticeSusceptibility(const std::array<KPoint, D>& q)
{
    return getStaticLatticeSusceptibility(q,_fGrid);
}

} // end of namespace FK
#endif // endif :: #ifndef __FK_DF_H__
