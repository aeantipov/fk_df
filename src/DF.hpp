#ifndef __FK_DF_HPP__
#define __FK_DF_HPP__

#include "DF.h"
#include "FFT.hpp"

namespace FK {

// Some hints
// GLocalType - g(w)
// GKType - g(w,k...)
// GLocalType::FunctionType g(w) - analytic
// GKType::FunctionType g(w,k...) - analytic
// __repeater<KMesh,D>::get_tuple(_kGrid) - returns a tuple of D kmeshes

template <typename LatticeT, size_t D>
template <typename ...LatticeParams> 
DFLadder<LatticeT,D>::DFLadder(const FKImpuritySolver &S, const FMatsubaraGrid& fGrid, KMesh kGrid, LatticeParams ... lattice_p):
    LatticeDMFTSC<LatticeT>(S,kGrid,lattice_p...),
    _fGrid(fGrid),
    GD0(std::tuple_cat(std::make_tuple(_fGrid),__repeater<KMesh,D>::get_tuple(_kGrid))),
    GD(GD0.getGrids()),
    SigmaD(GD0.getGrids()), 
    GLat(GD0.getGrids()),
    GLatLoc(_fGrid)
{
    _initialize();
    SigmaD = 0.0;
};


template <typename LatticeT, size_t D>
void DFLadder<LatticeT,D>::_initialize()
{
    GLat = LatticeDMFTSC<LatticeT>::getGLat(_fGrid);
    for (auto iw : _fGrid.getPoints()) {
        size_t iwn = size_t(iw);
        GD0[iwn] = GLat[iwn] - _S.gw(iw);
    };

    auto gd_f = [&](const wkArgTupleType &in)->ComplexType{
            ComplexType w = std::get<0>(in);
            auto ktuple = __tuple_tail(in);
            ComplexType e = _ek(ktuple);
            return (-e)/std::abs(w*w) -(e*e-lattice.disp_square_sum() - 2.0*e*(_S.mu - _S.Sigma._f(w)))/w/std::abs(w*w); 
            //return (-e)/std::abs(w*w) -(e*e-_t*_t*2*RealType(D) - 2.0*e*(_S.mu - _S.w_1*_S.U))/w/std::abs(w*w); 
            };
    GD0._f = __fun_traits<typename GKType::FunctionType>::getFromTupleF(gd_f);
    GD=1.0/(1.0/GD0 - SigmaD);
    GD._f = GD0._f;
};

template <typename LatticeT, size_t D>
inline typename DFLadder<LatticeT,D>::GKType DFLadder<LatticeT,D>::getGLat(const FMatsubaraGrid &gridF ) const
{
    GKType out(std::tuple_cat(std::forward_as_tuple(gridF), __repeater<KMesh,D>::get_tuple(_kGrid)));
    auto f1 = [&](const typename GKType::PointTupleType &in){return this->GLat(in);};
    auto f2 = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(f1);
    out.fill(f2);
    out._f = GLat._f;
    return out;
}

template <typename LatticeT, size_t D>
inline typename DFLadder<LatticeT,D>::GKType DFLadder<LatticeT,D>::getGLatDMFT(const FMatsubaraGrid& gridF) const 
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
    INFO("Using DF Ladder self-consistency in " << D << " dimensions on a cubic lattice of " << _kGrid.getSize() << "^" << D <<" atoms.");
    RealType beta = _fGrid._beta;
    RealType T = 1.0/beta;
    GLocalType gw(_fGrid); // non-const method. Better copy.
    gw = _S.gw;
    GLocalType Delta(_fGrid); // non-const method. Better copy. 
    GLocalType GDLoc(_fGrid); 
    Delta = _S.Delta;
    GLocalType Delta_out(_fGrid); Delta_out=0.0;
    auto wkgrids = std::tuple_cat(std::make_tuple(_fGrid),__repeater<KMesh,D>::get_tuple(_kGrid));
    //typedef typename ArgBackGenerator<D,KMeshPatch,GridObject,ComplexType,FMatsubaraGrid>::type SigmaPatchType;
    //auto reduced_grids = std::tuple_cat(FMatsubaraGrid(0,_fGrid._max,bea), __repeater<KMeshPatch,D>::get_tuple(
    _initialize();
    GKType GD_initial(GD); // copy dressed GD if non-zero Sigma was provided in the beginning

    // Put here operations with GD
    //DEBUG("GD0 = " << GD0);
    //DEBUG("GD  = " << GD);
    //DEBUG("SigmaD = " << SigmaD);
    SigmaD = 0.0;

    // Generate a list of unique q-points
    //size_t ksize = _kGrid.getSize();
    //const auto all_q_pts = lattice_traits::getAllBZPoints(_kGrid); 
    const auto unique_q_pts = lattice_traits::getUniqueBZPoints(_kGrid);
    size_t totalqpts = size_t(pow(_kGrid.getSize(),D)); 
    RealType knorm = RealType(totalqpts);

    // Prepare static vertex
    auto mult = _S.beta*_S.U*_S.U*_S.w_0*_S.w_1;
    GLocalType Lambda(_fGrid);
    Lambda.copyInterpolate(_S.getLambda());
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> StaticVertex4(std::forward_as_tuple(_fGrid,_fGrid)); 
    decltype(StaticVertex4)::PointFunctionType VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return _S.getVertex4(0.0, w1,w2);};

    #ifdef vertex_U_inf
        auto U = _S.U;
        typename GLocalType::FunctionType lambdaf = [mult,U](ComplexType w){return 1. - U*U/4./w/w;};
        Lambda.fill(lambdaf);
        VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return mult*Lambda(w1)*Lambda(w2) - RealType(w1._index == w2._index)*Lambda(w1)*Lambda(w1);};
    #elif vertex_Hubbard
        INFO("Using Hubbard (singlet) vertex");
        auto U = _S.U;
        typename GLocalType::FunctionType lambdaf = [mult,U](ComplexType w){return 1. - U*U/4./w/w;};
        Lambda.fill(lambdaf);
        VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2)->ComplexType{
            return  mult*Lambda(w1)*Lambda(w2)*(2. + RealType(w1._index == w2._index));
        };
    GLocalType triplet_vertex(_fGrid);
    triplet_vertex = mult*Lambda*Lambda*(-3.);
    #endif
    StaticVertex4.fill(VertexF2);
    auto StaticV4 = StaticVertex4.getData().getAsMatrix();

    GKType FullStaticVertex(wkgrids);
    
    // Prepare dynamic vertex
   
    std::unique_ptr<FullVertexType> full_dyn_vertex; 
    if (_EvaluateDynamicDiagrams) { 
        INFO2("Allocating memory for dynamic vertex");
        full_dyn_vertex.reset(new FullVertexType(std::tuple_cat(std::make_tuple(_fGrid),std::make_tuple(_fGrid),__repeater<KMesh,D>::get_tuple(_kGrid))));
    }
    // Miscelanneous - stream for checking convergence
    std::ofstream diffDF_stream("diffDF.dat",std::ios::out);
    diffDF_stream.close();

    INFO("Starting ladder dual fermion calculations")
    GLocalType GDsum(_fGrid);
    for (auto iw : _fGrid.getPoints()) { GDsum[iw._index] = std::abs(GD[iw._index].sum())/knorm; }; 
    INFO("Beginning with GD sum = " << std::abs(GDsum.sum())/RealType(_fGrid.getSize()));


    RealType diffGD = 1.0, diffGD_min = 1.0; 
    size_t diffGD_min_count = 0; // A counter to estimate, how many iterations have passed after the minimal achieved diff
    for (size_t nd_iter=0; nd_iter<_n_GD_iter && diffGD > _SC_cutoff; ++nd_iter) { 
        INFO("DF iteration " << nd_iter << ". Evaluating BS equation.");

        size_t nq = 1;
        // iterate over all unique kpoints (patch of BZ)
        for (auto pts_it = unique_q_pts.begin(); pts_it != unique_q_pts.end(); pts_it++) { 
            std::array<KMesh::point, D> q = pts_it->first; // point
            RealType q_weight = RealType(pts_it->second.size()); // it's weight
            auto other_pts = pts_it -> second; // other points, equivalent by symmetry
            INFO_NONEWLINE(nq << "/" << unique_q_pts.size() << ": [");
            for (size_t i=0; i<D; ++i) INFO_NONEWLINE(RealType(q[i]) << " "); INFO_NONEWLINE("]. Weight : " << q_weight << ". ");

            auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);
            auto GD_shift = GD.shift(Wq_args_static);
            INFO("");

            if (_EvaluateStaticDiagrams) {
                INFO_NONEWLINE("\tStatic contribution...");
                auto dual_bubble = Diagrams::getBubble(this->GD, GD_shift);
                #ifdef bs_matrix
                auto dual_bubble_matrix = dual_bubble.getData().getAsDiagonalMatrix();
                auto FullStaticV4_m = Diagrams::BS(dual_bubble_matrix, StaticV4, true, _eval_BS_SC, _n_BS_iter, _BSmix);
                auto FullStaticV4 = FullStaticV4_m.diagonal();
                auto triplet_bs = Diagrams::BS(dual_bubble, triplet_vertex, true, _eval_BS_SC, _n_BS_iter, _BSmix);
                for (FMatsubaraGrid::point iw1 : _fGrid.getPoints())  {
                    auto f_val = FullStaticV4(iw1._index);//, iw1._index);
                    auto t_val = triplet_bs(iw1);//, iw1._index);
                    for (auto q_pt : other_pts) { 
                        FullStaticVertex.get(std::tuple_cat(std::make_tuple(iw1),q_pt)) = 0.5*(f_val+t_val);
                        };
                    };
                #else
                auto m1 = mult*dual_bubble*Lambda*Lambda;
                ComplexType B_=(m1/(1.0+m1)).getData().sum();
                if (std::imag(B_)>1e-5) throw (exRuntimeError("B is imaginary."));
                RealType B = std::real(B_);
                INFO("\t\tB = "<<B);
                GLocalType B1=m1*Lambda/(1.0+m1);
                GLocalType FullVertex11(_fGrid);
                if (B < 1.0 && (!_eval_BS_SC)) {
                    FullVertex11 = mult*Lambda/(1.0+m1)*(Lambda*B - B1)/(1.0-B); // Diagonal part of vertex 
                    }
                else {
                    //size_t n_iter = 10;
                    auto dual_bubble_matrix = dual_bubble.getData().getAsDiagonalMatrix();
                    ERROR("DF iteration" << nd_iter << ". Evaluating BS equation using iterations.");
                    decltype(StaticV4) FullStaticV4 = Diagrams::BS(dual_bubble_matrix, StaticV4, true, true, _n_BS_iter, _BSmix).diagonal();
                    std::copy(FullStaticV4.data(), FullStaticV4.data()+FullStaticV4.size(), FullVertex11.getData()._data.data());
                    }
                for (FMatsubaraGrid::point iw1 : _fGrid.getPoints())  {
                    auto f_val = FullVertex11(iw1);
                    for (auto q_pt : other_pts) { 
                        FullStaticVertex.get(std::tuple_cat(std::make_tuple(iw1),q_pt)) = f_val;
                        };
                    };
                #endif
                };

            
            if (_EvaluateDynamicDiagrams) { 
                INFO_NONEWLINE("\tDynamic contribution...");
                decltype(StaticVertex4) DualBubbleDynamic(StaticVertex4.getGrids());
                decltype(StaticVertex4)::PointFunctionType dbfill = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){
                    return -T*(GD[size_t(w1)]*GD_shift[size_t(w2)]).sum()/RealType(totalqpts);
                    };
                decltype(StaticVertex4) DynamicFullVertex4(StaticVertex4.getGrids());
                DualBubbleDynamic.fill(dbfill);
                DynamicFullVertex4 = Diagrams::BS(DualBubbleDynamic, StaticVertex4*(-1.0), true, _eval_BS_SC,(_n_BS_iter>0?_n_BS_iter-1:0),_BSmix);
                DynamicFullVertex4*=(-1.0)*StaticVertex4*DualBubbleDynamic;
                for (FMatsubaraGrid::point iw1 : _fGrid.getPoints())  {
                    for (FMatsubaraGrid::point iw2 : _fGrid.getPoints())  {
                        auto f_val = DynamicFullVertex4(iw1,iw2);
                        for (auto q_pt : other_pts) { 
                            full_dyn_vertex->get(std::tuple_cat(std::make_tuple(iw1),std::make_tuple(iw2),q_pt)) = f_val;
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
            for (FMatsubaraGrid::point iw1 : _fGrid.getPoints())  {
                auto v4r = run_fft(FullStaticVertex[iw1._index], FFTW_FORWARD)/knorm;
                auto gdr = run_fft(GD[iw1._index], FFTW_BACKWARD);
                SigmaD[iw1._index]+= (1.0*T)*run_fft(v4r*gdr, FFTW_FORWARD); 
                };
        };
        if (_EvaluateDynamicDiagrams) {
            INFO2("Dynamic diagrams ... Running FFT");
            for (FMatsubaraGrid::point iw1 : _fGrid.getPoints())  {
                for (FMatsubaraGrid::point iw2 : _fGrid.getPoints())  {
                    auto v4r = run_fft((*full_dyn_vertex)[iw1._index][iw2._index], FFTW_FORWARD)/knorm;
                    auto gdr = run_fft(GD[iw2._index], FFTW_BACKWARD);
                    SigmaD[iw1._index]+= (1.0*T)*run_fft(v4r*gdr, FFTW_FORWARD); 
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
        GD._f = GD0._f; // assume DMFT asymptotics are good 
        SigmaD = 0.0;

        /*if (diffGD<= _SC_cutoff && _GDmix < 1.0) {
        ERROR("\n\tRestoring back DF mix");
        _SC_cutoff = diffGD;
        _GDmix=std::min(1.0,_GDmix*1.1);
        diffGD_min_count = 0;
        diffGD_min = diffGD*1.5;
        diffGD=1.0;
        }*/

       for (auto iw : _fGrid.getPoints()) { GDsum[iw._index] = std::abs(GD[iw._index].sum())/knorm; }; 
       INFO2("GD sum = " << std::abs(GDsum.sum())/RealType(_fGrid.getSize()));
    };
    INFO("Finished DF iterations");
        
    SigmaD = 1.0/GD0 - 1.0/GD;
    for (auto iw : _fGrid.getPoints()) {
        size_t iwn = size_t(iw);
        GLat[iwn] = 1.0/(Delta(iw) - _ek.getData()) + 1.0/(Delta(iw) - _ek.getData())/gw(iw)*GD[iwn]/gw(iw)/(Delta(iw) - _ek.getData());
        GDLoc[iwn] = GD[iwn].sum()/knorm; 
        GLatLoc[iwn] = GLat[iwn].sum()/knorm; 
        };

    full_dyn_vertex.release(); 
    // Finish - prepare all lattice quantities
    Delta_out = Delta + 1.0/gw * GDLoc / GLatLoc;
    // Assume DMFT asymptotics
    Delta_out._f = std::bind([&](ComplexType w)->ComplexType{return lattice.disp_square_sum()*((_S.mu-_S.w_1*_S.U)/std::abs(w*w) + 1.0/w);}, std::placeholders::_1);
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
std::tuple<typename DFLadder<LatticeT,D>::SuscType> DFLadder<LatticeT,D>::calculateLatticeData(const BMatsubaraGrid& gridB)
{
    KMeshPatch fullgrid(_kGrid);
    std::array<KMeshPatch, D> grids = __repeater<KMeshPatch,D>::get_array(fullgrid); 
    auto t1 = __repeater<KMeshPatch,D>::get_tuple(fullgrid); 
    return calculateLatticeData(gridB, grids);
}

template <typename LatticeT, size_t D>
std::tuple<typename DFLadder<LatticeT,D>::SuscType> DFLadder<LatticeT,D>::calculateLatticeData(const BMatsubaraGrid& gridB, const std::array<KMeshPatch, D>& qgrids)
{
    SuscType LatticeSusc(std::tuple_cat(std::forward_as_tuple(gridB),qgrids));
    //RealType T = 1.0/_fGrid._beta;
    GKType Lwk = this->getGLatDMFT(_fGrid)/GD0*(-1.0);
    Lwk._f = __fun_traits<typename GKType::FunctionType>::constant(-1.0);
    auto GDL = GD*Lwk;
    GLocalType DynVertex4(_fGrid), Chi0(_fGrid), FullDualDynVertex4(_fGrid), ChiLat(_fGrid), GDL_chi0(_fGrid);
        
    size_t totalqpts = 1;
    for (const auto& qmesh : qgrids) { 
            totalqpts*=qmesh.getSize();
        }

    bool calc_static = false; 
    auto find_r = gridB.find(0.0);
    if (std::get<0>(find_r)) calc_static = true;
    BMatsubaraGrid::point zero_point(std::get<1>(find_r), 0.0);

    INFO2("Calculating lattice susceptibility");
    for (auto iW : gridB.getPoints()) {
        INFO2("iW = " << iW);
        typename GLocalType::PointFunctionType VertexFillf = 
            std::bind(&FKImpuritySolver::getBVertex4<BMatsubaraGrid::point, FMatsubaraGrid::point>, std::cref(_S), iW, std::placeholders::_1); 
        typename GLocalType::FunctionType Vertexf = 
            std::bind(&FKImpuritySolver::getBVertex4<ComplexType, ComplexType>, std::cref(_S), ComplexType(iW), std::placeholders::_1); 
        DynVertex4.fill(VertexFillf);
        DynVertex4._f = Vertexf;
 
        std::array<KMesh::point, D> q;
        for (size_t nq=0; nq<totalqpts; ++nq) { // iterate over all kpoints
            size_t offset = 0;
            for (size_t i=0; i<D; ++i) { 
                q[D-1-i]=qgrids[i][(nq-offset)/(int(pow(qgrids[i].getSize(),i)))%qgrids[i].getSize()]; 
                offset+=(int(pow(qgrids[i].getSize(),i)))*size_t(q[D-1-i]); 
                };
            WQTupleType Wq_args = std::tuple_cat(std::forward_as_tuple(iW),q);
            
            INFO_NONEWLINE("\t\t" << nq << "/" << totalqpts<< ". "); 
            size_t count=0; for (auto qval : q) { INFO_NONEWLINE("q_" << count++ << "="<<RealType(qval) << " "); }; 
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
std::vector<ComplexType> DFLadder<LatticeT,D>::getStaticLatticeSusceptibility(const std::vector<std::array<KPoint, D>>& qpts, const FMatsubaraGrid& gridF)
{

    auto grids = std::tuple_cat(std::make_tuple(gridF),__repeater<KMesh,D>::get_tuple(_kGrid));

    
    // Prepate interpolated Green's functions
    GKType GD0_interp (grids), GD_interp(grids), GLat_interp(grids);
    if (gridF._w_max != _fGrid._w_max || gridF._w_min != _fGrid._w_min) {
        GD0_interp.copyInterpolate(GD0);
        GD_interp.copyInterpolate(GD);
        GLat_interp.copyInterpolate(GLat);
        }
    else { 
        GD0_interp = GD0;
        GD_interp = GD;
        GLat_interp=GLat;
        };
        
    auto mult = _S.beta*_S.U*_S.U*_S.w_0*_S.w_1;
    GLocalType Lambda(gridF);
    Lambda.copyInterpolate(_S.getLambda());
    #ifdef bs_matrix
        GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> StaticVertex4(std::forward_as_tuple(gridF,gridF)); 
        decltype(StaticVertex4)::PointFunctionType VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return _S.getVertex4(0.0, w1,w2);};
    #endif
    #ifdef vertex_U_inf
        auto U = _S.U;
        typename GLocalType::FunctionType lambdaf = [mult,U](ComplexType w){return 1. - U*U/4./w/w;};
        Lambda.fill(lambdaf);
        VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return mult*Lambda(w1)*Lambda(w2) - RealType(w1._index == w2._index)*Lambda(w1)*Lambda(w1);};
    #elif vertex_Hubbard
        auto U = _S.U;
        typename GLocalType::FunctionType lambdaf = [mult,U](ComplexType w){return 1. - U*U/4./w/w;};
        Lambda.fill(lambdaf);
        VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2)->ComplexType{
            return  mult*Lambda(w1)*Lambda(w2)*(2. + RealType(w1._index == w2._index));
        };
    #endif
    #ifdef bs_matrix
        StaticVertex4.fill(VertexF2);
        auto StaticV4 = StaticVertex4.getData().getAsMatrix();
    #endif
    GKType Lwk = this->getGLatDMFT(gridF)/GD0_interp*(-1.0);
    Lwk._f = __fun_traits<typename GKType::FunctionType>::constant(-1.0);
    auto GDL = GD_interp*Lwk;

    // Prepare output
    size_t nqpts = qpts.size();
    std::vector<ComplexType> out;
    out.reserve(nqpts);

    for (auto q : qpts) {

        ComplexType susc=0.0;
        INFO_NONEWLINE("Evaluation of static susceptibility for q=["); for (int i=0; i<D; ++i) INFO_NONEWLINE(RealType(q[i])<<" "); INFO("]");
        auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);
        auto GD_shift = GD_interp.shift(Wq_args_static);
        auto LatticeBubble = Diagrams::getBubble(GLat_interp, Wq_args_static);

        auto GDL_bubble = Diagrams::getBubble(GDL, Wq_args_static);
        auto GDL_bubble_vector = GDL_bubble.getData().getAsVector();

        auto dual_bubble = Diagrams::getBubble(GD_interp, Wq_args_static);
        auto dual_bubble_matrix = dual_bubble.getData().getAsDiagonalMatrix();

        #ifdef bs_matrix
        auto size = StaticV4.rows();
        auto V4Chi = MatrixType<ComplexType>::Identity(size,size) - StaticV4*dual_bubble_matrix;
        auto D1 = V4Chi.determinant();
        if (std::imag(D1)<1e-7 && std::real(D1)>0) { 
            auto FullStaticV4 = Diagrams::BS(dual_bubble_matrix, StaticV4, true, false);
            susc = (GDL_bubble_vector.transpose()*FullStaticV4*GDL_bubble_vector)(0,0)*0.5;
            }
        else susc = -1.;
        #else
        //auto mult = _S.beta*_S.U*_S.U*_S.w_0*_S.w_1;
        auto m1 = mult*dual_bubble*Lambda*Lambda;
        ComplexType B=(m1/(1.0+m1)).getData().sum();
        GLocalType B1=m1*Lambda/(1.0+m1);
    
        for (auto w1 : gridF.getPoints()) {
            auto F = mult/(1.0+m1(w1))*Lambda(w1)/(1.0-B);
            for (auto w2 : gridF.getPoints()) {
                RealType kronecker = RealType(w1._index == w2._index);
                susc+=GDL_bubble(w1)*F*(Lambda(w2)-Lambda(w1)*kronecker+B*Lambda(w1)*kronecker-B1(w2))*GDL_bubble(w2);
            }
        }
        #endif

        #ifndef vertex_Hubbard
        susc+=LatticeBubble.sum();
        #endif
        INFO("Static susceptibility at q=" << q[0] << " = " << susc);
        out.push_back(susc);
        };

    return out;
}



template <typename LatticeT, size_t D>
template <typename KPoint>
ComplexType DFLadder<LatticeT,D>::getStaticLatticeSusceptibility(const std::array<KPoint, D>& q, const FMatsubaraGrid& gridF)
{
    std::vector<std::array<KPoint,D>> qpts;
    qpts.push_back(q);
    return getStaticLatticeSusceptibility(qpts,gridF)[0];
}


template <typename LatticeT, size_t D>
template <typename KPoint>
ComplexType DFLadder<LatticeT,D>::getStaticLatticeSusceptibility(const std::array<KPoint, D>& q)
{
    return getStaticLatticeSusceptibility(q,_fGrid);
}

} // end of namespace FK
#endif // endif :: #ifndef __FK_DF_H__
