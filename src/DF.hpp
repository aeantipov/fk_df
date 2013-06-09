#ifndef __FK_DF_HPP__
#define __FK_DF_HPP__

#include "DF.h"

// Some hints
// GLocalType - g(w)
// GKType - g(w,k...)
// GLocalType::FunctionType g(w) - analytic
// GKType::FunctionType g(w,k...) - analytic
// __repeater<KMesh,D>::get_tuple(_kGrid) - returns a tuple of D kmeshes

namespace FK {

// Some hints
// GLocalType - g(w)
// GKType - g(w,k...)
// GLocalType::FunctionType g(w) - analytic
// GKType::FunctionType g(w,k...) - analytic
// __repeater<KMesh,D>::get_tuple(_kGrid) - returns a tuple of D kmeshes

template <size_t D>
DFLadderCubic<D>::DFLadderCubic(const FKImpuritySolver &S, const FMatsubaraGrid& fGrid, KMesh kGrid, RealType t):
    CubicDMFTSC<D>(S,t,kGrid),
    _fGrid(fGrid),
    GD0(std::tuple_cat(std::make_tuple(_fGrid),__repeater<KMesh,D>::get_tuple(_kGrid))),
    GD(GD0.getGrids()),
    SigmaD(GD0.getGrids()), 
    GLat(GD0.getGrids()),
    GLatLoc(_fGrid)
{
    _initialize();
};


template <size_t D>
void DFLadderCubic<D>::_initialize()
{
    GLat = CubicDMFTSC<D>::getGLat(_fGrid);
    for (auto iw : _fGrid.getPoints()) {
        size_t iwn = size_t(iw);
        GD0[iwn] = GLat[iwn] - _S.gw(iw);
    };

    auto gd_f = [&](const wkArgTupleType &in)->ComplexType{
            ComplexType w = std::get<0>(in);
            auto ktuple = __tuple_tail(in);
            ComplexType e = _ek(ktuple);
            return (-e)/std::abs(w*w) -(e*e-_t*_t*2*RealType(D) - 2.0*e*(_S.mu - _S.Sigma._f(w)))/w/std::abs(w*w); 
            //return (-e)/std::abs(w*w) -(e*e-_t*_t*2*RealType(D) - 2.0*e*(_S.mu - _S.w_1*_S.U))/w/std::abs(w*w); 
            };
    GD0._f = __fun_traits<typename GKType::FunctionType>::getFromTupleF(gd_f);
    GD=GD0;
};

/*
template <size_t D>
template <typename ...KP>
ComplexType DFLadderCubic<D>::getBubble2(BMatsubaraGrid::point W, KP... q, FMatsubaraGrid::point w1) const
{
    static auto pts = std::make_tuple(W,q...);
    static auto GD_shift = this->GD0.shift(W,q...);
    auto pts2 = std::make_tuple(W,q...);
    if (pts2!=pts) {pts=pts2; GD_shift = this->GD0.shift(W,q...); };
    typename EkStorage::PointFunctionType GD0Function = [&](KP...k)->ComplexType { return GD0(w1,k...); };
    //typename EkStorage::PointFunctionType ShiftFunction = [&](KP...k)->ComplexType { return GD0(w1,k...); };
    typename EkStorage::PointFunctionType ShiftFunction = [&](KP...k)->ComplexType { return GD_shift(w1,k...); };
 
    EkStorage G1(_ek.getGrids()), G2(_ek.getGrids());
    G1.fill(GD0Function);
    G2.fill(ShiftFunction);
    G1*=G2;
    return G1.sum()/RealType(__power<ksize,D>::value)/(-this->_fGrid._beta);
}*/

template <size_t D>
inline typename DFLadderCubic<D>::GKType DFLadderCubic<D>::getGLat(const FMatsubaraGrid &gridF ) const
{
    GKType out(std::tuple_cat(std::forward_as_tuple(gridF), __repeater<KMesh,D>::get_tuple(_kGrid)));
    auto f1 = [&](const typename GKType::PointTupleType &in){return this->GLat(in);};
    auto f2 = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(f1);
    out.fill(f2);
    out._f = GLat._f;
    return out;
}

template <size_t D>
inline typename DFLadderCubic<D>::GKType DFLadderCubic<D>::getGLatDMFT(const FMatsubaraGrid& gridF) const 
{ 
    return CubicDMFTSC<D>::getGLat(gridF); 
};

template <size_t D>
inline typename DFLadderCubic<D>::GKType DFLadderCubic<D>::getGLat() const 
{ 
    return GLat; 
};

template <size_t D>
typename DFLadderCubic<D>::GLocalType DFLadderCubic<D>::operator()()
{
    INFO("Using DF Ladder self-consistency in " << D << " dimensions on a cubic lattice of " << _kGrid.getSize() << "^" << D <<" atoms.");
    SigmaD = 0.0;
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

    // Put here operations with GD

    // Generate a list of unique q-points
    //size_t ksize = _kGrid.getSize();
    //const auto all_q_pts = CubicTraits<D>::getAllBZPoints(_kGrid); 
    const auto unique_q_pts = CubicTraits<D>::getUniqueBZPoints(_kGrid);
    size_t totalqpts = size_t(pow(_kGrid.getSize(),D)); 
    RealType knorm = RealType(totalqpts);

    // Prepare static vertex
    GLocalType Lambda(_fGrid);
    Lambda.copyInterpolate(_S.getLambda());
    GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> StaticVertex4(std::forward_as_tuple(_fGrid,_fGrid)); 
    decltype(StaticVertex4)::PointFunctionType VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return _S.getVertex4(0.0, w1,w2);};
    StaticVertex4.fill(VertexF2);
    auto StaticV4 = StaticVertex4.getData().getAsMatrix();

    // Prepare dynamic vertex
    GLocalType DynVertex4(_fGrid), FullDualDynVertex4(_fGrid), DualDynBubble(_fGrid);
    // Miscelanneous - stream for checking convergence
    std::ofstream diffDF_stream("diffDF.dat",std::ios::out);
    diffDF_stream.close();

    INFO("Starting ladder dual fermion calculations")
    INFO("Beginning with GD sum = " << std::abs(GD.sum())/RealType(_fGrid.getSize())/knorm);

    RealType diffGD = 1.0, diffGD_min = 1.0; 
    size_t diffGD_min_count = 0; // A counter to estimate, how many iterations have passed after the minimal achieved diff
    for (size_t nd_iter=0; nd_iter<_n_GD_iter && diffGD > _SC_cutoff; ++nd_iter) { 
        INFO("DF iteration " << nd_iter << ". Evaluating BS equation.");
        GKType addSigma(this->SigmaD.getGrids()), GD_shift(this->GD.getGrids()), GD_sum(GD.getGrids());

        size_t nq = 1;
        // iterate over all unique kpoints (patch of BZ)
        for (auto pts_it = unique_q_pts.begin(); pts_it != unique_q_pts.end(); pts_it++) { 
            std::array<KMesh::point, D> q = pts_it->first; // point
            RealType q_weight = RealType(pts_it->second.size()); // it's weight
            auto other_pts = pts_it -> second; // other points, equivalent by symmetry
            INFO_NONEWLINE(nq << "/" << unique_q_pts.size() << ": [");
            for (size_t i=0; i<D; ++i) INFO_NONEWLINE(RealType(q[i]) << " "); INFO_NONEWLINE("]. Weight : " << q_weight << ". ");

            auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);
            GD_shift = GD.shift(Wq_args_static);
            GD_sum = 0.0;
            for (auto q_pt : other_pts) { // Sum over different G(w,k+q) to obey symmetry
                INFO_NONEWLINE("+");
                auto Wq_args1 = std::tuple_cat(std::make_tuple(0.0),q_pt);
                GD_sum += GD.shift(Wq_args1);
                };
            INFO("");

            //DEBUG(GD);
            //DEBUG(GD_shift);

            if (_EvaluateStaticDiagrams) {
                INFO_NONEWLINE("\tStatic contribution...");
                auto dual_bubble = Diagrams::getBubble(this->GD, Wq_args_static);
                auto mult = _S.beta*_S.U*_S.U*_S.w_0*_S.w_1;
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
                    size_t n_iter = 10;
                    auto dual_bubble_matrix = dual_bubble.getData().getAsDiagonalMatrix();
                    ERROR("DF iteration" << nd_iter << ". Evaluating BS equation using iterations.");
                    decltype(StaticV4) FullStaticV4 = Diagrams::BS(dual_bubble_matrix, StaticV4, true, true, n_iter, _BSmix).diagonal();
                    std::copy(FullStaticV4.data(), FullStaticV4.data()+FullStaticV4.size(), FullVertex11.getData()._data.data());
                    //DEBUG(FullVertex11);
                    //exit(1);
                    }
                typename GKType::PointFunctionType SigmaF;
                auto SigmaF2 = [&](wkPointTupleType in)->ComplexType { 
                    FMatsubaraGrid::point w = std::get<0>(in);
                    auto kpts = __tuple_tail(in); 
                    auto out_static = (1.0*T)*( FullVertex11(w)*GD_sum(in)); // Here the fact that StaticVertex4(w,w) = 0 is used.
                    return out_static;
                    };
                SigmaF = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(SigmaF2);
                addSigma.fill(SigmaF);
                INFO2("Sigma static contribution diff = " << addSigma.diff(SigmaD*0)/q_weight);
                SigmaD+=addSigma/RealType(totalqpts);


                /*
                auto dual_bubble_matrix = dual_bubble.getData().getAsDiagonalMatrix();
                decltype(StaticV4) FullStaticV4;
                FullStaticV4 = Diagrams::BS(dual_bubble_matrix, StaticV4, true, _eval_BS_SC, _n_BS_iter, _BSmix);
                typename GKType::PointFunctionType SigmaF;
                auto SigmaF2 = [&](wkPointTupleType in)->ComplexType { 
                    FMatsubaraGrid::point w = std::get<0>(in);
                    auto kpts = __tuple_tail(in); 
                    auto out_static = (1.0*T)*( (FullStaticV4(size_t(w),size_t(w)))*GD_shift(in)); // Here the fact that StaticVertex4(w,w) = 0 is used.
                    return out_static;
                    };
                SigmaF = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(SigmaF2);

                addSigma.fill(SigmaF);
                */

                };

            
            if (_EvaluateDynamicDiagrams) { 
                INFO_NONEWLINE("\tDynamic contribution...");
                decltype(StaticVertex4) DualBubbleDynamic(StaticVertex4.getGrids());
                decltype(StaticVertex4)::PointFunctionType dbfill = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){
                    return -T*(GD[size_t(w1)]*GD_shift[size_t(w2)]).sum()/RealType(totalqpts);
                    };
                decltype(StaticVertex4) DynamicFullVertex4(StaticVertex4.getGrids());
                DualBubbleDynamic.fill(dbfill);
                //DEBUG(DualBubbleDynamic);
                DynamicFullVertex4 = Diagrams::BS(DualBubbleDynamic, StaticVertex4*(-1.0), true, _eval_BS_SC,(_n_BS_iter>0?_n_BS_iter-1:0),_BSmix);
                //DEBUG(DualBubbleDynamic);
                DualBubbleDynamic*=(-1.0)*StaticVertex4*DynamicFullVertex4;
                //DEBUG(DualBubbleDynamic);
                addSigma=0.0;
                for (auto w : _fGrid.getPoints()) { 
                    for (auto w1 : _fGrid.getPoints()) { 
                        //addSigma[size_t(w)] += T*GD_shift[size_t(w1)]*(-1.0)*StaticVertex4(w,w1)*DualBubbleDynamic(w,w1)*DynamicFullVertex4(w,w1);
                        addSigma[size_t(w)] += T*GD_sum[size_t(w1)]*DualBubbleDynamic(w,w1);
                        //DEBUG(GD_shift[size_t(w1)]);
                        //DEBUG(DualBubbleDynamic(w,w1));
                        //DEBUG(T*GD_shift[size_t(w1)]*DualBubbleDynamic(w,w1));
                        //DEBUG(addSigma[size_t(w)]);
                        //addSigma[size_t(w)] += T*GD_shift[size_t(w1)]*(-_S.getVertex4(0.0, w,w1))/(1.0 + _S.getVertex4(0.0, w,w1)*DualBubbleDynamic(w,w1));
                            }
                        };
                INFO2("Sigma dynamic contribution diff = " << addSigma.diff(SigmaD*0)/q_weight);
                SigmaD+=addSigma/RealType(totalqpts);
                };

            nq++;
            }; // end of q loop

        //GLocalType SigmaD0(_fGrid),__GD0(_fGrid), __GDnew0(_fGrid);
        //std::function<ComplexType(FMatsubaraGrid::point)> sf = [&](FMatsubaraGrid::point w){return SigmaD[size_t(w)][0][0];};
        //SigmaD0.fill(sf);        
        //DEBUG(SigmaD0);

        auto GD_new = _GDmix/(1.0/GD0 - SigmaD) + GD*(1.0-_GDmix); // Dyson eq;
        diffGD = GD_new.diff(GD);
        if (diffGD<diffGD_min-_SC_cutoff/10.) { diffGD_min = diffGD; diffGD_min_count = 0; }
        else diffGD_min_count++;
        INFO2("DF diff = " << diffGD);
        if (diffGD_min_count > 7 && std::abs(_GDmix-0.05)>1e-3 && _update_mixing) {
            ERROR("\n\tCaught loop cycle. Reducing DF mixing to " << _GDmix/2 << " .\n");
            _GDmix=std::max(_GDmix/2., 0.05);
            GD_new = GD0;
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

       INFO2("GD sum = " << std::abs(GD.sum())/RealType(_fGrid.getSize())/knorm);
    };
    INFO("Finished DF iterations");
        
    SigmaD = 1.0/GD - 1.0/GD0;
    for (auto iw : _fGrid.getPoints()) {
        size_t iwn = size_t(iw);
        GLat[iwn] = 1.0/(Delta(iw) - _ek.getData()) + 1.0/(Delta(iw) - _ek.getData())/gw(iw)*GD[iwn]/gw(iw)/(Delta(iw) - _ek.getData());
        GDLoc[iwn] = GD[iwn].sum()/knorm; 
        GLatLoc[iwn] = GLat[iwn].sum()/knorm; 
        };
 
    // Finish - prepare all lattice quantities
    Delta_out = Delta + 1.0/gw * GDLoc / GLatLoc;
    // Assume DMFT asymptotics
    Delta_out._f = std::bind([&](ComplexType w)->ComplexType{return _t*_t*2*RealType(D)*((_S.mu-_S.w_1*_S.U)/std::abs(w*w) + 1.0/w);}, std::placeholders::_1);
    return Delta_out;
}

template <size_t D>
inline typename DFLadderCubic<D>::GLocalType DFLadderCubic<D>::getGLoc()
{
    return GLatLoc; 
}

template <size_t D>
std::tuple<typename DFLadderCubic<D>::SuscType> DFLadderCubic<D>::calculateLatticeData(const BMatsubaraGrid& gridB)
{
    KMeshPatch fullgrid(_kGrid);
    std::array<KMeshPatch, D> grids = __repeater<KMeshPatch,D>::get_array(fullgrid); 
    auto t1 = __repeater<KMeshPatch,D>::get_tuple(fullgrid); 
    return calculateLatticeData(gridB, grids);
}

template <size_t D>
std::tuple<typename DFLadderCubic<D>::SuscType> DFLadderCubic<D>::calculateLatticeData(const BMatsubaraGrid& gridB, const std::array<KMeshPatch, D>& qgrids)
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


template <size_t D>
template <typename KPoint>
std::vector<ComplexType> DFLadderCubic<D>::getStaticLatticeSusceptibility(const std::vector<std::array<KPoint, D>>& qpts, const FMatsubaraGrid& gridF)
{

    auto grids = std::tuple_cat(std::make_tuple(gridF),__repeater<KMesh,D>::get_tuple(_kGrid));

    // Prepare static vertex (outdated)
    //GridObject<ComplexType,FMatsubaraGrid,FMatsubaraGrid> FullStaticV4(std::forward_as_tuple(gridF,gridF)); 
    //decltype(StaticVertex4)::PointFunctionType VertexF2 = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){return _S.getVertex4(0.0, w1,w2);};
    //StaticVertex4.fill(VertexF2);
    //auto StaticV4 = StaticVertex4.getData().getAsMatrix();

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
        
    GLocalType Lambda(gridF);
    Lambda.copyInterpolate(_S.getLambda());
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

        auto mult = _S.beta*_S.U*_S.U*_S.w_0*_S.w_1;
        auto m1 = mult*dual_bubble*Lambda*Lambda;
        ComplexType B=(m1/(1.0+m1)).getData().sum();
        GLocalType B1=m1*Lambda/(1.0+m1);
    
        ComplexType susc1;
        for (auto w1 : gridF.getPoints()) {
            auto F = mult/(1.0+m1(w1))*Lambda(w1)/(1.0-B);
            for (auto w2 : gridF.getPoints()) {
                RealType kronecker = RealType(w1._index == w2._index);
                susc1+=GDL_bubble(w1)*F*(Lambda(w2)-Lambda(w1)*kronecker+B*Lambda(w1)*kronecker-B1(w2))*GDL_bubble(w2);
            }
        }

        //auto FullStaticV4 = Diagrams::BS(dual_bubble_matrix, StaticV4, true);
        //auto susc1 = (GDL_bubble_vector.transpose()*FullStaticV4*GDL_bubble_vector);
        //susc = susc1(0,0);
        
        susc = susc1;

        susc+=LatticeBubble.sum();
        INFO("Static susceptibility at q=" << q[0] << " = " << susc);
        out.push_back(susc);
        };

    return out;
}



template <size_t D>
template <typename KPoint>
ComplexType DFLadderCubic<D>::getStaticLatticeSusceptibility(const std::array<KPoint, D>& q, const FMatsubaraGrid& gridF)
{
    std::vector<std::array<KPoint,D>> qpts;
    qpts.push_back(q);
    return getStaticLatticeSusceptibility(qpts,gridF)[0];
}


template <size_t D>
template <typename KPoint>
ComplexType DFLadderCubic<D>::getStaticLatticeSusceptibility(const std::array<KPoint, D>& q)
{
    return getStaticLatticeSusceptibility(q,_fGrid);
}

} // end of namespace FK
#endif // endif :: #ifndef __FK_DF_H__
