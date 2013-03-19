#include "DF.h"

namespace FK {

// Some hints
// GLocalType - g(w)
// GKType - g(w,k...)
// GLocalType::FunctionType g(w) - analytic
// GKType::FunctionType g(w,k...) - analytic
// __repeater<KMesh,D>::get_tuple(_kGrid) - returns a tuple of D kmeshes

template <size_t D>
DFLadder<D>::DFLadder(const FKImpuritySolver &S, const FMatsubaraGrid& fGrid, KMesh kGrid, const BMatsubaraGrid& bGrid, RealType t):
    CubicDMFTSC<D>(S,t,kGrid),
    _fGrid(fGrid),
    _bGrid(bGrid),
    GD0(std::tuple_cat(std::make_tuple(_fGrid),__repeater<KMesh,D>::get_tuple(_kGrid))),
    GD(GD0.getGrids()),
    SigmaD(GD0.getGrids()), 
    GLat(GD0.getGrids()),
    GLatLoc(_fGrid)
{
    _initialize();
};


template <size_t D>
void DFLadder<D>::_initialize()
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
ComplexType DFLadder<D>::getBubble2(BMatsubaraGrid::point W, KP... q, FMatsubaraGrid::point w1) const
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
inline typename DFLadder<D>::GKType DFLadder<D>::getGLat(const FMatsubaraGrid &gridF ) const
{
    GKType out(std::tuple_cat(std::forward_as_tuple(gridF), __repeater<KMesh,D>::get_tuple(_kGrid)));
    auto f1 = [&](const typename GKType::PointTupleType &in){return this->GLat(in);};
    auto f2 = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(f1);
    out.fill(f2);
    out._f = GLat._f;
    return out;
}

template <size_t D>
inline typename DFLadder<D>::GKType DFLadder<D>::getGLatDMFT(const FMatsubaraGrid& gridF) const 
{ 
    return CubicDMFTSC<D>::getGLat(gridF); 
};

template <size_t D>
inline typename DFLadder<D>::GKType DFLadder<D>::getGLat() const 
{ 
    return GLat; 
};

template <size_t D>
typename DFLadder<D>::GLocalType DFLadder<D>::operator()()
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
        GKType addSigma(this->SigmaD.getGrids()), GD_shift(this->GD.getGrids());

        size_t nq = 1;
        // iterate over all unique kpoints (patch of BZ)
        for (auto pts_it = unique_q_pts.begin(); pts_it != unique_q_pts.end(); pts_it++) { 
            std::array<KMesh::point, D> q = pts_it->first; // point
            RealType q_weight = RealType(pts_it->second); // it's weight
            INFO_NONEWLINE(nq << "/" << unique_q_pts.size() << ": [");
            for (size_t i=0; i<D; ++i) INFO_NONEWLINE(RealType(q[i]) << " "); INFO("]. Weight : " << q_weight);

            auto Wq_args_static = std::tuple_cat(std::make_tuple(0.0),q);
            GD_shift = GD.shift(Wq_args_static);
            auto dual_bubble = Diagrams::getBubble(this->GD, Wq_args_static);
            auto dual_bubble_matrix = dual_bubble.getData().getAsDiagonalMatrix();
            decltype(StaticV4) FullStaticV4;
            if (_EvaluateStaticDiagrams) {
                INFO_NONEWLINE("\tStatic contribution...");
                FullStaticV4 = Diagrams::BS(dual_bubble_matrix, StaticV4, true, _eval_BS_SC, _n_BS_iter, _BSmix);
            };
            
            decltype(StaticVertex4) DualBubbleDynamic(StaticVertex4.getGrids());
            decltype(StaticVertex4)::PointFunctionType dbfill = [&](FMatsubaraGrid::point w1, FMatsubaraGrid::point w2){
                return -T*(GD[size_t(w1)]*GD_shift[size_t(w2)]).sum()/RealType(totalqpts);
                };
            decltype(StaticVertex4) DynamicFullVertex4(StaticVertex4.getGrids());
            if (_EvaluateDynamicDiagrams) { 
                INFO_NONEWLINE("\tDynamic contribution...");
                DualBubbleDynamic.fill(dbfill);
                DynamicFullVertex4 = Diagrams::BS(DualBubbleDynamic, StaticVertex4*(-1.0), true, _eval_BS_SC,(_n_BS_iter>0?_n_BS_iter-1:0),_BSmix);
                DualBubbleDynamic*=(-1.0)*StaticVertex4*DynamicFullVertex4;
                addSigma=0.0;
                for (auto w : _fGrid.getPoints()) { 
                    for (auto w1 : _fGrid.getPoints()) { 
                        //addSigma[size_t(w)] += T*GD_shift[size_t(w1)]*(-1.0)*StaticVertex4(w,w1)*DualBubbleDynamic(w,w1)*DynamicFullVertex4(w,w1);
                        addSigma[size_t(w)] += T*GD_shift[size_t(w1)]*DualBubbleDynamic(w,w1);
                            }
                        };
                INFO2("Sigma dynamic contribution diff = " << addSigma.diff(SigmaD*0));
                SigmaD+=addSigma/RealType(totalqpts)*q_weight;
                };

            if (_EvaluateStaticDiagrams) { 
                typename GKType::PointFunctionType SigmaF;
                auto SigmaF2 = [&](wkPointTupleType in)->ComplexType { 
                    FMatsubaraGrid::point w = std::get<0>(in);
                    auto kpts = __tuple_tail(in); 
                    auto out_static = (1.0*T)*( (FullStaticV4(size_t(w),size_t(w)))*GD_shift(in)); // Here the fact that StaticVertex4(w,w) = 0 is used.
                    return out_static;
                    };
                SigmaF = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(SigmaF2);
                addSigma.fill(SigmaF);
                INFO2("Sigma static contribution diff = " << addSigma.diff(SigmaD*0));
                SigmaD+=addSigma/RealType(totalqpts)*q_weight;
                };
            nq++;
            }; // end of q loop

        //GLocalType SigmaD0(_fGrid),__GD0(_fGrid), __GDnew0(_fGrid);
        //std::function<ComplexType(FMatsubaraGrid::point)> sf = [&](FMatsubaraGrid::point w){return SigmaD[size_t(w)][0][0];};
        //SigmaD0.fill(sf);        
        //DEBUG(SigmaD0);

        auto GD_new = _GDmix/(1.0/GD0 - SigmaD) + GD*(1.0-_GDmix); // Dyson eq;
        //std::function<ComplexType(FMatsubaraGrid::point)> gdf = [&](FMatsubaraGrid::point w){return GD[size_t(w)][0][0];};
        //std::function<ComplexType(FMatsubaraGrid::point)> gdnewf = [&](FMatsubaraGrid::point w){return GD_new[size_t(w)][0][0];};
        //__GD0.fill(gdf);
        //__GDnew0.fill(gdnewf);
        //DEBUG(__GD0);
        //DEBUG(__GDnew0);
        diffGD = GD_new.diff(GD);
        if (diffGD<diffGD_min) { diffGD_min = diffGD; diffGD_min_count = 0; }
        else diffGD_min_count++;
        INFO2("DF diff = " << diffGD);
        if (diffGD_min_count > 6 ) {
            ERROR("\n\tCaught loop cycle. Reducing DF mixing to " << _GDmix/2 << " .\n");
            _GDmix/=2.;
            diffGD_min_count = 0;
            }

        diffDF_stream.open("diffDF.dat",std::ios::app);
        diffDF_stream << diffGD << "  " << _GDmix << std::endl;
        diffDF_stream.close();
        GD=GD_new;
        GD._f = GD0._f; // assume DMFT asymptotics are good 
        SigmaD = 0.0;

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
inline typename DFLadder<D>::GLocalType DFLadder<D>::getGLoc()
{
    return GLatLoc; 
}

template <size_t D>
std::tuple<typename DFLadder<D>::SuscType> DFLadder<D>::calculateLatticeData(const BMatsubaraGrid& gridB)
{
    KMeshPatch fullgrid(_kGrid);
    std::array<KMeshPatch, D> grids = __repeater<KMeshPatch,D>::get_array(fullgrid); 
    auto t1 = __repeater<KMeshPatch,D>::get_tuple(fullgrid); 
    return calculateLatticeData(gridB, grids);
}

template <size_t D>
std::tuple<typename DFLadder<D>::SuscType> DFLadder<D>::calculateLatticeData(const BMatsubaraGrid& gridB, const std::array<KMeshPatch, D>& qgrids)
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

// Implementations to compile

template struct DFLadder<1>;  
template struct DFLadder<2>;  
template struct DFLadder<3>;  
//template struct DFLadder<4>;  

/*
            for (auto iW : _bGrid.getPoints()) {
                INFO2("iW = " << iW);
                auto Wq_args = std::tuple_cat(std::make_tuple(iW),q);
                typename GLocalType::PointFunctionType VertexFillf = std::bind(
                    &FKImpuritySolver::getVertex4<BMatsubaraGrid::point, FMatsubaraGrid::point, FMatsubaraGrid::point>, std::cref(_S), iW, std::placeholders::_1, std::placeholders::_1); 
                typename GLocalType::FunctionType Vertexf = 
                    std::bind(&FKImpuritySolver::getVertex4<ComplexType, ComplexType, ComplexType>, std::cref(_S), ComplexType(iW), std::placeholders::_1, std::placeholders::_1); 
                DynVertex4.fill(VertexFillf);
                DynVertex4._f = Vertexf;
                // Bethe-Salpeter vertex
                DualDynBubble = Diagrams::getBubble(GD, Wq_args);
                FullDualDynVertex4 = Diagrams::BS(DualDynBubble, DynVertex4, true, _eval_BS_SC, (_n_BS_iter>0?_n_BS_iter-1:0), _BSmix);
                // Sigma
                auto GD_shift_dyn = GD.shift(Wq_args);
                typename GKType::PointFunctionType SigmaF;
                auto SigmaF2 = [&](wkPointTupleType in)->ComplexType { 
                    auto w = std::get<0>(in);
                    return DynVertex4(w)*DualDynBubble(w)*GD_shift(in)*FullDualDynVertex4(w);
                    };
                SigmaF = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(SigmaF2);
                addSigma.fill(SigmaF);
                addSigma*=T;
                INFO3("Sigma contribution diff = " << addSigma.diff(addSigma*0));
                SigmaD+=addSigma/totalqpts;
                };
            */

}; 
