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

template <class Solver, size_t D, size_t ksize>
DFLadder<Solver,D,ksize>::DFLadder(const Solver &S, const FMatsubaraGrid& fGrid, const BMatsubaraGrid& bGrid, RealType t):
    CubicDMFTSC<Solver,D,ksize>(S,t),
    _fGrid(fGrid),
    _bGrid(bGrid),
    Diagrams(DFDiagrams<D>(_fGrid,_kGrid)),
    GD0(std::tuple_cat(std::make_tuple(_fGrid),__repeater<KMesh,D>::get_tuple(_kGrid))),
    GD(GD0.getGrids()),
    SigmaD(GD0.getGrids()), 
    GLat(GD0.getGrids()),
    GLatLoc(_fGrid)
{
    _initialize();
};


template <class Solver, size_t D, size_t ksize>
void DFLadder<Solver,D,ksize>::_initialize()
{
    GLat = CubicDMFTSC<Solver, D, ksize>::getGLat(_fGrid);
    for (auto iw : _fGrid.getVals()) {
        size_t iwn = size_t(iw);
        GD0[iwn] = GLat[iwn] - _S.gw(iw);
    };

    auto gd_f = [&](const wkTupleType &in)->ComplexType{
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
template <class Solver, size_t D, size_t ksize>
template <typename ...KP>
ComplexType DFLadder<Solver,D,ksize>::getBubble2(BMatsubaraGrid::point W, KP... q, FMatsubaraGrid::point w1) const
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

template <class Solver, size_t D, size_t ksize>
inline typename DFLadder<Solver,D,ksize>::GKType DFLadder<Solver,D,ksize>::getGLat(const FMatsubaraGrid &gridF ) const
{
    GKType out(std::tuple_cat(std::forward_as_tuple(gridF), __repeater<KMesh,D>::get_tuple(_kGrid)));
    auto f1 = [&](const typename GKType::PointTupleType &in){return this->GLat(in);};
    auto f2 = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(f1);
    out.fill(f1);
    out._f = GLat._f;
    return out;
}

template <class Solver, size_t D, size_t ksize>
typename DFLadder<Solver,D,ksize>::GLocalType DFLadder<Solver,D,ksize>::operator()()
{
    INFO("Using DF Ladder self-consistency in " << D << " dimensions on a cubic lattice of " << ksize << "^" << D <<" atoms.");
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
    auto Wqgrids = std::tuple_cat(std::make_tuple(_bGrid),__repeater<KMesh,D>::get_tuple(_kGrid));
    _initialize();

    // Put here operations with GD
    GLocalType Vertex4(_fGrid);
    GLocalType Chi0(_fGrid);

    INFO2("Starting ladder dual fermion calculations")
    INFO2("Beginning with GD sum = " << std::abs(GD.sum())/RealType(_fGrid.getSize())/RealType(__power<ksize,D>::value));
    RealType diffGD = 1.0;
    for (size_t nd_iter=0; nd_iter<_n_GD_iter && diffGD > 1e-8; ++nd_iter) { 
        INFO2("DF iteration " << nd_iter);
        GKType addSigma(this->SigmaD.getGrids());
        INFO2("Evaluating Bethe-Salpeter equation");

        for (auto iW : _bGrid.getVals()) {
            INFO2("iW = " << iW);
            typename GLocalType::PointFunctionType VertexFillf = 
                std::bind(&FKImpuritySolver::getBVertex4<BMatsubaraGrid::point, FMatsubaraGrid::point>, std::cref(_S), iW, std::placeholders::_1); 
            typename GLocalType::FunctionType Vertexf = 
                std::bind(&FKImpuritySolver::getBVertex4<ComplexType, ComplexType>, std::cref(_S), ComplexType(iW), std::placeholders::_1); 
            Vertex4.fill(VertexFillf);
            Vertex4._f = Vertexf;
 
            size_t totalqpts = __power<ksize,D>::value;
            std::array<KMesh::point, D> q;
            for (size_t nq=0; nq<totalqpts; ++nq) { // iterate over all kpoints
                size_t offset = 0;
                for (size_t i=0; i<D; ++i) { q[D-1-i]=_kGrid[(nq-offset)/(int(pow(ksize,i)))%ksize]; offset+=(int(pow(ksize,i)))*size_t(q[D-1-i]); };
                WQTupleType Wq_args = std::tuple_cat(std::forward_as_tuple(iW),q);
                
                INFO_NONEWLINE("\t\t" << nq << "/" << totalqpts<< ". ");
                //DEBUG("qx = " << q[0] << ", qy = " << q[1]);

                // Bethe-Salpeter vertex
                Chi0 = Diagrams.getBubble(GD, Wq_args);
                auto FullDualVertex4 = Diagrams.BS(Chi0, Vertex4, _eval_BS_SC, _n_BS_iter, _BSmix);
                
                // Sigma
                auto GD_shift = GD.shift(Wq_args);
                typename GKType::PointFunctionType SigmaF;
                auto SigmaF2 = [&](wkTupleType in)->ComplexType { 
                    ComplexType w = std::get<0>(in);
                    return 0.5*Vertex4(w)*Chi0(w)*GD_shift(in)*FullDualVertex4(w);
                    };
                SigmaF = __fun_traits<typename GKType::PointFunctionType>::getFromTupleF(SigmaF2);
                addSigma.fill(SigmaF);
                SigmaD+=addSigma*T/2.0/totalqpts;
                } // end of q loop
        }; //end of iW loop

        auto GD_new = 1.0/(1.0/GD0 - SigmaD); // Dyson eq;
        diffGD = GD_new.diff(GD);
        INFO2("DF diff = " << diffGD);
        GD=GD_new*_GDmix + GD*(1.0-_GDmix);
        GD._f = GD0._f; // assume DMFT asymptotics are good 

       INFO2("GD sum = " << std::abs(GD.sum())/RealType(_fGrid.getSize())/RealType(__power<ksize,D>::value));
    };
    INFO("Finished DF iterations");
        
    for (auto iw : _fGrid.getVals()) {
        size_t iwn = size_t(iw);
        GLat[iwn] = 1.0/(Delta(iw) - _ek.getData()) + 1.0/(Delta(iw) - _ek.getData())/gw(iw)*GD[iwn]/gw(iw)/(Delta(iw) - _ek.getData());
        GDLoc[iwn] = GD[iwn].sum()/RealType(__power<ksize,D>::value);
        GLatLoc[iwn] = GLat[iwn].sum()/RealType(__power<ksize,D>::value);
        };
 
    // Finish - prepare all lattice quantities
    Delta_out = Delta + 1.0/gw * GDLoc / GLatLoc;
    // Assume DMFT asymptotics
    Delta_out._f = std::bind([&](ComplexType w)->ComplexType{return _t*_t*2*RealType(D)*((_S.mu-_S.w_1*_S.U)/std::abs(w*w) + 1.0/w);}, std::placeholders::_1);
    return Delta_out;
}

template <class Solver, size_t D, size_t ksize>
std::tuple<typename DFLadder<Solver,D,ksize>::SuscType> DFLadder<Solver,D,ksize>::calculateLatticeData(const BMatsubaraGrid& gridB)
{
    KMeshPatch fullgrid(_kGrid);
    std::array<KMeshPatch, D> grids = __repeater<KMeshPatch,D>::get_array(fullgrid); 
    auto t1 = __repeater<KMeshPatch,D>::get_tuple(fullgrid); 
    return calculateLatticeData(gridB, grids);
}

template <class Solver, size_t D, size_t ksize>
std::tuple<typename DFLadder<Solver,D,ksize>::SuscType> DFLadder<Solver,D,ksize>::calculateLatticeData(const BMatsubaraGrid& gridB, const std::array<KMeshPatch, D>& qgrids)
{
    SuscType LatticeSusc(std::tuple_cat(std::forward_as_tuple(gridB),qgrids));
    RealType T = 1.0/_fGrid._beta;
    GKType Lwk = GLat/GD0*(-1.0);
    Lwk._f = __fun_traits<typename GKType::FunctionType>::constant(-1.0);
    auto GDL = GD*Lwk;
    GLocalType Vertex4(_fGrid), Chi0(_fGrid);
        
    size_t totalqpts = 1;
    for (const auto& qmesh : qgrids) { 
            totalqpts*=qmesh.getSize();
        }

    INFO2("Calculating lattice susceptibility");
    for (auto iW : gridB.getVals()) {
        INFO2("iW = " << iW);
        typename GLocalType::PointFunctionType VertexFillf = 
            std::bind(&FKImpuritySolver::getBVertex4<BMatsubaraGrid::point, FMatsubaraGrid::point>, std::cref(_S), iW, std::placeholders::_1); 
        typename GLocalType::FunctionType Vertexf = 
            std::bind(&FKImpuritySolver::getBVertex4<ComplexType, ComplexType>, std::cref(_S), ComplexType(iW), std::placeholders::_1); 
        Vertex4.fill(VertexFillf);
        Vertex4._f = Vertexf;
 
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
            Chi0 = Diagrams.getBubble(GD, Wq_args);
            auto FullDualVertex4 = Diagrams.BS(Chi0, Vertex4, _eval_BS_SC, _n_BS_iter, _BSmix);
            
            // Lattice susceptibility
            GLocalType ChiLat0 = Diagrams.getBubble(GLat, Wq_args);
            auto GDL_chi0 = Diagrams.getBubble(GDL, Wq_args);
            auto ChiLat = GDL_chi0*FullDualVertex4*GDL_chi0 + ChiLat0;
            auto ChiVal = T*ChiLat.sum();
            LatticeSusc.get(Wq_args) = ChiVal;
        } // end of q loop
    }; //end of iW loop
    return std::forward_as_tuple(LatticeSusc);
}

} // end of namespace FK
#endif // endif :: #ifndef __FK_DF_H__
