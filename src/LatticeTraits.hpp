#ifndef __FK_LATTICE_TRAITS_HPP__
#define __FK_LATTICE_TRAITS_HPP__

namespace FK {

/** A typedef for a point in the Brillouin zone. */
template <size_t D>
using BZPoint = typename std::array<KMesh::point,D>;

/** Stream the BZPoint. */
template <size_t D>
std::ostream& operator<<(std::ostream& out, BZPoint<D> in)
{for (size_t i=0; i<D; ++i) out << RealType(in[i]) << " "; return out; } 

template <size_t D, class Derived>
struct LatticeTraitsBase {
    constexpr static size_t _D = D;
    typedef typename ArgFunGenerator<D,RealType,RealType>::type ArgFunType;
    typedef typename ArgBackGenerator<D,RealType,std::tuple>::type ArgTupleType;

    template <typename Arg1, typename ...Args> 
        typename std::enable_if<sizeof...(Args) == D-1 && std::is_convertible<Arg1,RealType>::value, RealType>::type 
            dispersion(Arg1 in1, Args... in) { return static_cast<Derived*>(this)->dispersion(in1, in...); };

    template <typename ...Args> 
        typename std::enable_if<sizeof...(Args) == D, RealType>::type 
        dispersion(std::tuple<Args...> in){
             __caller_tuple<RealType,std::tuple<Args...>> caller; 
             caller._f = [this](Args... in){return dispersion(in...);};
             caller._params = in;
             return caller.call();
        };
    /** Returns an analytic std::function of the dispersion. */
    ArgFunType get_dispersion(){ return __fun_traits<ArgFunType>::getFromTupleF ( [this](ArgTupleType in){return dispersion(in);});};
    /** Finds the equivalent point, which is used is calculations. */
    static BZPoint<D> findSymmetricBZPoint(const BZPoint<D>& in, const KMesh& kGrid)
        { return Derived::findSymmetricBZPoint(in, kGrid); };
    /** Returns a vector of D-dimensional arrays of points on the KMesh, which covers the Brillouin zone */
    static std::vector<BZPoint<D>> getAllBZPoints(const KMesh& in);
    /** Returns a vector of a pair of a D-dimensional arrays of points on the KMesh and the amount of points that can be obtained from a symmetry operation in the lattice. */
    static std::map<BZPoint<D>, std::vector<BZPoint<D>>> getUniqueBZPoints(const KMesh& kGrid);
    RealType disp_square_sum(){return static_cast<Derived*>(this)->disp_square_sum();}; 

};
    
template <size_t D> 
struct CubicTraits : LatticeTraitsBase<D,CubicTraits<D>>{ 
    RealType _t = 1.0;
    CubicTraits(RealType t):_t(t){};
    static BZPoint<D> findSymmetricBZPoint(const BZPoint<D>& in, const KMesh& kGrid);
    template <typename Arg1, typename ...Args> 
        typename std::enable_if<sizeof...(Args) == D-1 && std::is_convertible<Arg1,RealType>::value, RealType>::type 
        dispersion(Arg1 in1, Args... in) {
            return -2.0*_t*cos(RealType(in1)) + CubicTraits<D-1>(_t).dispersion(in...);
            };
    RealType disp_square_sum(){return 2.*_t*_t*D;}; 
};

template <>
struct CubicTraits<0>{ 
    RealType _t = 1.0;
    constexpr static size_t _D = 0;
    CubicTraits(RealType t):_t(t){};
    RealType dispersion(){return 0.0;};
    RealType dispersion(std::tuple<>){return 0.0;};
};

struct TriangularTraits : LatticeTraitsBase<2,TriangularTraits> {
    RealType _t = 1.0;
    RealType _tp = 1.0;

    RealType dispersion(RealType kx,RealType ky){return -2.*_t*(cos(kx)+cos(ky)) - 2.*_tp*cos(kx-ky);};
    RealType disp_square_sum(){return 4.*_t*_t + 2.*_tp*_tp;}; 
    static BZPoint<2> findSymmetricBZPoint(const BZPoint<2>& in, const KMesh& kGrid);
    };


//
// CubicTraits
//
template <size_t D>
inline BZPoint<D> CubicTraits<D>::findSymmetricBZPoint(const BZPoint<D>& in, const KMesh& kGrid)
{
    BZPoint<D> out(in);
    // Flip all pi+x to pi-x
    for (size_t i=0; i<D; ++i) {
        if (RealType(in[i])>PI) out[i]=kGrid.findClosest(2.0*PI-RealType(in[i]));
        }
    // Order x,y,z. Ensures x<=y<=z
    std::sort(out.begin(), out.end());
    return out;
}
//
// TriangularTraits
//
inline BZPoint<2> TriangularTraits::findSymmetricBZPoint(const BZPoint<2>& in, const KMesh& kGrid)
{
    RealType x = in[0]; RealType y = in[1];
    if (x>PI && y>PI) { x = PI-x; y = PI-y; };
    BZPoint<2> out(in); 
    out[0] = kGrid.findClosest(x); out[1] = kGrid.findClosest(y);
    std::sort(out.begin(), out.end());
    return out;
}


//
// LatticeTraitsBase
//
template <size_t D, class Derived>
std::vector<BZPoint<D>> LatticeTraitsBase<D,Derived>::getAllBZPoints(const KMesh& kGrid)
{
    size_t ksize = kGrid.getSize();
    size_t totalqpts = size_t(pow(ksize,D));
    std::vector<BZPoint<D>> out(totalqpts);

    std::array<KMesh::point, D> q;
    for (size_t nq=0; nq<totalqpts; ++nq) { // iterate over all kpoints
        size_t offset = 0;
        for (size_t i=0; i<D; ++i) { 
            q[D-1-i]=kGrid[(nq-offset)/(int(pow(ksize,i)))%ksize]; 
            offset+=(int(pow(ksize,i)))*size_t(q[D-1-i]); 
            };
        out[nq]=q;
        }
    return out;
}

template <size_t D, class Derived>
std::map<BZPoint<D>, std::vector<BZPoint<D>>> LatticeTraitsBase<D,Derived>::getUniqueBZPoints(const KMesh& kGrid)
{
    auto all_pts = getAllBZPoints(kGrid);
    auto totalqpts = all_pts.size();
    std::map<std::array<KMesh::point, D>, std::vector<BZPoint<D>>> unique_pts;
    for (size_t nq=0; nq<totalqpts; ++nq) {
        auto q = all_pts[nq];
//        INFO_NONEWLINE("Considering: " << nq << "/" << totalqpts << " " << q);
        BZPoint<D> q_unique = findSymmetricBZPoint(q, kGrid);

        if (unique_pts.find(q_unique)==unique_pts.end()) {
            unique_pts[q_unique]=std::vector<BZPoint<D>>();
            unique_pts[q_unique].push_back(q_unique);
            }
        else if (q_unique != q)  
            unique_pts[q_unique].push_back(q);
        };
    size_t count = 0;
    for (auto it = unique_pts.begin(); it!=unique_pts.end(); it++) { 
//        DEBUG(it->first << " : " << it->second.size()); 
        count+=it->second.size(); 
        };
//    DEBUG(totalqpts << " == " << count);
    assert(totalqpts == count);
    return unique_pts;
} 

/*
template <size_t D, class Derived>
inline std::array<RealType,D> LatticeTraitsBase<D,Derived>::findSymmetricBZPoint(const std::array<RealType,D>& in)
{
    std::array<RealType,D> out;
    for (size_t i=0; i<D; ++i) {
            in[i] = std::fmod(in[i],2.0*PI);
            if (RealType(in[i])>PI) out[i]=2.0*PI-in[i];
            }
        // Order x,y,z. Ensures x<=y<=z
        std::sort(out.begin(), out.end());
    return out;
}
*/

} // end of namespace FK

#endif // endif::#ifndef __FK_LATTICE_TRAITS_HPP__
