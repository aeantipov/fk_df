#ifndef __FK_LATTICE_TRAITS_HPP__
#define __FK_LATTICE_TRAITS_HPP__

namespace FK {

/** A typedef for a point in the Brillouin zone. */
template <size_t D>
using BZPoint = typename std::array<kmesh::point,D>;

/** Stream the BZPoint. */
template <size_t D>
std::ostream& operator<<(std::ostream& out, BZPoint<D> in)
{for (size_t i=0; i<D; ++i) out << real_type(in[i]) << " "; return out; } 

template <size_t D, class Derived>
struct LatticeTraitsBase {
    constexpr static size_t _D = D;
    typedef typename tools::ArgFunGenerator<D,real_type,real_type>::type ArgFunType;
    typedef typename tools::ArgBackGenerator<D,real_type,std::tuple>::type arg_tuple;

    template <typename Arg1, typename ...Args> 
        typename std::enable_if<sizeof...(Args) == D-1 && std::is_convertible<Arg1,real_type>::value, real_type>::type 
            dispersion(Arg1 in1, Args... in) { return static_cast<Derived*>(this)->dispersion(in1, in...); };

    template <typename ...Args> 
        typename std::enable_if<sizeof...(Args) == D, real_type>::type 
        dispersion(std::tuple<Args...> in){
             auto f = [this](Args... in){return dispersion(in...);};
             return tuple_tools::unfold_tuple(f,in);
        };
    /** Returns an analytic std::function of the dispersion. */
    ArgFunType get_dispersion(){ return tools::extract_tuple_f(std::function<real_type(arg_tuple)>([this](arg_tuple in){return dispersion(in);})); }
    /** Finds the equivalent point, which is used is calculations. */
    static BZPoint<D> findSymmetricBZPoint(const BZPoint<D>& in, const kmesh& kGrid)
        { return Derived::findSymmetricBZPoint(in, kGrid); };
    /** Returns a vector of D-dimensional arrays of points on the kmesh, which covers the Brillouin zone */
    static std::vector<BZPoint<D>> getAllBZPoints(const kmesh& in);
    /** Returns a vector of a pair of a D-dimensional arrays of points on the kmesh and the amount of points that can be obtained from a symmetry operation in the lattice. */
    static std::map<BZPoint<D>, std::vector<BZPoint<D>>> getUniqueBZPoints(const kmesh& kGrid);
    real_type disp_square_sum(){return static_cast<Derived*>(this)->disp_square_sum();}; 

};
    
template <size_t D> 
struct CubicTraits : LatticeTraitsBase<D,CubicTraits<D>>{ 
    real_type _t = 1.0;
    CubicTraits(real_type t):_t(t){};
    static BZPoint<D> findSymmetricBZPoint(const BZPoint<D>& in, const kmesh& kGrid);
    template <typename Arg1, typename ...Args> 
        typename std::enable_if<sizeof...(Args) == D-1 && std::is_convertible<Arg1,real_type>::value, real_type>::type 
        dispersion(Arg1 in1, Args... in) {
            return -2.0*_t*cos(real_type(in1)) + CubicTraits<D-1>(_t).dispersion(in...);
            };
    real_type disp_square_sum(){return 2.*_t*_t*D;}; 
};

template <>
struct CubicTraits<0>{ 
    real_type _t = 1.0;
    constexpr static size_t _D = 0;
    CubicTraits(real_type t):_t(t){};
    real_type dispersion(){return 0.0;};
    real_type dispersion(std::tuple<>){return 0.0;};
};

struct TriangularTraits : LatticeTraitsBase<2,TriangularTraits> {
    real_type _t = 1.0;
    real_type _tp = 1.0;
    TriangularTraits(real_type t, real_type tp):_t(t),_tp(tp){};

    real_type dispersion(real_type kx,real_type ky){return -2.*_t*(cos(kx)+cos(ky)) - 2.*_tp*cos(kx-ky);};
    real_type disp_square_sum(){return 4.*_t*_t + 2.*_tp*_tp;}; 
    static BZPoint<2> findSymmetricBZPoint(const BZPoint<2>& in, const kmesh& kGrid);
    };


//
// CubicTraits
//
template <size_t D>
inline BZPoint<D> CubicTraits<D>::findSymmetricBZPoint(const BZPoint<D>& in, const kmesh& kGrid)
{
    BZPoint<D> out(in);
    // Flip all pi+x to pi-x
    for (size_t i=0; i<D; ++i) {
        if (real_type(in[i])>PI) out[i]=kGrid.find_nearest(2.0*PI-real_type(in[i]));
        }
    // Order x,y,z. Ensures x<=y<=z
    std::sort(out.begin(), out.end());
    return out;
}
//
// TriangularTraits
//
inline BZPoint<2> TriangularTraits::findSymmetricBZPoint(const BZPoint<2>& in, const kmesh& kGrid)
{
    real_type x = in[0]; real_type y = in[1];
    if (x>PI && y>PI) { x = 2.0*PI-x; y = 2.0*PI-y; };
    BZPoint<2> out(in); 
    out[0] = kGrid.find_nearest(x); out[1] = kGrid.find_nearest(y);
    std::sort(out.begin(), out.end());
    return out;
}


//
// LatticeTraitsBase
//
template <size_t D, class Derived>
std::vector<BZPoint<D>> LatticeTraitsBase<D,Derived>::getAllBZPoints(const kmesh& kGrid)
{
    size_t ksize = kGrid.size();
    size_t totalqpts = size_t(pow(ksize,D));
    std::vector<BZPoint<D>> out(totalqpts,tuple_tools::repeater<typename kmesh::point, D>::get_array(kGrid[0]));

    std::array<kmesh::point, D> q = tuple_tools::repeater<typename kmesh::point, D>::get_array(kGrid[0]);
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
std::map<BZPoint<D>, std::vector<BZPoint<D>>> LatticeTraitsBase<D,Derived>::getUniqueBZPoints(const kmesh& kGrid)
{
    auto all_pts = getAllBZPoints(kGrid);
    auto totalqpts = all_pts.size();
    std::map<std::array<kmesh::point, D>, std::vector<BZPoint<D>>> unique_pts;
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
inline std::array<real_type,D> LatticeTraitsBase<D,Derived>::findSymmetricBZPoint(const std::array<real_type,D>& in)
{
    std::array<real_type,D> out;
    for (size_t i=0; i<D; ++i) {
            in[i] = std::fmod(in[i],2.0*PI);
            if (real_type(in[i])>PI) out[i]=2.0*PI-in[i];
            }
        // Order x,y,z. Ensures x<=y<=z
        std::sort(out.begin(), out.end());
    return out;
}
*/

} // end of namespace FK

#endif // endif::#ifndef __FK_LATTICE_TRAITS_HPP__
