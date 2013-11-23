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
    
template <size_t M> 
struct CubicTraits{ 
    /** Returns an analytic std::function of the dispersion. */
    template <typename FunctionType, typename ...ArgTypes> static FunctionType get_dispersion(RealType t);
    /** Finds the equivalent point, which is used is calculations. */
    static std::array<RealType,M> findSymmetricBZPoint(const std::array<RealType,M>& in);
    static BZPoint<M> findSymmetricBZPoint(const BZPoint<M>& in, const KMesh& kGrid);
    /** Returns a vector of D-dimensional arrays of points on the KMesh, which covers the Brillouin zone */
    static std::vector<BZPoint<M>> getAllBZPoints(const KMesh& in);
    /** Returns a vector of a pair of a D-dimensional arrays of points on the KMesh and the amount of points that can be obtained from a symmetry operation in the lattice. */
    static std::map<BZPoint<M>, std::vector<BZPoint<M>>> getUniqueBZPoints(const KMesh& kGrid);
};

template <>
struct CubicTraits<0>{ 
    /** Actual dispersion relation. */
    template <typename ArgType1, typename ...ArgTypes> static RealType ek(RealType t, ArgType1 kpoint1, ArgTypes... kpoints); 
    template <typename ArgType1> static RealType ek(RealType t, ArgType1 kpoint1); 
    /** Returns an analytic std::function of the dispersion. */
    template <typename FunctionType, typename ...ArgTypes> static FunctionType get_dispersion(RealType t) { return [t](ArgTypes ... in){return ek(t,in...);};};
};

//
// CubicTraits
//

template <size_t M> 
template <typename FunctionType, typename ...ArgTypes> 
inline FunctionType CubicTraits<M>::get_dispersion(RealType t)
{
    return CubicTraits<M-1>::template get_dispersion<FunctionType,ArgTypes...,RealType>(t); 
}

template <typename ArgType1, typename ...ArgTypes> 
inline RealType CubicTraits<0>::ek(RealType t, ArgType1 kpoint1, ArgTypes... kpoints) 
{
    static_assert(std::is_convertible<ArgType1, RealType>::value,"Wrong kpoint");
    assert (kpoint1>=0 && kpoint1 < 2*PI);
    return -2.0*t*cos(kpoint1)+ek(t, kpoints...);
}
 
template <typename ArgType1> 
inline RealType CubicTraits<0>::ek(RealType t, ArgType1 kpoint1)
{
    static_assert(std::is_convertible<ArgType1, RealType>::value,"Wrong kpoint");
    assert (kpoint1>=0 && kpoint1 < 2*PI);
    return -2*t*cos(kpoint1);
}

template <size_t M>
std::vector<BZPoint<M>> CubicTraits<M>::getAllBZPoints(const KMesh& kGrid)
{
    size_t ksize = kGrid.getSize();
    size_t totalqpts = size_t(pow(ksize,M));
    std::vector<BZPoint<M>> out(totalqpts);

    std::array<KMesh::point, M> q;
    for (size_t nq=0; nq<totalqpts; ++nq) { // iterate over all kpoints
        size_t offset = 0;
        for (size_t i=0; i<M; ++i) { 
            q[M-1-i]=kGrid[(nq-offset)/(int(pow(ksize,i)))%ksize]; 
            offset+=(int(pow(ksize,i)))*size_t(q[M-1-i]); 
            };
        out[nq]=q;
        }
    return out;
}

template <size_t M>
inline std::array<RealType,M> CubicTraits<M>::findSymmetricBZPoint(const std::array<RealType,M>& in)
{
    std::array<RealType,M> out;
    for (size_t i=0; i<M; ++i) {
            in[i] = std::fmod(in[i],2.0*PI);
            if (RealType(in[i])>PI) out[i]=2.0*PI-in[i];
            }
        // Order x,y,z. Ensures x<=y<=z
        std::sort(out.begin(), out.end());
    return out;
}

template <size_t M>
inline BZPoint<M> CubicTraits<M>::findSymmetricBZPoint(const BZPoint<M>& in, const KMesh& kGrid)
{
    BZPoint<M> out(in);
    // Flip all pi+x to pi-x
    for (size_t i=0; i<M; ++i) {
        if (RealType(in[i])>PI) out[i]=kGrid.findClosest(2.0*PI-RealType(in[i]));
        }
    // Order x,y,z. Ensures x<=y<=z
    std::sort(out.begin(), out.end());
    return out;
}

template <size_t M>
std::map<BZPoint<M>, std::vector<BZPoint<M>>> CubicTraits<M>::getUniqueBZPoints(const KMesh& kGrid)
{
    auto all_pts = getAllBZPoints(kGrid);
    auto totalqpts = all_pts.size();
    std::map<std::array<KMesh::point, M>, std::vector<BZPoint<M>>> unique_pts;
    for (size_t nq=0; nq<totalqpts; ++nq) {
        auto q = all_pts[nq];
//        INFO_NONEWLINE("Considering: " << nq << "/" << totalqpts << " " << q);
        BZPoint<M> q_unique = findSymmetricBZPoint(q, kGrid);

        if (unique_pts.find(q_unique)==unique_pts.end()) {
            unique_pts[q_unique]=std::vector<BZPoint<M>>();
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

} // end of namespace FK

#endif // endif::#ifndef __FK_LATTICE_TRAITS_HPP__
