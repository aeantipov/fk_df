#ifndef ___FK_GRID_HPP___
#define ___FK_GRID_HPP___

#include <numeric>
#include "Grid.h"

namespace FK { 

//
// Grid
//

template <typename ValueType, class Derived>
Grid<ValueType,Derived>::Grid()
{};

template <typename ValueType, class Derived>
Grid<ValueType,Derived>::Grid(const std::vector<point> &vals):_vals(vals)
{
};


template <typename ValueType, class Derived>
Grid<ValueType,Derived>::Grid(const std::vector<ValueType> &vals):_vals(vals)
{
    for (size_t i=0; i<_vals.size(); ++i) _vals[i]._index = i;
};

template <typename ValueType, class Derived>
Grid<ValueType,Derived>::Grid(int min, int max, std::function<ValueType (const int&)> f)
{
    if (max<min) std::swap(min,max);
    size_t n_points = max-min;
    _vals.resize(n_points); 
    for (int i=0; i<n_points; ++i) _vals[i]= point(f(min+i), i) ; 
}

template <typename ValueType, class Derived>
typename Grid<ValueType,Derived>::point Grid<ValueType,Derived>::operator[](size_t index) const
{
    if (index>_vals.size()) throw exWrongIndex();
    return _vals[index];
}

template <typename ValueType, class Derived>
const std::vector<typename Grid<ValueType,Derived>::point> & Grid<ValueType,Derived>::getVals() const
{
    return _vals;
}

template <typename ValueType, class Derived>
size_t Grid<ValueType,Derived>::getSize() const
{
    return _vals.size();
}

template <typename ValueType, class Derived>
template <class ArgType>
typename Grid<ValueType,Derived>::point Grid<ValueType,Derived>::shift(point in, ArgType shift_arg) const
{
    point out;
    out._val = in._val + ValueType(shift_arg);
    auto find_result = this->find(out._val);
    if (std::get<0>(find_result)) { out._index = std::get<1>(find_result); return (*this)[out._index]; }
    else { out._index = this->getSize(); 
           #ifndef NDEBUG
           ERROR("Returning point with an invalid index after shift.");
           #endif
           return out; };
}

template <typename ValueType, class Derived>
template <class ArgType>
ValueType Grid<ValueType,Derived>::shift(ValueType in, ArgType shift_arg) const
{
    return in+ValueType(shift_arg);
}


template <typename ValueType, class Derived>
std::ostream& operator<<(std::ostream& lhs, const Grid<ValueType,Derived> &gr)
{ 
    lhs << "{";
    //lhs << gr._vals;
    std::ostream_iterator<ValueType> out_it (lhs,", ");
    std::transform(gr._vals.begin(),gr._vals.end(), out_it, [](const typename Grid<ValueType,Derived>::point &x){return ValueType(x);});
    lhs << "}";
    return lhs;
}

template <typename ValueType, class Derived>
const char* Grid<ValueType,Derived>::exWrongIndex::what() const throw(){ 
     return "Index out of bounds";
};


//
// MatsubaraGrid
//

template <bool F>
inline MatsubaraGrid<F>::MatsubaraGrid(int min, int max, RealType beta):
    Grid<ComplexType, MatsubaraGrid<F>>(min,max,std::bind(Matsubara<F>, std::placeholders::_1, beta)),
    _beta(beta), 
    _spacing(PI/beta), 
    _w_min(min),
    _w_max(max)
{
}

template <bool F>
inline MatsubaraGrid<F>::MatsubaraGrid(const MatsubaraGrid<F> &rhs) : 
    Grid<ComplexType, MatsubaraGrid<F>>(rhs._vals),
    _beta(rhs._beta), 
    _spacing(rhs._spacing), 
    _w_min(rhs._w_min), 
    _w_max(rhs._w_max)
{
}

template <bool F>
inline MatsubaraGrid<F>::MatsubaraGrid(MatsubaraGrid<F>&& rhs):
    Grid<ComplexType, MatsubaraGrid>(rhs._vals), 
    _beta(rhs._beta), 
    _spacing(rhs._spacing), 
    _w_min(rhs._w_min), 
    _w_max(rhs._w_max)
{
}

template <bool F>
template <class Obj> 
inline auto MatsubaraGrid<F>::integrate(const Obj &in) const -> decltype(in(_vals[0]))
{
    decltype(in(this->_vals[0])) R = in(this->_vals[0]);
    R=std::accumulate(_vals.begin()+1, _vals.end(), R,[&](decltype(in(_vals[0]))& y,const decltype(_vals[0]) &x) {return y+in(x);}); 
    return R/_beta;
}

template <bool F>
template <class Obj, typename ...OtherArgTypes> 
inline auto MatsubaraGrid<F>::integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...))
{
    decltype(in(_vals[0],Args...)) R = in(_vals[0],Args...);
    R=std::accumulate(_vals.begin()+1, _vals.end(), R,[&](decltype(in(_vals[0]))& y,const decltype(_vals[0]) &x) {return y+in(x,Args...);}); 
    return R/_beta;
}

template <bool F>
template <class Obj> 
inline auto MatsubaraGrid<F>::prod(const Obj &in) const -> decltype(in(_vals[0]))
{
    decltype(in(_vals[0])) R = in(_vals[0]);
    R=std::accumulate(_vals.begin()+1, _vals.end(), R,[&](decltype(in(_vals[0]))& y,const decltype(_vals[0]) &x) {return y*in(x);}); 
    //decltype(in(_vals[0])) R = in(_vals[_vals.size()/2]);
    //R=std::accumulate(_vals.begin()+1+_vals.size()/2, _vals.end(), R,[&](decltype(in(_vals[0]))& y,decltype(_vals[0]) &x) {DEBUG(x << "|" << in(x) << "|" << y << "->" << y*in(x)); return y*in(x);}); 
    //R=std::accumulate(_vals.begin(), _vals.begin()+_vals.size()/2, R,[&](decltype(in(_vals[0]))& y,decltype(_vals[0]) &x) {DEBUG(x << "|" << in(x) << "|" << y << "->" << y*in(x)); return y*in(x);}); 
    return R;
}


/*
template <class Obj> 
auto MatsubaraGrid::gridIntegrate(const std::vector<Obj> &in) const -> Obj
{
    decltype(in[0]) R = in[0];
    R=std::accumulate(_vals.begin()+1, _vals.end(),R,[&](decltype(in[0])& y, decltype(in[0]) &x) {return y+x;}); 
    return R/_beta;
}
*/

template <bool F>
inline std::tuple <bool, size_t, RealType> MatsubaraGrid<F>::find (ComplexType in) const
{
    int n=getNumber(in);
    #ifndef NDEBUG
    DEBUG("Invoking matsubara find");
    #endif
    if (n<_w_min) { 
        #ifndef NDEBUG
        ERROR("Frequency to find is out of bounds, " << in << "<" << FMatsubara(_w_min,_beta)); 
        #endif
        return std::make_tuple(0,0,0); 
        };
    if (n>=_w_max) { 
        #ifndef NDEBUG
        ERROR("Frequency to find is out of bounds, " << in << ">" << FMatsubara(_w_max,_beta)); 
        #endif
        return std::make_tuple(0,_vals.size(),0); 
        };
    return std::make_tuple (1,n-_w_min,1);
}

template <bool F>
inline int MatsubaraGrid<F>::getNumber(ComplexType in) const
{
    assert (std::abs(real(in))<std::numeric_limits<RealType>::epsilon());
    return std::lround(imag(in)/_spacing-F)/2;
};

template <bool F>
template <class Obj>
inline auto MatsubaraGrid<F>::getValue(Obj &in, ComplexType x) const ->decltype(in[0]) 
{
    const auto find_result=this->find(x);
    if (!std::get<0>(find_result)) { throw (exWrongIndex()); } 
    return in[std::get<1>(find_result)];
}


template <bool F>
template <class Obj>
inline auto MatsubaraGrid<F>::getValue(Obj &in, MatsubaraGrid::point x) const ->decltype(in[0]) 
{
    if (x._index < _vals.size() && x == _vals[x._index])
    return in[x._index];
    else { 
        #ifndef NDEBUG
        ERROR ("Point not found"); 
        #endif
        return this->getValue(in, ComplexType(x)); 
         };
}

//
// RealGrid
//
inline RealGrid::RealGrid(RealType min, RealType max, size_t n_points):
    Grid(0,n_points,[n_points,max,min](size_t in){return (max-min)/n_points*in+min;}),
    _min(min),
    _max(max)
{
}

template <class Obj> 
inline auto RealGrid::integrate(const Obj &in) const -> decltype(in(_vals[0]))
{
    decltype(in(_vals[0])) R=0.0;
    for (int i=0; i<_vals.size()-1; ++i) {
        R+=0.5*(in(_vals[i])+in(_vals[i+1]))*(_vals[i+1]-_vals[i]);
        }
    return R;
}

template <class Obj, typename ...OtherArgTypes> 
inline auto RealGrid::integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...))
{
    decltype(in(_vals[0],Args...)) R=0.0;

    for (int i=0; i<_vals.size()-1; ++i) {
        R+=0.5*(in(_vals[i],Args...)+in(_vals[i+1],Args...))*(_vals[i+1]-_vals[i]);
        }
    return R;
}

inline std::tuple <bool, size_t, RealType> RealGrid::find (RealType in) const
{
    #ifndef NDEBUG
    DEBUG("Invoking find");
    #endif
    if (in<_min) { ERROR("Point to find is out of bounds, " << in << "<" << _min ); return std::make_tuple(0,0,0); };
    if (in>=_max) { ERROR("Point to find is out of bounds, " << in << ">" << _max ); return std::make_tuple(0,_vals.size(),0); };
    auto out = std::lower_bound (_vals.begin(), _vals.end(), in);
    size_t i = size_t(out-_vals.begin());
    RealType val_i = _vals[i];
    RealType weight=(in-val_i)/(_vals[i+1]/val_i);
    return std::make_tuple (1,i,weight);
}


template <class Obj>
inline auto RealGrid::getValue(Obj &in, RealType x) const ->decltype(in[0]) 
{
    const auto find_result=this->find(x);
    if (!std::get<0>(find_result)) throw (exWrongIndex()); 
    return in[std::get<1>(find_result)];
}

template <class Obj>
inline auto RealGrid::getValue(Obj &in, RealGrid::point x) const ->decltype(in[0]) 
{
    if (x._index < _vals.size() && x == _vals[x._index])
    return in[x._index];
    else { ERROR ("Point not found"); return this->getValue(in, RealType(x)); };
}


//
// KMesh
//

inline KMesh::KMesh(size_t n_points):
Grid(0,n_points,[n_points](size_t in){return 2.0*PI/n_points*in;}),
_points(n_points)
{
}

inline KMesh::KMesh(const KMesh& rhs):Grid(rhs._vals),_points(rhs._points)
{
}

inline KMesh::KMesh(KMesh &&rhs):Grid(rhs._vals),_points(rhs._points)
{
}

inline std::tuple <bool, size_t, RealType> KMesh::find (RealType in) const
{
    assert(in>=0 && in < 2.0*PI);
    int n = std::lround(in/2.0/PI*_points);
    if (n<0) { ERROR("KMesh point is out of bounds, " << in << "<" << 0); return std::make_tuple(0,0,0); };
    if (n>=_points) { ERROR("KMesh point is out of bounds, " << in << "> 2*PI"); return std::make_tuple(0,_points,0); };
    RealType weight=in/2.0/PI*_points-RealType(n);
    return std::make_tuple (1,n,weight);
}

template <class Obj>
inline auto KMesh::getValue(Obj &in, RealType x) const ->decltype(in[0]) 
{
    const auto find_result=this->find(x);
    if (!std::get<0>(find_result)) throw (exWrongIndex()); 
    return in[std::get<1>(find_result)];
}


template <class Obj>
inline auto KMesh::getValue(Obj &in, KMesh::point x) const ->decltype(in[0]) 
{
    if (x._index < _vals.size() && x == _vals[x._index])
    return in[x._index];
    else { 
        #ifndef NDEBUG
        ERROR ("Point not found"); 
        #endif
        return this->getValue(in, RealType(x)); 
         };
}

template <class Obj> 
inline auto KMesh::integrate(const Obj &in) const -> decltype(in(_vals[0]))
{
    decltype(in(_vals[0])) R = in(RealType(_vals[0]));
    R=std::accumulate(_vals.begin()+1, _vals.end(), R,[&](decltype(in(_vals[0]))& y,const decltype(_vals[0]) &x) {return y+in(x);}); 
    return R/_points;
}

template <class Obj, typename ...OtherArgTypes> 
inline auto KMesh::integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...))
{
    decltype(in(_vals[0],Args...)) R = in(_vals[0],Args...);
    R=std::accumulate(_vals.begin()+1, _vals.end(), R,[&](decltype(in(_vals[0]))& y,const decltype(_vals[0]) &x) {return y+in(x,Args...);}); 
    return R/_points;
}

template <class ArgType>
inline RealType KMesh::shift(RealType in, ArgType shift_arg) const
{
    assert (in>=0 && in < 2.0*PI);
    RealType out;
    out = in + RealType(shift_arg); 
    out-= std::floor(out/(2.0*PI))*2.0*PI;
    return out;
}


template <class ArgType>
inline typename KMesh::point KMesh::shift(point in, ArgType shift_arg) const
{
    point out;
    out._val = this->shift(in._val, shift_arg);
    auto find_result = this->find(out._val);
    if (!std::get<0>(find_result)) throw (exWrongIndex());
    out._index = std::get<1>(find_result);
    return (*this)[out._index];
}

//
// KMeshPatch
//


inline KMeshPatch::KMeshPatch(const KMesh& parent, std::vector<size_t> indices):
    _parent(parent),
    _npoints(indices.size())
{
    _vals.resize(_npoints); 
    for (size_t i=0; i<_npoints; ++i) {
        _vals[i]=_parent[indices[i]]; 
        }
}

inline KMeshPatch::KMeshPatch(const KMesh& parent):
    _parent(parent),
    _npoints(parent.getSize())
{
    _vals = parent.getVals();
}

} // end of namespace FK
#endif
