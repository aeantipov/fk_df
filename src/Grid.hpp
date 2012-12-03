#ifndef ___FK_GRID_HPP___
#define ___FK_GRID_HPP___

#include "Grid.h"

namespace FK { 

//
// Grid
//

template <typename ValueType, class Derived>
Grid<ValueType,Derived>::Grid()
{};

template <typename ValueType, class Derived>
Grid<ValueType,Derived>::Grid(const std::vector<ValueType> &vals):_vals(vals)
{
};

template <typename ValueType, class Derived>
Grid<ValueType,Derived>::Grid(int min, int max, std::function<ValueType (const int&)> f)
{
    if (max<min) std::swap(min,max);
    unsigned int n_points = max-min;
    _vals.resize(n_points); 
    for (int i=0; i<n_points; ++i) _vals[i]=f(min+i); 
}

template <typename ValueType, class Derived>
ValueType Grid<ValueType,Derived>::operator[](unsigned int index) const
{
    if (index>_vals.size()) throw exWrongIndex();
    return _vals[index];
}

template <typename ValueType, class Derived>
const std::vector<ValueType> & Grid<ValueType,Derived>::getVals() const
{
    return _vals;
}

template <typename ValueType, class Derived>
unsigned int Grid<ValueType,Derived>::getSize() const
{
    return _vals.size();
}

template <typename ValueType, class Derived>
std::ostream& operator<<(std::ostream& lhs, const Grid<ValueType,Derived> &gr)
{ 
    lhs << "{";
    //lhs << gr._vals;
    std::ostream_iterator<ValueType> out_it (lhs,", ");
    std::copy(index_begin<ValueType>(gr._vals),index_end<ValueType>(gr._vals) , out_it);
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

inline FMatsubaraGrid::FMatsubaraGrid(int min, int max, RealType beta):
    Grid(min,max,std::bind(FMatsubara, std::placeholders::_1, beta)),
    _beta(beta), _spacing(PI/beta), _w_min(min), _w_max(max)
{
}

template <class Obj> 
auto FMatsubaraGrid::integrate(const Obj &in) const -> decltype(in(_vals[0]))
{
    decltype(in(_vals[0])) R = in(_vals[0]);
    R=std::accumulate(_vals.begin()+1, _vals.end(), R,[&](decltype(in(_vals[0]))& y,decltype(_vals[0]) &x) {return y+in(x);}); 
    return R/_beta;
}

template <class Obj> 
auto FMatsubaraGrid::gridIntegrate(const std::vector<Obj> &in) const -> Obj
{
    decltype(in[0]) R = in[0];
    R=std::accumulate(_vals.begin()+1, _vals.end(),R,[&](decltype(in[0])& y, decltype(in[0]) &x) {return y+x;}); 
    return R/_beta;
}



inline std::tuple <bool, unsigned int, RealType> FMatsubaraGrid::find (ComplexType in) const
{
    assert (std::abs(real(in))<std::numeric_limits<RealType>::epsilon());
    int n=(imag(in)/_spacing-1)/2;
    if (n<_w_min) { ERROR("Frequency to find is out of bounds, " << in << "<" << FMatsubara(_w_min,_beta)); return std::make_tuple(0,0,0); };
    if (n>_w_max) { ERROR("Frequency to find is out of bounds, " << in << ">" << FMatsubara(_w_max,_beta)); return std::make_tuple(0,_vals.size(),0); };
    return std::make_tuple (1,n,1);
}


template <class Obj>
inline auto FMatsubaraGrid::getValue(Obj &in, ComplexType x) const ->decltype(in[0]) 
{
    const auto find_result=this->find(x);
    if (!std::get<0>(find_result)) throw (exWrongIndex()); 
    return in[std::get<1>(find_result)];
}
//
// RealGrid
//


template <class Obj> 
auto RealGrid::integrate(const Obj &in) -> decltype(in(_vals[0]))
{
    decltype(in(_vals[0])) R;
    for (int i=0; i<_vals.size()-1; ++i) {
        R+=0.5*(in(_vals[i])+in(_vals[i+1]))*(_vals[i+1]-_vals[i]);
        }
    return R;
}

} // end of namespace FK
#endif
