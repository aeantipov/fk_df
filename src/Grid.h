#ifndef ___FK_GRID_H___
#define ___FK_GRID_H___

#include "FKCommon.h"
#include "Logger.h"

namespace FK { 

/** A representatin of a one-dimensional grid, which stores an array of ValuType values. */
template <typename ValueType, class Derived>
class Grid {
protected:
    VectorType<ValueType> _vals;
public:
    /** Empty constructor. */
    Grid();
    /** Copy from vector. */
    Grid(const VectorType<ValueType> & vals);
    /** Initialize the values from a given function, that maps the integer values
     * to the ValueType values. 
     */
    Grid(int min, int max, std::function<ValueType (const int&)> f);
    /** Returns a value at given index. */
    ValueType operator[](unsigned int in) const;
    /** Returns all values. */
    const VectorType<ValueType> & getVals() const;
    /** Returns size of grid. */
    unsigned int getSize() const;

    /** Returns a tuple of left closest index, weight, right closest index and weight, which are the closest to input value. */
    std::tuple <bool, unsigned int, RealType> find (ValueType in) { return static_cast<Derived*>(this)->find(in); };

    /** A CRTP reference to one of the inherited objects. */
    template <class Obj> auto integrate(const Obj &in)->decltype(in[_vals[0]]) { return static_cast<Derived*>(this)->integrate(in); };
    /** Make the object printable. */
    template <typename ValType, class Derived2> friend std::ostream& operator<<(std::ostream& lhs, const Grid<ValType,Derived2> &gr);

    class exWrongIndex : public std::exception { virtual const char* what() const throw(); }; 
};

/** A particular Grid for Matsubara frequencies. */
class FMatsubaraGrid : public Grid<ComplexType, FMatsubaraGrid>
{
    /** Inverse temperature. */
    const RealType _beta;
    /** Spacing between values. */
    const RealType _spacing;
    /** Min and max numbers of freq. - useful for searching. */
    int _w_min, _w_max;
public:
    FMatsubaraGrid(int min, int max, RealType beta);
    std::tuple <bool, unsigned int, RealType> find (ComplexType in);
    template <class Obj> auto integrate(const Obj &in) -> decltype(in(_vals[0]));
    template <class Obj> auto gridIntegrate(VectorType<Obj> &in) -> Obj;
};

/** A grid of real frequencies. */
class RealGrid : public Grid<RealType, RealGrid>
{
    public:
    template <class Obj> auto integrate(const Obj &in)->decltype(in(_vals[0]));
    template <class Obj> auto gridIntegrate(VectorType<Obj> &in) -> Obj;
    //std::tuple <unsigned int, RealType, unsigned int, RealType> find (ValueType in);
};


/* ===========================================================================*/

//
// Grid
//

template <typename ValueType, class Derived>
Grid<ValueType,Derived>::Grid()
{};

template <typename ValueType, class Derived>
Grid<ValueType,Derived>::Grid(const VectorType<ValueType> &vals):_vals(vals)
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
const VectorType<ValueType> & Grid<ValueType,Derived>::getVals() const
{
    return _vals;
}

template <typename ValueType, class Derived>
unsigned int Grid<ValueType,Derived>::getSize() const
{
    return _vals.size();
}
/*
template <typename ValueType, class Derived>
void Grid<ValueType,Derived>::defineMap()
{
    for (unsigned int i=0; i<_vals.size(); ++i) _back_map[_vals[i]]=i;
}
*/

template <typename ValueType, class Derived>
std::ostream& operator<<(std::ostream& lhs, const Grid<ValueType,Derived> &gr)
{ 
    //lhs << gr._vals;
    std::ostream_iterator<ValueType> out_it (lhs,", ");
    std::copy(index_begin<ValueType>(gr._vals),index_end<ValueType>(gr._vals) , out_it);
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
auto FMatsubaraGrid::integrate(const Obj &in) -> decltype(in(_vals[0]))
{
    decltype(in(_vals[0])) R = in(_vals[0]);
    R=std::accumulate(++(index_begin<ComplexType>(_vals)), index_end<ComplexType>(_vals),R,[&](decltype(in(_vals[0]))& y,decltype(_vals[0]) &x) {return y+in(x);}); 
    return R/_beta;
}

template <class Obj> 
auto FMatsubaraGrid::gridIntegrate(VectorType<Obj> &in) -> Obj
{
    decltype(in[0]) R = in[0];
    R=std::accumulate(++(index_begin<ComplexType>(_vals)), index_end<ComplexType>(_vals),R,[&](decltype(in[0])& y, decltype(in[0]) &x) {return y+x;}); 
    return R/_beta;
    //decltype(in[0]) R = in[0];
    //R=std::accumulate(_vals.begin()+1, _vals.end(),R,[&](decltype(in[0])& y, decltype(_vals[0]) &x) {return y+in(x);}); 
    //return R/_beta;
}



inline std::tuple <bool, unsigned int, RealType> FMatsubaraGrid::find (ComplexType in)
{
    assert (std::abs(real(in))<std::numeric_limits<RealType>::epsilon());
    int n=(imag(in)/_spacing-1)/2;
    if (n<_w_min) return std::make_tuple(0,0,0);
    if (n>_w_max) return std::make_tuple(0,_vals.size(),0);
    else return std::make_tuple (1,n,1);
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

} // end :: namespace FK

#endif // endin :: ifndef ___FK_GRID_H___
