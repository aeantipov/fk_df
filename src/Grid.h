#ifndef ___FK_GRID_H___
#define ___FK_GRID_H___

#include "FKCommon.h"
#include "Logger.h"

namespace FK { 

/** A representatin of a one-dimensional grid, which stores an array of ValuType values. */
template <typename ValueType, class Derived>
class Grid {
protected:
    std::vector<ValueType> _vals;
public:
    /** Empty constructor. */
    Grid();
    /** Copy from vector. */
    Grid(const std::vector<ValueType> & vals);
    /** Initialize the values from a given function, that maps the integer values
     * to the ValueType values. 
     */
    Grid(int min, int max, std::function<ValueType (const int&)> f);
    /** Returns a value at given index. */
    ValueType operator[](unsigned int in) const;
    /** Returns all values. */
    const std::vector<ValueType> & getVals() const;
    /** Returns size of grid. */
    unsigned int getSize() const;

    /** Returns a tuple of left closest index, weight, right closest index and weight, which are the closest to input value. */
    std::tuple <bool, unsigned int, RealType> find (ValueType in) { return static_cast<Derived*>(this)->find(in); };
    /** A CRTP reference to one of the inherited objects. */
    template <class Obj> auto integrate(const Obj &in)->decltype(in[_vals[0]]) { return static_cast<Derived*>(this)->integrate(in); };
    template <class Obj> auto getValue(const Obj &in, ComplexType x)->decltype(in[_vals[0]]) { return static_cast<Derived*>(this)->get_val(in); };
    /** Make the object printable. */
    template <typename ValType, class Derived2> friend std::ostream& operator<<(std::ostream& lhs, const Grid<ValType,Derived2> &gr);

    class exWrongIndex : public std::exception { virtual const char* what() const throw(); }; 
};

/** A Grid of Matsubara frequencies. */
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
    template <class Obj> auto gridIntegrate(std::vector<Obj> &in) -> Obj;
    template <class Obj> auto getValue(const Obj &in, ComplexType x)->decltype(in[_vals[0]]);
};

/** A grid of real values. */
class RealGrid : public Grid<RealType, RealGrid>
{
    public:
    template <class Obj> auto integrate(const Obj &in)->decltype(in(_vals[0]));
    template <class Obj> auto gridIntegrate(std::vector<Obj> &in) -> Obj;
    //std::tuple <unsigned int, RealType, unsigned int, RealType> find (ValueType in);
};

} // end :: namespace FK

#include "Grid.hpp"

#endif // endin :: ifndef ___FK_GRID_H___
