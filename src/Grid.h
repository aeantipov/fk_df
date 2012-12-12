#ifndef ___FK_GRID_H___
#define ___FK_GRID_H___

#include "FKCommon.h"
#include "Logger.h"

namespace FK { 

/** A representatin of a one-dimensional grid, which stores an array of ValuType values. */
template <typename ValueType, class Derived>
class Grid {
public:
    struct point {
        ValueType _val;
        size_t _index;
        inline operator ValueType() const { return _val; }
        explicit inline operator size_t() const { return _index; }
        explicit inline operator int() const { return _index; }
        point(){};
        inline point(ValueType val, size_t index):_val(val),_index(index){};
        inline point(const point& rhs):_val(rhs._val),_index(rhs._index){};
        inline point(point&& rhs) { _val = rhs._val, _index = rhs._index; }
        inline point& operator=(point&& rhs) { _val = rhs._val, _index = rhs._index; return *this;}
        inline point operator=(const point& rhs) { _val = rhs._val, _index = rhs._index; return *this;}
        inline point(ValueType in):_val(in){_index=0;};
        bool operator==(const point &rhs) const {return (_val == rhs._val) && (_index == rhs._index);}
        friend std::ostream& operator<<(std::ostream& lhs, const point &gr)
            {lhs<<"{"<<gr._val<<"<-["<<gr._index<<"]}"; return lhs;};
    };
protected:
    std::vector<point> _vals;
    Grid(const std::vector<point> & vals);
public:
    typedef typename std::vector<ValueType>::iterator iterator;
    typedef typename std::vector<ValueType>::const_iterator const_iterator;
    typename Grid<ValueType, Derived>::iterator begin() const { return _vals.begin(); };
    typename Grid<ValueType, Derived>::iterator end() const { return _vals.end(); };
    /** Empty constructor. */
    Grid();
    /** Copy from vector. */
    Grid(const std::vector<ValueType> & vals);
    /** Initialize the values from a given function, that maps the integer values
     * to the ValueType values. 
     */
    Grid(int min, int max, std::function<ValueType (const int&)> f);
    /** Copy constructor. */
    Grid(const Grid& rhs):_vals(rhs._vals){};
    /** Move constructor. */
    Grid(Grid&& rhs){_vals.swap(rhs._vals);};
    /** Returns a value at given index. */
    ValueType operator[](size_t in) const;
    /** Returns all values. */
    const std::vector<point> & getVals() const;
    /** Returns size of grid. */
    size_t getSize() const;

    /** Returns a tuple of left closest index, weight, right closest index and weight, which are the closest to input value. */
    std::tuple <bool, size_t, RealType> find (ValueType in) const { return static_cast<Derived*>(this)->find(in); };
    /** A CRTP reference to one of the inherited objects. */
    template <class Obj> auto integrate(const Obj &in)->decltype(in[_vals[0]]) { return static_cast<Derived*>(this)->integrate(in); };
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...))
        { return static_cast<Derived*>(this)->integrate(in, Args...); };
    template <class Obj> auto getValue(const Obj &in, ValueType x)->decltype(in[0]) const { return static_cast<Derived*>(this)->get_val(in); };
    template <class Obj> auto getValue(const Obj &in, point x)->decltype(in[0]) const { return static_cast<Derived*>(this)->get_val(in); };
    /** Make the object printable. */
    template <typename ValType, class Derived2> friend std::ostream& operator<<(std::ostream& lhs, const Grid<ValType,Derived2> &gr);

    class exWrongIndex : public std::exception { virtual const char* what() const throw(); }; 
};

/** A Grid of Matsubara frequencies. */
class FMatsubaraGrid : public Grid<ComplexType, FMatsubaraGrid>
{
public:
    /** Inverse temperature. */
    const RealType _beta;
    /** Spacing between values. */
    const RealType _spacing;
    /** Min and max numbers of freq. - useful for searching. */
    const int _w_min, _w_max;
    FMatsubaraGrid(int min, int max, RealType beta);
    FMatsubaraGrid(const FMatsubaraGrid &rhs);
    FMatsubaraGrid(FMatsubaraGrid&& rhs);
    std::tuple <bool, size_t, RealType> find (ComplexType in) const ;
    template <class Obj> auto integrate(const Obj &in) const -> decltype(in(_vals[0]));
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...));
    template <class Obj> auto prod(const Obj &in) const -> decltype(in(_vals[0]));
    //template <class Obj> auto gridIntegrate(const std::vector<Obj> &in) const -> Obj;
    template <class Obj> auto getValue(Obj &in, ComplexType x) const ->decltype(in[0]);
    template <class Obj> auto getValue(Obj &in, point x) const ->decltype(in[0]);
};

/** A grid of real values. */
class RealGrid : public Grid<RealType, RealGrid>
{
    RealType _min;
    RealType _max;
public:
    template <class Obj> auto integrate(const Obj &in)->decltype(in(_vals[0]));
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...));
    RealGrid(RealType min, RealType max, size_t npoints);
    //template <class Obj> auto gridIntegrate(std::vector<Obj> &in) -> Obj;
    //std::tuple <size_t, RealType, size_t, RealType> find (ValueType in);
};

class KMesh : public Grid<RealType, KMesh>
{
public:
    const int _points;
    KMesh(size_t n_points);
    KMesh(const KMesh& rhs);
    KMesh(KMesh &&rhs);
    std::tuple <bool, size_t, RealType> find (RealType in) const ;
    template <class Obj> auto integrate(const Obj &in) const ->decltype(in(_vals[0]));
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...));
    //template <class Obj> auto gridIntegrate(std::vector<Obj> &in) const -> Obj;
    template <class Obj> auto getValue(Obj &in, RealType x) const ->decltype(in[0]);
    template <class Obj> auto getValue(Obj &in, point x) const ->decltype(in[0]);
};

/** A tool to recursiverly integrate over a grid. */
template <typename GridType, class Obj> struct RecursiveGridIntegrator;

/* Integrate a function over a grid. */
template <typename GridType, typename ValueType, typename ArgType1> 
struct RecursiveGridIntegrator<GridType, ValueType(ArgType1)>
{
    inline static ValueType integrate(const GridType& grid, const std::function<ValueType(ArgType1)>& in){ 
        //DEBUG("Using this method, N=1 with "<< sizeof...(OtherArgs) << " other args" );
        return grid.integrate(in);
    }
};

/** An alias for a std::function template type. */
template <typename GridType, typename ValueType, typename ArgType1> 
struct RecursiveGridIntegrator<GridType, std::function<ValueType(ArgType1)>>:RecursiveGridIntegrator<GridType, ValueType(ArgType1)>{};

/** Multi-dimensional integration over the same grid. */
template <typename GridType, typename ValueType, typename ArgType1, typename ...ArgTypes> 
struct RecursiveGridIntegrator<GridType, ValueType(ArgType1, ArgTypes...)>
{
    typedef std::function<ValueType(ArgTypes...)> type;
    inline static ValueType integrate(
        const GridType& grid, 
        const std::function<ValueType(ArgType1, ArgTypes...)>& in) 
     {
        //DEBUG("Using this method, N!=1");
        auto f1 = [&](const ArgType1& arg1)->ValueType { 
            //DEBUG(arg1);
            std::function<ValueType(ArgTypes...)> f0 = [&](const ArgTypes&... args){return in(arg1,args...);};
            return RecursiveGridIntegrator<GridType, ValueType(ArgTypes...)>::integrate(grid, f0);
            };
        //DEBUG(f1(0.0));
        return grid.integrate(f1);
    } 
};

/** An alias for a std::function template type. */
template <typename GridType, typename ValueType, typename ArgType1, typename ...ArgTypes> 
struct RecursiveGridIntegrator<GridType, std::function<ValueType(ArgType1, ArgTypes...)>>:
    RecursiveGridIntegrator<GridType, ValueType(ArgType1, ArgTypes...)> {};

} // end :: namespace FK
#include "Grid.hpp"

#endif // endin :: ifndef ___FK_GRID_H___
