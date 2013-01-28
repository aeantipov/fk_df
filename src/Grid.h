#ifndef ___FK_GRID_H___
#define ___FK_GRID_H___

#include "FKCommon.h"
#include "Logger.h"

namespace FK { 

template <typename ValueType, typename ... > struct GridPointTypeExtractor;
template <typename ValueType, typename ... > struct GridPointExtractor;

/** A representatin of a one-dimensional grid, which stores an array of ValuType values. */
template <typename ValueType, class Derived>
class Grid {
public:
    struct point {
        ValueType _val;
        size_t _index;
        operator ValueType() const { return _val; }
        explicit inline operator size_t() const { return _index; }
        explicit inline operator int() const { return _index; }
        point(){};
        inline point(ValueType val, size_t index):_val(val),_index(index){};
        inline point(const point& rhs):_val(rhs._val),_index(rhs._index){};
        inline point(point&& rhs) { _val = rhs._val, _index = rhs._index; }
        inline point& operator=(point&& rhs) { _val = rhs._val, _index = rhs._index; return *this;}
        inline point operator=(const point& rhs) { _val = rhs._val, _index = rhs._index; return *this;}
        //inline point(ValueType in):_val(in){_index=2;};
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
    //typename Grid<ValueType, Derived>::iterator begin() const { return _vals.begin(); };
    //typename Grid<ValueType, Derived>::iterator end() const { return _vals.end(); };
    typename Grid<ValueType, Derived>::const_iterator begin() const { return _vals.begin(); };
    typename Grid<ValueType, Derived>::const_iterator end() const { return _vals.end(); };
    /** Empty constructor. */
    Grid();
    /** Copy from vector. */
    Grid(const std::vector<ValueType> & vals);
    /** Initialize the values from a given function, that maps the integer values
     * to the ValueType values. 
     */
    Grid(int min, int max, std::function<ValueType (int)> f);
    /** Copy constructor. */
    Grid(const Grid& rhs):_vals(rhs._vals){};
    /** Move constructor. */
    Grid(Grid&& rhs){_vals.swap(rhs._vals);};
    /** Returns a value at given index. */
    point operator[](size_t in) const;
    /** Returns all values. */
    const std::vector<point> & getVals() const;
    /** Returns size of grid. */
    size_t getSize() const;

    template <class ArgType>
        point shift(point in, ArgType shift_arg) const;
    template <class ArgType>
        ValueType shift(ValueType in, ArgType shift_arg) const;

    // CFTP forwards
    /** Get a value of an object at the given point, which is defined on a grid. */
    template <class Obj> auto getValue(Obj &in, Grid<ValueType,Derived>::point x) const ->decltype(in[0])
        { return static_cast<const Derived*>(this)->getValue(in,x); };
    /** Get a value of an object at the given coordinate, which is defined on a grid. */
    template <class Obj> auto getValue(Obj &in, ValueType x) const ->decltype(in[0])
        { return static_cast<const Derived*>(this)->getValue(in,x); };
    /** Returns a tuple of left closest index, weight, right closest index and weight, which are the closest to input value. */
    std::tuple <bool, size_t, RealType> find (ValueType in) const 
        { return static_cast<const Derived*>(this)->find(in); };
    /** Returns the closest point to the given value. */
    point findClosest(ValueType in) const;
    /** Integrate over grid. */
    template <class Obj> auto integrate(const Obj &in) const ->decltype(in[_vals[0]]) 
        { return static_cast<const Derived*>(this)->integrate(in); };
    /** Integrate over grid with extra arguments provided. */
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...))
        { return static_cast<const Derived*>(this)->integrate(in, Args...); };
    /** Make the object printable. */
    template <typename ValType, class Derived2> friend std::ostream& operator<<(std::ostream& lhs, const Grid<ValType,Derived2> &gr);

    class exWrongIndex : public std::exception { virtual const char* what() const throw(); }; 
};

/** A Grid of fermionic Matsubara frequencies. */
template <bool Fermion>
class MatsubaraGrid : public Grid<ComplexType, MatsubaraGrid<Fermion>>
{
public:
    using Grid<ComplexType, MatsubaraGrid<Fermion>>::_vals;
    using typename Grid<ComplexType, MatsubaraGrid<Fermion>>::exWrongIndex;
    //typedef typename Grid<ComplexType, MatsubaraGrid<Fermion>>::point point;
    using typename Grid<ComplexType, MatsubaraGrid<Fermion>>::point;
    /** Inverse temperature. */
    const RealType _beta;
    /** Spacing between values. */
    const RealType _spacing;
    /** Min and max numbers of freq. - useful for searching. */
    const int _w_min, _w_max;
    MatsubaraGrid(int min, int max, RealType beta);
    MatsubaraGrid(const MatsubaraGrid &rhs);
    MatsubaraGrid(MatsubaraGrid&& rhs);
    int getNumber(ComplexType in) const;
    std::tuple <bool, size_t, RealType> find (ComplexType in) const ;
    template <class Obj> auto integrate(const Obj &in) const -> decltype(in(_vals[0]));
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...));
    template <class Obj> auto prod(const Obj &in) const -> decltype(in(_vals[0]));
    //template <class Obj> auto gridIntegrate(const std::vector<Obj> &in) const -> Obj;
    template <class Obj> auto getValue(Obj &in, ComplexType x) const ->decltype(in[0]);
    template <class Obj> auto getValue(Obj &in, point x) const ->decltype(in[0]);
};


typedef MatsubaraGrid<1> FMatsubaraGrid;
typedef MatsubaraGrid<0> BMatsubaraGrid;

inline std::ostream& operator<<(std::ostream& lhs, const __num_format< typename FMatsubaraGrid::point> &in){lhs << std::setprecision(in._prec) << imag(in._v._val); return lhs;};
inline std::istream& operator>>(std::istream& lhs, __num_format<typename FMatsubaraGrid::point> &out){RealType im; lhs >> im; out._v._val = I*im; return lhs;};
inline std::ostream& operator<<(std::ostream& lhs, const __num_format< typename BMatsubaraGrid::point> &in){lhs << std::setprecision(in._prec) << imag(in._v._val); return lhs;};
inline std::istream& operator>>(std::istream& lhs, __num_format<typename BMatsubaraGrid::point> &out){RealType im; lhs >> im; out._v._val = I*im; return lhs;};
/** A grid of real values. */
class RealGrid : public Grid<RealType, RealGrid>
{
    RealType _min;
    RealType _max;
public:
    template <class Obj> auto integrate(const Obj &in) const ->decltype(in(_vals[0]));
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...));
    RealGrid(RealType min, RealType max, size_t npoints);
    RealGrid(int min, int max, const std::function<RealType(int)> &f);
    std::tuple <bool, size_t, RealType> find (RealType in) const ;
    //template <class Obj> auto gridIntegrate(std::vector<Obj> &in) -> Obj;
    template <class Obj> auto getValue(Obj &in, RealType x) const ->decltype(in[0]);
    template <class Obj> auto getValue(Obj &in, RealGrid::point x) const ->decltype(in[0]);
};

class KMesh : public Grid<RealType, KMesh>
{
public:
    int _points;
    KMesh(size_t n_points);
    KMesh(const KMesh& rhs);
    KMesh(KMesh &&rhs);
    KMesh(){};
    KMesh& operator=(KMesh &&rhs){_points = rhs._points; _vals.swap(rhs._vals); return (*this);};
    KMesh& operator=(const KMesh &rhs){_points = rhs._points; _vals = rhs._vals; return (*this);};
    std::tuple <bool, size_t, RealType> find (RealType in) const ;
    template <class Obj> auto integrate(const Obj &in) const ->decltype(in(_vals[0]));
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...));
    //template <class Obj> auto gridIntegrate(std::vector<Obj> &in) const -> Obj;
    template <class Obj> auto getValue(Obj &in, RealType x) const ->decltype(in[0]);
    template <class Obj> auto getValue(Obj &in, point x) const ->decltype(in[0]);
    template <class ArgType> point shift(point in, ArgType shift_arg) const;
    template <class ArgType> RealType shift(RealType in, ArgType shift_arg) const;
};

struct KMeshPatch : public KMesh 
{
    std::map<size_t,size_t> _map_vals;
public:
    const KMesh& _parent;
    size_t _npoints;
    using KMesh::_vals;
    KMeshPatch(const KMesh& parent, std::vector<size_t> indices);
    KMeshPatch(const KMesh& parent);
    template <class Obj> auto getValue(Obj &in, RealType x) const ->decltype(in[0]);
    template <class Obj> auto getValue(Obj &in, KMesh::point x) const ->decltype(in[0]);
    size_t getIndex(KMesh::point x) const;
};

/** A tool to generate a function of argtypes of grids. */
template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...GridTypes, typename ...ArgTypes>
struct GridPointTypeExtractor<ValueType, T<GridType1, GridTypes...>, ArgTypes...> : 
GridPointTypeExtractor<ValueType, T<GridTypes...>, ArgTypes...,decltype(GridType1::point::_val)>
{
};

template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...ArgTypes>
struct GridPointTypeExtractor<ValueType, T<GridType1>, ArgTypes...> {
    typedef std::function<ValueType(ArgTypes...,decltype(GridType1::point::_val))> type; 
    typedef std::function<ValueType(ArgTypes...,decltype(GridType1::point::_val))> arg_type; 
    typedef std::tuple<ArgTypes...,decltype(GridType1::point::_val)> arg_tuple_type;
};

/** A tool to generate a function of argtypes of grids. */
template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...GridTypes, typename ...ArgTypes>
struct GridPointExtractor<ValueType, T<GridType1, GridTypes...>, ArgTypes...> : 
GridPointExtractor<ValueType, T<GridTypes...>, ArgTypes...,typename GridType1::point>
{
};

template <typename ValueType, template <typename ...> class T, typename GridType1, typename ...ArgTypes>
struct GridPointExtractor<ValueType, T<GridType1>, ArgTypes...> {
    typedef std::function<ValueType(ArgTypes...,typename GridType1::point)> point_type; 
    typedef std::tuple<ArgTypes...,typename GridType1::point> arg_tuple_type;
};


/* A tool to generate an array of grid sizes from a given tuple of grids. */
template <size_t N>
struct GetGridSizes {
    template <typename... GridType, size_t M>
    static inline void TupleSizeToArray( const std::tuple<GridType...>& in, std::array<size_t, M> &out ) {
        static_assert(N>1,"!");
        std::get<N-1>(out) = std::get<N-1>(in).getSize();
        GetGridSizes<N-1>::TupleSizeToArray( in, out );
    }
};

template <>
template <typename... GridType, size_t M>
inline void GetGridSizes<1>::TupleSizeToArray( const std::tuple<GridType...>& in, std::array<size_t, M> &out ) {
    std::get<0>(out) = std::get<0>(in).getSize();
}


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
