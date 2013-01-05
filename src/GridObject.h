#ifndef ___FK_GRID_OBJECT_H___
#define ___FK_GRID_OBJECT_H___

#include "FKCommon.h"
#include "Logger.h"
#include "Grid.h"
#include "Container.h"

namespace FK {

/** A GridObject is a wrapper over Container class, that stores data,
 * defined on multiple different grids.
 */
template< typename ValueType, typename ...GridTypes> 
class GridObject 
{
public:
    /** A typedef for a function that gives the analytical value of the object, when it's not stored. */
    typedef typename GridPointTypeExtractor<ValueType, std::tuple<GridTypes...> >::arg_type FunctionType;
    /** A typedef for a function that gives the analytical value of the object, when it's not stored. */
    typedef typename GridPointExtractor<ValueType, std::tuple<GridTypes...> >::point_type PointFunctionType;
    /** A typedef for a tuple of grids. */
    typedef std::tuple<GridTypes...> GridTupleType;
    /** A typedef for a tuple of grid points. */
    typedef typename GridPointTypeExtractor<ValueType, std::tuple<GridTypes...> >::arg_tuple_type PointTupleType;
protected:
    static const size_t N = sizeof...(GridTypes);
    /** Grids on which the data is defined. */
    const std::tuple<GridTypes...> _grids;
    /** The dimensions of the Container - deduced from grids. */
    std::array<size_t, N> _dims;
    /** A pointer to the Container. A pointer is used as there exist no default 
     * constructor for the Container.
     */
    std::unique_ptr<Container<N, ValueType>> _data;

    /** A helper recursive template utility to extract and set data from the container. */
    template <size_t Nc, typename ArgType1, typename ...ArgTypes> struct ContainerExtractor {
        /** Gets the data by values. */
        static ValueType& get(Container<Nc, ValueType> &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1, const ArgTypes&... args);
        /** Fills the container from function
         * \param[in] data Container to fill
         * \param[in] grids Grids, on which the data is defined. 
         * \param[in] f A function that defines the data. 
         */
        static void set(Container<Nc, ValueType> &data, const std::tuple<GridTypes...> &grids, const std::function<ValueType(ArgType1, ArgTypes...)> &f);
             };
    /** Specialization of ContainerExtractor for 1-dim container. */
    template <typename ArgType1> struct ContainerExtractor<1,ArgType1> {
        static ValueType& get(Container<1, ValueType> &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1); 
        static void set(Container<1, ValueType> &data, const std::tuple<GridTypes...> &grids, const std::function<ValueType(ArgType1)> &f);
        };
public:
    /** This function returns the value of the object when the point is not in container. */
    FunctionType _f;
    /** Constructs a grid object out of a tuple containing various grids. */
    GridObject( const std::tuple<GridTypes...> &grids);
    /** Constructor of grids and data. */
    GridObject( const std::tuple<GridTypes...> &grids, const Container<sizeof...(GridTypes), ValueType>& data):
        _grids(grids),
        _data(new Container<sizeof...(GridTypes), ValueType>>(data)),
        _f(__fun_traits<FunctionType>::constant(0.0)) {};
    /** Copy constructor. */
    GridObject( const GridObject<ValueType, GridTypes...>& rhs);
    /** Move constructor. */
    GridObject( GridObject<ValueType, GridTypes...>&& rhs);

    const std::tuple<GridTypes...> getGrids() const;
    /** Returns an Mth grid in _grids. */
    template<size_t M> auto getGrid() const -> const decltype(std::get<M>(_grids));
    /** Returns the top level grid. */
    auto getGrid() const -> const decltype(std::get<0>(_grids));
    /** Returns element number i, which corresponds to (*_grid)[i]. */
    auto operator[](size_t i)->decltype((*_data)[i]);
    /** Const operator[]. */
    auto operator[](size_t i) const ->decltype((*_data)[i]) const;
    //template <size_t M> ValueType& operator[](const std::array<size_t,M>& in);
    /** Returns the _data Container. */
    Container<sizeof...(GridTypes), ValueType>& getData(){return *_data;};
    /** Fills the Container with a provided function. */
    template <typename ...ArgTypes> void fill(const std::function<ValueType(ArgTypes...)> &);
    void fill(const FunctionType &in);
    template <typename ...ArgTypes> void fill_tuple(const std::function<ValueType(const std::tuple<ArgTypes...>)> &);
    /** Fills the Container with any proper class with call operator. Untested */
    //template <template <typename, class> class Filler, typename ...ArgTypes> void fill(const Filler<ValueType,ArgTypes...> &);

    /** Return the value by grid values. */
    template <typename ...ArgTypes> ValueType& get(const ArgTypes&... in);
    template <typename ...ArgTypes> ValueType operator()(const ArgTypes&... in) const;
    template <typename ...ArgTypes> ValueType operator()(const std::tuple<ArgTypes...>& in) const;
    //template <typename ...ArgTypes> auto operator()(const ArgType1& in)->decltype() const;

    /** A shortcut for fill method. */
    //template <typename ...ArgTypes> GridObject& operator= (const std::function<ValueType(ArgTypes...)> &);
    /** Algebraic operators. */
    GridObject& operator= (const GridObject & rhs);
    GridObject& operator= (const ValueType & rhs);
    GridObject& operator*= (const GridObject & rhs);
    GridObject& operator*= (const ValueType& rhs);
    GridObject operator* (const GridObject & rhs) const;
    GridObject operator* (const ValueType & rhs) const;
    GridObject& operator+= (const GridObject & rhs);
    GridObject& operator+= (const ValueType& rhs);
    GridObject operator+ (const GridObject & rhs) const;
    GridObject operator+ (const ValueType & rhs) const;
    GridObject& operator-= (const GridObject & rhs);
    GridObject& operator-= (const ValueType& rhs);
    GridObject operator- (const GridObject & rhs) const;
    GridObject operator- (const ValueType & rhs) const;
    GridObject& operator/= (const GridObject & rhs);
    GridObject& operator/= (const ValueType& rhs);
    GridObject operator/ (const GridObject & rhs) const;
    GridObject operator/ (const ValueType & rhs) const;
    friend inline GridObject operator* (const ValueType & lhs, const GridObject & rhs) {return rhs*lhs;};
    friend inline GridObject operator+ (const ValueType & lhs, const GridObject & rhs) {return rhs+lhs;};
    friend inline GridObject operator- (const ValueType & lhs, const GridObject & rhs) {return rhs*(-1.0)+lhs;};
    friend inline GridObject operator/ (const ValueType & lhs, const GridObject & rhs) {GridObject out(rhs); out=lhs; return out/rhs;};

    /** Returns the complex conjugate of this object, if it's complex valued. */
    template <typename U = ValueType, typename std::enable_if<std::is_same<U, ComplexType>::value, int>::type=0>
        GridObject conj();
    /** Returns the sum of all elements in the container. */
    ValueType sum();
    /** Save the data to the txt file. */
    void savetxt(const std::string& fname) const;
    /** Loads the data to the txt file. */
    void loadtxt(const std::string& fname);
    /** Dumps the object to the stream. */
    template <typename ValType, class ...GridTypes2> friend std::ostream& operator<<(std::ostream& lhs, const GridObject<ValType,GridTypes2...> &in);
    
    class exIOProblem : public std::exception { virtual const char* what() const throw(){return "IO problem.";} }; 
};

} // end of namespace FK
#endif // endif::ifndef ___FK_GRID_OBJECT_H___

#include "GridObject.hpp"

