#ifndef ___FK_CONTAINER_H___
#define ___FK_CONTAINER_H___

#include "FKCommon.h"
#include "Logger.h"
#include "Grid.h"

#include <type_traits>

namespace FK { 

/** This class is a N-dimensional container to store the data of
 * the ValueType type. */
template <size_t N, typename ValueType>
class Container {
static_assert(N>1, "N=1 has a specialized variant");
protected:
    /** The data is stored as a recursive vector<vector ... > > structure. */
    std::vector<Container<N-1, ValueType> >  _vals;
public:
    typedef typename std::vector<Container<N-1, ValueType>>::iterator iterator;
    typedef typename std::vector<Container<N-1, ValueType>>::const_iterator const_iterator;
    /** Constructor from the std::array of size_t. */
    template <size_t M> Container ( const std::array<size_t, M> &in);
    /** Copy constructor. */
    Container(const Container<N,ValueType> &rhs):_vals(rhs._vals){};
    /** Move constructor. */
    Container(Container<N,ValueType> &&rhs){std::swap(_vals, rhs._vals);};
    /** Returns the value at point i. */
    auto operator[](size_t i)->decltype(_vals[0]);
    auto operator[](size_t i) const ->decltype(_vals[0]) const;
    /** Begin iterator. */
    typename Container<N,ValueType>::iterator begin();
    /** End iterator. */
    typename Container<N,ValueType>::iterator end();
    /** Algebraic operators. */
    Container<N,ValueType>& operator=(Container<N,ValueType> &&rhs);
    Container<N,ValueType>& operator=(const Container<N,ValueType> &rhs);
    Container<N,ValueType>& operator=(const ValueType &rhs);
    template <typename RhsArg> Container<N,ValueType>& operator+=(const RhsArg &rhs); 
    template <typename RhsArg> Container<N,ValueType> operator+(const RhsArg &rhs) const; 
    template <typename RhsArg> Container<N,ValueType>& operator*=(const RhsArg &rhs); 
    template <typename RhsArg> Container<N,ValueType> operator*(const RhsArg &rhs) const; 
    template <typename RhsArg> Container<N,ValueType>& operator/=(const RhsArg &rhs); 
    template <typename RhsArg> Container<N,ValueType> operator/(const RhsArg &rhs) const; 
    template <typename RhsArg> Container<N,ValueType>& operator-=(const RhsArg &rhs); 
    template <typename RhsArg> Container<N,ValueType> operator-(const RhsArg &rhs) const; 
    Container<N,ValueType>& operator+=(const Container<N,ValueType> &rhs); 
    Container<N,ValueType>& operator*=(const Container<N,ValueType> &rhs); 
    Container<N,ValueType>& operator/=(const Container<N,ValueType> &rhs); 
    friend inline Container<N,ValueType> operator* (const ValueType & lhs, const Container<N,ValueType> & rhs) {return rhs*lhs;};
    friend inline Container<N,ValueType> operator+ (const ValueType & lhs, const Container<N,ValueType> & rhs) {return rhs+lhs;};
    friend inline Container<N,ValueType> operator- (const ValueType & lhs, const Container<N,ValueType> & rhs) {return rhs*(-1.0)+lhs;};
    friend inline Container<N,ValueType> operator/ (const ValueType & lhs, const Container<N,ValueType> & rhs) {Container<N,ValueType> out(rhs); out=lhs; return out/rhs;};


 
    /** Conjugate. */
    template <typename U = ValueType, typename std::enable_if<std::is_same<U, ComplexType>::value, int>::type=0> 
        Container<N, ValueType> conj();
    /** Recursively iterates and sums all values in the container. */
    ValueType sum();
    /** Make the object streamable. */
    template <size_t M, typename ValType> friend std::ostream& operator<<(std::ostream& lhs, const Container<M,ValType> &in);

    class exWrongIndex : public std::exception { virtual const char* what() const throw(){return "Index out of bounds";}}; 
};

/* A specialization for the 1d Container. */
template <typename ValueType>
class Container<1,ValueType> {
protected:
    std::vector<ValueType>  _vals;
public:
    typedef typename std::vector<ValueType>::iterator iterator;
    typedef typename std::vector<ValueType>::const_iterator const_iterator;
    //Container(){};
    /** Copy constructor. */
    Container(const Container<1,ValueType> &rhs):_vals(rhs._vals){};
    /** Move constructor. */
    Container(Container<1,ValueType> &&rhs){std::swap(_vals, rhs._vals);};
    /** Constructor from the std::array of size_t. The last value is used. */
    template <size_t M> Container ( const std::array<size_t, M> &in);
    /** Returns the value at index i. */
    ValueType& operator[](size_t i);
    /** Begin iterator. */
    typename Container<1,ValueType>::iterator begin();
    /** End iterator. */
    typename Container<1,ValueType>::iterator end();

    /** Algebraic operators. */
    Container<1,ValueType>& operator=(Container<1,ValueType> &&rhs);
    Container<1,ValueType>& operator=(const Container<1,ValueType> &rhs);
    Container<1,ValueType>& operator=(const ValueType &rhs);
    template <typename RhsArg> Container<1,ValueType>& operator+=(const RhsArg &rhs); 
    template <typename RhsArg> Container<1,ValueType> operator+(const RhsArg &rhs) const; 
    template <typename RhsArg> Container<1,ValueType>& operator*=(const RhsArg &rhs); 
    template <typename RhsArg> Container<1,ValueType> operator*(const RhsArg &rhs) const; 
    template <typename RhsArg> Container<1,ValueType>& operator/=(const RhsArg &rhs); 
    template <typename RhsArg> Container<1,ValueType> operator/(const RhsArg &rhs) const; 
    template <typename RhsArg> Container<1,ValueType>& operator-=(const RhsArg &rhs); 
    template <typename RhsArg> Container<1,ValueType> operator-(const RhsArg &rhs) const; 
    Container<1,ValueType>& operator+=(const Container<1,ValueType> &rhs); 
    Container<1,ValueType>& operator*=(const Container<1,ValueType> &rhs); 
    Container<1,ValueType>& operator/=(const Container<1,ValueType> &rhs); 
    friend inline Container<1,ValueType> operator* (const ValueType & lhs, const Container<1,ValueType> & rhs) {return rhs*lhs;};
    friend inline Container<1,ValueType> operator+ (const ValueType & lhs, const Container<1,ValueType> & rhs) {return rhs+lhs;};
    friend inline Container<1,ValueType> operator- (const ValueType & lhs, const Container<1,ValueType> & rhs) {return rhs*(-1.0)+lhs;};
    friend inline Container<1,ValueType> operator/ (const ValueType & lhs, const Container<1,ValueType> & rhs) {Container<1,ValueType> out(rhs); out=lhs; return out/rhs;};

    /** Conjugate. */
    template <typename U = ValueType, typename std::enable_if<std::is_same<U, ComplexType>::value, int>::type=0>
        Container<1,ValueType> conj();
    /** Returns the sum of all values in the container. */
    ValueType sum();
    /** Saves the data to a plain text file */
    void savetxt(const std::string& fname);
    /** Make the object streamable. */
    template <typename ValType> friend std::ostream& operator<<(std::ostream& lhs, const Container<1,ValType> &in);
};

}; // end of namespace FK
#endif

#include "Container.hpp"
