#ifndef ___FK_CONTAINER_HPP___
#define ___FK_CONTAINER_HPP___

#include "Container.h"
#include <fstream>

namespace FK {


//
// Container, N!=1
//

template <size_t N, typename ValueType>
template <size_t M>
inline Container<N,ValueType>::Container ( const std::array<size_t, M> &in):
    _vals(std::vector<Container<N-1,ValueType> >(std::get<M-N>(in), Container<N-1,ValueType>(in)))
{
}

template <size_t N, typename ValueType>
inline auto Container<N,ValueType>::operator[](size_t i)->decltype(_vals[0])
{
    if (i>=_vals.size()) throw exWrongIndex();
    return _vals[i];
}

template <size_t N, typename ValueType> 
std::ostream& operator<<(std::ostream& lhs, const Container<N,ValueType> &in)
{
    lhs << "[";
    std::ostream_iterator<Container<N-1,ValueType> > out_it (lhs,", ");
    std::copy(&in._vals.data()[0],&in._vals.data()[in._vals.size()], out_it);
    lhs << "]";
    return lhs;
}

template <size_t N, typename ValueType> 
typename Container<N,ValueType>::iterator Container<N,ValueType>::begin()
{
    return _vals.begin();
}

template <size_t N, typename ValueType> 
typename Container<N,ValueType>::iterator Container<N,ValueType>::end()
{
    return _vals.end();
}

template <size_t N, typename ValueType> 
template <typename U, typename std::enable_if<std::is_same<U, ComplexType>::value, int>::type>
Container <N,ValueType> Container<N,ValueType>::conj()
{
    Container <N,ValueType> out(*this);
    for (iterator it1 = this->begin(); it1!=this->end(); it1++) {  
        *it1=it1->conj();
    }
    return out;
}

template <size_t N, typename ValueType> 
ValueType Container<N,ValueType>::sum()
{
    ValueType out=0.0;
    out = std::accumulate(_vals.begin(), _vals.end(), out, [](ValueType x, Container<N-1,ValueType> &in){return x+in.sum();});
    return out;
}


//
// Container, N=1
//

template <typename ValueType>
template <size_t M>
inline Container<1,ValueType>::Container ( const std::array<size_t, M> &in):_vals(std::vector<ValueType>(std::get<M-1>(in)))
{
//    DEBUG("Constructing from size_t.");
}

template <typename ValueType>
inline ValueType& Container<1,ValueType>::operator[] (size_t i)
{
    return _vals[i];
}

template <typename ValueType> 
std::ostream& operator<<(std::ostream& lhs, const Container<1,ValueType> &in)
{
    lhs << "[";
    std::ostream_iterator<ValueType> out_it (lhs,", ");
    std::copy(in._vals.begin(),in._vals.end(), out_it);
    lhs << "]";
    return lhs;
}

template <typename ValueType> 
typename Container<1,ValueType>::iterator Container<1,ValueType>::begin()
{
    return _vals.begin();
}

template <typename ValueType> 
typename Container<1,ValueType>::iterator Container<1,ValueType>::end()
{
    return _vals.end();
}


template <typename ValueType> 
template <typename U, typename std::enable_if<std::is_same<U, ComplexType>::value, int>::type>
Container <1,ValueType> Container<1,ValueType>::conj()
{
    Container <1,ValueType> out(*this);
    for_each(out._vals.begin(), out._vals.end(), [&](ComplexType &x){x=std::conj(x);});
    return out;
}

template <typename ValueType> 
ValueType Container<1,ValueType>::sum()
{
    ValueType out=0.0;
    out = std::accumulate(_vals.begin(), _vals.end(), out, std::plus<ValueType>());
    return out;
}

template <typename ValueType> 
void Container<1,ValueType>::savetxt(const std::string& fname)
{
    std::ofstream out;
    out.open(fname.c_str());
    std::ostream_iterator<ValueType> out_it (out,"\n");
    std::copy(_vals.begin(),_vals.end(), out_it);
    out.close();
}


//
// Algebraic operators
//

// Operator=

template <size_t N, typename ValueType> 
Container<N,ValueType>& Container<N,ValueType>::operator=(const Container<N,ValueType> &rhs)
{
    _vals = rhs._vals; 
    return (*this);
}

template <typename ValueType> 
Container<1,ValueType>& Container<1,ValueType>::operator=(const Container<1,ValueType> &rhs)
{
    _vals = rhs._vals; 
    return (*this);
}

template <size_t N, typename ValueType> 
Container<N,ValueType>& Container<N,ValueType>::operator=(Container<N,ValueType> &&rhs)
{
    _vals.swap(rhs._vals); 
    return (*this);
}

template <typename ValueType> 
Container<1,ValueType>& Container<1,ValueType>::operator=(Container<1,ValueType> &&rhs)
{
    _vals.swap(rhs._vals); 
    return (*this);
}



template <size_t N, typename ValueType> 
Container<N,ValueType>& Container<N,ValueType>::operator=(const ValueType &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](Container<N-1,ValueType> &x){x=rhs;});
    return (*this);
}

template <typename ValueType> 
Container<1,ValueType>& Container<1,ValueType>::operator=(const ValueType &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](ValueType &x){x=rhs;});
    return (*this);
}



// Operator+=
template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType>& Container<N,ValueType>::operator+=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](Container<N-1,ValueType> &x){x+=rhs;});
    return *this;
}

template <size_t N, typename ValueType> 
inline Container<N,ValueType>& Container<N,ValueType>::operator+=(const Container<N,ValueType> &rhs)
{
    assert(this->_vals.size() == rhs._vals.size());
    //std::transform(_vals.begin(), _vals.end(), rhs._vals.begin(), _vals.begin(), [](Container<N-1,ValueType>&x, const Container<N-1,ValueType>&y){x+=y;} );
  //  std::transform(_vals.begin(), _vals.end(), rhs._vals.begin(), [](Container<N-1,ValueType>&x, const Container<N-1,ValueType>&y){x+=y;} );
    //std::transform(_vals.begin(), _vals.end(), rhs._vals.begin(), _vals.begin(), std::plus<Container<N-1,ValueType> >());
    for (int i=0; i<_vals.size(); ++i) _vals[i]+=rhs._vals[i];
    return *this;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType>& Container<1,ValueType>::operator+=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](ValueType &x){x+=rhs;});
    return *this;
}

template <typename ValueType> 
inline Container<1,ValueType>& Container<1,ValueType>::operator+=(const Container<1,ValueType> &rhs)
{
    assert(this->_vals.size() == rhs._vals.size());
    std::transform(_vals.begin(), _vals.end(), rhs._vals.begin(), _vals.begin(), std::plus<ValueType>());
    return *this;
}

//Operator +

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType> Container<N,ValueType>::operator+(const RhsArg &rhs) const
{
    Container<N,ValueType> out(*this);
    out+=rhs;
    return out;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType> Container<1,ValueType>::operator+(const RhsArg &rhs) const
{
    Container<1,ValueType> out(*this);
    out+=rhs;
    return out;
}

//
// Operator*=
//
template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType>& Container<N,ValueType>::operator*=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](Container<N-1,ValueType> &x){x*=rhs;});
    return *this;
}

template <size_t N, typename ValueType> 
inline Container<N,ValueType>& Container<N,ValueType>::operator*=(const Container<N,ValueType> &rhs)
{
    assert(this->_vals.size() == rhs._vals.size());
    for (int i=0; i<_vals.size(); ++i) _vals[i]*=rhs._vals[i];
    return *this;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType>& Container<1,ValueType>::operator*=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](ValueType &x){x*=rhs;});
    return *this;
}

template <typename ValueType> 
inline Container<1,ValueType>& Container<1,ValueType>::operator*=(const Container<1,ValueType> &rhs)
{
    assert(this->_vals.size() == rhs._vals.size());
    std::transform(_vals.begin(), _vals.end(), rhs._vals.begin(), _vals.begin(), std::multiplies<ValueType>());
    return *this;
}


//
// Operator*
//

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType> Container<1,ValueType>::operator*(const RhsArg &rhs) const
{
    Container<1,ValueType> out(*this);
    out*=rhs;
    return out;
}

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType> Container<N,ValueType>::operator*(const RhsArg &rhs) const
{
    Container<N,ValueType> out(*this);
    out*=rhs;
    return out;
}

//
// Operator-=
//

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType>& Container<1,ValueType>::operator-=(const RhsArg &rhs)
{
    (*this)+=rhs*ValueType(-1);
    return *this;
}

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType>& Container<N,ValueType>::operator-=(const RhsArg &rhs)
{
    (*this)+=rhs*ValueType(-1);
    return *this;
}

//
// Operator-
//

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType> Container<N,ValueType>::operator-(const RhsArg &rhs) const
{
    Container<N,ValueType> out(*this);
    out-=rhs;
    return out;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType> Container<1,ValueType>::operator-(const RhsArg &rhs) const
{
    Container<1,ValueType> out(*this);
    out-=rhs;
    return out;
}

//
// Operator/=
//

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType>& Container<N,ValueType>::operator/=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](Container<N-1,ValueType> &x){x/=rhs;});
    return *this;
}

template <size_t N, typename ValueType> 
inline Container<N,ValueType>& Container<N,ValueType>::operator/=(const Container<N,ValueType> &rhs)
{
    assert(this->_vals.size() == rhs._vals.size());
    for (int i=0; i<_vals.size(); ++i) _vals[i]/=rhs._vals[i];
    return *this;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType>& Container<1,ValueType>::operator/=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](ValueType &x){x/=rhs;});
    return *this;
}

template <typename ValueType> 
inline Container<1,ValueType>& Container<1,ValueType>::operator/=(const Container<1,ValueType> &rhs)
{
    assert(this->_vals.size() == rhs._vals.size());
    std::transform(_vals.begin(), _vals.end(), rhs._vals.begin(), _vals.begin(), std::divides<ValueType>());
    return *this;
}

//
// Operator /
//

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType> Container<1,ValueType>::operator/(const RhsArg &rhs) const
{
    Container<1,ValueType> out(*this);
    out/=rhs;
    return out;
}

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType> Container<N,ValueType>::operator/(const RhsArg &rhs) const
{
    Container<N,ValueType> out(*this);
    out/=rhs;
    return out;
}




} // end of namespace FK
#endif
