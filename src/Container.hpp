#ifndef ___FK_CONTAINER_HPP___
#define ___FK_CONTAINER_HPP___

#include "Container.h"

namespace FK {

//
// Container, N=1
//

template <typename ValueType>
template <size_t M>
inline Container<1,ValueType>::Container ( const std::array<size_t, M> &in):_vals(VectorType<ValueType>(std::get<M-1>(in)))
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
//    std::copy(index_begin<ValueType>(in._vals),index_end<ValueType>(in._vals) , out_it);
    std::copy(&in._vals.data()[0],&in._vals.data()[in._vals.size()] , out_it);
    lhs << "]";
    return lhs;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType>& Container<1,ValueType>::operator+=(const RhsArg &rhs)
{
//    assert(this->vals.size() == rhs.vals());
    //for (int i=0; i<_vals.size(); ++i) _vals[i]+=rhs;
    this->_vals+=VectorType<ValueType>::Constant(_vals.size(), rhs);
    return *this;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType> Container<1,ValueType>::operator+(const RhsArg &rhs) const
{
    Container<1,ValueType> out(*this);
    out+=rhs;
    return out;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType>& Container<1,ValueType>::operator*=(const RhsArg &rhs)
{
//    assert(this->vals.size() == rhs.vals());
    //this->_vals.noalias()*=rhs;
    this->_vals*=rhs;
    return *this;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType> Container<1,ValueType>::operator*(const RhsArg &rhs) const
{
    Container<1,ValueType> out(*this);
    out._vals*=rhs;
    return out;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType>& Container<1,ValueType>::operator-=(const RhsArg &rhs)
{
    (*this)*=rhs*(-1);
    return *this;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType> Container<1,ValueType>::operator-(const RhsArg &rhs) const
{
    Container<1,ValueType> out(*this);
    out._vals.noalias()-=rhs;
    return out;
}

template <typename ValueType> 
inline Container<1,ValueType>& Container<1,ValueType>::operator+=(const Container<1,ValueType> &rhs)
{
    assert(this->vals.size() == rhs.vals());
    this->_vals+=rhs._vals;
    return *this;
}

template <typename ValueType> 
inline Container<1,ValueType> Container<1,ValueType>::operator+(const Container<1,ValueType> &rhs) const
{
    Container<1,ValueType> out(*this);
    out._vals+=rhs._vals;
    return out;
}

template <typename ValueType> 
inline Container<1,ValueType>& Container<1,ValueType>::operator*=(const Container<1,ValueType> &rhs)
{
    assert(this->vals.size() == rhs.vals());
    this->_vals*=rhs._vals;
    return *this;
}

template <typename ValueType> 
inline Container<1,ValueType> Container<1,ValueType>::operator*(const Container<1,ValueType> &rhs) const
{
    Container<1,ValueType> out(*this);
    out._vals*=rhs._vals;
    return out;
}

template <typename ValueType> 
inline Container<1,ValueType>& Container<1,ValueType>::operator-=(const Container<1,ValueType> &rhs)
{
    assert(this->vals.size() == rhs.vals());
    (*this)*=rhs*(-1);
    return *this;
}

template <typename ValueType> 
inline Container<1,ValueType> Container<1,ValueType>::operator-(const Container<1,ValueType> &rhs) const
{
    Container<1,ValueType> out(*this);
    out._vals-=rhs._vals;
    return out;
}



//
// Container, N!=1
//

template <size_t N, typename ValueType>
template <size_t M>
inline Container<N,ValueType>::Container ( const std::array<size_t, M> &in):
    _vals(VectorType<Container<N-1,ValueType> >::Constant(std::get<M-N>(in), Container<N-1,ValueType>(in)))
{
}



template <size_t N, typename ValueType>
inline auto Container<N,ValueType>::operator[](size_t i)->decltype(_vals[0])
{
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
template <typename RhsArg>
inline Container<N,ValueType>& Container<N,ValueType>::operator+=(const RhsArg &rhs)
{
    //for (int i=0; i<_vals.size(); ++i) _vals[i]+=rhs;
    std::for_each(&_vals.data()[0], &_vals.data()[_vals.size()], [&](Container<N-1,ValueType> &x){x+=rhs;});
    return *this;
}

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType> Container<N,ValueType>::operator+(const RhsArg &rhs) const
{
    Container<N,ValueType> out(*this);
    out+=rhs;
    return out;
}

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType>& Container<N,ValueType>::operator*=(const RhsArg &rhs)
{
    //for (int i=0; i<_vals.size(); ++i) _vals[i]*=rhs;
    std::for_each(&_vals.data()[0], &_vals.data()[_vals.size()], [&](Container<N-1,ValueType> &x){x*=rhs;});
    //_vals=_vals*rhs;
    return *this;
}

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType> Container<N,ValueType>::operator*(const RhsArg &rhs) const
{
    Container<N,ValueType> out(*this);
    out*=rhs;
    return out;
}

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType>& Container<N,ValueType>::operator-=(const RhsArg &rhs)
{
    (*this)+=rhs*(-1);
    return *this;
}

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType> Container<N,ValueType>::operator-(const RhsArg &rhs) const
{
    Container<N,ValueType> out(*this);
    out-=rhs;
    return out;
}

template <size_t N, typename ValueType> 
inline Container<N,ValueType>& Container<N,ValueType>::operator+=(const Container<N,ValueType> &rhs)
{
    assert(this->vals.size() == rhs.vals());
    this->_vals+=rhs._vals;
    return *this;
}

template <size_t N, typename ValueType> 
inline Container<N,ValueType> Container<N,ValueType>::operator+(const Container<N,ValueType> &rhs) const
{
    Container<N,ValueType> out(*this);
    out._vals+=rhs._vals;
    return out;
}


template <size_t N, typename ValueType> 
inline Container<N,ValueType>& Container<N,ValueType>::operator*=(const Container<N,ValueType> &rhs)
{
    assert(this->vals.size() == rhs.vals());
    this->_vals*=rhs._vals;
    return *this;
}

template <size_t N, typename ValueType> 
inline Container<N,ValueType> Container<N,ValueType>::operator*(const Container<N,ValueType> &rhs) const
{
    Container<N,ValueType> out(*this);
    out._vals*=rhs._vals;
    return out;
}

template <size_t N, typename ValueType> 
inline Container<N,ValueType>& Container<N,ValueType>::operator-=(const Container<N,ValueType> &rhs)
{
    (*this)+=rhs*(-1);
    return *this;
}

template <size_t N, typename ValueType> 
inline Container<N,ValueType> Container<N,ValueType>::operator-(const Container<N,ValueType> &rhs) const
{
    Container<N,ValueType> out(*this);
    out._vals-=rhs._vals;
    return out;
}

} // end of namespace FK
#endif
