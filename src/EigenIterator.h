#ifndef ___INCLUDE_EIGENITERATOR_H___
#define ___INCLUDE_EIGENITERATOR_H___

#include<Eigen/Core>

namespace FK {
 
// An iterator help for Eigen from http://forum.kde.org/viewtopic.php?f=74&t=94120#p191501
template<typename Value_t, typename Container_t>
class index_iterator : public std::iterator<std::random_access_iterator_tag, Value_t>
{
protected:
   Container_t* container_;
   int          index_;

public:
   index_iterator() : container_(0), index_(0) { }
   index_iterator(Container_t& container, int index) : container_(&container), index_(index) { }

   bool operator==(const index_iterator& other) { return container_ == other.container_ && index_ == other.index_; }
   bool operator!=(const index_iterator& other) { return !(*this == other); }

   Value_t&       operator*()       { return (*container_)[index_]; }
   Value_t const& operator*() const { return (*container_)[index_]; }
   
   Value_t*       operator->()       { return &((*container_)[index_]); }
   Value_t const* operator->() const { return &((*container_)[index_]); }

   index_iterator& operator++() { ++index_; return *this;}
   index_iterator operator++(int) { index_iterator prev(*this); operator++(); return prev;}

   index_iterator& operator--() { --index_; return *this;}
   index_iterator operator--(int) { index_iterator prev(*this); operator--(); return prev;}
   
   friend index_iterator operator+(const index_iterator& a, int b) { index_iterator ret(a); ret += b; return ret; }
   friend index_iterator operator-(const index_iterator& a, int b) { index_iterator ret(a); ret -= b; return ret; }
   friend index_iterator operator+(int a, const index_iterator& b) { index_iterator ret(b); ret += a; return ret; }
   friend index_iterator operator-(int a, const index_iterator& b) { index_iterator ret(b); ret -= a; return ret; }

   int operator-(const index_iterator& other) const { return index_ - other.index_; }

   bool operator< (const index_iterator& other) { return container_ == other.container_ && index_ <  other.index_; }
   bool operator<=(const index_iterator& other) { return container_ == other.container_ && index_ <= other.index_; }
   bool operator> (const index_iterator& other) { return container_ == other.container_ && index_ >  other.index_; }
   bool operator>=(const index_iterator& other) { return container_ == other.container_ && index_ >= other.index_; }

   index_iterator& operator+=(int b) { index_ += b; }
   index_iterator& operator-=(int b) { index_ -= b; }

   Value_t&       operator[](int i)       { return (*container_)[i]; }
   Value_t const& operator[](int i) const { return (*container_)[i]; }
};
template<typename Value_t, typename Container_t>
inline index_iterator<Value_t, Container_t> index_begin(Container_t& container) 
{ 
   return index_iterator<Value_t, Container_t>(container, 0); 
}
template<typename Value_t, typename Container_t>
inline index_iterator<Value_t, Container_t> index_end(Container_t& container) 
{ 
   return index_iterator<Value_t, Container_t>(container, container.size()); 
}
template<typename Value_t, typename Container_t>
inline index_iterator<const Value_t, const Container_t> index_begin(const Container_t& container) 
{ 
   return index_iterator<const Value_t, const Container_t>(container, 0); 
}
template<typename Value_t, typename Container_t>
inline index_iterator<const Value_t, const Container_t> index_end(const Container_t& container) 
{ 
   return index_iterator<const Value_t, const Container_t>(container, container.size()); 
}

} // end namespace FK
#endif // endif :: ___INCLUDE_EIGENITERATOR_H___
