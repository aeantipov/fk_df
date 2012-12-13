#include "GFWrap.h"

#include <cmath>

namespace FK {

GFWrap::GFWrap(const FMatsubaraGrid& in, size_t tail_size):
    GridObject<ComplexType, FMatsubaraGrid>(in), 
    _tail_coeffs(std::vector<ComplexType>(tail_size,0.0))
{
    assert(tail_size>2);
    _tail_coeffs[2]=1.0; // 1/w tail
}

GFWrap::GFWrap(const std::tuple<FMatsubaraGrid> & in, size_t tail_size):GFWrap(std::get<0>(in), tail_size)
{
}

GFWrap::GFWrap(GFWrap&& in):GridObject<ComplexType,FMatsubaraGrid>(in._grids)
{
    _data.swap(in._data);
    _tail_coeffs.swap(in._tail_coeffs);
}

GFWrap::GFWrap(const GFWrap& in):GridObject<ComplexType,FMatsubaraGrid>(in),_tail_coeffs(in._tail_coeffs)
{
}

GFWrap::GFWrap(const GridObject<ComplexType,FMatsubaraGrid>& in):GridObject<ComplexType,FMatsubaraGrid>(in)
{
    _tail_coeffs = std::vector<ComplexType>(5,0.0);
}



GFWrap::GFWrap(GridObject<ComplexType,FMatsubaraGrid>&& in):GridObject<ComplexType,FMatsubaraGrid>(in)
{
    _tail_coeffs = std::vector<ComplexType>(5,0.0);
}

ComplexType GFWrap::operator()(const ComplexType &in) const
{
    assert(std::abs(real(in))<std::numeric_limits<RealType>::epsilon());
    if (imag(in)<0) return -(*this)(std::conj(in));
    if (imag(in)>=imag(FMatsubara(std::get<0>(_grids)._w_max, std::get<0>(_grids)._beta))) return interp(in);
    return GridObject<ComplexType, FMatsubaraGrid>::operator()(in);
}

ComplexType GFWrap::interp(ComplexType in) const
{
    ComplexType out = 0;
    for (int i=0; i<_tail_coeffs.size(); ++i) out+=_tail_coeffs[i]*pow(in,1-i);
    return out;
}

GFWrap& GFWrap::operator=(GFWrap&& in)
{
    _data.swap(in._data);
    assert(std::get<0>(_grids).getVals() == std::get<0>(in._grids).getVals());
    _tail_coeffs.swap(in._tail_coeffs);
    return (*this);
}

GFWrap& GFWrap::operator=(const GFWrap& in)
{
    (*_data) = *(in._data);
    assert(std::get<0>(_grids).getVals() == std::get<0>(in._grids).getVals());
    _tail_coeffs = in._tail_coeffs;
    return (*this);
}




} // end of namespace FK
