#include "GFWrap.h"

#include <cmath>

namespace FK {

GFWrap::GFWrap(const FMatsubaraGrid& in):
    GridObject<ComplexType, FMatsubaraGrid>(std::make_tuple(in))
{
}

GFWrap::GFWrap(const std::tuple<FMatsubaraGrid> & in):GridObject<ComplexType, FMatsubaraGrid>(in)
{
}

GFWrap::GFWrap(GFWrap&& rhs):GridObject<ComplexType,FMatsubaraGrid>(rhs._grids)
{
    _data.swap(rhs._data);
    _dims.swap(rhs._dims);
    _f.swap(rhs._f);
}

GFWrap::GFWrap(const GFWrap& in):GridObject<ComplexType,FMatsubaraGrid>(in)
{
}

GFWrap::GFWrap(const GridObject<ComplexType,FMatsubaraGrid>& in):GridObject<ComplexType,FMatsubaraGrid>(in)
{
}



GFWrap::GFWrap(GridObject<ComplexType,FMatsubaraGrid>&& in):GridObject<ComplexType,FMatsubaraGrid>(in)
{
}

ComplexType GFWrap::operator()(const ComplexType &in) const
{
    assert(std::abs(real(in))<std::numeric_limits<RealType>::epsilon());
    if (imag(in)<0) return std::conj((*this)(std::conj(in)));
    if (std::get<0>(_grids).getNumber(in)>=std::get<0>(_grids)._w_max) { return _f(in); };
    return GridObject<ComplexType, FMatsubaraGrid>::operator()(in);
}

ComplexType GFWrap::operator()(const typename FMatsubaraGrid::point &in) const
{ 
    if (size_t(in) < std::get<0>(_grids).getSize() && std::abs(ComplexType(in) - ComplexType(std::get<0>(_grids)[size_t(in)]))<std::numeric_limits<RealType>::epsilon())
        return (*_data)[in._index]; 
    else return (*this)(in._val);
};

GFWrap& GFWrap::operator=(GFWrap&& in)
{
    if (std::get<0>(_grids).getPoints() == std::get<0>(in._grids).getPoints()) {
        _data.swap(in._data);
        _f.swap(in._f);
        return (*this);
        }
    else return this->copyAndInterpolate(in);
}

GFWrap& GFWrap::copyAndInterpolate(const GFWrap &in)
{
        _f = in._f;
        const FMatsubaraGrid &g = std::get<0>(_grids);
        const FMatsubaraGrid &g_rhs = std::get<0>(in._grids);
        #ifndef NDEBUG
        DEBUG("Interpolating GFWrap");
        #endif
        int min_index = 0; int max_index = g._w_max-g._w_min;
        if (g._w_min < g_rhs._w_min) { 
            for (int i=g._w_min; i<g_rhs._w_min; ++i) {
                auto w = g[i-g._w_min];
                this->get(w) = in(ComplexType(w)); // interpolation 
                //this->get(ComplexType(w)) = in(w);
                };
            min_index = g_rhs._w_min - g._w_min;
            }
        //DEBUG(*this);
        if (g._w_max > g_rhs._w_max) { 
            for (int i=g_rhs._w_max; i<g._w_max; ++i) {
                auto w = g[i-g._w_min];
                this->get(w) = in(ComplexType(w)); // interpolation 
                //this->get(ComplexType(w)) = in(w);
                };
            max_index = g_rhs._w_max-g._w_min;
            };
        //DEBUG(*this);
        
        for (int i=min_index; i<max_index; ++i) { 
            //DEBUG(i);
            auto w=g[i];
            (*this)[i] = in(ComplexType(w));
            }
    return (*this);
}


GFWrap& GFWrap::operator=(const GFWrap& in)
{
    assert(std::get<0>(_grids)._beta == std::get<0>(in._grids)._beta);
    if (std::get<0>(_grids).getPoints() == std::get<0>(in._grids).getPoints()) { 
        (*_data) = *(in._data);
        _f = in._f;
        return (*this);
        }
    else return this->copyAndInterpolate(in);
}




} // end of namespace FK
