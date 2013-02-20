#ifndef ___FK_GF_WRAP_H__
#define ___FK_GF_WRAP_H__

#include "MatsubaraGrid.hpp"
#include "GridObject.hpp"

namespace FK {

using namespace GFTools;

class GFWrap : public GridObject<ComplexType, FMatsubaraGrid>
{
    GFWrap& copyAndInterpolate(const GFWrap &in);
public:
    typedef FMatsubaraGrid::point point;
    GFWrap(const FMatsubaraGrid &in);
    GFWrap(const std::tuple<FMatsubaraGrid> &in);
    GFWrap(const GFWrap& in);
    GFWrap(GFWrap&& in);
    GFWrap(const GridObject<ComplexType, FMatsubaraGrid>& in);
    GFWrap(GridObject<ComplexType, FMatsubaraGrid>&& in);
    GFWrap(const std::string& fname, RealType beta);

    ComplexType operator()(const point &in) const;
    ComplexType operator()(const ComplexType &in) const;
    
    //using GridObject<ComplexType, FMatsubaraGrid>::get;
    using GridObject<ComplexType, FMatsubaraGrid>::operator=;
    GFWrap& operator=(const GridObject<ComplexType, FMatsubaraGrid>& in);
    GFWrap& operator=(GridObject<ComplexType, FMatsubaraGrid>&& in);
    GFWrap& operator=(const GFWrap& in);
    GFWrap& operator=(GFWrap&& in);
};

} // end::namespace FK

#endif // endif :: #ifndef ___FK_GF_WRAP_H__

