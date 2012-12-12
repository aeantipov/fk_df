#include "SelfConsistency.h"

namespace FK { 

//
// Bethe SC
//

BetheSC::BetheSC(RealType t):_t(t)
{
}

BetheSC::GFType BetheSC::operator()(const BetheSC::GFType &gw) const
{
    return gw*(_t*_t);
}

BetheSC::GFType BetheSC::operator()(const BetheSC::GFType &gw, const BetheSC::GFType &Delta) const
{
    return (*this)(gw);
}


//
// CubicInfDMFTSC
//


CubicInfDMFTSC :: CubicInfDMFTSC(RealType t, const RealGrid realgrid):
    _t(t),
    _realgrid(realgrid),
    _nominator(RealW(realgrid)),
    _denominator(ComplW(realgrid))
{
    std::function<RealType(RealType)> f1 = [=](RealType w){return std::exp(-w*w/t/t);};
    _nominator = f1;
}


CubicInfDMFTSC::GFType CubicInfDMFTSC::operator()(const CubicInfDMFTSC::GFType &gw, const CubicInfDMFTSC::GFType &Delta) const
{
    for (ComplexType w : gw.getGrid().getVals()) { 
    //std::function<ComplexType(RealType)> f2 = [&](RealType w){return };
        DEBUG(w)
    }; 
}

} //end of namespace FK

