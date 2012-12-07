#include "SelfConsistency.h"

namespace FK { 

//
// Bethe SC
//

BetheSC::BetheSC(RealType t):_t(t)
{
}

BetheSC::GFType BetheSC::operator()(const GFType &gw)
{
    return gw*(_t*_t);
}

} //end of namespace FK

