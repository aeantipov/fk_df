#include "SelfConsistency.h"

namespace FK { 

//
// Bethe SC
//

BetheSC::BetheSC(RealType t):_t(t)
{
}

BetheSC::GFType BetheSC::operator()(const BetheSC::GFType &gw)
{
    return gw*(_t*_t);
}

} //end of namespace FK

