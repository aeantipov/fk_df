#include <numeric>

#include "FKCommon.h"
#include "Grid.h"
#include "Container.h"
#include "GridObject.h"
#include "Solver.h"
#include "SelfConsistency.h"

#include <iostream>
#include <ctime>
#include <array>

using namespace FK;
typedef GFWrap GF;

template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7)
{
    return (std::abs(x-y)<tolerance);
}


int main()
{
    Log.setDebugging(true);
    size_t n_freq = 5;
    RealType beta = 10;
    FMatsubaraGrid gridF(0, n_freq, beta);
    BMatsubaraGrid gridB(0, n_freq, beta);

    GF Delta(gridF);
    std::function<ComplexType(ComplexType)> f1;
    RealType t = 0.5;
    f1 = [t](ComplexType w) -> ComplexType {return t*t/w;};
    Delta.fill(f1);
    
    auto D2 = Delta.shift(gridB[1]);
    if (D2[0]!=Delta[1]) return EXIT_FAILURE;

    KMesh kGrid(10);
    GridObject<RealType,KMesh,KMesh> e_k(std::forward_as_tuple(kGrid,kGrid));
    auto ekf = [](RealType kx, RealType ky){return 2.0*(cos(kx)+cos(ky));};
    e_k.fill(ekf);

    bool b1 = is_equal((e_k.shift(PI,PI)+e_k).sum(),0);
    if (!b1) return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}
