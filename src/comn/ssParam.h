#pragma once


#include <unordered_set>

#include "ssHandles.h>

namespace SolveSpace
{
    class Param {
public:
    int         tag;
    hParam      h;

    double      val;
    bool        known;
    bool        free;

    // Used only in the solver
    Param       *substd;

    static const hParam NO_PARAM;

    void Clear() {}
};

// Use a forward declaration in order to avoid pulling dsc.h in for units that
// don't need to use `ParamList`
template<class T, class H> class IdList;

using ParamList = IdList<Param, hParam>;

using ParamSet = std::unordered_set<hParam, HandleHasher<hParam>>;

} // namespace SolveSpace
