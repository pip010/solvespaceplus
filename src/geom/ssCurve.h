#pragma once

#include "ssVector.h"
#include <ssHandles.h>
#include <ssSurface.h>
#include <ssShell.h>
#include <ssBezier.h>
#include <ssQuaternion.h>

namespace SolveSpace
{


// Stuff for the surface trim curves: piecewise linear
class SCurvePt {
public:
    int         tag;
    Vector      p;
    bool        vertex;
};

class SCurve {
public:
    hSCurve         h;

    // In a Boolean, C = A op B. The curves in A and B get copied into C, and
    // therefore must get new hSCurves assigned. For the curves in A and B,
    // we use newH to record their new handle in C.
    hSCurve         newH;
    enum class Source : uint32_t {
        A             = 100,
        B             = 200,
        INTERSECTION  = 300
    };
    Source          source;

    bool            isExact;
    SBezier         exact;

    List<SCurvePt>  pts;

    hSSurface       surfA;
    hSSurface       surfB;

    static SCurve FromTransformationOf(SCurve *a, Vector t,
                                       Quaternion q, double scale);
    SCurve MakeCopySplitAgainst(SShell *agnstA, SShell *agnstB,
                                SSurface *srfA, SSurface *srfB) const;
    void RemoveShortSegments(SSurface *srfA, SSurface *srfB);
    SSurface *GetSurfaceA(SShell *a, SShell *b) const;
    SSurface *GetSurfaceB(SShell *a, SShell *b) const;

    void Clear();
    void GetAxisAlignedBounding(Vector *ptMax, Vector *ptMin) const;
};

} // namespace SolveSpace
