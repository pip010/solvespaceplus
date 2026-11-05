
#pragma once

#include <ssVector.h>
#include <ssPoint.h>
#include <ssList.hpp>

namespace SolveSpace
{
    class SContour {
    public:
        int             tag;
        int             timesEnclosed;
        Vector          xminPt;
        List<SPoint>    l;

        void AddPoint(Vector p);
        void MakeEdgesInto(SEdgeList *el) const;
        void Reverse();
        Vector ComputeNormal() const;
        double SignedAreaProjdToNormal(Vector n) const;
        bool IsClockwiseProjdToNormal(Vector n) const;
        bool ContainsPointProjdToNormal(Vector n, Vector p) const;
        void OffsetInto(SContour *dest, double r) const;
        void CopyInto(SContour *dest) const;
        void FindPointWithMinX();
        Vector AnyEdgeMidpoint() const;

        bool IsEmptyTriangle(int ap, int bp, int cp, double scaledEPS) const;
        bool IsEar(int bp, double scaledEps) const;
        bool BridgeToContour(SContour *sc, SEdgeList *el, List<Vector> *vl);
        void ClipEarInto(SMesh *m, int bp, double scaledEps);
        void UvTriangulateInto(SMesh *m, SSurface *srf);
    };
}//ss
