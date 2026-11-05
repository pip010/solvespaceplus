#pragma once

#include "ssTriangle.h"
#include <ssVector.h>
#include <ssList.hpp>

namespace SolveSpace
{

    class SMesh {
    public:
        List<STriangle>     l;

        bool    flipNormal;
        bool    keepInsideOtherShell;
        bool    keepCoplanar;
        bool    atLeastOneDiscarded;
        bool    isTransparent;

        void Clear();
        void AddTriangle(const STriangle *st);
        void AddTriangle(STriMeta meta, Vector a, Vector b, Vector c);
        void AddTriangle(STriMeta meta, Vector n,
                        Vector a, Vector b, Vector c);
        void DoBounding(Vector v, Vector *vmax, Vector *vmin) const;
        void GetBounding(Vector *vmax, Vector *vmin) const;

        void Simplify(int start);

        void AddAgainstBsp(SMesh *srcm, SBsp3 *bsp3);
        void MakeFromUnionOf(SMesh *a, SMesh *b);
        void MakeFromDifferenceOf(SMesh *a, SMesh *b);
        void MakeFromIntersectionOf(SMesh *a, SMesh *b);

        void MakeFromCopyOf(SMesh *a);
        void MakeFromTransformationOf(SMesh *a, Vector trans,
                                    Quaternion q, double scale);
        void MakeFromAssemblyOf(SMesh *a, SMesh *b);

        void MakeEdgesInPlaneInto(SEdgeList *sel, Vector n, double d);
        void MakeOutlinesInto(SOutlineList *sol, EdgeKind type);

        void PrecomputeTransparency();
        void RemoveDegenerateTriangles();
        double CalculateVolume() const;
        double CalculateSurfaceArea(const std::vector<uint32_t> &faces) const;

        bool IsEmpty() const;
        void RemapFaces(Group *g, int remap);

        uint32_t FirstIntersectionWith(Point2d mp) const;

        Vector GetCenterOfMass() const;
    };
}//ss