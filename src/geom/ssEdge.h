#pragma once

#include <ssVector.h>
#include "ssPoint.h"

namespace SolveSpace
{
    class SEdge {
    public:
        int    tag;
        int    auxA, auxB;
        Vector a, b;

        static SEdge From(Vector a, Vector b);
        bool EdgeCrosses(Vector a, Vector b, Vector *pi=NULL, SPointList *spl=NULL) const;
    };

    class SEdgeLl {
    public:
        SEdge       *se;
        SEdgeLl     *next;

        static SEdgeLl *Alloc();
    };

    class SEdgeList {
    public:
        List<SEdge>     l;

        void Clear();
        void AddEdge(Vector a, Vector b, int auxA=0, int auxB=0, int tag=0);
        bool AssemblePolygon(SPolygon *dest, SEdge *errorAt, bool keepDir=false) const;
        bool AssembleContour(Vector first, Vector last, SContour *dest,
                                SEdge *errorAt, bool keepDir, int start) const;
        int AnyEdgeCrossings(Vector a, Vector b,
            Vector *pi=NULL, SPointList *spl=NULL) const;
        bool ContainsEdgeFrom(const SEdgeList *sel) const;
        bool ContainsEdge(const SEdge *se) const;
        void CullExtraneousEdges(bool both=true);
        void MergeCollinearSegments(Vector a, Vector b);
    };
}//ss