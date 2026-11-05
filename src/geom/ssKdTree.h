#pragma once

#include "ssEdge.h"
#include "ssTriangle.h"

namespace SolveSpace
{


    enum class EdgeKind : uint32_t {
        NAKED_OR_SELF_INTER  = 100,
        SELF_INTER           = 200,
        TURNING              = 300,
        EMPHASIZED           = 400,
        SHARP                = 500,
    };

    class SKdNodeEdges {
    public:
        int which; // whether c is x, y, or z
        double c;
        SKdNodeEdges    *gt;
        SKdNodeEdges    *lt;

        SEdgeLl         *edges;

        static SKdNodeEdges *From(SEdgeList *sel);
        static SKdNodeEdges *From(SEdgeLl *sell);
        static SKdNodeEdges *Alloc();
        int AnyEdgeCrossings(Vector a, Vector b, int cnt,
            Vector *pi=NULL, SPointList *spl=NULL) const;
    };


    class SKdNode {
    public:
        struct EdgeOnInfo {
            int        count;
            bool       frontFacing;
            bool       intersectsMesh;
            STriangle *tr;
            int        ai;
            int        bi;
        };

        int which;  // whether c is x, y, or z
        double c;

        SKdNode      *gt;
        SKdNode      *lt;

        STriangleLl  *tris;

        static SKdNode *Alloc();
        static SKdNode *From(SMesh *m);
        static SKdNode *From(STriangleLl *tll);

        void AddTriangle(STriangle *tr);
        void MakeMeshInto(SMesh *m) const;
        void ListTrianglesInto(std::vector<STriangle *> *tl) const;
        void ClearTags() const;

        void FindEdgeOn(Vector a, Vector b, int cnt, bool coplanarIsInter, EdgeOnInfo *info) const;
        void MakeCertainEdgesInto(SEdgeList *sel, EdgeKind how, bool coplanarIsInter,
                                bool *inter, bool *leaky, int auxA = 0) const;
        void MakeOutlinesInto(SOutlineList *sel, EdgeKind tagKind) const;

        void OcclusionTestLine(SEdge orig, SEdgeList *sel, int cnt) const;
        void SplitLinesAgainstTriangle(SEdgeList *sel, STriangle *tr) const;

        void SnapToMesh(SMesh *m);
        void SnapToVertex(Vector v, SMesh *extras);
    };
}//ss