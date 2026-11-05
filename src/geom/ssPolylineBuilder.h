#pragma once

#include <ssVector.h>
#include <ssEdge.h>
#include <ssOutline.h>


namespace SolveSpace
{

    class PolylineBuilder {
    public:
        struct Edge;

        struct Vertex {
            Vector              pos;
            std::vector<Edge *> edges;

            bool GetNext(uint32_t kind, Vertex **next, Edge **nextEdge);
            bool GetNext(uint32_t kind, Vector plane, double d, Vertex **next, Edge **nextEdge);
            size_t CountEdgesWithTagAndKind(int tag, uint32_t kind) const;
        };

        struct VertexPairHash {
            size_t operator()(const std::pair<Vertex *, Vertex *> &v) const;
        };

        struct Edge {
            Vertex   *a;
            Vertex   *b;
            uint32_t  kind;
            int       tag;

            union {
                uintptr_t  data;
                SOutline  *outline;
                SEdge     *edge;
            };

            Vertex *GetOtherVertex(Vertex *v) const;
            bool GetStartAndNext(Vertex **start, Vertex **next, bool loop) const;
        };

        std::unordered_map<Vector, Vertex *, VectorHash, VectorPred> vertices;
        std::unordered_map<std::pair<Vertex *, Vertex *>, Edge *, VertexPairHash> edgeMap;
        std::vector<Edge *> edges;

        ~PolylineBuilder();
        void Clear();

        Vertex *AddVertex(const Vector &pos);
        Edge *AddEdge(const Vector &p0, const Vector &p1, uint32_t kind, uintptr_t data = 0);
        void Generate(std::function<void(Vertex *start, Vertex *next, Edge *edge)> const &startFunc,
                    std::function<void(Vertex *next, Edge *edge)> const &nextFunc,
                    std::function<void(Edge *)> const &aloneFunc,
                    std::function<void()> const &endFunc = []() {});

        void MakeFromEdges(const SEdgeList &sel);
        void MakeFromOutlines(const SOutlineList &sol);
        void GenerateEdges(SEdgeList *sel);
        void GenerateOutlines(SOutlineList *sol);
    };

}//ss