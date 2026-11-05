#pragma once

namespace SolveSpace
{

    class StepFileWriter {
    public:
        bool HasCartesianPointAnAlias(int number, Vector v, int vertex,
                                    bool *vertex_has_alias = nullptr);
        int InsertPoint(int number);
        int InsertVertex(int number);
        bool HasBSplineCurveAnAlias(int number, std::vector<int> points);
        int InsertCurve(int number);
        bool HasEdgeCurveAnAlias(int number, int prevFinish, int thisFinish, int curveId,
                            bool *flip = nullptr);
        int InsertEdgeCurve(int number);
        bool HasOrientedEdgeAnAlias(int number, int edgeCurveId, bool flip);
        int InsertOrientedEdge(int number);
        void ExportSurfacesTo(const Platform::Path &filename);
        void WriteHeader();
        void WriteProductHeader();
        int ExportCurve(SBezier *sb);
        int ExportCurveLoop(SBezierLoop *loop, bool inner);
        void ExportSurface(SSurface *ss, SBezierList *sbl);
        void WriteWireframe();
        void WriteFooter();

        List<int> curves;
        List<int> advancedFaces;
        FILE *f;
        int id;

        // Structs to keep track of duplicated entities.
        // Basic alias.
        typedef struct {
            int reference;
            std::vector<int> aliases;
        } alias_t;

        // Cartesian points.
        typedef struct {
            alias_t alias;
            alias_t vertexAlias;
            Vector v;
        } pointAliases_t;

        // Curves.
        typedef struct {
            alias_t alias;
            std::vector<int> memberPoints;
            RgbaColor color;
        } curveAliases_t;

        // Edges.
        typedef struct {
            alias_t alias;
            int prevFinish;
            int thisFinish;
            int curveId;
            RgbaColor color;
        } edgeCurveAliases_t;

        // Edges.
        typedef struct {
            alias_t alias;
            int edgeCurveId;
            bool flip;
        } orientedEdgeAliases_t;

        std::vector<pointAliases_t> pointAliases;
        std::vector<edgeCurveAliases_t> edgeCurveAliases;
        std::vector<curveAliases_t> curveAliases;
        std::vector<orientedEdgeAliases_t> orientedEdgeAliases;
        bool exportParts = true;
        RgbaColor currentColor;
    };

    class Step2dFileWriter : public VectorFileWriter {
    StepFileWriter sfw;
    void StartPath( RgbaColor strokeRgb, double lineWidth,
                    bool filled, RgbaColor fillRgb, hStyle hs) override;
    void FinishPath(RgbaColor strokeRgb, double lineWidth,
                    bool filled, RgbaColor fillRgb, hStyle hs, bool closed) override;
    void Triangle(STriangle *tr) override;
    void Bezier(SBezier *sb) override;
    void Background(RgbaColor color) override;
    void StartFile() override;
    void FinishAndCloseFile() override;
    bool HasCanvasSize() const override { return false; }
    bool CanOutputMesh() const override { return false; }
};

}//ss