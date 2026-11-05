#pragma once

#include "ssVectorFileWriter.h"

namespace SolveSpace
{
            

    class GCodeFileWriter : public VectorFileWriter {
    public:
        SEdgeList sel;
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
    
} // namespace SolveSpace

