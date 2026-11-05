
#pragma once 

#include <ssVector.h>
#include <ssRgbaColor.h>

namespace SolveSpace
{

    typedef struct {
        uint32_t face;
        RgbaColor color;
    } STriMeta;

    class STriangle {
    public:
        int         tag;
        STriMeta    meta;

        union {
            struct { Vector a, b, c; };
            Vector vertices[3];
        };

        union {
            struct { Vector an, bn, cn; };
            Vector normals[3];
        };

        static STriangle From(STriMeta meta, Vector a, Vector b, Vector c);
        Vector Normal() const;
        void FlipNormal();
        double MinAltitude() const;
        bool ContainsPoint(Vector p) const;
        bool ContainsPointProjd(Vector n, Vector p) const;
        STriangle Transform(Vector o, Vector u, Vector v) const;
        bool Raytrace(const Vector &rayPoint, const Vector &rayDir,
                    double *t, Vector *inters) const;
        double SignedVolume() const;
        double Area() const;
        bool IsDegenerate() const;
    };


    // A linked list of triangles
    class STriangleLl {
    public:
        STriangle       *tri;

        STriangleLl     *next;

        static STriangleLl *Alloc();
    };
    
}//ss