#include "ssTriangle.h"

namespace SolveSpace
{
    STriangleLl *STriangleLl::Alloc()
    { return (STriangleLl *)Platform::AllocTemporary(sizeof(STriangleLl)); }

    Vector STriangle::Normal() const {
        Vector ab = b.Minus(a), bc = c.Minus(b);
        return ab.Cross(bc);
    }

    double STriangle::MinAltitude() const {
        double altA = a.DistanceToLine(b, c.Minus(b)),
            altB = b.DistanceToLine(c, a.Minus(c)),
            altC = c.DistanceToLine(a, b.Minus(a));

        return min(altA, min(altB, altC));
    }

    bool STriangle::ContainsPoint(Vector p) const {
        Vector n = Normal();
        if(MinAltitude() < LENGTH_EPS) {
            // shouldn't happen; zero-area triangle
            return false;
        }
        return ContainsPointProjd(n.WithMagnitude(1), p);
    }

    bool STriangle::ContainsPointProjd(Vector n, Vector p) const {
        Vector ab = b.Minus(a), bc = c.Minus(b), ca = a.Minus(c);

        Vector no_ab = n.Cross(ab);
        if(no_ab.Dot(p) < no_ab.Dot(a) - LENGTH_EPS) return false;

        Vector no_bc = n.Cross(bc);
        if(no_bc.Dot(p) < no_bc.Dot(b) - LENGTH_EPS) return false;

        Vector no_ca = n.Cross(ca);
        if(no_ca.Dot(p) < no_ca.Dot(c) - LENGTH_EPS) return false;

        return true;
    }

    bool STriangle::Raytrace(const Vector &rayPoint, const Vector &rayDir,
                            double *t, Vector *inters) const {
        // Algorithm from: "Fast, Minimum Storage Ray/Triangle Intersection" by
        // Tomas Moeller and Ben Trumbore.

        // Find vectors for two edges sharing vertex A.
        Vector edge1 = b.Minus(a);
        Vector edge2 = c.Minus(a);

        // Begin calculating determinant - also used to calculate U parameter.
        Vector pvec = rayDir.Cross(edge2);

        // If determinant is near zero, ray lies in plane of triangle.
        // Also, cull back facing triangles here.
        double det = edge1.Dot(pvec);
        if(-det < LENGTH_EPS) return false;
        double inv_det = 1.0f / det;

        // Calculate distance from vertex A to ray origin.
        Vector tvec = rayPoint.Minus(a);

        // Calculate U parameter and test bounds.
        double u = tvec.Dot(pvec) * inv_det;
        if (u < 0.0f || u > 1.0f) return false;

        // Prepare to test V parameter.
        Vector qvec = tvec.Cross(edge1);

        // Calculate V parameter and test bounds.
        double v = rayDir.Dot(qvec) * inv_det;
        if (v < 0.0f || u + v > 1.0f) return false;

        // Calculate t, ray intersects triangle.
        *t = edge2.Dot(qvec) * inv_det;

        // Calculate intersection point.
        if(inters != NULL) *inters = rayPoint.Plus(rayDir.ScaledBy(*t));

        return true;
    }

    double STriangle::SignedVolume() const {
        return a.Dot(b.Cross(c)) / 6.0;
    }

    double STriangle::Area() const {
        Vector ab = a.Minus(b);
        Vector cb = c.Minus(b);
        return ab.Cross(cb).Magnitude() / 2.0;
    }

    bool STriangle::IsDegenerate() const {
        return a.OnLineSegment(b, c) ||
            b.OnLineSegment(a, c) ||
            c.OnLineSegment(a, b);
    }

    void STriangle::FlipNormal() {
        swap(a, b);
        swap(an, bn);
    }

    STriangle STriangle::Transform(Vector u, Vector v, Vector n) const {
        STriangle tr = *this;
        tr.a  = tr.a.ScaleOutOfCsys(u, v, n);
        tr.an = tr.an.ScaleOutOfCsys(u, v, n);
        tr.b  = tr.b.ScaleOutOfCsys(u, v, n);
        tr.bn = tr.bn.ScaleOutOfCsys(u, v, n);
        tr.c  = tr.c.ScaleOutOfCsys(u, v, n);
        tr.cn = tr.cn.ScaleOutOfCsys(u, v, n);
        return tr;
    }

    STriangle STriangle::From(STriMeta meta, Vector a, Vector b, Vector c) {
        STriangle tr = {};
        tr.meta = meta;
        tr.a = a;
        tr.b = b;
        tr.c = c;
        return tr;
    }


}//ss