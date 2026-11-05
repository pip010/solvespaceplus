#include "ssContour.h"

namespace SolveSpace
{

    void SContour::AddPoint(Vector p) {
        SPoint sp;
        sp.tag = 0;
        sp.p = p;

        l.Add(&sp);
    }

    void SContour::MakeEdgesInto(SEdgeList *el) const {
        int i;
        for(i = 0; i < (l.n - 1); i++) {
            el->AddEdge(l[i].p, l[i+1].p);
        }
    }

    void SContour::CopyInto(SContour *dest) const {
        for(const SPoint *sp = l.First(); sp; sp = l.NextAfter(sp)) {
            dest->AddPoint(sp->p);
        }
    }

    void SContour::FindPointWithMinX() {
        xminPt = {VERY_POSITIVE, VERY_POSITIVE, VERY_POSITIVE};
        for(const SPoint *sp = l.First(); sp; sp = l.NextAfter(sp)) {
            if(sp->p.x < xminPt.x) {
                xminPt = sp->p;
            }
        }
    }

    Vector SContour::ComputeNormal() const {
        Vector n = {};

        for(int i = 0; i < l.n - 2; i++) {
            Vector u = (l[i+1].p).Minus(l[i+0].p).WithMagnitude(1);
            Vector v = (l[i+2].p).Minus(l[i+1].p).WithMagnitude(1);
            Vector nt = u.Cross(v);
            if(nt.Magnitude() > n.Magnitude()) {
                n = nt;
            }
        }
        return n.WithMagnitude(1);
    }

    Vector SContour::AnyEdgeMidpoint() const {
        ssassert(l.n >= 2, "Need two points to find a midpoint");
        return ((l[0].p).Plus(l[1].p)).ScaledBy(0.5);
    }

    bool SContour::IsClockwiseProjdToNormal(Vector n) const {
        // Degenerate things might happen as we draw; doesn't really matter
        // what we do then.
        if(n.Magnitude() < 0.01) return true;

        return (SignedAreaProjdToNormal(n) < 0);
    }

    double SContour::SignedAreaProjdToNormal(Vector n) const {
        // An arbitrary 2d coordinate system that has n as its normal
        Vector u = n.Normal(0);
        Vector v = n.Normal(1);

        double area = 0;
        for(int i = 0; i < (l.n - 1); i++) {
            double u0 = (l[i  ].p).Dot(u);
            double v0 = (l[i  ].p).Dot(v);
            double u1 = (l[i+1].p).Dot(u);
            double v1 = (l[i+1].p).Dot(v);

            area += ((v0 + v1)/2)*(u1 - u0);
        }
        return area;
    }

    bool SContour::ContainsPointProjdToNormal(Vector n, Vector p) const {
        Vector u = n.Normal(0);
        Vector v = n.Normal(1);

        double up = p.Dot(u);
        double vp = p.Dot(v);

        bool inside = false;
        for(int i = 0; i < (l.n - 1); i++) {
            double ua = (l[i  ].p).Dot(u);
            double va = (l[i  ].p).Dot(v);
            // The curve needs to be exactly closed; approximation is death.
            double ub = (l[(i+1)%(l.n-1)].p).Dot(u);
            double vb = (l[(i+1)%(l.n-1)].p).Dot(v);

            if ((((va <= vp) && (vp < vb)) ||
                ((vb <= vp) && (vp < va))) &&
                (up < (ub - ua) * (vp - va) / (vb - va) + ua))
            {
            inside = !inside;
            }
        }

        return inside;
    }

    void SContour::Reverse() {
        l.Reverse();
    }

    //-----------------------------------------------------------------------------
    // Calculate the intersection point of two lines. Returns true for success,
    // false if they're parallel.
    //-----------------------------------------------------------------------------
    static bool IntersectionOfLines(double x0A, double y0A, double dxA, double dyA,
                                    double x0B, double y0B, double dxB, double dyB,
                                    double *xi, double *yi)
    {
        double A[2][2];
        double b[2];

        // The line is given to us in the form:
        //    (x(t), y(t)) = (x0, y0) + t*(dx, dy)
        // so first rewrite it as
        //    (x - x0, y - y0) dot (dy, -dx) = 0
        //    x*dy - x0*dy - y*dx + y0*dx = 0
        //    x*dy - y*dx = (x0*dy - y0*dx)

        // So write the matrix, pre-pivoted.
        if(fabs(dyA) > fabs(dyB)) {
            A[0][0] = dyA;  A[0][1] = -dxA;  b[0] = x0A*dyA - y0A*dxA;
            A[1][0] = dyB;  A[1][1] = -dxB;  b[1] = x0B*dyB - y0B*dxB;
        } else {
            A[1][0] = dyA;  A[1][1] = -dxA;  b[1] = x0A*dyA - y0A*dxA;
            A[0][0] = dyB;  A[0][1] = -dxB;  b[0] = x0B*dyB - y0B*dxB;
        }

        // Check the determinant; if it's zero then no solution.
        if(fabs(A[0][0]*A[1][1] - A[0][1]*A[1][0]) < LENGTH_EPS) {
            return false;
        }

        // Solve
        double v = A[1][0] / A[0][0];
        A[1][0] -= A[0][0]*v;
        A[1][1] -= A[0][1]*v;
        b[1] -= b[0]*v;

        // Back-substitute.
        *yi = b[1] / A[1][1];
        *xi = (b[0] - A[0][1]*(*yi)) / A[0][0];

        return true;
    }
    void SContour::OffsetInto(SContour *dest, double r) const {
        int i;

        for(i = 0; i < l.n; i++) {
            Vector a, b, c;
            Vector dp, dn;
            double thetan, thetap;

            a = l[WRAP(i-1, (l.n-1))].p;
            b = l[WRAP(i,   (l.n-1))].p;
            c = l[WRAP(i+1, (l.n-1))].p;

            dp = a.Minus(b);
            thetap = atan2(dp.y, dp.x);

            dn = b.Minus(c);
            thetan = atan2(dn.y, dn.x);

            // A short line segment in a badly-generated polygon might look
            // okay but screw up our sense of direction.
            if(dp.Magnitude() < LENGTH_EPS || dn.Magnitude() < LENGTH_EPS) {
                continue;
            }

            if(thetan > thetap && (thetan - thetap) > PI) {
                thetap += 2*PI;
            }
            if(thetan < thetap && (thetap - thetan) > PI) {
                thetan += 2*PI;
            }

            if(fabs(thetan - thetap) < (1*PI)/180) {
                Vector p = { b.x - r*sin(thetap), b.y + r*cos(thetap), 0 };
                dest->AddPoint(p);
            } else if(thetan < thetap) {
                // This is an inside corner. We have two edges, Ep and En. Move
                // out from their intersection by radius, normal to En, and
                // then draw a line parallel to En. Move out from their
                // intersection by radius, normal to Ep, and then draw a second
                // line parallel to Ep. The point that we want to generate is
                // the intersection of these two lines--it removes as much
                // material as we can without removing any that we shouldn't.
                double px0, py0, pdx, pdy;
                double nx0, ny0, ndx, ndy;
                double x = 0.0, y = 0.0;

                px0 = b.x - r*sin(thetap);
                py0 = b.y + r*cos(thetap);
                pdx = cos(thetap);
                pdy = sin(thetap);

                nx0 = b.x - r*sin(thetan);
                ny0 = b.y + r*cos(thetan);
                ndx = cos(thetan);
                ndy = sin(thetan);

                IntersectionOfLines(px0, py0, pdx, pdy,
                                    nx0, ny0, ndx, ndy,
                                    &x, &y);

                dest->AddPoint({x, y, 0});
            } else {
                if(fabs(thetap - thetan) < (6*PI)/180) {
                    Vector pp = { b.x - r*sin(thetap),
                                b.y + r*cos(thetap), 0 };
                    dest->AddPoint(pp);

                    Vector pn = { b.x - r*sin(thetan),
                                b.y + r*cos(thetan), 0 };
                    dest->AddPoint(pn);
                } else {
                    double theta;
                    for(theta = thetap; theta <= thetan; theta += (6*PI)/180) {
                        Vector p = { b.x - r*sin(theta),
                                    b.y + r*cos(theta), 0 };
                        dest->AddPoint(p);
                    }
                }
            }
        }
    }


    bool SContour::BridgeToContour(SContour *sc, SEdgeList *avoidEdges, List<Vector> *avoidPts)
    {
        int i, j;
        bool withbridge = true;

        // Start looking for a bridge on our new hole near its leftmost (min x)
        // point.
        int sco = 0;
        for(i = 0; i < (sc->l.n - 1); i++) {
            if((sc->l[i].p).EqualsExactly(sc->xminPt)) {
                sco = i;
            }
        }

        // And start looking on our merged contour at whichever point is nearest
        // to the leftmost point of the new segment.
        int thiso = 0;
        double dmin = 1e10;
        for(i = 0; i < l.n-1; i++) {
            Vector p = l[i].p;
            double d = (p.Minus(sc->xminPt)).MagSquared();
            if(d < dmin) {
                dmin = d;
                thiso = i;
            }
        }

        int thisp, scp;

        Vector a, b, *f;

        // First check if the contours share a point; in that case we should
        // merge them there, without a bridge.
        for(i = 0; i < l.n; i++) {
            thisp = WRAP(i+thiso, l.n-1);
            a = l[thisp].p;

            for(f = avoidPts->First(); f; f = avoidPts->NextAfter(f)) {
                if(f->Equals(a)) break;
            }
            if(f) continue;

            for(j = 0; j < (sc->l.n - 1); j++) {
                scp = WRAP(j+sco, (sc->l.n - 1));
                b = sc->l[scp].p;

                if(a.Equals(b)) {
                    withbridge = false;
                    goto haveEdge;
                }
            }
        }

        // If that fails, look for a bridge that does not intersect any edges.
        for(i = 0; i < l.n; i++) {
            thisp = WRAP(i+thiso, l.n);
            a = l[thisp].p;

            for(f = avoidPts->First(); f; f = avoidPts->NextAfter(f)) {
                if(f->Equals(a)) break;
            }
            if(f) continue;

            for(j = 0; j < (sc->l.n - 1); j++) {
                scp = WRAP(j+sco, (sc->l.n - 1));
                b = sc->l[scp].p;

                for(f = avoidPts->First(); f; f = avoidPts->NextAfter(f)) {
                    if(f->Equals(b)) break;
                }
                if(f) continue;

                if(avoidEdges->AnyEdgeCrossings(a, b) > 0) {
                    // doesn't work, bridge crosses an existing edge
                } else {
                    goto haveEdge;
                }
            }
        }

        // Tried all the possibilities, didn't find an edge
        return false;

    haveEdge:
        SContour merged = {};
        for(i = 0; i < l.n; i++) {
            if(withbridge || (i != thisp)) {
                merged.AddPoint(l[i].p);
            }
            if(i == thisp) {
                // less than or equal; need to duplicate the join point
                for(j = 0; j <= (sc->l.n - 1); j++) {
                    int jp = WRAP(j + scp, (sc->l.n - 1));
                    merged.AddPoint((sc->l[jp]).p);
                }
                // and likewise duplicate join point for the outer curve
                if(withbridge) {
                    merged.AddPoint(l[i].p);
                }
            }
        }

        // and future bridges mustn't cross our bridge, and it's tricky to get
        // things right if two bridges come from the same point
        if(withbridge) {
            avoidEdges->AddEdge(a, b);
            avoidPts->Add(&a);
        }
        avoidPts->Add(&b);

        l.Clear();
        l = merged.l;
        return true;
    }

    bool SContour::IsEmptyTriangle(int ap, int bp, int cp, double scaledEPS) const {

        STriangle tr = {};
        tr.a = l[ap].p;
        tr.b = l[bp].p;
        tr.c = l[cp].p;

        // Accelerate with an axis-aligned bounding box test
        Vector maxv = tr.a, minv = tr.a;
        (tr.b).MakeMaxMin(&maxv, &minv);
        (tr.c).MakeMaxMin(&maxv, &minv);

        Vector n = {0, 0, -1};

        int i;
        for(i = 0; i < l.n; i++) {
            if(i == ap || i == bp || i == cp) continue;

            Vector p = l[i].p;
            if(p.OutsideAndNotOn(maxv, minv)) continue;

            // A point on the edge of the triangle is considered to be inside,
            // and therefore makes it a non-ear; but a point on the vertex is
            // "outside", since that's necessary to make bridges work.
            if(p.EqualsExactly(tr.a)) continue;
            if(p.EqualsExactly(tr.b)) continue;
            if(p.EqualsExactly(tr.c)) continue;

            if(tr.ContainsPointProjd(n, p)) {
                return false;
            }
        }
        return true;
    }

    // Test if ray b->d passes through triangle a,b,c
    static bool RayIsInside(Vector a, Vector c, Vector b, Vector d) {
        // coincident edges are not considered to intersect the triangle
        if (d.Equals(a)) return false;
        if (d.Equals(c)) return false;
        // if d and c are on opposite sides of ba, we are ok
        // likewise if d and a are on opposite sides of bc
        Vector ba = a.Minus(b);
        Vector bc = c.Minus(b);
        Vector bd = d.Minus(b);

        // perpendicular to (x,y) is (x,-y) so dot that with the two points. If they
        // have opposite signs their product will be negative. If bd and bc are on
        // opposite sides of ba the ray does not intersect. Likewise for bd,ba and bc.
        if ( (bd.x*(ba.y) + (bd.y * (-ba.x))) * ( bc.x*(ba.y) + (bc.y * (-ba.x))) < LENGTH_EPS)
            return false;
        if ( (bd.x*(bc.y) + (bd.y * (-bc.x))) * ( ba.x*(bc.y) + (ba.y * (-bc.x))) < LENGTH_EPS)
            return false;

        return true;
    }

    bool SContour::IsEar(int bp, double scaledEps) const {
        int ap = WRAP(bp-1, l.n),
            cp = WRAP(bp+1, l.n);

        STriangle tr = {};
        tr.a = l[ap].p;
        tr.b = l[bp].p;
        tr.c = l[cp].p;

        if((tr.a).Equals(tr.c)) {
            // This is two coincident and anti-parallel edges. Zero-area, so
            // won't generate a real triangle, but we certainly can clip it.
            return true;
        }

        Vector n = {0, 0, -1};
        if((tr.Normal()).Dot(n) < scaledEps) {
            // This vertex is reflex, or between two collinear edges; either way,
            // it's not an ear.
            return false;
        }

        // Accelerate with an axis-aligned bounding box test
        Vector maxv = tr.a, minv = tr.a;
        (tr.b).MakeMaxMin(&maxv, &minv);
        (tr.c).MakeMaxMin(&maxv, &minv);

        int i;
        for(i = 0; i < l.n; i++) {
            if(i == ap || i == bp || i == cp) continue;

            Vector p = l[i].p;
            if(p.OutsideAndNotOn(maxv, minv)) continue;

            // A point on the edge of the triangle is considered to be inside,
            // and therefore makes it a non-ear; but a point on the vertex is
            // "outside", since that's necessary to make bridges work.
            if(p.EqualsExactly(tr.a)) continue;
            if(p.EqualsExactly(tr.c)) continue;
            // points coincident with bp have to be allowed for bridges but edges
            // from that other point must not cross through our triangle.
            if(p.EqualsExactly(tr.b)) {
                int j = WRAP(i-1, l.n);
                int k = WRAP(i+1, l.n);
                Vector jp = l[j].p;
                Vector kp = l[k].p;

                // two consecutive bridges (A,B,C) and later (C,B,A) are not an ear
                if (jp.Equals(tr.c) && kp.Equals(tr.a)) return false;
                // check both edges from the point in question
                if (!RayIsInside(tr.a, tr.c, p,jp) && !RayIsInside(tr.a, tr.c, p,kp))
                    continue;
            }

            if(tr.ContainsPointProjd(n, p)) {
                return false;
            }
        }
        return true;
    }

    void SContour::ClipEarInto(SMesh *m, int bp, double scaledEps) {
        int ap = WRAP(bp-1, l.n),
            cp = WRAP(bp+1, l.n);

        STriangle tr = {};
        tr.a = l[ap].p;
        tr.b = l[bp].p;
        tr.c = l[cp].p;
        if(tr.Normal().MagSquared() < scaledEps*scaledEps) {
            // A vertex with more than two edges will cause us to generate
            // zero-area triangles, which must be culled.
        } else {
            m->AddTriangle(&tr);
        }

        // By deleting the point at bp, we may change the ear-ness of the points
        // on either side.
        l[ap].ear = EarType::UNKNOWN;
        l[cp].ear = EarType::UNKNOWN;

        l.ClearTags();
        l[bp].tag = 1;
        l.RemoveTagged();
    }

    void SContour::UvTriangulateInto(SMesh *m, SSurface *srf) {
        Vector tu, tv;
        srf->TangentsAt(0.5, 0.5, &tu, &tv);
        double s = sqrt(tu.MagSquared() + tv.MagSquared());
        // We would like to apply our tolerances in xyz; but that would be a lot
        // of work, so at least scale the epsilon semi-reasonably. That's
        // perfect for square planes, less perfect for anything else.
        double scaledEps = LENGTH_EPS / s;

        int i;
        // Clean the original contour by removing any zero-length edges.
        // initialize eartypes to unknown while we're going over them.
        l.ClearTags();
        l[0].ear = EarType::UNKNOWN;
        for(i = 1; i < l.n; i++) {
        l[i].ear = EarType::UNKNOWN;
        if((l[i].p).Equals(l[i-1].p)) {
                l[i].tag = 1;
            }
        }
        if( (l[0].p).Equals(l[l.n-1].p) ) {
            l[l.n-1].tag = 1;
        }
        l.RemoveTagged();

        // Handle simple triangle fans all at once. This pass is optional.
        if(srf->degm == 1 && srf->degn == 1) {
            l.ClearTags();
            int j=0;
            int pstart = 0;
            double elen = -1.0;
            double oldspan = 0.0;
            for(i = 1; i < l.n; i++) {
                Vector ab = l[i].p.Minus(l[i-1].p);
                // first time just measure the segment
                if (elen < 0.0) {
                    elen = ab.Dot(ab);
                    oldspan = elen;
                    j = 1;
                    continue;
                }
                // check for consecutive segments of similar size which are also
                // ears and where the group forms a convex ear
                bool end = false;
                double ratio = ab.Dot(ab) / elen;
                if ((ratio < 0.25) || (ratio > 4.0)) end = true;

                double slen = l[pstart].p.Minus(l[i].p).MagSquared();
                if (slen < oldspan) end = true;

                if (!IsEar(i-1, scaledEps) ) end = true;
    //            if ((j>0) && !IsEar(pstart, i-1, i, scaledEps)) end = true;
                if ((j>0) && !IsEmptyTriangle(pstart, i-1, i, scaledEps)) end = true;
                // the new segment is valid so add to the fan
                if (!end) {
                    j++;
                    oldspan = slen;
                }
                // we need to stop at the end of polygon but may still
                if (i == l.n-1) {
                    end = true;
                }
                if (end) {  // triangulate the fan and tag the vertices
                    if (j > 3) {
                        Vector center = l[pstart+1].p.Plus(l[pstart+j-1].p).ScaledBy(0.5);
                        for (int x=0; x<j; x++) {
                            STriangle tr = {};
                            tr.a = center;
                            tr.b = l[pstart+x].p;
                            tr.c = l[pstart+x+1].p;
                            m->AddTriangle(&tr);
                        }
                        for (int x=1; x<j; x++) {
                            l[pstart+x].tag = 1;
                        }
                        STriangle tr = {};
                        tr.a = center;
                        tr.b = l[pstart+j].p;
                        tr.c = l[pstart].p;
                        m->AddTriangle(&tr);    
                    }
                    pstart = i-1;
                    elen = ab.Dot(ab);
                    oldspan = elen;
                    j = 1;
                }
            }
            l.RemoveTagged();
        }  // end optional fan creation pass

        bool toggle = false;
        while(l.n > 3) {
            int bestEar = -1;
            double bestChordTol = VERY_POSITIVE;
            // Alternate the starting position so we generate strip-like
            // triangulations instead of fan-like
            toggle = !toggle;
            int offset = toggle ? -1 : 0;
            for(i = 0; i < l.n; i++) {
                int ear = WRAP(i+offset, l.n);
                if(l[ear].ear == EarType::UNKNOWN) {
                    (l[ear]).ear = IsEar(ear, scaledEps) ? EarType::EAR : EarType::NOT_EAR;
                }
                if(l[ear].ear == EarType::EAR) {
                    if(srf->degm == 1 && srf->degn == 1) {
                        // This is a plane; any ear is a good ear.
                        bestEar = ear;
                        break;
                    }
                    // If we are triangulating a curved surface, then try to
                    // clip ears that have a small chord tolerance from the
                    // surface.
                    Vector prev = l[WRAP((i+offset-1), l.n)].p,
                        next = l[WRAP((i+offset+1), l.n)].p;
                    double tol = srf->ChordToleranceForEdge(prev, next);
                    if(tol < bestChordTol - scaledEps) {
                        bestEar = ear;
                        bestChordTol = tol;
                    }
                    if(bestChordTol < 0.1*SS.ChordTolMm()) {
                        break;
                    }
                }
            }
            if(bestEar < 0) {
                dbp("couldn't find an ear! fail");
                return;
            }
            ClipEarInto(m, bestEar, scaledEps);
        }

        ClipEarInto(m, 0, scaledEps); // add the last triangle
    }

}//ss