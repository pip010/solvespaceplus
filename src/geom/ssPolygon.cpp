#include "ssPolygon.h"


namespace SolveSpace
{

    void SPolygon::Clear() {
        int i;
        for(i = 0; i < l.n; i++) {
            (l[i]).l.Clear();
        }
        l.Clear();
    }

    void SPolygon::AddEmptyContour() {
        SContour c = {};
        l.Add(&c);
    }

    void SPolygon::MakeEdgesInto(SEdgeList *el) const {
        int i;
        for(i = 0; i < l.n; i++) {
            (l[i]).MakeEdgesInto(el);
        }
    }

    Vector SPolygon::ComputeNormal() const {
        if(l.IsEmpty())
            return {};
        return (l[0]).ComputeNormal();
    }

    double SPolygon::SignedArea() const {
        double area = 0;
        // This returns the true area only if the contours are all oriented
        // correctly, with the holes backwards from the outer contours.
        for(const SContour *sc = l.First(); sc; sc = l.NextAfter(sc)) {
            area += sc->SignedAreaProjdToNormal(normal);
        }
        return area;
    }

    bool SPolygon::ContainsPoint(Vector p) const {
        return (WindingNumberForPoint(p) % 2) == 1;
    }

    size_t SPolygon::WindingNumberForPoint(Vector p) const {
        auto winding = std::count_if(l.begin(), l.end(), [&](const SContour &sc) {
            return sc.ContainsPointProjdToNormal(normal, p);
        });
        return winding;
    }

    void SPolygon::FixContourDirections() {
        // At output, the contour's tag will be 1 if we reversed it, else 0.
        l.ClearTags();

        // Outside curve looks counterclockwise, projected against our normal.
        int i, j;
        for(i = 0; i < l.n; i++) {
            SContour *sc = &(l[i]);
            if(sc->l.n < 2) continue;
            // The contours may not intersect, but they may share vertices; so
            // testing a vertex for point-in-polygon may fail, but the midpoint
            // of an edge is okay.
            Vector pt = (((sc->l[0]).p).Plus(sc->l[1].p)).ScaledBy(0.5);

            sc->timesEnclosed = 0;
            bool outer = true;
            for(j = 0; j < l.n; j++) {
                if(i == j) continue;
                SContour *sct = &(l[j]);
                if(sct->ContainsPointProjdToNormal(normal, pt)) {
                    outer = !outer;
                    (sc->timesEnclosed)++;
                }
            }

            bool clockwise = sc->IsClockwiseProjdToNormal(normal);
            if((clockwise && outer) || (!clockwise && !outer)) {
                sc->Reverse();
                sc->tag = 1;
            }
        }
    }

    bool SPolygon::IsEmpty() const {
        if(l.IsEmpty() || l[0].l.IsEmpty())
            return true;
        return false;
    }

    Vector SPolygon::AnyPoint() const {
        ssassert(!IsEmpty(), "Need at least one point");
        return l[0].l[0].p;
    }

    bool SPolygon::SelfIntersecting(Vector *intersectsAt) const {
        SEdgeList el = {};
        MakeEdgesInto(&el);
        SKdNodeEdges *kdtree = SKdNodeEdges::From(&el);

        int cnt = 1;
        el.l.ClearTags();

        bool ret = false;
        SEdge *se;
        for(se = el.l.First(); se; se = el.l.NextAfter(se)) {
            int inters = kdtree->AnyEdgeCrossings(se->a, se->b, cnt, intersectsAt);
            if(inters != 1) {
                ret = true;
                break;
            }
            cnt++;
        }

        el.Clear();
        return ret;
    }

    void SPolygon::InverseTransformInto(SPolygon *sp, Vector u, Vector v, Vector n) const {
        for(const SContour &sc : l) {
            SContour tsc = {};
            tsc.timesEnclosed = sc.timesEnclosed;
            for(const SPoint &sp : sc.l) {
                tsc.AddPoint(sp.p.DotInToCsys(u, v, n));
            }
            sp->l.Add(&tsc);
        }
    }

    //-----------------------------------------------------------------------------
    // Low-quality routines to cutter radius compensate a polygon. Assumes the
    // polygon is in the xy plane, and the contours all go in the right direction
    // with respect to normal (0, 0, -1).
    //-----------------------------------------------------------------------------
    void SPolygon::OffsetInto(SPolygon *dest, double r) const {
        int i;
        dest->Clear();
        for(i = 0; i < l.n; i++) {
            const SContour *sc = &(l[i]);
            dest->AddEmptyContour();
            sc->OffsetInto(&(dest->l[dest->l.n-1]), r);
        }
    }


    void SPolygon::UvTriangulateInto(SMesh *m, SSurface *srf) {
        if(l.n <= 0) return;

        //int64_t in = GetMilliseconds();

        normal = {0, 0, 1};

        while(l.n > 0) {
            FixContourDirections();
            l.ClearTags();

            // Find a top-level contour, and start with that. Then build bridges
            // in order to merge all its islands into a single contour.
            SContour *top;
            for(top = l.First(); top; top = l.NextAfter(top)) {
                if(top->timesEnclosed == 0) {
                    break;
                }
            }
            if(!top) {
                dbp("polygon has no top-level contours?");
                return;
            }

            // Start with the outer contour
            SContour merged = {};
            top->tag = 1;
            top->CopyInto(&merged);
            merged.l.RemoveLast(1);

            // List all of the edges, for testing whether bridges work.
            SEdgeList el = {};
            top->MakeEdgesInto(&el);
            List<Vector> vl = {};

            // And now find all of its holes. Note that we will also find any
            // outer contours that lie entirely within this contour, and any
            // holes for those contours. But that's okay, because we can merge
            // those too.
            SContour *sc;
            for(sc = l.First(); sc; sc = l.NextAfter(sc)) {
                if(sc->timesEnclosed != 1) continue;
                if(sc->l.n < 2) continue;

                // Test the midpoint of an edge. Our polygon may not be self-
                // intersecting, but two contours may share a vertex; so a
                // vertex could be on the edge of another polygon, in which
                // case ContainsPointProjdToNormal returns indeterminate.
                Vector tp = sc->AnyEdgeMidpoint();
                if(top->ContainsPointProjdToNormal(normal, tp)) {
                    sc->tag = 2;
                    sc->MakeEdgesInto(&el);
                    sc->FindPointWithMinX();
                }
            }

    //        dbp("finished finding holes: %d ms", (int)(GetMilliseconds() - in));
            for(;;) {
                double xmin = 1e10;
                SContour *scmin = NULL;

                for(sc = l.First(); sc; sc = l.NextAfter(sc)) {
                    if(sc->tag != 2) continue;

                    if(sc->xminPt.x < xmin) {
                        xmin = sc->xminPt.x;
                        scmin = sc;
                    }
                }
                if(!scmin) break;

                if(!merged.BridgeToContour(scmin, &el, &vl)) {
                    dbp("couldn't merge our hole");
                    return;
                }
    //            dbp("   bridged to contour: %d ms", (int)(GetMilliseconds() - in));
                scmin->tag = 3;
            }
    //        dbp("finished merging holes: %d ms", (int)(GetMilliseconds() - in));

            merged.UvTriangulateInto(m, srf);
    //        dbp("finished ear clippping: %d ms", (int)(GetMilliseconds() - in));
            merged.l.Clear();
            el.Clear();
            vl.Clear();

            // Careful, need to free the points within the contours, and not just
            // the contours themselves. This was a tricky memory leak.
            for(sc = l.First(); sc; sc = l.NextAfter(sc)) {
                if(sc->tag) {
                    sc->l.Clear();
                }
            }
            l.RemoveTagged();
        }
    }

    void SPolygon::UvGridTriangulateInto(SMesh *mesh, SSurface *srf) {
        SEdgeList orig = {};
        MakeEdgesInto(&orig);

        SEdgeList holes = {};

        normal = {0, 0, 1};
        FixContourDirections();

        // Build a rectangular grid, with horizontal and vertical lines in the
        // uv plane. The spacing of these lines is adaptive, so calculate that.
        List<double> li, lj;
        li = {};
        lj = {};
        double v[5] = {0.0, 0.25, 0.5, 0.75, 1.0};
        li.Add(&v[0]);
        srf->MakeTriangulationGridInto(&li, 0, 1, /*swapped=*/true, 0);
        lj.Add(&v[0]);
        srf->MakeTriangulationGridInto(&lj, 0, 1, /*swapped=*/false, 0);

        // force 2nd order grid to have at least 4 segments in each direction
        if ((li.n < 5) && (srf->degm>1)) { // 4 segments minimum
            li.Clear();
            li.Add(&v[0]);li.Add(&v[1]);li.Add(&v[2]);li.Add(&v[3]);li.Add(&v[4]);
        }
        if ((lj.n < 5) && (srf->degn>1)) { // 4 segments minimum
            lj.Clear();
            lj.Add(&v[0]);lj.Add(&v[1]);lj.Add(&v[2]);lj.Add(&v[3]);lj.Add(&v[4]);
        }

        if ((li.n > 3) && (lj.n > 3)) {
            // Now iterate over each quad in the grid. If it's outside the polygon,
            // or if it intersects the polygon, then we discard it. Otherwise we
            // generate two triangles in the mesh, and cut it out of our polygon.
            // Quads around the perimeter would be rejected by AnyEdgeCrossings.
            std::vector<bool> bottom(lj.n, false); // did we use this quad?
            Vector tu = {}, tv = {}; 
            int i, j;
            for(i = 1; i < (li.n-1); i++) {
                bool prev_flag = false;
                for(j = 1; j < (lj.n-1); j++) {
                    bool this_flag = true;
                    double us = li[i], uf = li[i+1],
                        vs = lj[j], vf = lj[j+1];

                    Vector a = {us, vs, 0},
                        b = {us, vf, 0},
                        c = {uf, vf, 0},
                        d = {uf, vs, 0};

                    //  |   d-----c
                    //  |   |     |
                    //  |   |     |
                    //  |   a-----b
                    //  |
                    //  +-------------> j/v axis

                    if( (i==(li.n-2)) || (j==(lj.n-2)) ||
                    orig.AnyEdgeCrossings(a, b, NULL) ||
                    orig.AnyEdgeCrossings(b, c, NULL) ||
                    orig.AnyEdgeCrossings(c, d, NULL) ||
                    orig.AnyEdgeCrossings(d, a, NULL))
                    {
                        this_flag = false;
                    }

                    // There's no intersections, so it doesn't matter which point
                    // we decide to test.
                    if(!this->ContainsPoint(a)) {
                        this_flag = false;
                    }
                    
                    if (this_flag) {
                        // Add the quad to our mesh
                        srf->TangentsAt(us,vs, &tu,&tv);
                        if(tu.Dot(tv) < LENGTH_EPS) {
                            /* Split "the other way" if angle>90
                            compare to LENGTH_EPS instead of zero to avoid alternating triangle
                            "orientations" when the tangents are orthogonal (revolve, lathe etc.)
                            this results in a higher quality mesh. */
                            STriangle tr = {};
                            tr.a = a;
                            tr.b = b;
                            tr.c = c;
                            mesh->AddTriangle(&tr);
                            tr.a = a;
                            tr.b = c;
                            tr.c = d;
                            mesh->AddTriangle(&tr);
                        } else{
                            STriangle tr = {};
                            tr.a = a;
                            tr.b = b;
                            tr.c = d;
                            mesh->AddTriangle(&tr);
                            tr.a = b;
                            tr.b = c;
                            tr.c = d;
                            mesh->AddTriangle(&tr);
                        }
                        if (!prev_flag) // add our own left edge
                            holes.AddEdge(d, a);
                        if (!bottom[j]) // add our own bottom edge
                            holes.AddEdge(a, b);
                    } else {
                        if (prev_flag) // add our left neighbors right edge
                            holes.AddEdge(a, d);
                        if (bottom[j]) // add our bottom neighbors top edge
                            holes.AddEdge(b, a);
                    }
                    prev_flag = this_flag;
                    bottom[j] = this_flag;
                }
            }

            // Because no duplicate edges were created we do not need to cull them.
            SPolygon hp = {};
            holes.AssemblePolygon(&hp, NULL, /*keepDir=*/true);

            SContour *sc;
            for(sc = hp.l.First(); sc; sc = hp.l.NextAfter(sc)) {
                l.Add(sc);
            }
            hp.l.Clear();
        }
        orig.Clear();
        holes.Clear();
        li.Clear();
        lj.Clear();
        UvTriangulateInto(mesh, srf);
    }

    void SPolygon::TriangulateInto(SMesh *m) const {
        Vector n = normal;
        if(n.Equals({0, 0, 0})) {
        n = ComputeNormal();
        }
        Vector u = n.Normal(0);
        Vector v = n.Normal(1);

        SPolygon p = {};
        this->InverseTransformInto(&p, u, v, n);

        SSurface srf = SSurface::FromPlane({0.0, 0.0, 0.0},
                                        {1.0, 0.0, 0.0},
                                        {0.0, 1.0, 0.0});
        SMesh pm = {};
        p.UvTriangulateInto(&pm, &srf);
        for(STriangle st : pm.l) {
            st = st.Transform(u, v, n);
            m->AddTriangle(&st);
        }

        p.Clear();
        pm.Clear();
    }

}//ss