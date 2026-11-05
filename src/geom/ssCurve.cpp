//-----------------------------------------------------------------------------
// Anything involving curves and sets of curves (except for the real math,
// which is in ratpoly.cpp).
//
// Copyright 2008-2013 Jonathan Westhues.
//-----------------------------------------------------------------------------
#include "ssCurve.h"

namespace SolveSpace {

void SCurve::GetAxisAlignedBounding(Vector *ptMax, Vector *ptMin) const {
    *ptMax = {VERY_NEGATIVE, VERY_NEGATIVE, VERY_NEGATIVE};
    *ptMin = {VERY_POSITIVE, VERY_POSITIVE, VERY_POSITIVE};

    for(int i = 0; i <= exact.deg; i++) {
        exact.ctrl[i].MakeMaxMin(ptMax, ptMin);
    }
}

SCurve SCurve::FromTransformationOf(SCurve *a, Vector t,
                                    Quaternion q, double scale)
{
    bool needRotate    = !EXACT(q.vx == 0.0 && q.vy == 0.0 && q.vz == 0.0 && q.w == 1.0);
    bool needTranslate = !EXACT(t.x  == 0.0 && t.y  == 0.0 && t.z  == 0.0);
    bool needScale     = !EXACT(scale == 1.0);

    SCurve ret = {};
    ret.h = a->h;
    ret.isExact = a->isExact;
    ret.exact = (a->exact).TransformedBy(t, q, scale);
    ret.surfA = a->surfA;
    ret.surfB = a->surfB;

    SCurvePt *p;
    ret.pts.ReserveMore(a->pts.n);
    for(p = a->pts.First(); p; p = a->pts.NextAfter(p)) {
        SCurvePt pp = *p;
        if(needScale) {
            pp.p = (pp.p).ScaledBy(scale);
        }
        if(needRotate) {
            pp.p = q.Rotate(pp.p);
        }
        if(needTranslate) {
            pp.p = pp.p.Plus(t);
        }
        ret.pts.Add(&pp);
    }
    return ret;
}

void SCurve::Clear() {
    pts.Clear();
}

SSurface *SCurve::GetSurfaceA(SShell *a, SShell *b) const {
    if(source == Source::A) {
        return a->surface.FindById(surfA);
    } else if(source == Source::B) {
        return b->surface.FindById(surfA);
    } else if(source == Source::INTERSECTION) {
        return a->surface.FindById(surfA);
    } else ssassert(false, "Unexpected curve source");
}

SSurface *SCurve::GetSurfaceB(SShell *a, SShell *b) const {
    if(source == Source::A) {
        return a->surface.FindById(surfB);
    } else if(source == Source::B) {
        return b->surface.FindById(surfB);
    } else if(source == Source::INTERSECTION) {
        return b->surface.FindById(surfB);
    } else ssassert(false, "Unexpected curve source");
}

//-----------------------------------------------------------------------------
// When we split line segments wherever they intersect a surface, we introduce
// extra pwl points. This may create very short edges that could be removed
// without violating the chord tolerance. Those are ugly, and also break
// stuff in the Booleans. So remove them.
//-----------------------------------------------------------------------------
void SCurve::RemoveShortSegments(SSurface *srfA, SSurface *srfB) {
    if(pts.n <= 2) return;
    pts.ClearTags();

    Vector prev = pts[0].p;
    double tprev = 0;
    double t = 0;
    double tnext = 0;

    int i, a;
    for(i = 1; i < pts.n - 1; i++) {
        SCurvePt *sct = &(pts[i]),
                 *scn = &(pts[i+1]);

        if(sct->vertex) {
            prev = sct->p;
            continue;
        }

        // if the curve is exact and points are >0.05 apart wrt t, point is there
        // deliberately regardless of chord tolerance (ex: small circles)
        tprev = t = tnext = 0;
        if (isExact) {
            exact.ClosestPointTo(prev, &tprev, /*mustconverge=*/ true);
            exact.ClosestPointTo(sct->p, &t, /*mustconverge=*/ true);
            exact.ClosestPointTo(scn->p, &tnext, /*mustconverge=*/ true);
        }
        if ( (t - tprev > 0.05) && (tnext - t > 0.05) ) {
            prev = sct->p;
            continue;
        }

        bool mustKeep = false;

        // We must check against both surfaces; the piecewise linear edge
        // may have a different chord tolerance in the two surfaces. (For
        // example, a circle in the surface of a cylinder is just a straight
        // line, so it always has perfect chord tol, but that circle in
        // a plane is a circle so it doesn't).
        for(a = 0; a < 2; a++) {
            SSurface *srf = (a == 0) ? srfA : srfB;
            Vector puv, nuv;
            srf->ClosestPointTo(prev,   &(puv.x), &(puv.y));
            srf->ClosestPointTo(scn->p, &(nuv.x), &(nuv.y));

            if(srf->ChordToleranceForEdge(nuv, puv) > SS.ChordTolMm() ) {
                mustKeep = true;
            }
        }

        if(mustKeep) {
            prev = sct->p;
        } else {
            sct->tag = 1;
            // and prev is unchanged, since there's no longer any point
            // in between
        }
    }

    pts.RemoveTagged();
}

STrimBy STrimBy::EntireCurve(SShell *shell, hSCurve hsc, bool backwards) {
    STrimBy stb = {};
    stb.curve = hsc;
    SCurve *sc = shell->curve.FindById(hsc);

    if(backwards) {
        stb.finish = sc->pts[0].p;
        stb.start = sc->pts.Last()->p;
        stb.backwards = true;
    } else {
        stb.start = sc->pts[0].p;
        stb.finish = sc->pts.Last()->p;
        stb.backwards = false;
    }

    return stb;
}


// We will be inserting other curve vertices into our curves to split them.
// This is helpful when curved surfaces become tangent along a trim and the
// usual tests for curve-surface intersection don't split the curve at a vertex.
// This is faster than the previous version that split at surface corners and
// handles more buggy cases. It's not clear this is the best way but it works ok.
static void FindVertsOnCurve(List<SInter> *l, const SCurve *curve, SShell *sh) {

    Vector amax, amin;
    curve->GetAxisAlignedBounding(&amax, &amin);

    // The vertices of the shell are the endpoints of its trims, not the
    // endpoints of its curves: an exact intersection curve may extend past
    // the real geometry on both sides (e.g. out to the padded bounds of a
    // surface that got enlarged by SShell::MergeCoincidentSurfaces()), and
    // splitting some other curve at such a phantom point--which is not a
    // vertex of the trims adjacent across that other curve--produces
    // T-junctions and naked edges in the triangulated shell (issue #1452).
    for(const SSurface &ss : sh->surface) {
        for(const STrimBy &stb : ss.trim) {
            for(int i = 0; i < 2; i++) {
                Vector pt = (i == 0) ? stb.start : stb.finish;
                if(pt.OutsideAndNotOn(amax, amin)) continue;

                double t;
                curve->exact.ClosestPointTo(pt, &t, /*mustConverge=*/false);
                double d = pt.Minus(curve->exact.PointAt(t)).Magnitude();
                if((t > LENGTH_EPS) && (t < (1.0 - LENGTH_EPS)) && (d < LENGTH_EPS)) {
                    SInter inter = {};
                    inter.p = pt;
                    l->Add(&inter);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
// Take our original pwl curve. Wherever an edge intersects a surface within
// either agnstA or agnstB, split the piecewise linear element. Then refine
// the intersection so that it lies on all three relevant surfaces: the
// intersecting surface, srfA, and srfB. (So the pwl curve should lie at
// the intersection of srfA and srfB.) Return a new pwl curve with everything
// split.
//-----------------------------------------------------------------------------
SCurve SCurve::MakeCopySplitAgainst(SShell *agnstA, SShell *agnstB,
                                    SSurface *srfA, SSurface *srfB) const
{
    SCurve ret;
    ret = *this;
    ret.pts = {};

    // First find any vertex that lies on our curve.
    List<SInter> vertpts = {};
    if(isExact) {
        if(agnstA)
            FindVertsOnCurve(&vertpts, this, agnstA);
        if(agnstB)
            FindVertsOnCurve(&vertpts, this, agnstB);
    }
    
    const SCurvePt *p = pts.First();
    ssassert(p != NULL, "Cannot split an empty curve");
    SCurvePt prev = *p;
    ret.pts.Add(p);
    p = pts.NextAfter(p);
            
    for(; p; p = pts.NextAfter(p)) {
        List<SInter> il = {};

        // Find all the intersections with the two passed shells
        if(agnstA)
            agnstA->AllPointsIntersecting(prev.p, p->p, &il,
                /*asSegment=*/true, /*trimmed=*/false, /*inclTangent=*/true);
        if(agnstB)
            agnstB->AllPointsIntersecting(prev.p, p->p, &il,
                /*asSegment=*/true, /*trimmed=*/false, /*inclTangent=*/true);

        if(!il.IsEmpty()) {
            // The intersections were generated by intersecting the pwl
            // edge against a surface; so they must be refined to lie
            // exactly on the original curve.
            il.ClearTags();
            SInter *pi;
            for(pi = il.First(); pi; pi = il.NextAfter(pi)) {
                if(pi->srf == srfA || pi->srf == srfB) {
                    // The edge certainly intersects the surfaces that it
                    // trims (at its endpoints), but those ones don't count.
                    // They are culled later, but no sense calculating them
                    // and they will cause numerical problems (since two
                    // of the three surfaces they're refined to lie on will
                    // be identical, so the matrix will be singular).
                    pi->tag = 1;
                    continue;
                }

                Point2d puv;
                (pi->srf)->ClosestPointTo(pi->p, &puv, /*mustConverge=*/false);

                // Split the edge if the intersection lies within the surface's
                // trim curves, or within the chord tol of the trim curve; want
                // some slop if points are close to edge and pwl is too coarse,
                // and it doesn't hurt to split unnecessarily.
                Point2d dummy = { 0, 0 };
                SBspUv::Class c = (pi->srf->bsp) ? pi->srf->bsp->ClassifyPoint(puv, dummy, pi->srf) : SBspUv::Class::OUTSIDE;
                if(c == SBspUv::Class::OUTSIDE) {
                    double d = VERY_POSITIVE;
                    if(pi->srf->bsp) d = pi->srf->bsp->MinimumDistanceToEdge(puv, pi->srf);
                    if(d > SS.ChordTolMm()) {
                        pi->tag = 1;
                        continue;
                    }
                }

                // We're keeping the intersection, so actually refine it. Finding the intersection
                // to within EPS is important to match the ends of different chopped trim curves.
                // The general 3-surface intersection fails to refine for trims where surfaces
                // are tangent at the curve, but those trims are usually exact, so…
                if(isExact) {
                    (pi->srf)->PointOnCurve(&exact, &(puv.x), &(puv.y));
                } else {
                    (pi->srf)->PointOnSurfaces(srfA, srfB, &(puv.x), &(puv.y));
                }
                pi->p = (pi->srf)->PointAt(puv);
            }
            il.RemoveTagged();
        }
        // Now add any vertex that is on this segment
        const Vector lineStart     = prev.p;
        const Vector lineDirection = (p->p).Minus(prev.p);
        for(const auto &vtx : vertpts) {
            double t = (vtx.p.Minus(lineStart)).DivProjected(lineDirection);
            if((0.0 < t) && (t < 1.0)) {
                il.Add(&vtx);
            }
        }
        if(!il.IsEmpty()) {
            SInter *pi;

            // And now sort them in order along the line. Note that we must
            // do that after refining, in case the refining would make two
            // points switch places.
            std::sort(il.begin(), il.end(), [&](const SInter &a, const SInter &b) {
                double ta = (a.p.Minus(lineStart)).DivProjected(lineDirection);
                double tb = (b.p.Minus(lineStart)).DivProjected(lineDirection);

                return (ta < tb);
            });

            // And now uses the intersections to generate our split pwl edge(s)
            Vector prev = {VERY_POSITIVE, 0, 0};
            for(pi = il.First(); pi; pi = il.NextAfter(pi)) {
                // On-edge intersection will generate same split point for
                // both surfaces, so don't create zero-length edge.
                if(!prev.Equals(pi->p)) {
                    SCurvePt scpt;
                    scpt.tag    = 0;
                    scpt.p      = pi->p;
                    scpt.vertex = true;
                    ret.pts.Add(&scpt);
                }
                prev = pi->p;
            }
        }

        il.Clear();
        ret.pts.Add(p);
        prev = *p;
    }
    vertpts.Clear();
    return ret;
}


} // namespace SolveSpace
