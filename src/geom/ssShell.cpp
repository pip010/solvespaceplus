//-----------------------------------------------------------------------------
// Routines for ray-casting: intersecting a line segment or an infinite line
// with a surface or shell. Ray-casting against a shell is used for point-in-
// shell testing, and the intersection of edge line segments against surfaces
// is used to get rough surface-curve intersections, which are later refined
// numerically.
//
// Copyright 2008-2013 Jonathan Westhues.
//-----------------------------------------------------------------------------
#include "ssShell.h"

namespace SolveSpace {

    void SShell::CopyCurvesSplitAgainst(bool opA, SShell *agnst, SShell *into) {
#pragma omp parallel for
    for(int i=0; i<curve.n; i++) {
        SCurve *sc = &curve[i];
        SCurve scn = sc->MakeCopySplitAgainst(agnst, NULL,
                                surface.FindById(sc->surfA),
                                surface.FindById(sc->surfB));
        scn.source = opA ? SCurve::Source::A : SCurve::Source::B;
#pragma omp critical
        {
            hSCurve hsc = into->curve.AddAndAssignId(&scn);
            // And note the new ID so that we can rewrite the trims appropriately
            sc->newH = hsc;
        }
    }
}

//from: raycast.cpp


// Dot product tolerance for perpendicular; this is on the direction cosine,
// so it's about 0.001 degrees.
const double SShell::DOTP_TOL = 1e-5;

// TODO wtf?
extern int FLAG;

void SShell::MakeFromUnionOf(SShell *a, SShell *b) {
    MakeFromBoolean(a, b, SSurface::CombineAs::UNION);
}

void SShell::MakeFromDifferenceOf(SShell *a, SShell *b) {
    MakeFromBoolean(a, b, SSurface::CombineAs::DIFFERENCE);
}

void SShell::MakeFromIntersectionOf(SShell *a, SShell *b) {
    MakeFromBoolean(a, b, SSurface::CombineAs::INTERSECTION);
}

void SShell::AllPointsIntersecting(Vector a, Vector b,
                                   List<SInter> *il,
                                   bool asSegment, bool trimmed, bool inclTangent)
{
    for(SSurface &ss : surface) {
        ss.AllPointsIntersecting(a, b, il,
            asSegment, trimmed, inclTangent);
    }
}



SShell::Class SShell::ClassifyRegion(Vector edge_n, Vector inter_surf_n,
                           Vector edge_surf_n) const
{
    double dot = inter_surf_n.DirectionCosineWith(edge_n);
    if(fabs(dot) < DOTP_TOL) {
        // The edge's surface and the edge-on-face surface
        // are coincident. Test the edge's surface normal
        // to see if it's with same or opposite normals.
        if(inter_surf_n.Dot(edge_surf_n) > 0) {
            return Class::SURF_COINC_SAME;
        } else {
            return Class::SURF_COINC_OPP;
        }
    } else if(dot > 0) {
        return Class::SURF_OUTSIDE;
    } else {
        return Class::SURF_INSIDE;
    }
}

//-----------------------------------------------------------------------------
// Does the given point lie on our shell? There are many cases; inside and
// outside are obvious, but then there's all the edge-on-edge and edge-on-face
// possibilities.
//
// To calculate, we intersect a ray through p with our shell, and classify
// using the closest intersection point. If the ray hits a surface on edge,
// then just reattempt in a different random direction.
//-----------------------------------------------------------------------------

// table of vectors in 6 arbitrary directions covering 4 of the 8 octants.
// use overlapping sets of 3 to reduce memory usage.
static const double Random[8] = {1.278, 5.0103, 9.427, -2.331, 7.13, 2.954, 5.034, -4.777};
 
//-----------------------------------------------------------------------------
// A face of our shell meets the surface being classified edge-on, with its
// normal parallel to that surface's normal at p, so the face lies in that
// surface's tangent plane at p. It might be coincident with the surface there,
// or it might merely be tangent to it at the shell's edge and curve away from
// it past that edge. The normals at p can't tell those apart, so probe the
// face's actual geometry a little way in from the shell's edge, and report how
// far it then deviates from the tangent plane (measured along nu, a unit
// normal to the surface at p). The direction we probed in, which points into
// the face perpendicular to the shell's edge, is returned in f.
//-----------------------------------------------------------------------------
static double ProbeTangentFace(SSurface *srf, Vector p, Vector inter_edge_n,
                               Vector nu, Vector *f, bool *flat)
{
    // Probe at a small fraction of the face's own size, so that the deviation
    // of a curved face is well above numerical noise no matter the model's
    // scale.
    Vector c00 = srf->ctrl[0][0],
           cm0 = srf->ctrl[srf->degm][0],
           c0n = srf->ctrl[0][srf->degn],
           cmn = srf->ctrl[srf->degm][srf->degn];
    double size = max(max((cm0.Minus(c00)).Magnitude(),
                          (c0n.Minus(c00)).Magnitude()),
                      (cmn.Minus(c00)).Magnitude());
    // inter_edge_n points out of the face, so negate it to probe into it.
    *f = (inter_edge_n.ScaledBy(-1)).WithMagnitude(max(size / 20, LENGTH_EPS * 100));

    Point2d fuv;
    srf->ClosestPointTo(p.Plus(*f), &fuv, /*mustConverge=*/false);
    double dev = ((srf->PointAt(fuv)).Minus(p)).Dot(nu);
    *flat = fabs(dev) < max(LENGTH_EPS, 1e-3 * f->Magnitude());
    return dev;
}

bool SShell::ClassifyEdge(Class *indir, Class *outdir,
                          Vector ea, Vector eb,
                          Vector p,
                          Vector edge_n_in, Vector edge_n_out, Vector surf_n)
{
    List<SInter> l = {};

    // First, check for edge-on-edge
    int edge_inters = 0;
    Vector inter_surf_n[2], inter_edge_n[2];
    SSurface *inter_srf[2] = {};
    for(SSurface &srf : surface) {
        if(srf.LineEntirelyOutsideBbox(ea, eb, /*asSegment=*/true)) continue;

        SEdgeList *sel = &(srf.edges);
        SEdge *se;
        for(se = sel->l.First(); se; se = sel->l.NextAfter(se)) {
            if((ea.Equals(se->a) && eb.Equals(se->b)) ||
               (eb.Equals(se->a) && ea.Equals(se->b)) ||
                p.OnLineSegment(se->a, se->b))
            {
                if(edge_inters < 2) {
                    // Edge-on-edge case
                    Point2d pm;
                    srf.ClosestPointTo(p,  &pm, /*mustConverge=*/false);
                    // A vector normal to the surface, at the intersection point
                    inter_surf_n[edge_inters] = srf.NormalAt(pm);
                    // A vector normal to the intersecting edge (but within the
                    // intersecting surface) at the intersection point, pointing
                    // out.
                    inter_edge_n[edge_inters] =
                      (inter_surf_n[edge_inters]).Cross((se->b).Minus((se->a)));
                    inter_srf[edge_inters] = &srf;
                }

                edge_inters++;
            }
        }
    }

    if(edge_inters == 2) {
        //! @todo make this use the appropriate curved normals
        double dotp[2];
        for(int i = 0; i < 2; i++) {
            dotp[i] = edge_n_out.DirectionCosineWith(inter_surf_n[i]);
        }

        if(fabs(dotp[1]) < DOTP_TOL) {
            swap(dotp[0],         dotp[1]);
            swap(inter_surf_n[0], inter_surf_n[1]);
            swap(inter_edge_n[0], inter_edge_n[1]);
            swap(inter_srf[0],    inter_srf[1]);
        }

        // The class of the side of our edge that face 0 extends into, when
        // face 0 lies in our tangent plane at p and face 1 does not. That's
        // coincident with the shell if face 0 really is coincident with our
        // surface, but see below for when it's only tangent to it.
        Class towards_0 =
            (surf_n.Dot(inter_surf_n[0])) > 0 ? Class::SURF_COINC_SAME : Class::SURF_COINC_OPP;

        if(fabs(dotp[0]) < DOTP_TOL && fabs(dotp[1]) >= DOTP_TOL) {
            // Face 0 might be merely tangent to our surface at the shell's
            // edge, curving away from it past that edge, e.g. where the
            // spline face of a cut is tangent to the face it cuts through
            // and the cut's other face crosses that face there (issue
            // #1291). Then the side of our edge towards face 0 is not
            // coincident with the shell at all; it's inside or outside it
            // depending on which way face 0 curves away.
            Vector f;
            bool   flat;
            double dev = ProbeTangentFace(inter_srf[0], p, inter_edge_n[0],
                                          surf_n.WithMagnitude(1), &f, &flat);
            if(!flat) {
                // As in the both-tangent case below: our edge is on the
                // shell's material side of face 0 iff face 0 curves away
                // opposite its own outward normal.
                towards_0 = (dev * surf_n.Dot(inter_surf_n[0]) > 0) ?
                                Class::SURF_INSIDE : Class::SURF_OUTSIDE;
            }
        }

        if(fabs(dotp[0]) < DOTP_TOL && fabs(dotp[1]) < DOTP_TOL) {
            // Both faces meeting at the shell's edge lie in our surface's
            // tangent plane at p. Usually the shell's surface is just
            // split into two pieces joining at our edge, and both sides
            // of our edge are coincident with the shell. But a face may
            // instead be *tangent* to our surface at the shell's edge,
            // curving away from it past that edge (e.g. a fillet cylinder
            // meeting the flat cap it's tangent to, issue #1291); the
            // side of our edge towards such a face is not coincident with
            // the shell. The normals at p can't distinguish these, so
            // probe each face's actual geometry a little way in from the
            // shell's edge.
            Vector f[2];    // points into face i, perpendicular to edge
            double dev[2];  // face i's deviation from our tangent plane
            bool   flat[2]; // does face i lie in our tangent plane?
            Vector nu = surf_n.WithMagnitude(1);
            for(int i = 0; i < 2; i++) {
                dev[i] = ProbeTangentFace(inter_srf[i], p, inter_edge_n[i], nu,
                                          &f[i], &flat[i]);
            }
            Class *dir[2] = { indir, outdir };
            Vector  en[2] = { edge_n_in, edge_n_out };
            for(int s = 0; s < 2; s++) {
                // Of the faces extending into this side of our edge, find
                // the one lying nearest our tangent plane.
                int best = -1;
                for(int i = 0; i < 2; i++) {
                    if(f[i].Dot(en[s]) <= 0) continue;
                    if(best < 0 || fabs(dev[i]) < fabs(dev[best])) best = i;
                }
                if(best >= 0 && flat[best]) {
                    // This side of our edge lies on that face.
                    *dir[s] = (surf_n.Dot(inter_surf_n[best]) > 0) ?
                                Class::SURF_COINC_SAME : Class::SURF_COINC_OPP;
                } else if(best >= 0) {
                    // The nearest face on this side curves away from our
                    // tangent plane, so this side of our edge is not on
                    // the shell; it's on the shell's material side of
                    // that face iff the face curves away opposite the
                    // face's outward normal.
                    *dir[s] = (dev[best] * surf_n.Dot(inter_surf_n[best]) > 0) ?
                                Class::SURF_INSIDE : Class::SURF_OUTSIDE;
                } else {
                    // Both faces fold back to the other side of our edge,
                    // meeting our surface tangentially at the shell's
                    // edge. They enclose a zero-angle wedge, of material
                    // if the farther face lies on the material side of
                    // the nearer one (then this side of our edge is
                    // beyond that material, outside the shell), and of
                    // void cut into material otherwise.
                    int nr = (fabs(dev[0]) < fabs(dev[1])) ? 0 : 1;
                    bool material = (dev[1-nr] - dev[nr]) *
                                        surf_n.Dot(inter_surf_n[nr]) < 0;
                    *dir[s] = material ? Class::SURF_OUTSIDE
                                       : Class::SURF_INSIDE;
                }
            }
        } else if(fabs(dotp[0]) < DOTP_TOL && dotp[1] > DOTP_TOL) {
            if(edge_n_out.Dot(inter_edge_n[0]) > 0) {
                *indir  = towards_0;
                *outdir = Class::SURF_OUTSIDE;
            } else {
                *indir  = Class::SURF_INSIDE;
                *outdir = towards_0;
            }
        } else if(fabs(dotp[0]) < DOTP_TOL && dotp[1] < -DOTP_TOL) {
            if(edge_n_out.Dot(inter_edge_n[0]) > 0) {
                *indir  = towards_0;
                *outdir = Class::SURF_INSIDE;
            } else {
                *indir  = Class::SURF_OUTSIDE;
                *outdir = towards_0;
            }
        } else if(dotp[0] > DOTP_TOL && dotp[1] > DOTP_TOL) {
            *indir  = Class::SURF_INSIDE;
            *outdir = Class::SURF_OUTSIDE;
        } else if(dotp[0] < -DOTP_TOL && dotp[1] < -DOTP_TOL) {
            *indir  = Class::SURF_OUTSIDE;
            *outdir = Class::SURF_INSIDE;
        } else {
            // Our edge is tangent to the shell at the shell's edge: the
            // two faces meeting at the shell's edge lie on the same side
            // of our surface, so the regions on both sides of our edge
            // get the same class. That's outside the shell if the shell's
            // edge is convex, and inside if it's concave (reflex).
            // inter_edge_n[i] points away from face i's material, so the
            // shell's edge is convex iff face 1's material direction lies
            // behind face 0's plane.
            Class c = (inter_edge_n[1].Dot(inter_surf_n[0]) > 0) ?
                            Class::SURF_OUTSIDE : Class::SURF_INSIDE;
            *indir  = c;
            *outdir = c;
        }
        return true;
    }

    if(edge_inters >= 4 && edge_inters % 2 == 0) {
        // More than two surfaces of our shell meet edge-on-edge along the
        // given edge: regions of the shell join along a knife edge there
        // (e.g. two prisms sharing a single edge, issue #1091). The faces
        // along the knife edge divide the directions perpendicular to the
        // edge at p into sectors that are alternately inside and outside
        // the shell, so classify each side of our edge against the face
        // that lies angularly closest to it: that face bounds the sector
        // containing it, and its normal points out of that sector iff the
        // sector is outside the shell.
        Vector edge_d = (eb.Minus(ea)).WithMagnitude(1);
        Vector en[2] = { edge_n_in, edge_n_out };
        Vector ep[2];
        double bestdot[2] = { VERY_NEGATIVE, VERY_NEGATIVE };
        Vector bestn[2];
        for(int i = 0; i < 2; i++) {
            ep[i] = (en[i].Minus(edge_d.ScaledBy(en[i].Dot(edge_d))))
                        .WithMagnitude(1);
        }

        for(SSurface &srf : surface) {
            if(srf.LineEntirelyOutsideBbox(ea, eb, /*asSegment=*/true)) continue;

            SEdgeList *sel = &(srf.edges);
            SEdge *se;
            for(se = sel->l.First(); se; se = sel->l.NextAfter(se)) {
                if((ea.Equals(se->a) && eb.Equals(se->b)) ||
                   (eb.Equals(se->a) && ea.Equals(se->b)) ||
                    p.OnLineSegment(se->a, se->b))
                {
                    Point2d pm;
                    srf.ClosestPointTo(p, &pm, /*mustConverge=*/false);
                    Vector n = srf.NormalAt(pm);
                    // A vector tangent to the intersecting surface at p,
                    // pointing into that face, perpendicular to our edge.
                    Vector f = ((se->b).Minus(se->a)).Cross(n);
                    f = f.Minus(edge_d.ScaledBy(f.Dot(edge_d)));
                    if(f.Magnitude() < LENGTH_EPS) continue;
                    f = f.WithMagnitude(1);
                    for(int i = 0; i < 2; i++) {
                        double dot = ep[i].Dot(f);
                        if(dot > bestdot[i]) {
                            bestdot[i] = dot;
                            bestn[i]   = n;
                        }
                    }
                }
            }
        }

        Class *dir[2] = { indir, outdir };
        for(int i = 0; i < 2; i++) {
            double dotn = en[i].DirectionCosineWith(bestn[i]);
            if(bestdot[i] > 0 && fabs(dotn) < DOTP_TOL) {
                // This side of our edge lies along the nearest face, so
                // our surface is locally coincident with it there.
                *dir[i] = (surf_n.Dot(bestn[i]) > 0) ? Class::SURF_COINC_SAME
                                                     : Class::SURF_COINC_OPP;
            } else {
                *dir[i] = (dotn > 0) ? Class::SURF_OUTSIDE
                                     : Class::SURF_INSIDE;
            }
        }
        return true;
    }

    if(edge_inters != 0) dbp("bad, edge_inters=%d", edge_inters);

    // Next, check for edge-on-surface. The ray-casting for edge-inside-shell
    // would catch this too, but test separately, for speed (since many edges
    // are on surface) and for numerical stability, so we don't pick up
    // the additional error from the line intersection.

    for(SSurface &srf : surface) {
        if(srf.LineEntirelyOutsideBbox(ea, eb, /*asSegment=*/true)) continue;

        Point2d puv;
        srf.ClosestPointTo(p, &(puv.x), &(puv.y), /*mustConverge=*/false);
        Vector pp = srf.PointAt(puv);

        if((pp.Minus(p)).Magnitude() > LENGTH_EPS) continue;
        Point2d dummy = { 0, 0 };
        SBspUv::Class c = (srf.bsp) ? srf.bsp->ClassifyPoint(puv, dummy, &srf) : SBspUv::Class::OUTSIDE;
        if(c == SBspUv::Class::OUTSIDE) continue;

        // Edge-on-face (unless edge-on-edge above superseded)
        Point2d pin, pout;
        srf.ClosestPointTo(p.Plus(edge_n_in),  &pin,  /*mustConverge=*/false);
        srf.ClosestPointTo(p.Plus(edge_n_out), &pout, /*mustConverge=*/false);

        Vector surf_n_in  = srf.NormalAt(pin),
               surf_n_out = srf.NormalAt(pout);

        *indir  = ClassifyRegion(edge_n_in,  surf_n_in,  surf_n);
        *outdir = ClassifyRegion(edge_n_out, surf_n_out, surf_n);
        return true;
    }

    // Edge is not on face or on edge; so it's either inside or outside
    // the shell, and we'll determine which by raycasting.
    int cnt = 0;
    for(;;) {
        // Cast a ray in a random direction (two-sided so that we test if
        // the point lies on a surface, but use only one side for in/out
        // testing)
        Vector ray = {Random[cnt], Random[cnt+1], Random[cnt+2]};

        AllPointsIntersecting(
            p.Minus(ray), p.Plus(ray), &l,
                /*asSegment=*/false, /*trimmed=*/true, /*inclTangent=*/false);

        // no intersections means it's outside
        *indir  = Class::SURF_OUTSIDE;
        *outdir = Class::SURF_OUTSIDE;
        double dmin = VERY_POSITIVE;
        bool onEdge = false;
        edge_inters = 0;

        SInter *si;
        for(si = l.First(); si; si = l.NextAfter(si)) {
            double t = ((si->p).Minus(p)).DivProjected(ray);
            if(t*ray.Magnitude() < -LENGTH_EPS) {
                // wrong side, doesn't count
                continue;
            }

            double d = ((si->p).Minus(p)).Magnitude();

            // We actually should never hit this case; it should have been
            // handled above.
            if(d < LENGTH_EPS && si->onEdge) {
                edge_inters++;
            }

            if(d < dmin) {
                dmin = d;
                // Edge does not lie on surface; either strictly inside
                // or strictly outside
                if((si->surfNormal).Dot(ray) > 0) {
                    *indir  = Class::SURF_INSIDE;
                    *outdir = Class::SURF_INSIDE;
                } else {
                    *indir  = Class::SURF_OUTSIDE;
                    *outdir = Class::SURF_OUTSIDE;
                }
                onEdge = si->onEdge;
            }
        }
        l.Clear();

        // If the point being tested lies exactly on an edge of the shell,
        // then our ray always lies on edge, and that's okay. Otherwise
        // try again in a different random direction.
        if(!onEdge) break;
        cnt++;
        if(cnt > 5) {
            dbp("can't find a ray that doesn't hit on edge!");
            dbp("on edge = %d, edge_inters = %d", onEdge, edge_inters);
            SS.nakedEdges.AddEdge(ea, eb);
            break;
        }
    }

    return true;
}





  // ------------------- MERGE ---------------------
    //-----------------------------------------------------------------------------
// Routines to merge multiple coincident surfaces (each with their own trim
// curves) into a single surface, with all of the trim curves.
//
// Copyright 2008-2013 Jonathan Westhues.
//-----------------------------------------------------------------------------

void SShell::MergeCoincidentSurfaces() {
    surface.ClearTags();

    int i, j;
    SSurface *si, *sj;

    for(i = 0; i < surface.n; i++) {
        si = &(surface[i]);
        if(si->tag) continue;
        // Let someone else clean up the empty surfaces; we can certainly merge
        // them, but we don't know how to calculate a reasonable bounding box.
        if(si->trim.IsEmpty())
            continue;
        // And for now we handle only coincident planes, so no sense wasting
        // time on other surfaces.
        if(si->degm != 1 || si->degn != 1) continue;

        SEdgeList sel = {};
        si->MakeEdgesInto(this, &sel, SSurface::MakeAs::XYZ);

        bool mergedThisTime, merged = false;
        do {
            mergedThisTime = false;

            for(j = i + 1; j < surface.n; j++) {
                sj = &(surface[j]);
                if(sj->tag) continue;
                if(!sj->CoincidentWith(si, /*sameNormal=*/true)) continue;
                if(!sj->color.Equals(si->color)) continue;
                // But we do merge surfaces with different face entities, since
                // otherwise we'd hardly ever merge anything.

                // This surface is coincident. But let's not merge coincident
                // surfaces if they contain disjoint contours; that just makes
                // the bounding box tests less effective, and possibly things
                // less robust.
                SEdgeList tel = {};
                sj->MakeEdgesInto(this, &tel, SSurface::MakeAs::XYZ);
                if(!sel.ContainsEdgeFrom(&tel)) {
                    tel.Clear();
                    continue;
                }
                tel.Clear();

                sj->tag = 1;
                merged = true;
                mergedThisTime = true;
                sj->MakeEdgesInto(this, &sel, SSurface::MakeAs::XYZ);
                sj->trim.Clear();

                // All the references to this surface get replaced with the
                // new srf
                for(SCurve &sc : curve) {
                    if(sc.surfA == sj->h) sc.surfA = si->h;
                    if(sc.surfB == sj->h) sc.surfB = si->h;
                }
            }

            // If this iteration merged a contour onto ours, then we have to
            // go through the surfaces again; that might have made a new
            // surface touch us.
        } while(mergedThisTime);

        if(merged) {
            sel.CullExtraneousEdges();
            si->trim.Clear();
            si->TrimFromEdgeList(&sel, /*asUv=*/false);

            // And we must choose control points such that all the trims lie
            // with u and v in [0, 1], so that the bbox tests work.
            Vector u, v, n;
            si->TangentsAt(0.5, 0.5, &u, &v);
            u = u.WithMagnitude(1);
            v = v.WithMagnitude(1);
            n = si->NormalAt(0.5, 0.5).WithMagnitude(1);
            v = (n.Cross(u)).WithMagnitude(1);

            double umax = VERY_NEGATIVE, umin = VERY_POSITIVE,
                   vmax = VERY_NEGATIVE, vmin = VERY_POSITIVE;
            SEdge *se;
            for(se = sel.l.First(); se; se = sel.l.NextAfter(se)) {
                double ut = (se->a).Dot(u), vt = (se->a).Dot(v);
                umax = max(umax, ut);
                vmax = max(vmax, vt);
                umin = min(umin, ut);
                vmin = min(vmin, vt);
            }

            // An interesting problem here; the real curve could extend
            // slightly beyond the bounding box of the piecewise linear
            // bits. Not a problem for us, but some apps won't import STEP
            // in that case. So give a bit of extra room; in theory just
            // a chord tolerance, but more can't hurt.
            double muv = max((umax - umin), (vmax - vmin));
            double tol = muv/50 + 3*SS.ChordTolMm();
            umax += tol;
            vmax += tol;
            umin -= tol;
            vmin -= tol;

            // We move in the +v direction as v goes from 0 to 1, and in the
            // +u direction as u goes from 0 to 1. So our normal ends up
            // pointed the same direction.
            double nt = (si->ctrl[0][0]).Dot(n);
            si->ctrl[0][0] =
                Vector::From(umin, vmin, nt).ScaleOutOfCsys(u, v, n);
            si->ctrl[0][1] =
                Vector::From(umin, vmax, nt).ScaleOutOfCsys(u, v, n);
            si->ctrl[1][1] =
                Vector::From(umax, vmax, nt).ScaleOutOfCsys(u, v, n);
            si->ctrl[1][0] =
                Vector::From(umax, vmin, nt).ScaleOutOfCsys(u, v, n);
        }
        sel.Clear();
    }

    surface.RemoveTagged();
}


//-----------------------------------------------------------------------------
// Anything involving NURBS shells (i.e., shells); except
// for the real math, which is in ratpoly.cpp.
//
// Copyright 2008-2013 Jonathan Westhues.
//-----------------------------------------------------------------------------

typedef struct {
    hSCurve     hc;
    hSSurface   hs;
} TrimLine;


void SShell::MakeFromExtrusionOf(SBezierLoopSet *sbls, Vector t0, Vector t1, RgbaColor color)
{
    // Make the extrusion direction consistent with respect to the normal
    // of the sketch we're extruding.
    if((t0.Minus(t1)).Dot(sbls->normal) < 0) {
        swap(t0, t1);
    }

    // Define a coordinate system to contain the original sketch, and get
    // a bounding box in that csys
    Vector n = sbls->normal.ScaledBy(-1);
    Vector u = n.Normal(0), v = n.Normal(1);
    Vector orig = sbls->point;
    double umax = VERY_NEGATIVE, umin = VERY_POSITIVE;
    sbls->GetBoundingProjd(u, orig, &umin, &umax);
    double vmax = VERY_NEGATIVE, vmin = VERY_POSITIVE;
    sbls->GetBoundingProjd(v, orig, &vmin, &vmax);
    // and now fix things up so that all u and v lie between 0 and 1
    orig = orig.Plus(u.ScaledBy(umin));
    orig = orig.Plus(v.ScaledBy(vmin));
    u = u.ScaledBy(umax - umin);
    v = v.ScaledBy(vmax - vmin);

    // So we can now generate the top and bottom surfaces of the extrusion,
    // planes within a translated (and maybe mirrored) version of that csys.
    SSurface s0, s1;
    s0 = SSurface::FromPlane(orig.Plus(t0), u, v);
    s0.color = color;
    s1 = SSurface::FromPlane(orig.Plus(t1).Plus(u), u.ScaledBy(-1), v);
    s1.color = color;
    hSSurface hs0 = surface.AddAndAssignId(&s0),
              hs1 = surface.AddAndAssignId(&s1);

    // Now go through the input curves. For each one, generate its surface
    // of extrusion, its two translated trim curves, and one trim line. We
    // go through by loops so that we can assign the lines correctly.
    SBezierLoop *sbl;
    for(sbl = sbls->l.First(); sbl; sbl = sbls->l.NextAfter(sbl)) {
        SBezier *sb;
        List<TrimLine> trimLines = {};

        for(sb = sbl->l.First(); sb; sb = sbl->l.NextAfter(sb)) {
            // Generate the surface of extrusion of this curve, and add
            // it to the list
            SSurface ss = SSurface::FromExtrusionOf(sb, t0, t1);
            ss.color = color;
            hSSurface hsext = surface.AddAndAssignId(&ss);

            // Translate the curve by t0 and t1 to produce two trim curves
            SCurve sc = {};
            sc.isExact = true;
            sc.exact = sb->TransformedBy(t0, Quaternion::IDENTITY, 1.0);
            (sc.exact).MakePwlInto(&(sc.pts));
            sc.surfA = hs0;
            sc.surfB = hsext;
            hSCurve hc0 = curve.AddAndAssignId(&sc);

            sc = {};
            sc.isExact = true;
            sc.exact = sb->TransformedBy(t1, Quaternion::IDENTITY, 1.0);
            (sc.exact).MakePwlInto(&(sc.pts));
            sc.surfA = hs1;
            sc.surfB = hsext;
            hSCurve hc1 = curve.AddAndAssignId(&sc);

            STrimBy stb0, stb1;
            // The translated curves trim the flat top and bottom surfaces.
            stb0 = STrimBy::EntireCurve(this, hc0, /*backwards=*/false);
            stb1 = STrimBy::EntireCurve(this, hc1, /*backwards=*/true);
            (surface.FindById(hs0))->trim.Add(&stb0);
            (surface.FindById(hs1))->trim.Add(&stb1);

            // The translated curves also trim the surface of extrusion.
            stb0 = STrimBy::EntireCurve(this, hc0, /*backwards=*/true);
            stb1 = STrimBy::EntireCurve(this, hc1, /*backwards=*/false);
            (surface.FindById(hsext))->trim.Add(&stb0);
            (surface.FindById(hsext))->trim.Add(&stb1);

            // And form the trim line
            Vector pt = sb->Finish();
            sc = {};
            sc.isExact = true;
            sc.exact = SBezier::From(pt.Plus(t0), pt.Plus(t1));
            (sc.exact).MakePwlInto(&(sc.pts));
            hSCurve hl = curve.AddAndAssignId(&sc);
            // save this for later
            TrimLine tl;
            tl.hc = hl;
            tl.hs = hsext;
            trimLines.Add(&tl);
        }

        int i;
        for(i = 0; i < trimLines.n; i++) {
            TrimLine *tl = &(trimLines[i]);
            SSurface *ss = surface.FindById(tl->hs);

            TrimLine *tlp = &(trimLines[WRAP(i-1, trimLines.n)]);

            STrimBy stb;
            stb = STrimBy::EntireCurve(this, tl->hc, /*backwards=*/true);
            ss->trim.Add(&stb);
            stb = STrimBy::EntireCurve(this, tlp->hc, /*backwards=*/false);
            ss->trim.Add(&stb);

            (curve.FindById(tl->hc))->surfA = ss->h;
            (curve.FindById(tlp->hc))->surfB = ss->h;
        }
        trimLines.Clear();
    }
}

bool SShell::CheckNormalAxisRelationship(SBezierLoopSet *sbls, Vector pt, Vector axis, double da, double dx)
// Check that the direction of revolution/extrusion ends up parallel to the normal of
// the sketch, on the side of the axis where the sketch is.
{
    SBezierLoop *sbl;
    Vector pto;
    double md = VERY_NEGATIVE;
    for(sbl = sbls->l.First(); sbl; sbl = sbls->l.NextAfter(sbl)) {
        SBezier *sb;
        for(sb = sbl->l.First(); sb; sb = sbl->l.NextAfter(sb)) {
            // Choose the point farthest from the axis; we'll get garbage
            // if we choose a point that lies on the axis, for example.
            // (And our surface will be self-intersecting if the sketch
            // spans the axis, so don't worry about that.)
            for(int i = 0; i <= sb->deg; i++) {
                Vector p = sb->ctrl[i];
                double d = p.DistanceToLine(pt, axis);
                if(d > md) {
                    md  = d;
                    pto = p;
                }
            }
        }
    }
    Vector ptc = pto.ClosestPointOnLine(pt, axis),
           up = axis.Cross(pto.Minus(ptc)).ScaledBy(da),
           vp = up.Plus(axis.ScaledBy(dx));
   
    return (vp.Dot(sbls->normal) > 0);
}

// sketch must not contain the axis of revolution as a non-construction line for helix
void SShell::MakeFromHelicalRevolutionOf(SBezierLoopSet *sbls, Vector pt, Vector axis,
                                         RgbaColor color, Group *group, double angles,
                                         double anglef, double dists, double distf) {
    int i0 = surface.n; // number of pre-existing surfaces
    SBezierLoop *sbl;
    // for testing - hard code the axial distance, and number of sections.
    // distance will need to be parameters in the future.
    double dist  = distf - dists;
    int sections = (int)(fabs(anglef - angles) / (PI / 2) + 1);
    double wedge = (anglef - angles) / sections;
    int startMapping = Group::REMAP_LATHE_START, endMapping = Group::REMAP_LATHE_END;

    if(CheckNormalAxisRelationship(sbls, pt, axis, anglef-angles, distf-dists)) {
        swap(angles, anglef);
        swap(dists, distf);
        dist  = -dist;
        wedge = -wedge;
        swap(startMapping, endMapping);
    }

    // Define a coordinate system to contain the original sketch, and get
    // a bounding box in that csys
    Vector n = sbls->normal.ScaledBy(-1);
    Vector u = n.Normal(0), v = n.Normal(1);
    Vector orig = sbls->point;
    double umax = VERY_NEGATIVE, umin = VERY_POSITIVE;
    sbls->GetBoundingProjd(u, orig, &umin, &umax);
    double vmax = VERY_NEGATIVE, vmin = VERY_POSITIVE;
    sbls->GetBoundingProjd(v, orig, &vmin, &vmax);
    // and now fix things up so that all u and v lie between 0 and 1
    orig = orig.Plus(u.ScaledBy(umin));
    orig = orig.Plus(v.ScaledBy(vmin));
    u    = u.ScaledBy(umax - umin);
    v    = v.ScaledBy(vmax - vmin);

    // So we can now generate the end caps of the extrusion within
    // a translated and rotated (and maybe mirrored) version of that csys.
    SSurface s0, s1;
    s0 = SSurface::FromPlane(orig.RotatedAbout(pt, axis, angles).Plus(axis.ScaledBy(dists)),
                             u.RotatedAbout(axis, angles), v.RotatedAbout(axis, angles));
    s0.color = color;

    hEntity face0 = group->Remap(Entity::NO_ENTITY, startMapping);
    s0.face = face0.v;

    s1 = SSurface::FromPlane(
        orig.Plus(u).RotatedAbout(pt, axis, anglef).Plus(axis.ScaledBy(distf)),
        u.ScaledBy(-1).RotatedAbout(axis, anglef), v.RotatedAbout(axis, anglef));
    s1.color = color;

    hEntity face1 = group->Remap(Entity::NO_ENTITY, endMapping);
    s1.face = face1.v;

    hSSurface hs0 = surface.AddAndAssignId(&s0);
    hSSurface hs1 = surface.AddAndAssignId(&s1);

    // Now we actually build and trim the swept surfaces. One loop at a time.
    for(sbl = sbls->l.First(); sbl; sbl = sbls->l.NextAfter(sbl)) {
        int i, j;
        SBezier *sb;
        List<std::vector<hSSurface>> hsl = {};

        // This is where all the NURBS are created and Remapped to the generating curve
        for(sb = sbl->l.First(); sb; sb = sbl->l.NextAfter(sb)) {
            std::vector<hSSurface> revs(sections);
            for(j = 0; j < sections; j++) {
                if((dist == 0) && sb->deg == 1 &&
                   (sb->ctrl[0]).DistanceToLine(pt, axis) < LENGTH_EPS &&
                   (sb->ctrl[1]).DistanceToLine(pt, axis) < LENGTH_EPS) {
                    // This is a line on the axis of revolution; it does
                    // not contribute a surface.
                    revs[j].v = 0;
                } else {
                    SSurface ss = SSurface::FromRevolutionOf(
                        sb, pt, axis, angles + (wedge)*j, angles + (wedge) * (j + 1),
                        dists + j * dist / sections, dists + (j + 1) * dist / sections);
                    ss.color = color;
                    if(sb->entity != 0) {
                        hEntity he;
                        he.v          = sb->entity;
                        hEntity hface = group->Remap(he, Group::REMAP_LINE_TO_FACE);
                        if(SK.entity.FindByIdNoOops(hface) != NULL) {
                            ss.face = hface.v;
                        }
                    }
                    revs[j] = surface.AddAndAssignId(&ss);
                }
            }
            hsl.Add(&revs);
        }
        // Still the same loop. Need to create trim curves
        for(i = 0; i < sbl->l.n; i++) {
            std::vector<hSSurface> revs = hsl[i], revsp = hsl[WRAP(i - 1, sbl->l.n)];

            sb = &(sbl->l[i]);

            // we will need the grid t-values for this entire row of surfaces
            List<double> t_values;
            t_values = {};
            if (revs[0].v) { 
                double ps = 0.0;
                t_values.Add(&ps);
                (surface.FindById(revs[0]))->MakeTriangulationGridInto(
                        &t_values, 0.0, 1.0, true, 0);
            }
            // we generate one more curve than we did surfaces
            for(j = 0; j <= sections; j++) {
                SCurve sc;
                Quaternion qs = Quaternion::From(axis, angles + wedge * j);
                // we want Q*(x - p) + p = Q*x + (p - Q*p)
                Vector ts =
                    pt.Minus(qs.Rotate(pt)).Plus(axis.ScaledBy(dists + j * dist / sections));

                // If this input curve generated a surface, then trim that
                // surface with the rotated version of the input curve.
                if(revs[0].v) { // not d[j] because crash on j==sections
                    sc         = {};
                    sc.isExact = true;
                    sc.exact   = sb->TransformedBy(ts, qs, 1.0);
                    // make the PWL for the curve based on t value list
                    for(int x = 0; x < t_values.n; x++) {
                        SCurvePt scpt;
                        scpt.tag    = 0;
                        scpt.p      = sc.exact.PointAt(t_values[x]);
                        scpt.vertex = (x == 0) || (x == (t_values.n - 1));
                        sc.pts.Add(&scpt);
                    }

                    // the surfaces already exists so trim with this curve
                    if(j < sections) {
                        sc.surfA = revs[j];
                    } else {
                        sc.surfA = hs1; // end cap
                    }

                    if(j > 0) {
                        sc.surfB = revs[j - 1];
                    } else {
                        sc.surfB = hs0; // staring cap
                    }

                    hSCurve hcb = curve.AddAndAssignId(&sc);

                    STrimBy stb;
                    stb = STrimBy::EntireCurve(this, hcb, /*backwards=*/true);
                    (surface.FindById(sc.surfA))->trim.Add(&stb);
                    stb = STrimBy::EntireCurve(this, hcb, /*backwards=*/false);
                    (surface.FindById(sc.surfB))->trim.Add(&stb);
                } else if(j == 0) { // curve was on the rotation axis and is shared by the end caps.
                    sc         = {};
                    sc.isExact = true;
                    sc.exact   = sb->TransformedBy(ts, qs, 1.0);
                    (sc.exact).MakePwlInto(&(sc.pts));
                    sc.surfA    = hs1; // end cap
                    sc.surfB    = hs0; // staring cap
                    hSCurve hcb = curve.AddAndAssignId(&sc);

                    STrimBy stb;
                    stb = STrimBy::EntireCurve(this, hcb, /*backwards=*/true);
                    (surface.FindById(sc.surfA))->trim.Add(&stb);
                    stb = STrimBy::EntireCurve(this, hcb, /*backwards=*/false);
                    (surface.FindById(sc.surfB))->trim.Add(&stb);
                }

                // And if this input curve and the one after it both generated
                // surfaces, then trim both of those by the appropriate
                // curve based on the control points.
                if((j < sections) && revs[j].v && revsp[j].v) {
                    SSurface *ss = surface.FindById(revs[j]);

                    sc         = {};
                    sc.isExact = true;
                    sc.exact   = SBezier::From(ss->ctrl[0][0], ss->ctrl[0][1], ss->ctrl[0][2]);
                    sc.exact.weight[1] = ss->weight[0][1];
                    double max_dt = 0.5;
                    if (sc.exact.deg > 1) max_dt = 0.125;
                    (sc.exact).MakePwlInto(&(sc.pts), 0.0, max_dt);
                    sc.surfA = revs[j];
                    sc.surfB = revsp[j];

                    hSCurve hcc = curve.AddAndAssignId(&sc);

                    STrimBy stb;
                    stb = STrimBy::EntireCurve(this, hcc, /*backwards=*/false);
                    (surface.FindById(sc.surfA))->trim.Add(&stb);
                    stb = STrimBy::EntireCurve(this, hcc, /*backwards=*/true);
                    (surface.FindById(sc.surfB))->trim.Add(&stb);
                }
            }
            t_values.Clear();
        }

        hsl.Clear();
    }

    if(dist == 0) {
        MakeFirstOrderRevolvedSurfaces(pt, axis, i0);
    }
}

void SShell::MakeFromRevolutionOf(SBezierLoopSet *sbls, Vector pt, Vector axis, RgbaColor color,
                                  Group *group) {
    int i0 = surface.n; // number of pre-existing surfaces
    SBezierLoop *sbl;

    if(CheckNormalAxisRelationship(sbls, pt, axis, 1.0, 0.0)) {
        axis = axis.ScaledBy(-1);
    }

    // Now we actually build and trim the surfaces.
    for(sbl = sbls->l.First(); sbl; sbl = sbls->l.NextAfter(sbl)) {
        int i, j;
        SBezier *sb;
        List<std::vector<hSSurface>> hsl = {};

        for(sb = sbl->l.First(); sb; sb = sbl->l.NextAfter(sb)) {
            std::vector<hSSurface> revs(4);
            for(j = 0; j < 4; j++) {
                if(sb->deg == 1 &&
                    (sb->ctrl[0]).DistanceToLine(pt, axis) < LENGTH_EPS &&
                    (sb->ctrl[1]).DistanceToLine(pt, axis) < LENGTH_EPS)
                {
                    // This is a line on the axis of revolution; it does
                    // not contribute a surface.
                    revs[j].v = 0;
                } else {
                    SSurface ss = SSurface::FromRevolutionOf(sb, pt, axis, (PI / 2) * j,
                                                             (PI / 2) * (j + 1), 0.0, 0.0);
                    ss.color = color;
                    if(sb->entity != 0) {
                        hEntity he;
                        he.v = sb->entity;
                        hEntity hface = group->Remap(he, Group::REMAP_LINE_TO_FACE);
                        if(SK.entity.FindByIdNoOops(hface) != NULL) {
                            ss.face = hface.v;
                        }
                    }
                    revs[j] = surface.AddAndAssignId(&ss);
                }
            }
            hsl.Add(&revs);
        }

        for(i = 0; i < sbl->l.n; i++) {
            std::vector<hSSurface> revs  = hsl[i],
                     revsp = hsl[WRAP(i-1, sbl->l.n)];

            sb   = &(sbl->l[i]);

            for(j = 0; j < 4; j++) {
                SCurve sc;
                Quaternion qs = Quaternion::From(axis, (PI/2)*j);
                // we want Q*(x - p) + p = Q*x + (p - Q*p)
                Vector ts = pt.Minus(qs.Rotate(pt));

                // If this input curve generate a surface, then trim that
                // surface with the rotated version of the input curve.
                if(revs[j].v) {
                    sc = {};
                    sc.isExact = true;
                    sc.exact = sb->TransformedBy(ts, qs, 1.0);
                    (sc.exact).MakePwlInto(&(sc.pts));
                    sc.surfA = revs[j];
                    sc.surfB = revs[WRAP(j-1, 4)];

                    hSCurve hcb = curve.AddAndAssignId(&sc);

                    STrimBy stb;
                    stb = STrimBy::EntireCurve(this, hcb, /*backwards=*/true);
                    (surface.FindById(sc.surfA))->trim.Add(&stb);
                    stb = STrimBy::EntireCurve(this, hcb, /*backwards=*/false);
                    (surface.FindById(sc.surfB))->trim.Add(&stb);
                }

                // And if this input curve and the one after it both generated
                // surfaces, then trim both of those by the appropriate
                // circle.
                if(revs[j].v && revsp[j].v) {
                    SSurface *ss = surface.FindById(revs[j]);

                    sc = {};
                    sc.isExact = true;
                    sc.exact = SBezier::From(ss->ctrl[0][0],
                                             ss->ctrl[0][1],
                                             ss->ctrl[0][2]);
                    sc.exact.weight[1] = ss->weight[0][1];
                    (sc.exact).MakePwlInto(&(sc.pts));
                    sc.surfA = revs[j];
                    sc.surfB = revsp[j];

                    hSCurve hcc = curve.AddAndAssignId(&sc);

                    STrimBy stb;
                    stb = STrimBy::EntireCurve(this, hcc, /*backwards=*/false);
                    (surface.FindById(sc.surfA))->trim.Add(&stb);
                    stb = STrimBy::EntireCurve(this, hcc, /*backwards=*/true);
                    (surface.FindById(sc.surfB))->trim.Add(&stb);
                }
            }
        }

        hsl.Clear();
    }

    MakeFirstOrderRevolvedSurfaces(pt, axis, i0);
}

void SShell::MakeFirstOrderRevolvedSurfaces(Vector pt, Vector axis, int i0) {
    int i;

    for(i = i0; i < surface.n; i++) {
        SSurface *srf = &(surface[i]);

        // Revolution of a line; this is potentially a plane, which we can
        // rewrite to have degree (1, 1).
        if(srf->degm == 1 && srf->degn == 2) {
            // close start, far start, far finish
            Vector cs, fs, ff;
            double d0, d1;
            d0 = (srf->ctrl[0][0]).DistanceToLine(pt, axis);
            d1 = (srf->ctrl[1][0]).DistanceToLine(pt, axis);

            if(d0 > d1) {
                cs = srf->ctrl[1][0];
                fs = srf->ctrl[0][0];
                ff = srf->ctrl[0][2];
            } else {
                cs = srf->ctrl[0][0];
                fs = srf->ctrl[1][0];
                ff = srf->ctrl[1][2];
            }

            // origin close, origin far
            Vector oc = cs.ClosestPointOnLine(pt, axis),
                   of = fs.ClosestPointOnLine(pt, axis);

            if(oc.Equals(of)) {
                // This is a plane, not a (non-degenerate) cone.
                Vector oldn = srf->NormalAt(0.5, 0.5);

                Vector u = fs.Minus(of), v;

                v = (axis.Cross(u)).WithMagnitude(1);

                double vm = (ff.Minus(of)).Dot(v);
                v = v.ScaledBy(vm);

                srf->degm = 1;
                srf->degn = 1;
                srf->ctrl[0][0] = of;
                srf->ctrl[0][1] = of.Plus(u);
                srf->ctrl[1][0] = of.Plus(v);
                srf->ctrl[1][1] = of.Plus(u).Plus(v);
                srf->weight[0][0] = 1;
                srf->weight[0][1] = 1;
                srf->weight[1][0] = 1;
                srf->weight[1][1] = 1;

                if(oldn.Dot(srf->NormalAt(0.5, 0.5)) < 0) {
                    swap(srf->ctrl[0][0], srf->ctrl[1][0]);
                    swap(srf->ctrl[0][1], srf->ctrl[1][1]);
                }
                continue;
            }

            if(fabs(d0 - d1) < LENGTH_EPS) {
                // This is a cylinder; so transpose it so that we'll recognize
                // it as a surface of extrusion.
                SSurface sn = *srf;

                // Transposing u and v flips the normal, so reverse u to
                // flip it again and put it back where we started.
                sn.degm = 2;
                sn.degn = 1;
                int dm, dn;
                for(dm = 0; dm <= 1; dm++) {
                    for(dn = 0; dn <= 2; dn++) {
                        sn.ctrl  [dn][dm] = srf->ctrl  [1-dm][dn];
                        sn.weight[dn][dm] = srf->weight[1-dm][dn];
                    }
                }

                *srf = sn;
                continue;
            }
        }
    }
}

void SShell::MakeFromCopyOf(SShell *a) {
    ssassert(this != a, "Can't make from copy of self");
    MakeFromTransformationOf(a, {0, 0, 0}, Quaternion::IDENTITY, 1.0);
}

void SShell::MakeFromTransformationOf(SShell *a,
                                      Vector t, Quaternion q, double scale)
{
    booleanFailed = false;
    surface.ReserveMore(a->surface.n);
    for(SSurface &s : a->surface) {
        SSurface n;
        n = SSurface::FromTransformationOf(&s, t, q, scale, /*includingTrims=*/true);
        surface.Add(&n); // keeping the old ID
    }

    curve.ReserveMore(a->curve.n);
    for(SCurve &c : a->curve) {
        SCurve n;
        n = SCurve::FromTransformationOf(&c, t, q, scale);
        curve.Add(&n); // keeping the old ID
    }
}

void SShell::MakeEdgesInto(SEdgeList *sel) {
    for(SSurface &s : surface) {
        s.MakeEdgesInto(this, sel, SSurface::MakeAs::XYZ);
    }
}

void SShell::MakeSectionEdgesInto(Vector n, double d, SEdgeList *sel, SBezierList *sbl)
{
    for(SSurface &s : surface) {
        if(s.CoincidentWithPlane(n, d)) {
            s.MakeSectionEdgesInto(this, sel, sbl);
        }
    }
}

void SShell::TriangulateInto(SMesh *sm) {
#pragma omp parallel for
    for(int i=0; i<surface.n; i++) {
        SSurface *s = &surface[i];
        SMesh m;
        s->TriangulateInto(this, &m);
        #pragma omp critical
        sm->MakeFromCopyOf(&m);
        m.Clear();
    }
}

bool SShell::IsEmpty() const {
    return surface.IsEmpty();
}

void SShell::Clear() {
    for(SSurface &s : surface) {
        s.Clear();
    }
    surface.Clear();

    for(SCurve &c : curve) {
        c.Clear();
    }
    curve.Clear();
}



void SShell::CopySurfacesTrimAgainst(SShell *sha, SShell *shb, SShell *into, SSurface::CombineAs type) {
    std::vector <SSurface> ssn(surface.n);
#pragma omp parallel for
    for (int i = 0; i < surface.n; i++)
    {
        SSurface *ss = &surface[i];
        ssn[i] = ss->MakeCopyTrimAgainst(this, sha, shb, into, type, i);
    }

    for (int i = 0; i < surface.n; i++)
    {
        surface[i].newH = into->surface.AddAndAssignId(&ssn[i]);
    }
    I += surface.n;
}

void SShell::MakeIntersectionCurvesAgainst(SShell *agnst, SShell *into) {
#pragma omp parallel for
    for(int i = 0; i< surface.n; i++) {
        SSurface *sa = &surface[i];

        for(SSurface &sb : agnst->surface){
            // Intersect every surface from our shell against every surface
            // from agnst; this will add zero or more curves to the curve
            // list for into.
            sa->IntersectAgainst(&sb, this, agnst, into);
        }
    }
}

void SShell::CleanupAfterBoolean() {
    for(SSurface &ss : surface) {
        ss.edges.Clear();
    }
}

//-----------------------------------------------------------------------------
// All curves contain handles to the two surfaces that they trim. After a
// Boolean or assembly, we must rewrite those handles to refer to the curves
// by their new IDs.
//-----------------------------------------------------------------------------
void SShell::RewriteSurfaceHandlesForCurves(SShell *a, SShell *b) {
    for(SCurve &sc : curve) {
        sc.surfA = sc.GetSurfaceA(a, b)->newH,
        sc.surfB = sc.GetSurfaceB(a, b)->newH;
    }
}

//-----------------------------------------------------------------------------
// Copy all the surfaces and curves from two different shells into a single
// shell. The only difficulty is to rewrite all of their handles; we don't
// look for any surface intersections, so if two objects interfere then the
// result is just self-intersecting. This is used for assembly, since it's
// much faster than merging as union.
//-----------------------------------------------------------------------------
void SShell::MakeFromAssemblyOf(SShell *a, SShell *b) {
    booleanFailed = false;

    Vector t = {};
    Quaternion q = Quaternion::IDENTITY;
    int i = 0;
    SShell *ab;

    // First, copy over all the curves. Note which shell (a or b) each curve
    // came from, but assign it a new ID.
    curve.ReserveMore(a->curve.n + b->curve.n);
    SCurve cn;
    for(i = 0; i < 2; i++) {
        ab = (i == 0) ? a : b;
        for(SCurve &c : ab->curve) {
            cn = SCurve::FromTransformationOf(&c, t, q, 1.0);
            cn.source = (i == 0) ? SCurve::Source::A : SCurve::Source::B;
            // surfA and surfB are wrong now, and we can't fix them until
            // we've assigned IDs to the surfaces. So we'll get that later.
            c.newH = curve.AddAndAssignId(&cn);
        }
    }

    // Likewise copy over all the surfaces.
    surface.ReserveMore(a->surface.n + b->surface.n);
    SSurface sn;
    for(i = 0; i < 2; i++) {
        ab = (i == 0) ? a : b;
        for(SSurface &s : ab->surface) {
            sn = SSurface::FromTransformationOf(&s, t, q, 1.0, /*includingTrims=*/true);
            // All the trim curve IDs get rewritten; we know the new handles
            // to the curves since we recorded them in the previous step.
            STrimBy *stb;
            for(stb = sn.trim.First(); stb; stb = sn.trim.NextAfter(stb)) {
                stb->curve = ab->curve.FindById(stb->curve)->newH;
            }
            s.newH = surface.AddAndAssignId(&sn);
        }
    }

    // Finally, rewrite the surfaces associated with each curve to use the
    // new handles.
    RewriteSurfaceHandlesForCurves(a, b);
}

void SShell::MakeFromBoolean(SShell *a, SShell *b, SSurface::CombineAs type) {
    booleanFailed = false;

    a->MakeClassifyingBsps(NULL);
    b->MakeClassifyingBsps(NULL);

    // Copy over all the original curves, splitting them so that a
    // piecewise linear segment never crosses a surface from the other
    // shell.
    a->CopyCurvesSplitAgainst(/*opA=*/true,  b, this);
    b->CopyCurvesSplitAgainst(/*opA=*/false, a, this);

    // Generate the intersection curves for each surface in A against all
    // the surfaces in B (which is all of the intersection curves).
    a->MakeIntersectionCurvesAgainst(b, this);

    for(SCurve &sc : curve) {
        SSurface *srfA = sc.GetSurfaceA(a, b),
                 *srfB = sc.GetSurfaceB(a, b);

        sc.RemoveShortSegments(srfA, srfB);
    }

    // And clean up the piecewise linear things we made as a calculation aid
    a->CleanupAfterBoolean();
    b->CleanupAfterBoolean();
    // Remake the classifying BSPs with the split (and short-segment-removed)
    // curves
    a->MakeClassifyingBsps(this);
    b->MakeClassifyingBsps(this);

    if(b->surface.IsEmpty() || a->surface.IsEmpty()) {
        I = 1000000;
    } else {
        I = 0;
    }
    // Then trim and copy the surfaces
    a->CopySurfacesTrimAgainst(a, b, this, type);
    b->CopySurfacesTrimAgainst(a, b, this, type);

    // Now that we've copied the surfaces, we know their new hSurfaces, so
    // rewrite the curves to refer to the surfaces by their handles in the
    // result.
    RewriteSurfaceHandlesForCurves(a, b);

    // And clean up the piecewise linear things we made as a calculation aid
    a->CleanupAfterBoolean();
    b->CleanupAfterBoolean();
}

//-----------------------------------------------------------------------------
// All of the BSP routines that we use to perform and accelerate polygon ops.
//-----------------------------------------------------------------------------
void SShell::MakeClassifyingBsps(SShell *useCurvesFrom) {
#pragma omp parallel for
    for(int i = 0; i<surface.n; i++) {
        surface[i].MakeClassifyingBsp(this, useCurvesFrom);
    }
}

} // namespace SolveSpace
