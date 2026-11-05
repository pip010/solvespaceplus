
#include "ssSurface.h"

namespace SolveSpace
{
        
    double SSurface::ChordToleranceForEdge(Vector a, Vector b) const {
        Vector as = PointAt(a.x, a.y), bs = PointAt(b.x, b.y);

        double worst = VERY_NEGATIVE;
        int i;
        for(i = 1; i <= 3; i++) {
            Vector p  = a. Plus((b. Minus(a )).ScaledBy(i/4.0)),
                ps = as.Plus((bs.Minus(as)).ScaledBy(i/4.0));

            Vector pps = PointAt(p.x, p.y);
            worst = max(worst, (pps.Minus(ps)).MagSquared());
        }
        return sqrt(worst);
    }

    Vector SSurface::PointAtMaybeSwapped(double u, double v, bool swapped) const {
        if(swapped) {
            return PointAt(v, u);
        } else {
            return PointAt(u, v);
        }
    }

    Vector SSurface::NormalAtMaybeSwapped(double u, double v, bool swapped) const {
        Vector du, dv;
        if(swapped) {
            TangentsAt(v, u, &dv, &du);
        } else {
            TangentsAt(u, v, &du, &dv);
        }
        return du.Cross(dv).WithMagnitude(1.0);
    }

    void SSurface::MakeTriangulationGridInto(List<double> *l, double vs, double vf,
                                            bool swapped, int depth) const
    {
        double worst = 0;

        // Try piecewise linearizing four curves, at u = 0, 1/3, 2/3, 1; choose
        // the worst chord tolerance of any of those.
        double worst_twist = 1.0;
        int i;
        for(i = 0; i <= 3; i++) {
            double u = i/3.0;

            // This chord test should be identical to the one in SBezier::MakePwl
            // to make the piecewise linear edges line up with the grid more or
            // less.
            Vector ps = PointAtMaybeSwapped(u, vs, swapped),
                pf = PointAtMaybeSwapped(u, vf, swapped);

            double vm1 = (2*vs + vf) / 3,
                vm2 = (vs + 2*vf) / 3;

            Vector pm1 = PointAtMaybeSwapped(u, vm1, swapped),
                pm2 = PointAtMaybeSwapped(u, vm2, swapped);

            // 0.999 is about 2.5 degrees of twist over the middle 1/3 V-span.
            // we don't check at the ends because the derivative may not be valid there.
            double twist = 1.0;
            if (degm == 1) twist = NormalAtMaybeSwapped(u, vm1, swapped).Dot(
                                NormalAtMaybeSwapped(u, vm2, swapped) );
            if (twist < worst_twist) worst_twist = twist;

            worst = max(worst, pm1.DistanceToLine(ps, pf.Minus(ps)));
            worst = max(worst, pm2.DistanceToLine(ps, pf.Minus(ps)));
        }

        double step = 1.0/SS.GetMaxSegments();
        if( ((vf - vs) < step || worst < SS.ChordTolMm())
            && ((worst_twist > 0.999) || (depth > 3)) ) {
            l->Add(&vf);
        } else {
            MakeTriangulationGridInto(l, vs, (vs+vf)/2, swapped, depth+1);
            MakeTriangulationGridInto(l, (vs+vf)/2, vf, swapped, depth+1);
        }
    }

    SSurface SSurface::FromExtrusionOf(SBezier *sb, Vector t0, Vector t1) {
        SSurface ret = {};

        ret.degm = sb->deg;
        ret.degn = 1;

        int i;
        for(i = 0; i <= ret.degm; i++) {
            ret.ctrl[i][0] = (sb->ctrl[i]).Plus(t0);
            ret.weight[i][0] = sb->weight[i];

            ret.ctrl[i][1] = (sb->ctrl[i]).Plus(t1);
            ret.weight[i][1] = sb->weight[i];
        }

        return ret;
    }

    bool SSurface::IsExtrusion(SBezier *of, Vector *alongp) const {
        int i;

        if(degn != 1) return false;

        Vector along = (ctrl[0][1]).Minus(ctrl[0][0]);
        for(i = 0; i <= degm; i++) {
            if((fabs(weight[i][1] - weight[i][0]) < LENGTH_EPS) &&
            ((ctrl[i][1]).Minus(ctrl[i][0])).Equals(along))
            {
                continue;
            }
            return false;
        }

        // yes, we are a surface of extrusion; copy the original curve and return
        if(of) {
            for(i = 0; i <= degm; i++) {
                of->weight[i] = weight[i][0];
                of->ctrl[i] = ctrl[i][0];
            }
            of->deg = degm;
            *alongp = along;
        }
        return true;
    }

    bool SSurface::IsCylinder(Vector *axis, Vector *center, double *r,
                                Vector *start, Vector *finish) const
    {
        SBezier sb;
        if(!IsExtrusion(&sb, axis)) return false;
        if(!sb.IsCircle(*axis, center, r)) return false;

        *start = sb.ctrl[0];
        *finish = sb.ctrl[2];
        return true;
    }

    // Create a surface patch by revolving and possibly translating a curve.
    // Works for sections up to but not including 180 degrees.
    SSurface SSurface::FromRevolutionOf(SBezier *sb, Vector pt, Vector axis, double thetas,
                                        double thetaf, double dists,
                                        double distf) { // s is start, f is finish
        SSurface ret = {};
        ret.degm = sb->deg;
        ret.degn = 2;

        double dtheta = fabs(WRAP_SYMMETRIC(thetaf - thetas, 2*PI));
        double w      = cos(dtheta / 2);

        // Revolve the curve about the z axis
        int i;
        for(i = 0; i <= ret.degm; i++) {
            Vector p = sb->ctrl[i];

            Vector ps = p.RotatedAbout(pt, axis, thetas),
                pf = p.RotatedAbout(pt, axis, thetaf);

            // The middle control point should be at the intersection of the tangents at ps and pf.
            // This is equivalent but works for 0 <= angle < 180 degrees.
            Vector mid = ps.Plus(pf).ScaledBy(0.5);
            Vector c   = ps.ClosestPointOnLine(pt, axis);
            Vector ct  = mid.Minus(c).ScaledBy(1 / (w * w)).Plus(c);

            // not sure this is needed
            if(ps.Equals(pf)) {
                ps = c;
                ct = c;
                pf = c;
            }
            // moving along the axis can create hilical surfaces (or straight extrusion if
            // thetas==thetaf)
            ret.ctrl[i][0] = ps.Plus(axis.ScaledBy(dists));
            ret.ctrl[i][1] = ct.Plus(axis.ScaledBy((dists + distf) / 2));
            ret.ctrl[i][2] = pf.Plus(axis.ScaledBy(distf));

            ret.weight[i][0] = sb->weight[i];
            ret.weight[i][1] = sb->weight[i] * w;
            ret.weight[i][2] = sb->weight[i];
        }

        return ret;
    }

    SSurface SSurface::FromPlane(Vector pt, Vector u, Vector v) {
        SSurface ret = {};

        ret.degm = 1;
        ret.degn = 1;

        ret.weight[0][0] = ret.weight[0][1] = 1;
        ret.weight[1][0] = ret.weight[1][1] = 1;

        ret.ctrl[0][0] = pt;
        ret.ctrl[0][1] = pt.Plus(u);
        ret.ctrl[1][0] = pt.Plus(v);
        ret.ctrl[1][1] = pt.Plus(v).Plus(u);

        return ret;
    }

    SSurface SSurface::FromTransformationOf(SSurface *a, Vector t, Quaternion q, double scale,
                                            bool includingTrims)
    {
        bool needRotate    = !EXACT(q.vx == 0.0 && q.vy == 0.0 && q.vz == 0.0 && q.w == 1.0);
        bool needTranslate = !EXACT(t.x  == 0.0 && t.y  == 0.0 && t.z  == 0.0);
        bool needScale     = !EXACT(scale == 1.0);

        SSurface ret = {};
        ret.h = a->h;
        ret.color = a->color;
        ret.face = a->face;

        ret.degm = a->degm;
        ret.degn = a->degn;
        int i, j;
        for(i = 0; i <= 3; i++) {
            for(j = 0; j <= 3; j++) {
                Vector ctrl = a->ctrl[i][j];
                if(needScale) {
                    ctrl = ctrl.ScaledBy(scale);
                }
                if(needRotate) {
                    ctrl = q.Rotate(ctrl);
                }
                if(needTranslate) {
                    ctrl = ctrl.Plus(t);
                }
                ret.ctrl[i][j] = ctrl;
                ret.weight[i][j] = a->weight[i][j];
            }
        }

        if(includingTrims) {
            STrimBy *stb;
            ret.trim.ReserveMore(a->trim.n);
            for(stb = a->trim.First(); stb; stb = a->trim.NextAfter(stb)) {
                STrimBy n = *stb;
                if(needScale) {
                    n.start  = n.start.ScaledBy(scale);
                    n.finish = n.finish.ScaledBy(scale);
                }
                if(needRotate) {
                    n.start  = q.Rotate(n.start);
                    n.finish = q.Rotate(n.finish);
                }
                if(needTranslate) {
                    n.start  = n.start.Plus(t);
                    n.finish = n.finish.Plus(t);
                }
                ret.trim.Add(&n);
            }
        }

        if(scale < 0) {
            // If we mirror every surface of a shell, then it will end up inside
            // out. So fix that here.
            ret.Reverse();
        }

        return ret;
    }

    void SSurface::GetAxisAlignedBounding(Vector *ptMax, Vector *ptMin) const {
        *ptMax = {VERY_NEGATIVE, VERY_NEGATIVE, VERY_NEGATIVE};
        *ptMin = {VERY_POSITIVE, VERY_POSITIVE, VERY_POSITIVE};

        int i, j;
        for(i = 0; i <= degm; i++) {
            for(j = 0; j <= degn; j++) {
                (ctrl[i][j]).MakeMaxMin(ptMax, ptMin);
            }
        }
    }

    bool SSurface::LineEntirelyOutsideBbox(Vector a, Vector b, bool asSegment) const {
        Vector amax, amin;
        GetAxisAlignedBounding(&amax, &amin);
        if(!Vector::BoundingBoxIntersectsLine(amax, amin, a, b, asSegment)) {
            // The line segment could fail to intersect the bbox, but lie entirely
            // within it and intersect the surface.
            if(a.OutsideAndNotOn(amax, amin) && b.OutsideAndNotOn(amax, amin)) {
                return true;
            }
        }
        return false;
    }

    //-----------------------------------------------------------------------------
    // Generate the piecewise linear approximation of the trim stb, which applies
    // to the curve sc.
    //-----------------------------------------------------------------------------
    void SSurface::MakeTrimEdgesInto(SEdgeList *sel, MakeAs flags,
                                    SCurve *sc, STrimBy *stb)
    {
        Vector prev = {};
        bool inCurve = false, empty = true;
        double u = 0, v = 0;

        int i, first, last, increment;
        if(stb->backwards) {
            first = sc->pts.n - 1;
            last = 0;
            increment = -1;
        } else {
            first = 0;
            last = sc->pts.n - 1;
            increment = 1;
        }
        for(i = first; i != (last + increment); i += increment) {
            Vector tpt, *pt = &(sc->pts[i].p);

            if(flags == MakeAs::UV) {
                ClosestPointTo(*pt, &u, &v);
                tpt = {u, v, 0};
            } else {
                tpt = *pt;
            }

            if(inCurve) {
                sel->AddEdge(prev, tpt, sc->h.v, stb->backwards);
                empty = false;
            }

            prev = tpt;     // either uv or xyz, depending on flags

            if(pt->Equals(stb->start)) inCurve = true;
            if(pt->Equals(stb->finish)) inCurve = false;
        }
        if(inCurve) dbp("trim was unterminated");
        if(empty)   dbp("trim was empty");
    }

    //-----------------------------------------------------------------------------
    // Generate all of our trim curves, in piecewise linear form. We can do
    // so in either uv or xyz coordinates. And if requested, then we can use
    // the split curves from useCurvesFrom instead of the curves in our own
    // shell.
    //-----------------------------------------------------------------------------
    void SSurface::MakeEdgesInto(SShell *shell, SEdgeList *sel, MakeAs flags,
                                SShell *useCurvesFrom)
    {
        STrimBy *stb;
        for(stb = trim.First(); stb; stb = trim.NextAfter(stb)) {
            SCurve *sc = shell->curve.FindById(stb->curve);

            // We have the option to use the curves from another shell; this
            // is relevant when generating the coincident edges while doing the
            // Booleans, since the curves from the output shell will be split
            // against any intersecting surfaces (and the originals aren't).
            if(useCurvesFrom) {
                sc = useCurvesFrom->curve.FindById(sc->newH);
            }

            MakeTrimEdgesInto(sel, flags, sc, stb);
        }
    }

    //-----------------------------------------------------------------------------
    // Compute the exact tangent to the intersection curve between two surfaces,
    // by taking the cross product of the surface normals. We choose the direction
    // of this tangent so that its dot product with dir is positive.
    //-----------------------------------------------------------------------------
    Vector SSurface::ExactSurfaceTangentAt(Vector p, SSurface *srfA, SSurface *srfB, Vector dir)
    {
        Point2d puva, puvb;
        srfA->ClosestPointTo(p, &puva);
        srfB->ClosestPointTo(p, &puvb);
        Vector ts = (srfA->NormalAt(puva)).Cross(
                    (srfB->NormalAt(puvb)));
        ts = ts.WithMagnitude(1);
        if(ts.Dot(dir) < 0) {
            ts = ts.ScaledBy(-1);
        }
        return ts;
    }

    //-----------------------------------------------------------------------------
    // Report our trim curves. If a trim curve is exact and sbl is not null, then
    // add its exact form to sbl. Otherwise, add its piecewise linearization to
    // sel.
    //-----------------------------------------------------------------------------
    void SSurface::MakeSectionEdgesInto(SShell *shell, SEdgeList *sel, SBezierList *sbl)
    {
        STrimBy *stb;
        for(stb = trim.First(); stb; stb = trim.NextAfter(stb)) {
            SCurve *sc = shell->curve.FindById(stb->curve);
            SBezier *sb = &(sc->exact);

            if(sbl && sc->isExact && (sb->deg != 1 || !sel)) {
                double ts, tf;
                if(stb->backwards) {
                    sb->ClosestPointTo(stb->start,  &tf);
                    sb->ClosestPointTo(stb->finish, &ts);
                } else {
                    sb->ClosestPointTo(stb->start,  &ts);
                    sb->ClosestPointTo(stb->finish, &tf);
                }
                SBezier junk_bef, keep_aft;
                sb->SplitAt(ts, &junk_bef, &keep_aft);
                // In the kept piece, the range that used to go from ts to 1
                // now goes from 0 to 1; so rescale tf appropriately.
                tf = (tf - ts)/(1 - ts);

                SBezier keep_bef, junk_aft;
                keep_aft.SplitAt(tf, &keep_bef, &junk_aft);

                sbl->l.Add(&keep_bef);
            } else if(sbl && !sel && !sc->isExact) {
                // We must approximate this trim curve, as piecewise cubic sections.
                SSurface *srfA = shell->surface.FindById(sc->surfA);
                SSurface *srfB = shell->surface.FindById(sc->surfB);

                Vector s = stb->backwards ? stb->finish : stb->start,
                    f = stb->backwards ? stb->start : stb->finish;

                int sp, fp;
                for(sp = 0; sp < sc->pts.n; sp++) {
                    if(s.Equals(sc->pts[sp].p)) break;
                }
                if(sp >= sc->pts.n) return;
                for(fp = sp; fp < sc->pts.n; fp++) {
                    if(f.Equals(sc->pts[fp].p)) break;
                }
                if(fp >= sc->pts.n) return;
                // So now the curve we want goes from elem[sp] to elem[fp]

                while(sp < fp) {
                    // Initially, we'll try approximating the entire trim curve
                    // as a single Bezier segment
                    int fpt = fp;

                    for(;;) {
                        // So construct a cubic Bezier with the correct endpoints
                        // and tangents for the current span.
                        Vector st = sc->pts[sp].p,
                            ft = sc->pts[fpt].p,
                            sf = ft.Minus(st);
                        double m = sf.Magnitude() / 3;

                        Vector stan = ExactSurfaceTangentAt(st, srfA, srfB, sf),
                            ftan = ExactSurfaceTangentAt(ft, srfA, srfB, sf);

                        SBezier sb = SBezier::From(st,
                                                st.Plus (stan.WithMagnitude(m)),
                                                ft.Minus(ftan.WithMagnitude(m)),
                                                ft);

                        // And test how much this curve deviates from the
                        // intermediate points (if any).
                        int i;
                        bool tooFar = false;
                        for(i = sp + 1; i <= (fpt - 1); i++) {
                            Vector p = sc->pts[i].p;
                            double t;
                            sb.ClosestPointTo(p, &t, /*mustConverge=*/false);
                            Vector pp = sb.PointAt(t);
                            if((pp.Minus(p)).Magnitude() > SS.ChordTolMm()/2) {
                                tooFar = true;
                                break;
                            }
                        }

                        if(tooFar) {
                            // Deviates by too much, so try a shorter span
                            fpt--;
                            continue;
                        } else {
                            // Okay, so use this piece and break.
                            sbl->l.Add(&sb);
                            break;
                        }
                    }

                    // And continue interpolating, starting wherever the curve
                    // we just generated finishes.
                    sp = fpt;
                }
            } else {
                if(sel) MakeTrimEdgesInto(sel, MakeAs::XYZ, sc, stb);
            }
        }
    }

    void SSurface::TriangulateInto(SShell *shell, SMesh *sm) {
        SEdgeList el = {};

        MakeEdgesInto(shell, &el, MakeAs::UV);

        SPolygon poly = {};
        if(el.AssemblePolygon(&poly, NULL, /*keepDir=*/true)) {
            int i, start = sm->l.n;
            if(degm == 1 && degn == 1) {
                // A surface with curvature along one direction only; so
                // choose the triangulation with chords that lie as much
                // as possible within the surface. And since the trim curves
                // have been pwl'd to within the desired chord tol, that will
                // produce a surface good to within roughly that tol.
                //
                // If this is just a plane (degree (1, 1)) then the triangulation
                // code will notice that, and not bother checking chord tols.
                poly.UvTriangulateInto(sm, this);
            } else {
                // A surface with compound curvature. So we must overlay a
                // two-dimensional grid, and triangulate around that.
                poly.UvGridTriangulateInto(sm, this);
            }

            STriMeta meta = { face, color };
            for(i = start; i < sm->l.n; i++) {
                STriangle *st = &(sm->l[i]);
                st->meta = meta;
                st->an = NormalAt(st->a.x, st->a.y);
                st->bn = NormalAt(st->b.x, st->b.y);
                st->cn = NormalAt(st->c.x, st->c.y);
                st->a = PointAt(st->a.x, st->a.y);
                st->b = PointAt(st->b.x, st->b.y);
                st->c = PointAt(st->c.x, st->c.y);
                // Works out that my chosen contour direction is inconsistent with
                // the triangle direction, sigh.
                st->FlipNormal();
            }
        } else {
            dbp("failed to assemble polygon to trim nurbs surface in uv space");
        }

        el.Clear();
        poly.Clear();
    }

    //-----------------------------------------------------------------------------
    // Reverse the parametrisation of one of our dimensions, which flips the
    // normal. We therefore must reverse all our trim curves too. The uv
    // coordinates change, but trim curves are stored as xyz so nothing happens
    //-----------------------------------------------------------------------------
    void SSurface::Reverse() {
        int i, j;
        for(i = 0; i < (degm+1)/2; i++) {
            for(j = 0; j <= degn; j++) {
                swap(ctrl[i][j], ctrl[degm-i][j]);
                swap(weight[i][j], weight[degm-i][j]);
            }
        }

        STrimBy *stb;
        for(stb = trim.First(); stb; stb = trim.NextAfter(stb)) {
            stb->backwards = !stb->backwards;
            swap(stb->start, stb->finish);
        }
    }

    void SSurface::ScaleSelfBy(double s) {
        int i, j;
        for(i = 0; i <= degm; i++) {
            for(j = 0; j <= degn; j++) {
                ctrl[i][j] = ctrl[i][j].ScaledBy(s);
            }
        }
    }

    void SSurface::Clear() {
        trim.Clear();
    }

    //------------------------------ RAY-CAST --------------------

    double SSurface::DepartureFromCoplanar() const {
        int i, j;
        int ia, ja, ib = 0, jb = 0, ic = 0, jc = 0;
        double best;

        // Grab three points to define a plane; first choose (0, 0) arbitrarily.
        ia = ja = 0;
        // Then the point farthest from pt a.
        best = VERY_NEGATIVE;
        for(i = 0; i <= degm; i++) {
            for(j = 0; j <= degn; j++) {
                if(i == ia && j == ja) continue;

                double dist = (ctrl[i][j]).Minus(ctrl[ia][ja]).Magnitude();
                if(dist > best) {
                    best = dist;
                    ib = i;
                    jb = j;
                }
            }
        }
        // Then biggest magnitude of ab cross ac.
        best = VERY_NEGATIVE;
        for(i = 0; i <= degm; i++) {
            for(j = 0; j <= degn; j++) {
                if(i == ia && j == ja) continue;
                if(i == ib && j == jb) continue;

                double mag =
                    ((ctrl[ia][ja].Minus(ctrl[ib][jb]))).Cross(
                    (ctrl[ia][ja].Minus(ctrl[i ][j ]))).Magnitude();
                if(mag > best) {
                    best = mag;
                    ic = i;
                    jc = j;
                }
            }
        }

        Vector n = ((ctrl[ia][ja].Minus(ctrl[ib][jb]))).Cross(
                    (ctrl[ia][ja].Minus(ctrl[ic][jc])));
        n = n.WithMagnitude(1);
        double d = (ctrl[ia][ja]).Dot(n);

        // Finally, calculate the deviation from each point to the plane.
        double farthest = VERY_NEGATIVE;
        for(i = 0; i <= degm; i++) {
            for(j = 0; j <= degn; j++) {
                double dist = fabs(n.Dot(ctrl[i][j]) - d);
                if(dist > farthest) {
                    farthest = dist;
                }
            }
        }
        return farthest;
    }

    void SSurface::WeightControlPoints() {
        int i, j;
        for(i = 0; i <= degm; i++) {
            for(j = 0; j <= degn; j++) {
                ctrl[i][j] = (ctrl[i][j]).ScaledBy(weight[i][j]);
            }
        }
    }
    void SSurface::UnWeightControlPoints() {
        int i, j;
        for(i = 0; i <= degm; i++) {
            for(j = 0; j <= degn; j++) {
                ctrl[i][j] = (ctrl[i][j]).ScaledBy(1.0/weight[i][j]);
            }
        }
    }
    void SSurface::CopyRowOrCol(bool row, int this_ij, SSurface *src, int src_ij) {
        if(row) {
            int j;
            for(j = 0; j <= degn; j++) {
                ctrl  [this_ij][j] = src->ctrl  [src_ij][j];
                weight[this_ij][j] = src->weight[src_ij][j];
            }
        } else {
            int i;
            for(i = 0; i <= degm; i++) {
                ctrl  [i][this_ij] = src->ctrl  [i][src_ij];
                weight[i][this_ij] = src->weight[i][src_ij];
            }
        }
    }
    void SSurface::BlendRowOrCol(bool row, int this_ij, SSurface *a, int a_ij,
                                                        SSurface *b, int b_ij)
    {
        if(row) {
            int j;
            for(j = 0; j <= degn; j++) {
                Vector c = (a->ctrl  [a_ij][j]).Plus(b->ctrl  [b_ij][j]);
                double w = (a->weight[a_ij][j]   +   b->weight[b_ij][j]);
                ctrl  [this_ij][j] = c.ScaledBy(0.5);
                weight[this_ij][j] = w / 2;
            }
        } else {
            int i;
            for(i = 0; i <= degm; i++) {
                Vector c = (a->ctrl  [i][a_ij]).Plus(b->ctrl  [i][b_ij]);
                double w = (a->weight[i][a_ij]   +   b->weight[i][b_ij]);
                ctrl  [i][this_ij] = c.ScaledBy(0.5);
                weight[i][this_ij] = w / 2;
            }
        }
    }
    void SSurface::SplitInHalf(bool byU, SSurface *sa, SSurface *sb) {
        sa->degm = sb->degm = degm;
        sa->degn = sb->degn = degn;

        // by de Casteljau's algorithm in a projective space; so we must work
        // on points (w*x, w*y, w*z, w) so create a temporary copy
        SSurface st;
        st = *this;
        st.WeightControlPoints();

        switch(byU ? degm : degn) {
            case 1:
                sa->CopyRowOrCol (byU, 0, &st, 0);
                sb->CopyRowOrCol (byU, 1, &st, 1);

                sa->BlendRowOrCol(byU, 1, &st, 0, &st, 1);
                sb->BlendRowOrCol(byU, 0, &st, 0, &st, 1);
                break;

            case 2:
                sa->CopyRowOrCol (byU, 0, &st, 0);
                sb->CopyRowOrCol (byU, 2, &st, 2);

                sa->BlendRowOrCol(byU, 1, &st, 0, &st, 1);
                sb->BlendRowOrCol(byU, 1, &st, 1, &st, 2);

                sa->BlendRowOrCol(byU, 2, sa,   1, sb,   1);
                sb->BlendRowOrCol(byU, 0, sa,   1, sb,   1);
                break;

            case 3: {
                sa->CopyRowOrCol (byU, 0, &st, 0);
                sb->CopyRowOrCol (byU, 3, &st, 3);

                sa->BlendRowOrCol(byU, 1, &st, 0, &st, 1);
                sb->BlendRowOrCol(byU, 2, &st, 2, &st, 3);
                st. BlendRowOrCol(byU, 0, &st, 1, &st, 2); // use row/col 0 as scratch

                sa->BlendRowOrCol(byU, 2, sa,   1, &st,  0);
                sb->BlendRowOrCol(byU, 1, sb,   2, &st,  0);

                sa->BlendRowOrCol(byU, 3, sa,   2, sb,   1);
                sb->BlendRowOrCol(byU, 0, sa,   2, sb,   1);
                break;
            }

            default: ssassert(false, "Unexpected degree of spline");
        }

        sa->UnWeightControlPoints();
        sb->UnWeightControlPoints();
    }

    //-----------------------------------------------------------------------------
    // Find all points where the indicated finite (if segment) or infinite (if not
    // segment) line intersects our surface. Report them in uv space in the list.
    // We first do a bounding box check; if the line doesn't intersect, then we're
    // done. If it does, then we check how small our surface is. If it's big,
    // then we subdivide into quarters and recurse. If it's small, then we refine
    // by Newton's method and record the point.
    //-----------------------------------------------------------------------------
    void SSurface::AllPointsIntersectingUntrimmed(Vector a, Vector b,
                                                int *cnt, int *level,
                                                List<Inter> *l, bool asSegment,
                                                SSurface *sorig)
    {
        // Test if the line intersects our axis-aligned bounding box; if no, then
        // no possibility of an intersection
        if(LineEntirelyOutsideBbox(a, b, asSegment)) return;

        if(*cnt > 2000) {
            dbp("!!! too many subdivisions (level=%d)!", *level);
            dbp("degm = %d degn = %d", degm, degn);
            return;
        }
        (*cnt)++;

        // If we might intersect, and the surface is small, then switch to Newton
        // iterations.
        if(DepartureFromCoplanar() < 0.2*SS.ChordTolMm()) {
            Vector p = (ctrl[0   ][0   ]).Plus(
                        ctrl[0   ][degn]).Plus(
                        ctrl[degm][0   ]).Plus(
                        ctrl[degm][degn]).ScaledBy(0.25);
            Inter inter;
            sorig->ClosestPointTo(p, &(inter.p.x), &(inter.p.y), /*mustConverge=*/false);
            if(sorig->PointIntersectingLine(a, b, &(inter.p.x), &(inter.p.y))) {
                Vector p = sorig->PointAt(inter.p.x, inter.p.y);
                // Debug check, verify that the point lies in both surfaces
                // (which it ought to, since the surfaces should be coincident)
                double u, v;
                ClosestPointTo(p, &u, &v);
                l->Add(&inter);
            } else {
                // Might not converge if line is almost tangent to surface...
                dbp("AllPointsIntersectingUntrimmed failed");
            }
            return;
        }

        // But the surface is big, so split it, alternating by u and v
        SSurface surf0, surf1;
        SplitInHalf((*level & 1) == 0, &surf0, &surf1);

        int nextLevel = (*level) + 1;
        (*level) = nextLevel;
        surf0.AllPointsIntersectingUntrimmed(a, b, cnt, level, l, asSegment, sorig);
        (*level) = nextLevel;
        surf1.AllPointsIntersectingUntrimmed(a, b, cnt, level, l, asSegment, sorig);
    }

    //-----------------------------------------------------------------------------
    // Find all points where a line through a and b intersects our surface, and
    // add them to the list. If seg is true then report only intersections that
    // lie within the finite line segment (not including the endpoints); otherwise
    // we work along the infinite line. And we report either just intersections
    // inside the trim curve, or any intersection with u, v in [0, 1]. And we
    // either disregard or report tangent points.
    //-----------------------------------------------------------------------------
    void SSurface::AllPointsIntersecting(Vector a, Vector b,
                                        List<SInter> *l,
                                        bool asSegment, bool trimmed, bool inclTangent)
    {
        if(LineEntirelyOutsideBbox(a, b, asSegment)) return;

        Vector ba = b.Minus(a);
        double bam = ba.Magnitude();

        List<Inter> inters = {};

        // All the intersections between the line and the surface; either special
        // cases that we can quickly solve in closed form, or general numerical.
        Vector center, axis, start, finish;
        double radius;
        if(degm == 1 && degn == 1) {
            // Against a plane, easy.
            Vector n = NormalAt(0, 0).WithMagnitude(1);
            double d = n.Dot(PointAt(0, 0));
            // Trim to line segment now if requested, don't generate points that
            // would just get discarded later.
            if(!asSegment ||
            (n.Dot(a) > d + LENGTH_EPS && n.Dot(b) < d - LENGTH_EPS) ||
            (n.Dot(b) > d + LENGTH_EPS && n.Dot(a) < d - LENGTH_EPS))
            {
                Vector p = Vector::AtIntersectionOfPlaneAndLine(n, d, a, b, NULL);
                Inter inter;
                ClosestPointTo(p, &(inter.p.x), &(inter.p.y));
                inters.Add(&inter);
            }
        } else if(IsCylinder(&axis, &center, &radius, &start, &finish)) {
            // This one can be solved in closed form too.
            Vector ab = b.Minus(a);
            if(axis.Cross(ab).Magnitude() < LENGTH_EPS) {
                // edge is parallel to axis of cylinder, no intersection points
                return;
            }
            // A coordinate system centered at the center of the circle, with
            // the edge under test horizontal
            Vector u, v, n = axis.WithMagnitude(1);
            u = (ab.Minus(n.ScaledBy(ab.Dot(n)))).WithMagnitude(1);
            v = n.Cross(u);
            Point2d ap = (a.Minus(center)).DotInToCsys(u, v, n).ProjectXy(),
                    bp = (b.Minus(center)).DotInToCsys(u, v, n).ProjectXy(),
                    sp = (start. Minus(center)).DotInToCsys(u, v, n).ProjectXy(),
                    fp = (finish.Minus(center)).DotInToCsys(u, v, n).ProjectXy();

            double thetas = atan2(sp.y, sp.x), thetaf = atan2(fp.y, fp.x);

            Point2d ip[2];
            int ip_n = 0;
            if(fabs(fabs(ap.y) - radius) < LENGTH_EPS) {
                // tangent
                if(inclTangent) {
                    ip[0] = Point2d::From(0, ap.y);
                    ip_n = 1;
                }
            } else if(fabs(ap.y) < radius) {
                // two intersections
                double xint = sqrt(radius*radius - ap.y*ap.y);
                ip[0] = Point2d::From(-xint, ap.y);
                ip[1] = Point2d::From( xint, ap.y);
                ip_n = 2;
            }
            int i;
            for(i = 0; i < ip_n; i++) {
                double t = (ip[i].Minus(ap)).DivProjected(bp.Minus(ap));
                // This is a point on the circle; but is it on the arc?
                Point2d pp = ap.Plus((bp.Minus(ap)).ScaledBy(t));
                double theta = atan2(pp.y, pp.x);
                double dp = WRAP_SYMMETRIC(theta  - thetas, 2*PI),
                    df = WRAP_SYMMETRIC(thetaf - thetas, 2*PI);
                double tol = LENGTH_EPS/radius;

                if((df > 0 && ((dp < -tol) || (dp > df + tol))) ||
                (df < 0 && ((dp >  tol) || (dp < df - tol))))
                {
                    continue;
                }

                Vector p = a.Plus((b.Minus(a)).ScaledBy(t));

                Inter inter;
                ClosestPointTo(p, &(inter.p.x), &(inter.p.y));
                inters.Add(&inter);
            }
        } else {
            // General numerical solution by subdivision, fallback
            int cnt = 0, level = 0;
            AllPointsIntersectingUntrimmed(a, b, &cnt, &level, &inters, asSegment, this);
        }

        // Remove duplicate intersection points
        inters.ClearTags();
        int i, j;
        for(i = 0; i < inters.n; i++) {
            for(j = i + 1; j < inters.n; j++) {
                if(inters[i].p.Equals(inters[j].p)) {
                    inters[j].tag = 1;
                }
            }
        }
        inters.RemoveTagged();

        for(i = 0; i < inters.n; i++) {
            Point2d puv = inters[i].p;

            // Make sure the point lies within the finite line segment
            Vector pxyz = PointAt(puv.x, puv.y);
            double t = (pxyz.Minus(a)).DivProjected(ba);
            if(asSegment && (t > 1 - LENGTH_EPS/bam || t < LENGTH_EPS/bam)) {
                continue;
            }

            // And that it lies inside our trim region
            Point2d dummy = { 0, 0 };
            SBspUv::Class c = (bsp) ? bsp->ClassifyPoint(puv, dummy, this) : SBspUv::Class::OUTSIDE;
            if(trimmed && c == SBspUv::Class::OUTSIDE) {
                continue;
            }

            // It does, so generate the intersection
            SInter si;
            si.p = pxyz;
            si.surfNormal = NormalAt(puv.x, puv.y);
            si.pinter = puv;
            si.srf = this;
            si.onEdge = (c != SBspUv::Class::INSIDE);
            l->Add(&si);
        }

        inters.Clear();
    }

  //----------------------- POLY -------------


Vector SSurface::PointAt(Point2d puv) const {
    return PointAt(puv.x, puv.y);
}
Vector SSurface::PointAt(double u, double v) const {
    Vector num = {};
    double den = 0;

    int i, j;
    for(i = 0; i <= degm; i++) {
        for(j = 0; j <= degn; j++) {
            double Bi = Bernstein(i, degm, u),
                   Bj = Bernstein(j, degn, v);

            num = num.Plus(ctrl[i][j].ScaledBy(Bi*Bj*weight[i][j]));
            den += weight[i][j]*Bi*Bj;
        }
    }
    num = num.ScaledBy(1.0/den);
    return num;
}

void SSurface::TangentsAt(double u, double v, Vector *tu, Vector *tv, bool retry) const {
    Vector num   = {},
           num_u = {},
           num_v = {};
    double den   = 0,
           den_u = 0,
           den_v = 0;

    int i, j;
    for(i = 0; i <= degm; i++) {
        for(j = 0; j <= degn; j++) {
            double Bi  = Bernstein(i, degm, u),
                   Bj  = Bernstein(j, degn, v),
                   Bip = BernsteinDerivative(i, degm, u),
                   Bjp = BernsteinDerivative(j, degn, v);

            num = num.Plus(ctrl[i][j].ScaledBy(Bi*Bj*weight[i][j]));
            den += weight[i][j]*Bi*Bj;

            num_u = num_u.Plus(ctrl[i][j].ScaledBy(Bip*Bj*weight[i][j]));
            den_u += weight[i][j]*Bip*Bj;

            num_v = num_v.Plus(ctrl[i][j].ScaledBy(Bi*Bjp*weight[i][j]));
            den_v += weight[i][j]*Bi*Bjp;
        }
    }
    // quotient rule; f(t) = n(t)/d(t), so f' = (n'*d - n*d')/(d^2)
    *tu = ((num_u.ScaledBy(den)).Minus(num.ScaledBy(den_u)));
    *tu = tu->ScaledBy(1.0/(den*den));

    *tv = ((num_v.ScaledBy(den)).Minus(num.ScaledBy(den_v)));
    *tv = tv->ScaledBy(1.0/(den*den));
    
    // Tangent is zero at sungularities like the north pole. Move away a bit and retry. 
    if(tv->Equals({0, 0, 0}) && retry)
        TangentsAt(u+(0.5-u)*0.00001, v, tu, tv, false);
    if(tu->Equals({0, 0, 0}) && retry)
        TangentsAt(u, v+(0.5-v)*0.00001, tu, tv, false);
}

Vector SSurface::NormalAt(Point2d puv) const {
    return NormalAt(puv.x, puv.y);
}

Vector SSurface::NormalAt(double u, double v) const {
    Vector tu, tv;
    TangentsAt(u, v, &tu, &tv);
    return tu.Cross(tv);
}

void SSurface::ClosestPointTo(Vector p, Point2d *puv, bool mustConverge) {
    ClosestPointTo(p, &(puv->x), &(puv->y), mustConverge);
}

void SSurface::ClosestPointTo(Vector p, double *u, double *v, bool mustConverge) {
    // A few special cases first; when control points are coincident the
    // derivative goes to zero at the control points, and would result in
    // nonconvergence. We avoid that here, and also guarantee a consistent
    // (u, v) (of the infinitely many possible in one parameter).
    if(p.Equals(ctrl[0]   [0]   )) { *u = 0; *v = 0; return; }
    if(p.Equals(ctrl[degm][0]   )) { *u = 1; *v = 0; return; }
    if(p.Equals(ctrl[degm][degn])) { *u = 1; *v = 1; return; }
    if(p.Equals(ctrl[0]   [degn])) { *u = 0; *v = 1; return; }

    // And planes are trivial, so don't waste time iterating over those.
    if(degm == 1 && degn == 1) {
        Vector orig =  ctrl[0][0],
               bu   = (ctrl[1][0]).Minus(orig),
               bv   = (ctrl[0][1]).Minus(orig);
        if((ctrl[1][1]).Equals(orig.Plus(bu).Plus(bv))) {

            Vector n = bu.Cross(bv);
            Vector ty = n.Cross(bu).ScaledBy(1.0/bu.MagSquared());
            Vector tx = bv.Cross(n).ScaledBy(1.0/bv.MagSquared());

            Vector dp = p.Minus(orig);
            *u = dp.Dot(bu) / tx.MagSquared();
            *v = dp.Dot(bv) / ty.MagSquared();
            return;
        }
    }

    // Try whatever the previous guess was. This is likely to do something
    // good if we're working our way along a curve or something else where
    // we project successive points that are close to each other; something
    // like a 20% speedup empirically.
    if(mustConverge) {
        double ut = cached.x, vt = cached.y;
        if(ClosestPointNewton(p, &ut, &vt, mustConverge)) {
            cached.x = *u = ut;
            cached.y = *v = vt;
            return;
        }
    }

    // Search for a reasonable initial guess
    int i, j;
    double minDist = VERY_POSITIVE;
    int res = (max(degm, degn) == 2) ? 7 : 20;
    for(i = 0; i < res; i++) {
        for(j = 0; j < res; j++) {
            double tryu = (i + 0.5)/res, tryv = (j + 0.5)/res;

            Vector tryp = PointAt(tryu, tryv);
            double d = (tryp.Minus(p)).Magnitude();
            if(d < minDist) {
                *u = tryu;
                *v = tryv;
                minDist = d;
            }
        }
    }

    if(ClosestPointNewton(p, u, v, mustConverge)) {
        cached.x = *u;
        cached.y = *v;
        return;
    }

    // If we failed to converge, then at least don't return NaN.
    if(mustConverge) {
// This is expected not to converge when the target point is not on the surface but nearby.
// let's not pollute the output window for normal use.
//        Vector p0 = PointAt(*u, *v);
//        dbp("didn't converge");
//        dbp("have %.3f %.3f %.3f", CO(p0));
//        dbp("want %.3f %.3f %.3f", CO(p));
//        dbp("distance = %g", (p.Minus(p0)).Magnitude());
    }
    if(IsReasonable(*u) || IsReasonable(*v)) {
        *u = *v = 0;
    }
}

bool SSurface::ClosestPointNewton(Vector p, double *u, double *v, bool mustConverge) const
{
    // Initial guess is in u, v; refine by Newton iteration.
    Vector p0 = {};
    for(int i = 0; i < (mustConverge ? 25 : 5); i++) {
        p0 = PointAt(*u, *v);
        if(mustConverge) {
            if(p0.Equals(p, RATPOLY_EPS)) {
                return true;
            }
        }

        Vector tu, tv, tx, ty;
        TangentsAt(*u, *v, &tu, &tv);
        Vector n = tu.Cross(tv);
        // since tu and tv may not be orthogonal, use y in place of v.
        // |y| = |v|sin(theta) where theta is the angle between tu and tv.
        ty = n.Cross(tu).ScaledBy(1.0/tu.MagSquared());
        tx = tv.Cross(n).ScaledBy(1.0/tv.MagSquared());

        // Project the point into a plane through p0, with basis tu, tv; a
        // second-order thing would converge faster but needs second
        // derivatives.
        Vector dp = p.Minus(p0);
        double du = dp.Dot(tx),
               dv = dp.Dot(ty);
        *u += du / (tx.MagSquared());
        *v += dv / (ty.MagSquared());

        if (*u < 0.0) *u = 0.0;
        else if (*u > 1.0) *u = 1.0;
        if (*v < 0.0) *v = 0.0;
        else if (*v > 1.0) *v = 1.0;

    }

    return false;
}

// This checks the nearest edge of our surface for intersection with curve.
// It is common for curve/surface intersection to fail because the curve
// is tangent to the surface at the intersection point, which is often on
// the edge of the surface - that's just how people design things...
bool SSurface::EdgeCurveIntersection(double *u, double *v, SBezier *curve) const
{
    SBezier edge = {};
    // we already have a starting point in space
    Vector p = PointAt(*u, *v);
    if(fabs(*u - 0.5) > fabs(*v - 0.5)) {
    // u is closer to 0 or 1 than v
      if(*u < 0.5) {
        for(int n=0;n<4;n++) { edge.ctrl[n]=ctrl[0][n]; edge.weight[n]=weight[0][n]; }
        edge.deg = degn;
        if (edge.PointOnNonparallelCurve(curve, &p)) {
            edge.ClosestPointTo(p, v); *u=0.0;
            dbp("found U=0.0, V=%g point=(%.3f %.3f %.3f)", *v, CO(p) );
            return true;
        }
      }
      else { // u > 0.5
        for(int n=0;n<4;n++) { edge.ctrl[n]=ctrl[degm][n]; edge.weight[n]=weight[degm][n]; }
        edge.deg = degn;              
        if (edge.PointOnNonparallelCurve(curve, &p)) {
            edge.ClosestPointTo(p, v); *u=1.0;
            dbp("found U=1.0, V=%g point=(%.3f %.3f %.3f)", *v, CO(p) );
            return true;
        }
      }
    } else {
    // v is closer to 0 or 1
      if(*v < 0.5) {
        for(int n=0;n<4;n++) { edge.ctrl[n]=ctrl[n][0]; edge.weight[n]=weight[n][0]; }
        edge.deg = degm;
        if (edge.PointOnNonparallelCurve(curve, &p)) {
            edge.ClosestPointTo(p, u); *v=0.0;
            dbp("found U=%g, V=0.0 point=(%.3f %.3f %.3f)", *u, CO(p));
            return true;
        }      
      } else { // v > 0.5
        for(int n=0;n<4;n++) { edge.ctrl[n]=ctrl[n][degn]; edge.weight[n]=weight[n][degn]; }
        edge.deg = degm;
        if (edge.PointOnNonparallelCurve(curve, &p)) {
            edge.ClosestPointTo(p, u); *v=1.0;
            dbp("found U=%g, V=1.0 point=(%.3f %.3f %.3f)", *u, CO(p));
            return true;
        }            
      }
    }
    return false;
}
bool SSurface::PointIntersectingLine(Vector p0, Vector p1, double *u, double *v) const
{
    int i;
    for(i = 0; i < 20; i++) {
        Vector pi, p, tu, tv;
        p = PointAt(*u, *v);
        TangentsAt(*u, *v, &tu, &tv);

        Vector n = (tu.Cross(tv)).WithMagnitude(1);
        double d = p.Dot(n);

        bool parallel;
        pi = Vector::AtIntersectionOfPlaneAndLine(n, d, p0, p1, &parallel);
        if(parallel) {
            dbp("parallel (surface intersecting line)");
            break;
        }

        // Check for convergence - success exits
        if(pi.Equals(p, RATPOLY_EPS)) return true;

        n = tu.Cross(tv);
        Vector ty = n.Cross(tu).ScaledBy(1.0/tu.MagSquared());
        Vector tx = tv.Cross(n).ScaledBy(1.0/tv.MagSquared());

        // Adjust our guess and iterate
        Vector dp = pi.Minus(p);
        double du = dp.Dot(tx), dv = dp.Dot(ty);
        *u += du / tx.MagSquared();
        *v += dv / ty.MagSquared();
    }
    // Lets try the edge of the surface against the line before giving up
    SBezier curve = {};
    curve.ctrl[0] = p0; curve.weight[0] = 1.0;
    curve.ctrl[1] = p1; curve.weight[1] = 1.0;
    curve.deg = 1;
    if (EdgeCurveIntersection(u,v, &curve)) {
      return true;
    }
    dbp("didn't converge (surface intersecting line)");
    return false;
}

Vector SSurface::ClosestPointOnThisAndSurface(SSurface *srf2, Vector p) {
    // This is untested.
    int i, j;
    Point2d puv[2];
    SSurface *srf[2] = { this, srf2 };

    for(j = 0; j < 2; j++) {
        (srf[j])->ClosestPointTo(p, &(puv[j]), /*mustConverge=*/false);
    }

    for(i = 0; i < 10; i++) {
        Vector tu[2], tv[2], cp[2], n[2];
        double d[2];

        for(j = 0; j < 2; j++) {
            (srf[j])->TangentsAt(puv[j].x, puv[j].y, &(tu[j]), &(tv[j]));

            cp[j] = (srf[j])->PointAt(puv[j]);

            n[j] = ((tu[j]).Cross(tv[j])).WithMagnitude(1);
            d[j] = (n[j]).Dot(cp[j]);
        }

        if((cp[0]).Equals(cp[1], RATPOLY_EPS)) break;

        Vector p0 = Vector::AtIntersectionOfPlanes(n[0], d[0], n[1], d[1]),
               dp = (n[0]).Cross(n[1]);

        Vector pc = p.ClosestPointOnLine(p0, dp);

        // Adjust our guess and iterate
        for(j = 0; j < 2; j++) {
            Vector n = tu[j].Cross(tv[j]);
            Vector ty = n.Cross(tu[j]).ScaledBy(1.0/tu[j].MagSquared());
            Vector tx = tv[j].Cross(n).ScaledBy(1.0/tv[j].MagSquared());

            Vector dc = pc.Minus(cp[j]);
            double du = dc.Dot(tx), dv = dc.Dot(ty);
            puv[j].x += du / tx.MagSquared();
            puv[j].y += dv / ty.MagSquared();
        }
    }
    if(i >= 10) {
        dbp("this and srf, didn't converge, d=%g",
            (puv[0].Minus(puv[1])).Magnitude());
    }

    // If this converged, then the two points are actually equal.
    return ((srf[0])->PointAt(puv[0])).Plus(
           ((srf[1])->PointAt(puv[1]))).ScaledBy(0.5);
}

void SSurface::PointOnSurfaces(SSurface *s1, SSurface *s2, double *up, double *vp)
{
    double u[3] = { *up, 0, 0 }, v[3] = { *vp, 0, 0 };
    SSurface *srf[3] = { this, s1, s2 };

    // Get initial guesses for (u, v) in the other surfaces
    Vector p = PointAt(*u, *v);
    (srf[1])->ClosestPointTo(p, &(u[1]), &(v[1]), /*mustConverge=*/false);
    (srf[2])->ClosestPointTo(p, &(u[2]), &(v[2]), /*mustConverge=*/false);

    int i, j;
    for(i = 0; i < 20; i++) {
        // Approximate each surface by a plane
        Vector p[3], tu[3], tv[3], n[3];
        double d[3];
        for(j = 0; j < 3; j++) {
            p[j] = (srf[j])->PointAt(u[j], v[j]);
            (srf[j])->TangentsAt(u[j], v[j], &(tu[j]), &(tv[j]));
            n[j] = ((tu[j]).Cross(tv[j])).WithMagnitude(1);
            d[j] = (n[j]).Dot(p[j]);
        }

        // If a = b and b = c, then does a = c? No, it doesn't.
        if((p[0]).Equals(p[1], RATPOLY_EPS) &&
           (p[1]).Equals(p[2], RATPOLY_EPS) &&
           (p[2]).Equals(p[0], RATPOLY_EPS))
        {
            *up = u[0];
            *vp = v[0];
            return;
        }

        bool parallel;
        Vector pi = Vector::AtIntersectionOfPlanes(n[0], d[0],
                                                   n[1], d[1],
                                                   n[2], d[2], &parallel);

        if(parallel) { // lets try something else for parallel planes
            pi = p[0].Plus(p[1]).Plus(p[2]).ScaledBy(1.0/3.0);
        }

        for(j = 0; j < 3; j++) {
            Vector n = tu[j].Cross(tv[j]);
            Vector ty = n.Cross(tu[j]).ScaledBy(1.0/tu[j].MagSquared());
            Vector tx = tv[j].Cross(n).ScaledBy(1.0/tv[j].MagSquared());

            Vector dp = pi.Minus(p[j]);
            double du = dp.Dot(tx), dv = dp.Dot(ty);

            u[j] += du / tx.MagSquared();
            v[j] += dv / ty.MagSquared();
        }
    }
    dbp("didn't converge (three surfaces intersecting)");
}

void SSurface::PointOnCurve(const SBezier *curve, double *up, double *vp)
{
    Vector tu,tv,n;
    double u = *up, v = *vp;
    Vector ps = PointAt(u, v);
    // Get initial guesses for t on the curve
    double tCurve = 0.5;
    curve->ClosestPointTo(ps, &tCurve, /*mustConverge=*/false);
    if(tCurve < 0.0) tCurve = 0.0;
    if(tCurve > 1.0) tCurve = 1.0;

    for(int i = 0; i < 30; i++) {
        // Approximate the surface by a plane
        Vector ps = PointAt(u, v);
        TangentsAt(u, v, &tu, &tv);
        n = tu.Cross(tv).WithMagnitude(1);

        // point on curve and tangent line direction
        Vector pc = curve->PointAt(tCurve);
        Vector tc = curve->TangentAt(tCurve);

        if(ps.Equals(pc, RATPOLY_EPS)) {
            *up = u;
            *vp = v;
            return;
        }

        //pi is where the curve tangent line intersects the surface tangent plane
        Vector pi;
        double d = tc.Dot(n);
        if (fabs(d) < 1e-10) { // parallel line and plane, guess the average rather than fail
            pi = pc.Plus(ps).ScaledBy(0.5);
        } else {
            pi = pc.Minus(tc.ScaledBy(pc.Minus(ps).Dot(n)/d));
        }

        // project the point onto the tangent plane and line
        {
            Vector n = tu.Cross(tv);
            Vector ty = n.Cross(tu).ScaledBy(1.0/tu.MagSquared());
            Vector tx = tv.Cross(n).ScaledBy(1.0/tv.MagSquared());

            Vector dp = pi.Minus(ps);
            double du = dp.Dot(tx), dv = dp.Dot(ty);

            u += du / tx.MagSquared();
            v += dv / ty.MagSquared();
        }
        tCurve += pi.Minus(pc).Dot(tc) / tc.MagSquared();
        if(tCurve < 0.0) tCurve = 0.0;
        if(tCurve > 1.0) tCurve = 1.0;
    }
    dbp("didn't converge (surface and curve intersecting)");
}


void SSurface::TrimFromEdgeList(SEdgeList *el, bool asUv) {
    el->l.ClearTags();

    STrimBy stb = {};
    for(;;) {
        // Find an edge, any edge; we'll start from there.
        SEdge *se;
        for(se = el->l.First(); se; se = el->l.NextAfter(se)) {
            if(se->tag) continue;
            break;
        }
        if(!se) break;
        se->tag = 1;
        stb.start = se->a;
        stb.finish = se->b;
        stb.curve.v = se->auxA;
        stb.backwards = se->auxB ? true : false;

        // Find adjoining edges from the same curve; those should be
        // merged into a single trim.
        bool merged;
        do {
            merged = false;
            for(se = el->l.First(); se; se = el->l.NextAfter(se)) {
                if(se->tag)                         continue;
                if(se->auxA != (int)stb.curve.v)    continue;
                if(( se->auxB && !stb.backwards) ||
                   (!se->auxB &&  stb.backwards))   continue;

                if((se->a).Equals(stb.finish)) {
                    stb.finish = se->b;
                    se->tag = 1;
                    merged = true;
                } else if((se->b).Equals(stb.start)) {
                    stb.start = se->a;
                    se->tag = 1;
                    merged = true;
                }
            }
        } while(merged);

        if(asUv) {
            stb.start  = PointAt(stb.start.x,  stb.start.y);
            stb.finish = PointAt(stb.finish.x, stb.finish.y);
        }

        // And add the merged trim, with xyz (not uv like the polygon) pts
        trim.Add(&stb);
    }
}

static bool KeepRegion(SSurface::CombineAs type, bool opA, SShell::Class shell, SShell::Class orig)
{
    bool inShell = (shell == SShell::Class::SURF_INSIDE),
         outSide = (shell == SShell::Class::SURF_OUTSIDE),
         coincSame  = (shell == SShell::Class::SURF_COINC_SAME),
         coincOpp  = (shell == SShell::Class::SURF_COINC_OPP),
         inOrig  = (orig == SShell::Class::SURF_INSIDE);

    // This one line is not really part of this functions logic
    if(!inOrig) return false;

    switch(type) {
        case SSurface::CombineAs::UNION:
            if(opA) {
                return outSide;
            } else {
                return outSide || coincSame;
            }

        case SSurface::CombineAs::DIFFERENCE:
            if(opA) {
                return outSide || coincOpp;
            } else {
                return inShell;
            }

        case SSurface::CombineAs::INTERSECTION:
            if(opA) {
                return inShell;
            } else {
                return inShell || coincSame;
            }

        default: ssassert(false, "Unexpected combine type");
    }
}
static bool KeepEdge(SSurface::CombineAs type, bool opA,
                     SShell::Class indir_shell, SShell::Class outdir_shell,
                     SShell::Class indir_orig, SShell::Class outdir_orig)
{
    bool keepIn  = KeepRegion(type, opA, indir_shell,  indir_orig),
         keepOut = KeepRegion(type, opA, outdir_shell, outdir_orig);

    // If the regions to the left and right of this edge are both in or both
    // out, then this edge is not useful and should be discarded.
    if(keepIn && !keepOut) return true;
    return false;
}

static void TagByClassifiedEdge(SBspUv::Class bspclass, SShell::Class *indir, SShell::Class *outdir)
{
    switch(bspclass) {
        case SBspUv::Class::INSIDE:
            *indir  = SShell::Class::SURF_INSIDE;
            *outdir = SShell::Class::SURF_INSIDE;
            break;

        case SBspUv::Class::OUTSIDE:
            *indir  = SShell::Class::SURF_OUTSIDE;
            *outdir = SShell::Class::SURF_OUTSIDE;
            break;

        case SBspUv::Class::EDGE_PARALLEL:
            *indir  = SShell::Class::SURF_INSIDE;
            *outdir = SShell::Class::SURF_OUTSIDE;
            break;

        case SBspUv::Class::EDGE_ANTIPARALLEL:
            *indir  = SShell::Class::SURF_OUTSIDE;
            *outdir = SShell::Class::SURF_INSIDE;
            break;

        default:
            dbp("TagByClassifiedEdge: fail!");
            *indir  = SShell::Class::SURF_OUTSIDE;
            *outdir = SShell::Class::SURF_OUTSIDE;
            break;
    }
}

static void DEBUGEDGELIST(SEdgeList *sel, SSurface *surf) {
    dbp("print %d edges", sel->l.n);
    SEdge *se;
    for(se = sel->l.First(); se; se = sel->l.NextAfter(se)) {
        Vector mid = (se->a).Plus(se->b).ScaledBy(0.5);
        Vector arrow = (se->b).Minus(se->a);
        swap(arrow.x, arrow.y);
        arrow.x *= -1;
        arrow = arrow.WithMagnitude(0.01);
        arrow = arrow.Plus(mid);

        SS.nakedEdges.AddEdge(surf->PointAt(se->a.x, se->a.y),
                              surf->PointAt(se->b.x, se->b.y));
        SS.nakedEdges.AddEdge(surf->PointAt(mid.x, mid.y),
                              surf->PointAt(arrow.x, arrow.y));
    }
}

//-----------------------------------------------------------------------------
// We are given src, with at least one edge, and avoid, a list of points to
// avoid. We return a chain of edges (that share endpoints), such that no
// point within the avoid list ever occurs in the middle of a chain. And we
// delete the edges in that chain from our source list.
//-----------------------------------------------------------------------------
void SSurface::FindChainAvoiding(SEdgeList *src, SEdgeList *dest,
                                 SPointList *avoid)
{
    ssassert(!src->l.IsEmpty(), "Need at least one edge");
    // Start with an arbitrary edge.
    dest->l.Add(src->l.First());
    src->l.ClearTags();
    src->l.First()->tag = 1;

    bool added;
    do {
        added = false;
        // The start and finish of the current edge chain
        Vector s = dest->l.First()->a,
               f = dest->l.Last()->b;

        // We can attach a new edge at the start or finish, as long as that
        // start or finish point isn't in the list of points to avoid.
        bool startOkay  = !avoid->ContainsPoint(s),
             finishOkay = !avoid->ContainsPoint(f);

        // Now look for an unused edge that joins at the start or finish of
        // our chain (if permitted by the avoid list).
        SEdge *se;
        for(se = src->l.First(); se; se = src->l.NextAfter(se)) {
            if(se->tag) continue;
            if(startOkay && s.Equals(se->b)) {
                dest->l.AddToBeginning(se);
                s = se->a;
                se->tag = 1;
                startOkay = !avoid->ContainsPoint(s);
            } else if(finishOkay && f.Equals(se->a)) {
                dest->l.Add(se);
                f = se->b;
                se->tag = 1;
                finishOkay = !avoid->ContainsPoint(f);
            } else {
                continue;
            }

            added = true;
        }
    } while(added);

    src->l.RemoveTagged();
}

void SSurface::EdgeNormalsWithinSurface(Point2d auv, Point2d buv,
                                        Vector *pt,
                                        Vector *enin, Vector *enout,
                                        Vector *surfn,
                                        uint32_t auxA,
                                        SShell *shell, SShell *sha, SShell *shb)
{
    // the midpoint of the edge
    Point2d muv  = (auv.Plus(buv)).ScaledBy(0.5);

    *pt    = PointAt(muv);

    // If this edge just approximates a curve, then refine our midpoint so
    // so that it actually lies on that curve too. Otherwise stuff like
    // point-on-face tests will fail, since the point won't actually lie
    // on the other face.
    hSCurve hc = { auxA };
    SCurve *sc = shell->curve.FindById(hc);
    if(sc->isExact && sc->exact.deg != 1) {
        double t;
        sc->exact.ClosestPointTo(*pt, &t, /*mustConverge=*/false);
        *pt = sc->exact.PointAt(t);
    } else if(!sc->isExact) {
        SSurface *trimmedA = sc->GetSurfaceA(sha, shb),
                       *trimmedB = sc->GetSurfaceB(sha, shb);
        *pt = trimmedA->ClosestPointOnThisAndSurface(trimmedB, *pt);
    }
    ClosestPointTo(*pt, &muv);

    *surfn = NormalAt(muv.x, muv.y);

    // Compute the edge's inner normal in xyz space. Use the configured
    // chord tolerance as the probe distance, but never more than a
    // fraction of this edge's own length: the classification must probe
    // the surfaces close to the edge, and a coarse tolerance on a small
    // model would otherwise evaluate them far away, or even extrapolate
    // them outside their domain, and misclassify.
    Vector ab    = (PointAt(auv)).Minus(PointAt(buv)),
           enxyz = (ab.Cross(*surfn)).WithMagnitude(
                        min(SS.ChordTolMm(), ab.Magnitude() / 10));
    // And based on that, compute the edge's inner normal in uv space. This
    // vector is perpendicular to the edge in xyz, but not necessarily in uv.
    Vector tu, tv, tx, ty;
    TangentsAt(muv.x, muv.y, &tu, &tv);
    Vector n = tu.Cross(tv);
    // since tu and tv may not be orthogonal, use y in place of v, x in place of u.
    // |y| = |v|sin(theta) where theta is the angle between tu and tv.
    ty = n.Cross(tu).ScaledBy(1.0/tu.MagSquared());
    tx = tv.Cross(n).ScaledBy(1.0/tv.MagSquared());

    Point2d enuv;
    enuv.x = enxyz.Dot(tx) / tx.MagSquared();
    enuv.y = enxyz.Dot(ty) / ty.MagSquared();

    // Compute the inner and outer normals of this edge (within the srf),
    // in xyz space. These are not necessarily antiparallel, if the
    // surface is curved.
    Vector pin   = PointAt(muv.Minus(enuv)),
           pout  = PointAt(muv.Plus(enuv));
    *enin  = pin.Minus(*pt),
    *enout = pout.Minus(*pt);
// ideally this should work (fail screwdriver file)    
//    *enin = enxyz.ScaledBy(-1.0);
//    *enout = enxyz;
}

//-----------------------------------------------------------------------------
// Trim this surface against the specified shell, in the way that's appropriate
// for the specified Boolean operation type (and which operand we are). We
// also need a pointer to the shell that contains our own surface, since that
// contains our original trim curves.
//-----------------------------------------------------------------------------
SSurface SSurface::MakeCopyTrimAgainst(SShell *parent,
                                       SShell *sha, SShell *shb,
                                       SShell *into,
                                       SSurface::CombineAs type,
                                       int dbg_index)
{
    bool opA = (parent == sha);
    SShell *agnst = opA ? shb : sha;

    SSurface ret;
    // The returned surface is identical, just the trim curves change
    ret = *this;
    ret.trim = {};

    // First, build a list of the existing trim curves; update them to use
    // the split curves.
    STrimBy *stb;
    for(stb = trim.First(); stb; stb = trim.NextAfter(stb)) {
        STrimBy stn = *stb;
        stn.curve = (parent->curve.FindById(stn.curve))->newH;
        ret.trim.Add(&stn);
    }

    if(type == SSurface::CombineAs::DIFFERENCE && !opA) {
        // The second operand of a Boolean difference gets turned inside out
        ret.Reverse();
    }

    // Build up our original trim polygon; remember the coordinates could
    // be changed if we just flipped the surface normal, and we are using
    // the split curves (not the original curves).
    SEdgeList orig = {};
    ret.MakeEdgesInto(into, &orig, MakeAs::UV);
    ret.trim.Clear();
    // which means that we can't necessarily use the old BSP...
    SBspUv *origBsp = SBspUv::From(&orig, &ret);

    // And now intersect the other shell against us
    SEdgeList inter = {};

    SSurface *ss;
    for(SCurve &sc : into->curve) {
        if(sc.source != SCurve::Source::INTERSECTION) continue;
        if(opA) {
            if(sc.surfA != h) continue;
            ss = shb->surface.FindById(sc.surfB);
        } else {
            if(sc.surfB != h) continue;
            ss = sha->surface.FindById(sc.surfA);
        }
        int i;
        for(i = 1; i < sc.pts.n; i++) {
            Vector a = sc.pts[i-1].p,
                   b = sc.pts[i].p;

            Point2d auv, buv;
            ss->ClosestPointTo(a, &(auv.x), &(auv.y));
            ss->ClosestPointTo(b, &(buv.x), &(buv.y));

            SBspUv::Class c = (ss->bsp) ? ss->bsp->ClassifyEdge(auv, buv, ss) : SBspUv::Class::OUTSIDE;
            if(c != SBspUv::Class::OUTSIDE) {
                Vector ta = {};
                Vector tb = {};
                ret.ClosestPointTo(a, &(ta.x), &(ta.y));
                ret.ClosestPointTo(b, &(tb.x), &(tb.y));

                Vector tn = ret.NormalAt(ta.x, ta.y);
                Vector sn = ss->NormalAt(auv.x, auv.y);

                // We are subtracting the portion of our surface that
                // lies in the shell, so the in-plane edge normal should
                // point opposite to the surface normal.
                Vector tnxd = tn.Cross(b.Minus(a));
                double dot = tnxd.Dot(sn);
                if(fabs(dot) < SShell::DOTP_TOL*tnxd.Magnitude()*sn.Magnitude()) {
                    // The surfaces are tangent along this curve, so their
                    // normals give no orientation for it. Add the edge in
                    // both directions; the classification against the
                    // shells keeps the correctly oriented one and
                    // discards the other.
                    inter.AddEdge(ta, tb, sc.h.v, 0);
                    inter.AddEdge(tb, ta, sc.h.v, 1);
                    continue;
                }
                bool bkwds = true;
                if(dot < 0) bkwds = !bkwds;
                if((type == SSurface::CombineAs::DIFFERENCE && !opA) ||
                   (type == SSurface::CombineAs::INTERSECTION)) { // Invert all newly created edges for intersection
                    bkwds = !bkwds;
                }
                if(bkwds) {
                    inter.AddEdge(tb, ta, sc.h.v, 1);
                } else {
                    inter.AddEdge(ta, tb, sc.h.v, 0);
                }
            }
        }
    }

    // Record all the points where more than two edges join, which I will call
    // the choosing points. If two edges join at a non-choosing point, then
    // they must either both be kept or both be discarded (since that would
    // otherwise create an open contour).
    SPointList choosing = {};
    SEdge *se;
    for(se = orig.l.First(); se; se = orig.l.NextAfter(se)) {
        choosing.IncrementTagFor(se->a);
        choosing.IncrementTagFor(se->b);
    }
    for(se = inter.l.First(); se; se = inter.l.NextAfter(se)) {
        choosing.IncrementTagFor(se->a);
        choosing.IncrementTagFor(se->b);
    }
    SPoint *sp;
    for(sp = choosing.l.First(); sp; sp = choosing.l.NextAfter(sp)) {
        if(sp->tag == 2) {
            sp->tag = 1;
        } else {
            sp->tag = 0;
        }
    }
    choosing.l.RemoveTagged();

    // The list of edges to trim our new surface, a combination of edges from
    // our original and intersecting edge lists.
    SEdgeList final = {};

    while(!orig.l.IsEmpty()) {
        SEdgeList chain = {};
        FindChainAvoiding(&orig, &chain, &choosing);

        // Arbitrarily choose an edge within the chain to classify; they
        // should all be the same, though.
        se = &(chain.l[chain.l.n/2]);

        Point2d auv  = (se->a).ProjectXy(),
                buv  = (se->b).ProjectXy();

        Vector pt, enin, enout, surfn;
        ret.EdgeNormalsWithinSurface(auv, buv, &pt, &enin, &enout, &surfn,
                                        se->auxA, into, sha, shb);

        // Initialize to a deterministic (if arbitrary) guess, so that if the
        // classification fails we don't use values from uninitialized stack.
        SShell::Class indir_shell  = SShell::Class::SURF_OUTSIDE,
                      outdir_shell = SShell::Class::SURF_OUTSIDE,
                      indir_orig, outdir_orig;

        indir_orig  = SShell::Class::SURF_INSIDE;
        outdir_orig = SShell::Class::SURF_OUTSIDE;

        if(!agnst->ClassifyEdge(&indir_shell, &outdir_shell,
                                ret.PointAt(auv), ret.PointAt(buv), pt,
                                enin, enout, surfn))
        {
            dbp("MakeCopyTrimAgainst: failed to classify orig edge (I=%d)", I+dbg_index);
        }

        if(KeepEdge(type, opA, indir_shell, outdir_shell,
                               indir_orig,  outdir_orig))
        {
            for(se = chain.l.First(); se; se = chain.l.NextAfter(se)) {
                final.AddEdge(se->a, se->b, se->auxA, se->auxB);
            }
        }
        chain.Clear();
    }

    while(!inter.l.IsEmpty()) {
        SEdgeList chain = {};
        FindChainAvoiding(&inter, &chain, &choosing);

        // Any edge in the chain, same as above.
        se = &(chain.l[chain.l.n/2]);

        Point2d auv = (se->a).ProjectXy(),
                buv = (se->b).ProjectXy();

        Vector pt, enin, enout, surfn;
        ret.EdgeNormalsWithinSurface(auv, buv, &pt, &enin, &enout, &surfn,
                                        se->auxA, into, sha, shb);

        // Initialize to a deterministic (if arbitrary) guess, so that if the
        // classification fails we don't use values from uninitialized stack.
        SShell::Class indir_shell  = SShell::Class::SURF_OUTSIDE,
                      outdir_shell = SShell::Class::SURF_OUTSIDE,
                      indir_orig, outdir_orig;

        SBspUv::Class c_this = (origBsp) ? origBsp->ClassifyEdge(auv, buv, &ret) : SBspUv::Class::OUTSIDE;

        if(c_this == SBspUv::Class::EDGE_PARALLEL) {
            // The intersection edge lies exactly along an edge of our
            // original trim polygon, in the same direction. Whatever trim
            // is required there, the original edge (which the loop above
            // classifies identically, since for both edges the decision
            // reduces to keeping iff the region on the in-side of the shared
            // line is kept) already provides it, so this edge is redundant.
            // Keeping it would at best duplicate the original edge exactly
            // (and get culled below), but if the two copies are split at
            // different interior points the duplicates survive the cull and
            // the trim polygon fails to assemble, which shows as a missing
            // face when several coplanar faces join edge-on-edge across
            // multiple union steps (issue #1452). Discard it instead.
            chain.Clear();
            continue;
        }

        TagByClassifiedEdge(c_this, &indir_orig, &outdir_orig);

        if(!agnst->ClassifyEdge(&indir_shell, &outdir_shell,
                                ret.PointAt(auv), ret.PointAt(buv), pt,
                                enin, enout, surfn))
        {
            dbp("MakeCopyTrimAgainst: failed to classify inter edge (I=%d)", I+dbg_index);
        }

        if(KeepEdge(type, opA, indir_shell, outdir_shell,
                               indir_orig,  outdir_orig))
        {
            for(se = chain.l.First(); se; se = chain.l.NextAfter(se)) {
                final.AddEdge(se->a, se->b, se->auxA, se->auxB);
            }
        }
        chain.Clear();
    }

    // Cull extraneous edges; duplicates or anti-parallel pairs. In particular,
    // we can get duplicate edges if our surface intersects the other shell
    // at an edge, so that both surfaces intersect coincident (and both
    // generate an intersection edge).
    final.CullExtraneousEdges(/*both=*/true);

    // Use our reassembled edges to trim the new surface.
    ret.TrimFromEdgeList(&final, /*asUv=*/true);

    SPolygon poly = {};
    final.l.ClearTags();
    if(!final.AssemblePolygon(&poly, NULL, /*keepDir=*/true))
#pragma omp critical
    {
        into->booleanFailed = true;
        dbp("failed: I=%d, avoid=%d", I+dbg_index, choosing.l.n);
        DEBUGEDGELIST(&final, &ret);
    }
    poly.Clear();

    choosing.Clear();
    final.Clear();
    inter.Clear();
    orig.Clear();
    return ret;
}


//-----------------------------------------------------------------------------
// All of the BSP routines that we use to perform and accelerate polygon ops.
//-----------------------------------------------------------------------------
void SSurface::MakeClassifyingBsp(SShell *shell, SShell *useCurvesFrom) {
    SEdgeList el = {};

    MakeEdgesInto(shell, &el, MakeAs::UV, useCurvesFrom);
    bsp = SBspUv::From(&el, this);
    el.Clear();

    edges = {};
    MakeEdgesInto(shell, &edges, MakeAs::XYZ, useCurvesFrom);
}

} // namespace SolveSpace
