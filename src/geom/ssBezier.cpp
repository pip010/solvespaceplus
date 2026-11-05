#include "ssBezier.h"

namespace SolveSpace
{
    static inline double Bernstein(int k, int deg, double t) {
    // indexed by [degree][k][exponent]
    static const double bernstein_coeff[4][4][4] = {
    { { 1.0,0.0,0.0,0.0 },  { 1.0,0.0,0.0,0.0 }, { 1.0,0.0,0.0,0.0 }, { 1.0,0.0,0.0,0.0 } },
    { { 1.0,-1.0,0.0,0.0 }, { 0.0,1.0,0.0,0.0 }, { 0.0,0.0,0.0,0.0 }, { 0.0,0.0,0.0,0.0 } },
    { { 1.0,-2.0,1.0,0.0 }, { 0.0,2.0,-2.0,0.0 },{ 0.0,0.0,1.0,0.0 }, { 0.0,0.0,0.0,0.0 } },
    { { 1.0,-3.0,3.0,-1.0 },{ 0.0,3.0,-6.0,3.0 },{ 0.0,0.0,3.0,-3.0}, { 0.0,0.0,0.0,1.0 } } };

    const double *c = bernstein_coeff[deg][k];
    return (((c[3]*t+c[2])*t)+c[1])*t+c[0];
}

static inline double BernsteinDerivative(int k, int deg, double t) {
    static const double bernstein_derivative_coeff[4][4][3] = {
    { { 0.0,0.0,0.0 },  { 0.0,0.0,0.0 }, { 0.0,0.0,0.0 }, { 0.0,0.0,0.0 } },
    { { -1.0,0.0,0.0 }, { 1.0,0.0,0.0 }, { 0.0,0.0,0.0 }, { 0.0,0.0,0.0 } },
    { { -2.0,2.0,0.0 }, { 2.0,-4.0,0.0 },{ 0.0,2.0,0.0 }, { 0.0,0.0,0.0 } },
    { { -3.0,6.0,-3.0 },{ 3.0,-12.0,9.0 },{ 0.0,6.0,-9.0}, { 0.0,0.0,3.0 } } };

    const double *c = bernstein_derivative_coeff[deg][k];
    return ((c[2]*t)+c[1])*t+c[0];
}

    
Vector SBezier::PointAt(double t) const {
    Vector pt = {};
    double d = 0;

    int i;
    for(i = 0; i <= deg; i++) {
        double B = Bernstein(i, deg, t);
        pt = pt.Plus(ctrl[i].ScaledBy(B*weight[i]));
        d += weight[i]*B;
    }
    pt = pt.ScaledBy(1.0/d);
    return pt;
}

Vector SBezier::TangentAt(double t) const {
    Vector pt = {}, pt_p = {};
    double d = 0, d_p = 0;

    int i;
    for(i = 0; i <= deg; i++) {
        double B  = Bernstein(i, deg, t),
               Bp = BernsteinDerivative(i, deg, t);

        pt = pt.Plus(ctrl[i].ScaledBy(B*weight[i]));
        d += weight[i]*B;

        pt_p = pt_p.Plus(ctrl[i].ScaledBy(Bp*weight[i]));
        d_p += weight[i]*Bp;
    }

    // quotient rule; f(t) = n(t)/d(t), so f' = (n'*d - n*d')/(d^2)
    Vector ret;
    ret = (pt_p.ScaledBy(d)).Minus(pt.ScaledBy(d_p));
    ret = ret.ScaledBy(1.0/(d*d));
    return ret;
}

void SBezier::ClosestPointTo(Vector p, double *t, bool mustConverge) const {
    int i;
    double minDist = VERY_POSITIVE;
    *t = 0;
    double res = (deg <= 2) ? 7.0 : 20.0;
    for(i = 0; i < (int)res; i++) {
        double tryt = (i/res);

        Vector tryp = PointAt(tryt);
        double d = (tryp.Minus(p)).Magnitude();
        if(d < minDist) {
            *t = tryt;
            minDist = d;
        }
    }

    Vector p0;
    for(i = 0; i < (mustConverge ? 15 : 5); i++) {
        p0 = PointAt(*t);
        if(p0.Equals(p, RATPOLY_EPS)) {
            return;
        }

        Vector dp = TangentAt(*t);
        Vector pc = p.ClosestPointOnLine(p0, dp);
        *t += (pc.Minus(p0)).DivProjected(dp);
    }
    if(mustConverge) {
        dbp("didn't converge (closest point on bezier curve)");
    }
}

bool SBezier::PointOnNonparallelCurve(const SBezier *curve, Vector *p) const {
    if(deg + curve->deg == 2) {  // check for parallel lines
        Vector d1 = ctrl[1].Minus(ctrl[0]).WithMagnitude(1.0);
        Vector d2 = curve->ctrl[1].Minus(curve->ctrl[0]).WithMagnitude(1.0);
        // I'm not sure what the angle tollerance should be here.
        if(d1.Cross(d2).Magnitude() < LENGTH_EPS) {
            return false;
        }
    }
    // we will not tail call here because of one more check afterward.
    if(!PointOnThisAndCurve(curve, p)) {
        return false;
    }
    // if it did converge we need to verify that it is not a parallel situation
    // since it will not have converged all the way and would produce bad results.
    double ta, tb;
    this->ClosestPointTo(*p, &ta, /*mustConverge=*/false);
    curve->ClosestPointTo(*p, &tb, /*mustConverge=*/false);
    Vector da = this->TangentAt(ta).WithMagnitude(1),
           db = curve->TangentAt(tb).WithMagnitude(1);
    // for some reason using RATPOLY_EPS here causes leaks in one of the test cases?
    if(da.Cross(db).Magnitude() < LENGTH_EPS) {
        return false;
    }
    return true;    
}

bool SBezier::PointOnThisAndCurve(const SBezier *sbb, Vector *p) const {
    double ta, tb;
    this->ClosestPointTo(*p, &ta, /*mustConverge=*/false);
    sbb ->ClosestPointTo(*p, &tb, /*mustConverge=*/false);

    int i;
    for(i = 0; i < 20; i++) {
        Vector pa = this->PointAt(ta),
               pb = sbb ->PointAt(tb),
               da = this->TangentAt(ta),
               db = sbb ->TangentAt(tb);

        if(pa.Equals(pb, RATPOLY_EPS)) {
            *p = pa;
            return true;
        }

        double tta, ttb;
        Vector::ClosestPointBetweenLines(pa, da, pb, db, &tta, &ttb);
        ta += tta;
        tb += ttb;
    }
    return false;
}

void SBezier::SplitAt(double t, SBezier *bef, SBezier *aft) const {
    Vector4 ct[4];
    int i;
    for(i = 0; i <= deg; i++) {
        ct[i] = Vector4::From(weight[i], ctrl[i]);
    }

    switch(deg) {
        case 1: {
            Vector4 cts = Vector4::Blend(ct[0], ct[1], t);
            *bef = SBezier::From(ct[0], cts);
            *aft = SBezier::From(cts, ct[1]);
            break;
        }
        case 2: {
            Vector4 ct01 = Vector4::Blend(ct[0], ct[1], t),
                    ct12 = Vector4::Blend(ct[1], ct[2], t),
                    cts  = Vector4::Blend(ct01,  ct12,  t);

            *bef = SBezier::From(ct[0], ct01, cts);
            *aft = SBezier::From(cts, ct12, ct[2]);
            break;
        }
        case 3: {
            Vector4 ct01    = Vector4::Blend(ct[0], ct[1], t),
                    ct12    = Vector4::Blend(ct[1], ct[2], t),
                    ct23    = Vector4::Blend(ct[2], ct[3], t),
                    ct01_12 = Vector4::Blend(ct01,  ct12,  t),
                    ct12_23 = Vector4::Blend(ct12,  ct23,  t),
                    cts     = Vector4::Blend(ct01_12, ct12_23, t);

            *bef = SBezier::From(ct[0], ct01, ct01_12, cts);
            *aft = SBezier::From(cts, ct12_23, ct23, ct[3]);
            break;
        }
        default: ssassert(false, "Unexpected degree of spline");
    }
}

void SBezier::MakePwlInto(SEdgeList *sel, double chordTol, double max_dt) const {
    List<Vector> lv = {};
    MakePwlInto(&lv, chordTol, max_dt);
    int i;
    for(i = 1; i < lv.n; i++) {
        sel->AddEdge(lv[i-1], lv[i]);
    }
    lv.Clear();
}
void SBezier::MakePwlInto(List<SCurvePt> *l, double chordTol, double max_dt) const {
    List<Vector> lv = {};
    MakePwlInto(&lv, chordTol, max_dt);
    int i;
    for(i = 0; i < lv.n; i++) {
        SCurvePt scpt;
        scpt.tag    = 0;
        scpt.p      = lv[i];
        scpt.vertex = (i == 0) || (i == (lv.n - 1));
        l->Add(&scpt);
    }
    lv.Clear();
}
void SBezier::MakePwlInto(SContour *sc, double chordTol, double max_dt) const {
    List<Vector> lv = {};
    MakePwlInto(&lv, chordTol, max_dt);
    int i;
    for(i = 0; i < lv.n; i++) {
        sc->AddPoint(lv[i]);
    }
    lv.Clear();
}
//--------------------------------------------------------------------------------------
// all variants of MakePwlInto come here. Split a rational Bezier into Piecewise Linear
// segments that don't deviate from the actual curve by more than the chordTol distance.
// max_dt allows to force curves to be split into spans of no more than a certain
// length based on t-parameter. RemoveShortSegments() may delete points when dt <= 0.1
//--------------------------------------------------------------------------------------
void SBezier::MakePwlInto(List<Vector> *l, double chordTol, double max_dt) const {
    if(EXACT(chordTol == 0)) {
        // Use the default chord tolerance.
        chordTol = SS.ChordTolMm();
    }
    // Never do fewer than three intermediate points for curves; people seem to get
    // unhappy when their circles turn into squares, but maybe less
    // unhappy with octagons. Now 16-gons.
    if (EXACT(max_dt == 0.0)) {
        max_dt = (deg == 1) ? 1.0 : 0.25;
    }
    l->Add(&(ctrl[0]));
    // don't split first degee (lines) unless asked to by the caller via max_dt
    if((deg == 1) && (max_dt >= 1.0)) {
        l->Add(&(ctrl[1]));
    } else {
        MakePwlInitialWorker(l, 0.0, 0.5, chordTol, max_dt);
        MakePwlInitialWorker(l, 0.5, 1.0, chordTol, max_dt);
    }
}
void SBezier::MakePwlWorker(List<Vector> *l, double ta, double tb, double chordTol, double max_dt) const
{
    Vector pa = PointAt(ta);
    Vector pb = PointAt(tb);

    Vector pm = PointAt((ta + tb) / 2.0);
    double d = pm.DistanceToLine(pa, pb.Minus(pa));

    double step = 1.0/SS.GetMaxSegments();
    if(((tb - ta) < step || d < chordTol) && ((tb-ta) <= max_dt) ) {
        // A previous call has already added the beginning of our interval.
        l->Add(&pb);
    } else {
        double tm = (ta + tb) / 2;
        MakePwlWorker(l, ta, tm, chordTol, max_dt);
        MakePwlWorker(l, tm, tb, chordTol, max_dt);
    }
}
void SBezier::MakePwlInitialWorker(List<Vector> *l, double ta, double tb, double chordTol, double max_dt) const
{
    Vector pa = PointAt(ta);
    Vector pb = PointAt(tb);

    double tm1 = ta + (tb - ta) * 0.25;
    double tm2 = ta + (tb - ta) * 0.5;
    double tm3 = ta + (tb - ta) * 0.75;

    Vector pm1 = PointAt(tm1);
    Vector pm2 = PointAt(tm2);
    Vector pm3 = PointAt(tm3);
    Vector dir = pb.Minus(pa);

    double d = max({
                   pm1.DistanceToLine(pa, dir),
                   pm2.DistanceToLine(pa, dir),
                   pm3.DistanceToLine(pa, dir)
                });

    double step = 1.0/SS.GetMaxSegments();
    if( ((tb - ta) < step || d < chordTol) && ((tb-ta) <= max_dt) ) {
        // A previous call has already added the beginning of our interval.
        l->Add(&pb);
    } else {
        double tm = (ta + tb) / 2;
        MakePwlWorker(l, ta, tm, chordTol, max_dt);
        MakePwlWorker(l, tm, tb, chordTol, max_dt);
    }
}

void SBezier::MakeNonrationalCubicInto(SBezierList *bl, double tolerance, int depth) const {
    Vector t0 = TangentAt(0), t1 = TangentAt(1);
    // The curve is correct, and the first derivatives are correct, at the
    // endpoints.
    SBezier bnr = SBezier::From(
                        Start(),
                        Start().Plus(t0.ScaledBy(1.0/3)),
                        Finish().Minus(t1.ScaledBy(1.0/3)),
                        Finish());

    bool closeEnough = true;
    int i;
    for(i = 1; i <= 3; i++) {
        double t = i/4.0;
        Vector p0 = PointAt(t),
               pn = bnr.PointAt(t);
        double d = (p0.Minus(pn)).Magnitude();
        if(d > tolerance) {
            closeEnough = false;
        }
    }

    if(closeEnough || depth > 3) {
        bl->l.Add(this);
    } else {
        SBezier bef, aft;
        SplitAt(0.5, &bef, &aft);
        bef.MakeNonrationalCubicInto(bl, tolerance, depth+1);
        aft.MakeNonrationalCubicInto(bl, tolerance, depth+1);
    }
}


SBezier SBezier::From(Vector4 p0, Vector4 p1) {
    SBezier ret = {};
    ret.deg = 1;
    ret.weight[0] = p0.w;
    ret.ctrl  [0] = p0.PerspectiveProject();
    ret.weight[1] = p1.w;
    ret.ctrl  [1] = p1.PerspectiveProject();
    return ret;
}

SBezier SBezier::From(Vector4 p0, Vector4 p1, Vector4 p2) {
    SBezier ret = {};
    ret.deg = 2;
    ret.weight[0] = p0.w;
    ret.ctrl  [0] = p0.PerspectiveProject();
    ret.weight[1] = p1.w;
    ret.ctrl  [1] = p1.PerspectiveProject();
    ret.weight[2] = p2.w;
    ret.ctrl  [2] = p2.PerspectiveProject();
    return ret;
}

SBezier SBezier::From(Vector4 p0, Vector4 p1, Vector4 p2, Vector4 p3) {
    SBezier ret = {};
    ret.deg = 3;
    ret.weight[0] = p0.w;
    ret.ctrl  [0] = p0.PerspectiveProject();
    ret.weight[1] = p1.w;
    ret.ctrl  [1] = p1.PerspectiveProject();
    ret.weight[2] = p2.w;
    ret.ctrl  [2] = p2.PerspectiveProject();
    ret.weight[3] = p3.w;
    ret.ctrl  [3] = p3.PerspectiveProject();
    return ret;
}

SBezier SBezier::From(Vector p0, Vector p1) {
    return SBezier::From(p0.Project4d(),
                         p1.Project4d());
}

SBezier SBezier::From(Vector p0, Vector p1, Vector p2) {
    return SBezier::From(p0.Project4d(),
                         p1.Project4d(),
                         p2.Project4d());
}

SBezier SBezier::From(Vector p0, Vector p1, Vector p2, Vector p3) {
    return SBezier::From(p0.Project4d(),
                         p1.Project4d(),
                         p2.Project4d(),
                         p3.Project4d());
}

Vector SBezier::Start() const {
    return ctrl[0];
}

Vector SBezier::Finish() const {
    return ctrl[deg];
}

void SBezier::Reverse() {
    int i;
    for(i = 0; i < (deg+1)/2; i++) {
        swap(ctrl[i], ctrl[deg-i]);
        swap(weight[i], weight[deg-i]);
    }
}

void SBezier::ScaleSelfBy(double s) {
    int i;
    for(i = 0; i <= deg; i++) {
        ctrl[i] = ctrl[i].ScaledBy(s);
    }
}

void SBezier::GetBoundingProjd(Vector u, Vector orig,
                               double *umin, double *umax) const
{
    int i;
    for(i = 0; i <= deg; i++) {
        double ut = ((ctrl[i]).Minus(orig)).Dot(u);
        if(ut < *umin) *umin = ut;
        if(ut > *umax) *umax = ut;
    }
}

SBezier SBezier::TransformedBy(Vector t, Quaternion q, double scale) const {
    SBezier ret = *this;
    int i;
    for(i = 0; i <= deg; i++) {
        ret.ctrl[i] = (ret.ctrl[i]).ScaledBy(scale);
        ret.ctrl[i] = (q.Rotate(ret.ctrl[i])).Plus(t);
    }
    return ret;
}

//-----------------------------------------------------------------------------
// Does this curve lie entirely within the specified plane? It does if all
// the control points lie in that plane.
//-----------------------------------------------------------------------------
bool SBezier::IsInPlane(Vector n, double d) const {
    int i;
    for(i = 0; i <= deg; i++) {
        if(fabs((ctrl[i]).Dot(n) - d) > LENGTH_EPS) {
            return false;
        }
    }
    return true;
}

//-----------------------------------------------------------------------------
// Is this Bezier exactly the arc of a circle, projected along the specified
// axis? If yes, return that circle's center and radius.
//-----------------------------------------------------------------------------
bool SBezier::IsCircle(Vector axis, Vector *center, double *r) const {
    if(deg != 2) return false;

    if(ctrl[1].DistanceToLine(ctrl[0], ctrl[2].Minus(ctrl[0])) < LENGTH_EPS) {
        // This is almost a line segment. So it's a circle with very large
        // radius, which is likely to make code that tries to handle circles
        // blow up. So return false.
        return false;
    }

    Vector t0 = (ctrl[0]).Minus(ctrl[1]),
           t2 = (ctrl[2]).Minus(ctrl[1]),
           r0 = axis.Cross(t0),
           r2 = axis.Cross(t2);

    *center = Vector::AtIntersectionOfLines(ctrl[0], (ctrl[0]).Plus(r0),
                                            ctrl[2], (ctrl[2]).Plus(r2),
                                            NULL, NULL, NULL);

    double rd0 = center->Minus(ctrl[0]).Magnitude(),
           rd2 = center->Minus(ctrl[2]).Magnitude();
    if(fabs(rd0 - rd2) > LENGTH_EPS) {
        return false;
    }
    *r = rd0;

    Vector u = r0.WithMagnitude(1),
           v = (axis.Cross(u)).WithMagnitude(1);
    Point2d c2  = center->Project2d(u, v),
            pa2 = (ctrl[0]).Project2d(u, v).Minus(c2),
            pb2 = (ctrl[2]).Project2d(u, v).Minus(c2);

    double thetaa = atan2(pa2.y, pa2.x), // in fact always zero due to csys
           thetab = atan2(pb2.y, pb2.x),
           dtheta = WRAP_NOT_0(thetab - thetaa, 2*PI);
    if(dtheta > PI) {
        // Not possible with a second order Bezier arc; so we must have
        // the points backwards.
        dtheta = 2*PI - dtheta;
    }

    if(fabs(weight[1] - cos(dtheta/2)) > LENGTH_EPS) {
        return false;
    }

    return true;
}

bool SBezier::IsRational() const {
    int i;
    for(i = 0; i <= deg; i++) {
        if(fabs(weight[i] - 1) > LENGTH_EPS) return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
// Apply a perspective transformation to a rational Bezier curve, calculating
// the new weights as required.
//-----------------------------------------------------------------------------
SBezier SBezier::InPerspective(Vector u, Vector v, Vector n,
                               Vector origin, double cameraTan) const
{
    Quaternion q = Quaternion::From(u, v);
    q = q.Inverse();
    // we want Q*(p - o) = Q*p - Q*o
    SBezier ret = this->TransformedBy(q.Rotate(origin).ScaledBy(-1), q, 1.0);
    int i;
    for(i = 0; i <= deg; i++) {
        Vector4 ct = Vector4::From(ret.weight[i], ret.ctrl[i]);
        // so the desired curve, before perspective, is
        //    (x/w, y/w, z/w)
        // and after perspective is
        //    ((x/w)/(1 - (z/w)*cameraTan, ...
        //  = (x/(w - z*cameraTan), ...
        // so we want to let w' = w - z*cameraTan
        ct.w = ct.w - ct.z*cameraTan;

        ret.ctrl[i] = ct.PerspectiveProject();
        ret.weight[i] = ct.w;
    }
    return ret;
}

bool SBezier::Equals(SBezier *b) const {
    // We just test of identical degree and control points, even though two
    // curves could still be coincident (even sharing endpoints).
    if(deg != b->deg) return false;
    int i;
    for(i = 0; i <= deg; i++) {
        if(!(ctrl[i]).Equals(b->ctrl[i])) return false;
        if(fabs(weight[i] - b->weight[i]) > LENGTH_EPS) return false;
    }
    return true;
}

void SBezierList::Clear() {
    l.Clear();
}

void SBezierList::ScaleSelfBy(double s) {
    SBezier *sb;
    for(sb = l.First(); sb; sb = l.NextAfter(sb)) {
        sb->ScaleSelfBy(s);
    }
}

//-----------------------------------------------------------------------------
// If our list contains multiple identical Beziers (in either forward or
// reverse order), then cull them. If both is true, both beziers are removed.
// Otherwise only one of them is removed.
//-----------------------------------------------------------------------------
void SBezierList::CullIdenticalBeziers(bool both) {
    int i, j;

    l.ClearTags();
    for(i = 0; i < l.n; i++) {
        SBezier *bi = &(l[i]), bir;
        bir = *bi;
        bir.Reverse();

        for(j = i + 1; j < l.n; j++) {
            SBezier *bj = &(l[j]);
            if(bj->Equals(bi) ||
               bj->Equals(&bir))
            {
                if (both) bi->tag = 1;
                bj->tag = 1;
            }
        }
    }
    l.RemoveTagged();
}

//-----------------------------------------------------------------------------
// Find all the points where a list of Bezier curves intersects another list
// of Bezier curves. We do this by intersecting their piecewise linearizations,
// and then refining any intersections that we find to lie exactly on the
// curves. So this will screw up on tangencies and stuff, but otherwise should
// be fine.
//-----------------------------------------------------------------------------
void SBezierList::AllIntersectionsWith(SBezierList *sblb, SPointList *spl) const {
    for(const SBezier *sba = l.First(); sba; sba = l.NextAfter(sba)) {
        for(const SBezier *sbb = sblb->l.First(); sbb; sbb = sblb->l.NextAfter(sbb)) {
            sbb->AllIntersectionsWith(sba, spl);
        }
    }
}
void SBezier::AllIntersectionsWith(const SBezier *sbb, SPointList *spl) const {
    SPointList splRaw = {};
    SEdgeList sea, seb;
    sea = {};
    seb = {};
    this->MakePwlInto(&sea);
    sbb ->MakePwlInto(&seb);
    SEdge *se;
    for(se = sea.l.First(); se; se = sea.l.NextAfter(se)) {
        // This isn't quite correct, since AnyEdgeCrossings doesn't count
        // the case where two pairs of line segments intersect at their
        // vertices. So this isn't robust, although that case isn't very
        // likely.
        seb.AnyEdgeCrossings(se->a, se->b, NULL, &splRaw);
    }
    SPoint *sp;
    for(sp = splRaw.l.First(); sp; sp = splRaw.l.NextAfter(sp)) {
        Vector p = sp->p;
        if(PointOnThisAndCurve(sbb, &p)) {
            if(!spl->ContainsPoint(p)) spl->Add(p);
        }
    }
    sea.Clear();
    seb.Clear();
    splRaw.Clear();
}

//-----------------------------------------------------------------------------
// Find a plane that contains all of the curves in this list. If the curves
// are all colinear (or coincident, or empty), then that plane is not exactly
// determined but we choose the additional degree(s) of freedom arbitrarily.
// Returns true if all the curves are coplanar, otherwise false.
//-----------------------------------------------------------------------------
bool SBezierList::GetPlaneContainingBeziers(Vector *p, Vector *u, Vector *v,
                        Vector *notCoplanarAt) const
{
    Vector pt, ptFar, ptOffLine, dp, n;
    double farMax, offLineMax;
    int i;

    // Get any point on any Bezier; or an arbitrary point if list is empty.
    if(!l.IsEmpty()) {
        pt = l[0].Start();
    } else {
        pt = {};
    }
    ptFar = ptOffLine = pt;

    // Get the point farthest from our arbitrary point.
    farMax = VERY_NEGATIVE;
    for(const SBezier *sb = l.First(); sb; sb = l.NextAfter(sb)) {
        for(i = 0; i <= sb->deg; i++) {
            double m = (pt.Minus(sb->ctrl[i])).Magnitude();
            if(m > farMax) {
                ptFar = sb->ctrl[i];
                farMax = m;
            }
        }
    }
    if(ptFar.Equals(pt)) {
        // The points are all coincident. So neither basis vector matters.
        *p = pt;
        *u = {1, 0, 0};
        *v = {0, 1, 0};
        return true;
    }

    // Get the point farthest from the line between pt and ptFar
    dp = ptFar.Minus(pt);
    offLineMax = VERY_NEGATIVE;
    for(const SBezier *sb = l.First(); sb; sb = l.NextAfter(sb)) {
        for(i = 0; i <= sb->deg; i++) {
            double m = (sb->ctrl[i]).DistanceToLine(pt, dp);
            if(m > offLineMax) {
                ptOffLine = sb->ctrl[i];
                offLineMax = m;
            }
        }
    }

    *p = pt;
    if(offLineMax < LENGTH_EPS) {
        // The points are all colinear; so choose the second basis vector
        // arbitrarily.
        *u = (ptFar.Minus(pt)).WithMagnitude(1);
        *v = (u->Normal(0)).WithMagnitude(1);
    } else {
        // The points actually define a plane.
        n = (ptFar.Minus(pt)).Cross(ptOffLine.Minus(pt));
        *u = (n.Normal(0)).WithMagnitude(1);
        *v = (n.Normal(1)).WithMagnitude(1);
    }

    // So we have a plane; but check whether all of the points lie in that
    // plane.
    n = u->Cross(*v);
    n = n.WithMagnitude(1);
    double d = p->Dot(n);
    for(const SBezier *sb = l.First(); sb; sb = l.NextAfter(sb)) {
        for(i = 0; i <= sb->deg; i++) {
            if(fabs(n.Dot(sb->ctrl[i]) - d) > LENGTH_EPS) {
                if(notCoplanarAt) *notCoplanarAt = sb->ctrl[i];
                return false;
            }
        }
    }
    return true;
}

//-----------------------------------------------------------------------------
// Assemble curves in sbl into a single loop. The curves may appear in any
// direction (start to finish, or finish to start), and will be reversed if
// necessary. The curves in the returned loop are removed from sbl, even if
// the loop cannot be closed.
//-----------------------------------------------------------------------------
SBezierLoop SBezierLoop::FromCurves(SBezierList *sbl,
                                    bool *allClosed, SEdge *errorAt)
{
    SBezierLoop loop = {};

    if(sbl->l.n < 1) {
        *allClosed = false;
        return loop;
    }
    sbl->l.ClearTags();

    SBezier *first = &(sbl->l[0]);
    first->tag = 1;
    loop.l.Add(first);
    Vector start = first->Start();
    Vector hanging = first->Finish();
    int auxA = first->auxA;

    sbl->l.RemoveTagged();

    while(!sbl->l.IsEmpty() && !hanging.Equals(start)) {
        int i;
        bool foundNext = false;
        for(i = 0; i < sbl->l.n; i++) {
            SBezier *test = &(sbl->l[i]);

            if((test->Finish()).Equals(hanging) && test->auxA == auxA) {
                test->Reverse();
                // and let the next test catch it
            }
            if((test->Start()).Equals(hanging) && test->auxA == auxA) {
                test->tag = 1;
                loop.l.Add(test);
                hanging = test->Finish();
                sbl->l.RemoveTagged();
                foundNext = true;
                break;
            }
        }
        if(!foundNext) {
            // The loop completed without finding the hanging edge, so
            // it's an open loop
            errorAt->a = hanging;
            errorAt->b = start;
            *allClosed = false;
            return loop;
        }
    }
    if(hanging.Equals(start)) {
        *allClosed = true;
    } else {
        // We ran out of edges without forming a closed loop.
        errorAt->a = hanging;
        errorAt->b = start;
        *allClosed = false;
    }

    return loop;
}

void SBezierLoop::Reverse() {
    l.Reverse();
    SBezier *sb;
    for(sb = l.First(); sb; sb = l.NextAfter(sb)) {
        // If we didn't reverse each curve, then the next curve in list would
        // share your start, not your finish.
        sb->Reverse();
    }
}

void SBezierLoop::GetBoundingProjd(Vector u, Vector orig,
                                   double *umin, double *umax) const
{
    for(const SBezier *sb = l.First(); sb; sb = l.NextAfter(sb)) {
        sb->GetBoundingProjd(u, orig, umin, umax);
    }
}

void SBezierLoop::MakePwlInto(SContour *sc, double chordTol) const {
    for(const SBezier *sb = l.First(); sb; sb = l.NextAfter(sb)) {
        sb->MakePwlInto(sc, chordTol);
        // Avoid double points at join between Beziers; except that
        // first and last points should be identical.
        if(l.NextAfter(sb) != NULL) {
            sc->l.RemoveLast(1);
        }
    }
    // Ensure that it's exactly closed, not just within a numerical tolerance.
    if((sc->l.Last()->p).Equals(sc->l.First()->p)) {
        *sc->l.Last() = *sc->l.First();
    }
}

bool SBezierLoop::IsClosed() const {
    if(l.IsEmpty()) return false;
    Vector s = l.First()->Start(),
           f = l.Last()->Finish();
    return s.Equals(f);
}


//-----------------------------------------------------------------------------
// Assemble the curves in sbl into multiple loops, and piecewise linearize the
// curves into poly. If we can't close a contour, then we add it to
// openContours (if that isn't NULL) and keep going; so this works even if the
// input contains a mix of open and closed curves.
//-----------------------------------------------------------------------------
SBezierLoopSet SBezierLoopSet::From(SBezierList *sbl, SPolygon *poly,
                                    double chordTol,
                                    bool *allClosed, SEdge *errorAt,
                                    SBezierLoopSet *openContours)
{
    SBezierLoopSet ret = {};

    *allClosed = true;
    while(!sbl->l.IsEmpty()) {
        bool thisClosed;
        SBezierLoop loop;
        loop = SBezierLoop::FromCurves(sbl, &thisClosed, errorAt);
        if(!thisClosed) {
            // Record open loops in a separate list, if requested.
            *allClosed = false;
            if(openContours) {
                openContours->l.Add(&loop);
            } else {
                loop.Clear();
            }
        } else {
            ret.l.Add(&loop);
            poly->AddEmptyContour();
            loop.MakePwlInto(poly->l.Last(), chordTol);
        }
    }

    poly->normal = poly->ComputeNormal();
    ret.normal = poly->normal;
    if(poly->l.n > 0) {
        ret.point = poly->AnyPoint();
    } else {
        ret.point = {};
    }

    return ret;
}

void SBezierLoopSet::GetBoundingProjd(Vector u, Vector orig,
                                      double *umin, double *umax) const
{
    for(const SBezierLoop *sbl = l.First(); sbl; sbl = l.NextAfter(sbl)) {
        sbl->GetBoundingProjd(u, orig, umin, umax);
    }
}

double SBezierLoopSet::SignedArea() {
    if(EXACT(area == 0.0)) {
        SPolygon sp = {};
        MakePwlInto(&sp);
        sp.normal = sp.ComputeNormal();
        area = sp.SignedArea();
        sp.Clear();
    }
    return area;
}

//-----------------------------------------------------------------------------
// Convert all the Beziers into piecewise linear form, and assemble that into
// a polygon, one contour per loop.
//-----------------------------------------------------------------------------
void SBezierLoopSet::MakePwlInto(SPolygon *sp) const {
    for(const SBezierLoop *sbl = l.First(); sbl; sbl = l.NextAfter(sbl)) {
        sp->AddEmptyContour();
        sbl->MakePwlInto(sp->l.Last());
    }
}

void SBezierLoopSet::Clear() {
    int i;
    for(i = 0; i < l.n; i++) {
        (l[i]).Clear();
    }
    l.Clear();
}

//-----------------------------------------------------------------------------
// An export helper function. We start with a list of Bezier curves, and
// assemble them into loops. We find the outer loops, and find the outer loops'
// inner loops, and group them accordingly.
//-----------------------------------------------------------------------------
void SBezierLoopSetSet::FindOuterFacesFrom(SBezierList *sbl, SPolygon *spxyz,
                                   SSurface *srfuv,
                                   double chordTol,
                                   bool *allClosed, SEdge *notClosedAt,
                                   bool *allCoplanar, Vector *notCoplanarAt,
                                   SBezierLoopSet *openContours)
{
    SSurface srfPlane;
    if(!srfuv) {
        Vector p, u, v;
        *allCoplanar =
            sbl->GetPlaneContainingBeziers(&p, &u, &v, notCoplanarAt);
        if(!*allCoplanar) {
            // Don't even try to assemble them into loops if they're not
            // all coplanar.
            if(openContours) {
                SBezier *sb;
                for(sb = sbl->l.First(); sb; sb = sbl->l.NextAfter(sb)) {
                    SBezierLoop sbl={};
                    sbl.l.Add(sb);
                    openContours->l.Add(&sbl);
                }
            }
            return;
        }
        // All the curves lie in a plane through p with basis vectors u and v.
        srfPlane = SSurface::FromPlane(p, u, v);
        srfuv = &srfPlane;
    }

    int i, j;
    // Assemble the Bezier trim curves into closed loops; we also get the
    // piecewise linearization of the curves (in the SPolygon spxyz), as a
    // calculation aid for the loop direction.
    SBezierLoopSet sbls = SBezierLoopSet::From(sbl, spxyz, chordTol,
                                               allClosed, notClosedAt,
                                               openContours);
    if(sbls.l.n != spxyz->l.n) return;

    // Convert the xyz piecewise linear to uv piecewise linear.
    SPolygon spuv = {};
    SContour *sc;
    for(sc = spxyz->l.First(); sc; sc = spxyz->l.NextAfter(sc)) {
        spuv.AddEmptyContour();
        SPoint *pt;
        for(pt = sc->l.First(); pt; pt = sc->l.NextAfter(pt)) {
            double u, v;
            srfuv->ClosestPointTo(pt->p, &u, &v);
            spuv.l.Last()->AddPoint({u, v, 0});
        }
    }
    spuv.normal = {0, 0, 1}; // must be, since it's in xy plane now

    static const int OUTER_LOOP = 10;
    static const int INNER_LOOP = 20;
    static const int USED_LOOP  = 30;
    // Fix the contour directions; we do this properly, in uv space, so it
    // works for curved surfaces too (important for STEP export).
    spuv.FixContourDirections();
    for(i = 0; i < spuv.l.n; i++) {
        SContour    *contour = &(spuv.l[i]);
        SBezierLoop *bl = &(sbls.l[i]);
        if(contour->tag) {
            // This contour got reversed in the polygon to make the directions
            // consistent, so the same must be necessary for the Bezier loop.
            bl->Reverse();
        }
        if(contour->IsClockwiseProjdToNormal(spuv.normal)) {
            bl->tag = INNER_LOOP;
        } else {
            bl->tag = OUTER_LOOP;
        }
    }

    bool loopsRemaining = true;
    while(loopsRemaining) {
        loopsRemaining = false;
        for(i = 0; i < sbls.l.n; i++) {
            SBezierLoop *loop = &(sbls.l[i]);
            if(loop->tag != OUTER_LOOP) continue;

            // Check if this contour contains any outer loops; if it does, then
            // we should do those "inner outer loops" first; otherwise we
            // will steal their holes, since their holes also lie inside this
            // contour.
            for(j = 0; j < sbls.l.n; j++) {
                SBezierLoop *outer = &(sbls.l[j]);
                if(i == j) continue;
                if(outer->tag != OUTER_LOOP) continue;

                Vector p = spuv.l[j].AnyEdgeMidpoint();
                if(spuv.l[i].ContainsPointProjdToNormal(spuv.normal, p)) {
                    break;
                }
            }
            if(j < sbls.l.n) {
                // It does, can't do this one yet.
                continue;
            }

            SBezierLoopSet outerAndInners = {};
            loopsRemaining = true;
            loop->tag = USED_LOOP;
            outerAndInners.l.Add(loop);
            int auxA = 0;
            if(loop->l.n > 0) auxA = loop->l[0].auxA;

            for(j = 0; j < sbls.l.n; j++) {
                SBezierLoop *inner = &(sbls.l[j]);
                if(inner->tag != INNER_LOOP) continue;
                if(inner->l.n < 1) continue;
                if(inner->l[0].auxA != auxA) continue;

                Vector p = spuv.l[j].AnyEdgeMidpoint();
                if(spuv.l[i].ContainsPointProjdToNormal(spuv.normal, p)) {
                    outerAndInners.l.Add(inner);
                    inner->tag = USED_LOOP;
                }
            }

            outerAndInners.point  = srfuv->PointAt(0, 0);
            outerAndInners.normal = srfuv->NormalAt(0, 0);
            l.Add(&outerAndInners);
        }
    }

    // If we have poorly-formed loops--for example, overlapping zero-area
    // stuff--then we can end up with leftovers. We use this function to
    // group stuff into closed paths for export when possible, so it's bad
    // to screw up on that stuff. So just add them onto the open curve list.
    // Very ugly, but better than losing curves.
    for(i = 0; i < sbls.l.n; i++) {
        SBezierLoop *loop = &(sbls.l[i]);
        if(loop->tag == USED_LOOP) continue;

        if(openContours) {
            openContours->l.Add(loop);
        } else {
            loop->Clear();
        }
        // but don't free the used loops, since we shallow-copied them to
        // ourself
    }

    sbls.l.Clear(); // not sbls.Clear(), since that would deep-clear
    spuv.Clear();
}

void SBezierLoopSetSet::AddOpenPath(SBezier *sb) {
    SBezierLoop sbl = {};
    sbl.l.Add(sb);

    SBezierLoopSet sbls = {};
    sbls.l.Add(&sbl);

    l.Add(&sbls);
}

void SBezierLoopSetSet::Clear() {
    SBezierLoopSet *sbls;
    for(sbls = l.First(); sbls; sbls = l.NextAfter(sbls)) {
        sbls->Clear();
    }
    l.Clear();
}


} // namespace SolveSpace
