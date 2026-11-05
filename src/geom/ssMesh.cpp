#include "ssMesh.h"

namespace SolveSpace
{
    
void SMesh::Clear() {
    l.Clear();
}

void SMesh::AddTriangle(STriMeta meta, Vector n, Vector a, Vector b, Vector c) {
    Vector ab = b.Minus(a), bc = c.Minus(b);
    Vector np = ab.Cross(bc);
    if(np.Magnitude() < 1e-10) {
        // ugh; gl sometimes tessellates to collinear triangles
        return;
    }
    if(np.Dot(n) > 0) {
        AddTriangle(meta, a, b, c);
    } else {
        AddTriangle(meta, c, b, a);
    }
}
void SMesh::AddTriangle(STriMeta meta, Vector a, Vector b, Vector c) {
    STriangle t = {};
    t.meta = meta;
    t.a = a;
    t.b = b;
    t.c = c;
    AddTriangle(&t);
}
void SMesh::AddTriangle(const STriangle *st) {
    l.Add(st);
}

void SMesh::DoBounding(Vector v, Vector *vmax, Vector *vmin) const {
    vmax->x = max(vmax->x, v.x);
    vmax->y = max(vmax->y, v.y);
    vmax->z = max(vmax->z, v.z);

    vmin->x = min(vmin->x, v.x);
    vmin->y = min(vmin->y, v.y);
    vmin->z = min(vmin->z, v.z);
}
void SMesh::GetBounding(Vector *vmax, Vector *vmin) const {
    int i;
    *vmin = {VERY_POSITIVE, VERY_POSITIVE, VERY_POSITIVE};
    *vmax = {VERY_NEGATIVE, VERY_NEGATIVE, VERY_NEGATIVE};
    for(i = 0; i < l.n; i++) {
        const STriangle *st = &(l[i]);
        DoBounding(st->a, vmax, vmin);
        DoBounding(st->b, vmax, vmin);
        DoBounding(st->c, vmax, vmin);
    }
}

//----------------------------------------------------------------------------
// Report the edges of the boundary of the region(s) of our mesh that lie
// within the plane n dot p = d.
//----------------------------------------------------------------------------
void SMesh::MakeEdgesInPlaneInto(SEdgeList *sel, Vector n, double d) {
    SMesh m = {};
    m.MakeFromCopyOf(this);

    // Delete all triangles in the mesh that do not lie in our export plane.
    m.l.ClearTags();
    int i;
    for(i = 0; i < m.l.n; i++) {
        STriangle *tr = &(m.l[i]);

        if((fabs(n.Dot(tr->a) - d) >= LENGTH_EPS) ||
           (fabs(n.Dot(tr->b) - d) >= LENGTH_EPS) ||
           (fabs(n.Dot(tr->c) - d) >= LENGTH_EPS))
        {
            tr->tag  = 1;
        }
    }
    m.l.RemoveTagged();

    // Select the naked edges in our resulting open mesh.
    SKdNode *root = SKdNode::From(&m);
    root->SnapToMesh(&m);
    root->MakeCertainEdgesInto(sel, EdgeKind::NAKED_OR_SELF_INTER,
                               /*coplanarIsInter=*/false, NULL, NULL);

    m.Clear();
}

void SMesh::MakeOutlinesInto(SOutlineList *sol, EdgeKind edgeKind) {
    SKdNode *root = SKdNode::From(this);
    root->MakeOutlinesInto(sol, edgeKind);
}

//-----------------------------------------------------------------------------
// When we are called, all of the triangles from l[start] to the end must
// be coplanar. So we try to find a set of fewer triangles that covers the
// exact same area, in order to reduce the number of triangles in the mesh.
// We use this after a triangle has been split against the BSP.
//
// This is really ugly code; basically it just pastes things together to
// form convex polygons, merging collinear edges when possible, then
// triangulates the convex poly.
//-----------------------------------------------------------------------------
void SMesh::Simplify(int start) {
    int maxTriangles = (l.n - start) + 10;

    STriMeta meta = l[start].meta;

    STriangle *tout = new STriangle[maxTriangles];
    int toutc = 0;

    Vector n = {};
    Vector *conv = new Vector[maxTriangles * 3];
    int convc = 0;

    int start0 = start;

    int i, j;
    for(i = start; i < l.n; i++) {
        STriangle *tr = &(l[i]);
        if(tr->MinAltitude() < LENGTH_EPS) {
            tr->tag = 1;
        } else {
            tr->tag = 0;
        }
    }

    for(;;) {
        bool didAdd;
        convc = 0;
        for(i = start; i < l.n; i++) {
            STriangle *tr = &(l[i]);
            if(tr->tag) continue;

            tr->tag = 1;
            n = (tr->Normal()).WithMagnitude(1);
            conv[convc++] = tr->a;
            conv[convc++] = tr->b;
            conv[convc++] = tr->c;

            start = i+1;
            break;
        }
        if(i >= l.n) break;

        do {
            didAdd = false;

            for(j = 0; j < convc; j++) {
                Vector a = conv[WRAP((j-1), convc)],
                       b = conv[j],
                       d = conv[WRAP((j+1), convc)],
                       e = conv[WRAP((j+2), convc)];

                Vector c;
                for(i = start; i < l.n; i++) {
                    STriangle *tr = &(l[i]);
                    if(tr->tag) continue;

                    if((tr->a).Equals(d) && (tr->b).Equals(b)) {
                        c = tr->c;
                    } else if((tr->b).Equals(d) && (tr->c).Equals(b)) {
                        c = tr->a;
                    } else if((tr->c).Equals(d) && (tr->a).Equals(b)) {
                        c = tr->b;
                    } else {
                        continue;
                    }
                    // The vertex at C must be convex; but the others must
                    // be tested
                    Vector ab = b.Minus(a);
                    Vector bc = c.Minus(b);
                    Vector cd = d.Minus(c);
                    Vector de = e.Minus(d);

                    double bDot = (ab.Cross(bc)).Dot(n);
                    double dDot = (cd.Cross(de)).Dot(n);

                    bDot /= min(ab.Magnitude(), bc.Magnitude());
                    dDot /= min(cd.Magnitude(), de.Magnitude());

                    if(fabs(bDot) < LENGTH_EPS && fabs(dDot) < LENGTH_EPS) {
                        conv[WRAP((j+1), convc)] = c;
                        // and remove the vertex at j, which is a dup
                        std::move(conv+j+1, conv+convc, conv+j);
                        convc--;
                    } else if(fabs(bDot) < LENGTH_EPS && dDot > 0) {
                        conv[j] = c;
                    } else if(fabs(dDot) < LENGTH_EPS && bDot > 0) {
                        conv[WRAP((j+1), convc)] = c;
                    } else if(bDot > 0 && dDot > 0) {
                        // conv[j] is unchanged, conv[j+1] goes to [j+2]
                        std::move_backward(conv+j+1, conv+convc, conv+convc+1);
                        conv[j+1] = c;
                        convc++;
                    } else {
                        continue;
                    }

                    didAdd = true;
                    tr->tag = 1;
                    break;
                }
            }
        } while(didAdd);

        // I need to debug why this is required; sometimes the above code
        // still generates a convex polygon
        for(i = 0; i < convc; i++) {
            Vector a = conv[WRAP((i-1), convc)],
                   b = conv[i],
                   c = conv[WRAP((i+1), convc)];
            Vector ab = b.Minus(a);
            Vector bc = c.Minus(b);
            double bDot = (ab.Cross(bc)).Dot(n);
            bDot /= min(ab.Magnitude(), bc.Magnitude());

            if(bDot < 0) return; // XXX, shouldn't happen
        }

        for(i = 0; i < convc - 2; i++) {
            STriangle tr = STriangle::From(meta, conv[0], conv[i+1], conv[i+2]);
            if(tr.MinAltitude() > LENGTH_EPS) {
                tout[toutc++] = tr;
            }
        }
    }

    l.n = start0;
    for(i = 0; i < toutc; i++) {
        AddTriangle(&(tout[i]));
    }
    delete[] tout;
    delete[] conv;
}

void SMesh::AddAgainstBsp(SMesh *srcm, SBsp3 *bsp3) {
    int i;

    for(i = 0; i < srcm->l.n; i++) {
        STriangle *st = &(srcm->l[i]);
        int pn = l.n;
        atLeastOneDiscarded = false;
        SBsp3::InsertOrCreate(bsp3, st, this);
        if(!atLeastOneDiscarded && (l.n != (pn+1))) {
            l.n = pn;
            if(flipNormal) {
                AddTriangle(st->meta, st->c, st->b, st->a);
            } else {
                AddTriangle(st->meta, st->a, st->b, st->c);
            }
        }
        if(l.n - pn > 1) {
            Simplify(pn);
        }
    }
}

void SMesh::MakeFromUnionOf(SMesh *a, SMesh *b) {
    SBsp3 *bspa = SBsp3::FromMesh(a);
    SBsp3 *bspb = SBsp3::FromMesh(b);

    flipNormal = false;
    keepInsideOtherShell = false;

    keepCoplanar = true;
    AddAgainstBsp(b, bspa);

    keepCoplanar = false;
    AddAgainstBsp(a, bspb);
}

void SMesh::MakeFromDifferenceOf(SMesh *a, SMesh *b) {
    SBsp3 *bspa = SBsp3::FromMesh(a);
    SBsp3 *bspb = SBsp3::FromMesh(b);

    flipNormal = true;
    keepCoplanar = true;
    keepInsideOtherShell = true;
    AddAgainstBsp(b, bspa);

    flipNormal = false;
    keepCoplanar = false;
    keepInsideOtherShell = false;
    AddAgainstBsp(a, bspb);
}

void SMesh::MakeFromIntersectionOf(SMesh *a, SMesh *b) {
    SBsp3 *bspa = SBsp3::FromMesh(a);
    SBsp3 *bspb = SBsp3::FromMesh(b);

    keepInsideOtherShell = true;
    flipNormal = false;

    keepCoplanar = false;
    AddAgainstBsp(a, bspb);

    keepCoplanar = true;
    AddAgainstBsp(b, bspa);
}

void SMesh::MakeFromCopyOf(SMesh *a) {
    ssassert(this != a, "Can't make from copy of self");
    for(int i = 0; i < a->l.n; i++) {
        AddTriangle(&(a->l[i]));
    }
}

void SMesh::MakeFromAssemblyOf(SMesh *a, SMesh *b) {
    MakeFromCopyOf(a);
    MakeFromCopyOf(b);
}

void SMesh::MakeFromTransformationOf(SMesh *a, Vector trans,
                                     Quaternion q, double scale)
{
    STriangle *tr;
    for(tr = a->l.First(); tr; tr = a->l.NextAfter(tr)) {
        STriangle tt = *tr;
        tt.a = (tt.a).ScaledBy(scale);
        tt.b = (tt.b).ScaledBy(scale);
        tt.c = (tt.c).ScaledBy(scale);
        if(scale < 0) {
            // The mirroring would otherwise turn a closed mesh inside out.
            swap(tt.a, tt.b);
        }
        tt.a = (q.Rotate(tt.a)).Plus(trans);
        tt.b = (q.Rotate(tt.b)).Plus(trans);
        tt.c = (q.Rotate(tt.c)).Plus(trans);
        AddTriangle(&tt);
    }
}

bool SMesh::IsEmpty() const { return (l.IsEmpty()); }

uint32_t SMesh::FirstIntersectionWith(Point2d mp) const {
    Vector rayPoint = SS.GW.UnProjectPoint3({mp.x, mp.y, 0.0});
    Vector rayDir = SS.GW.UnProjectPoint3({mp.x, mp.y, 1.0}).Minus(rayPoint);

    uint32_t face = 0;
    double faceT = VERY_NEGATIVE;
    for(int i = 0; i < l.n; i++) {
        const STriangle &tr = l[i];
        if(tr.meta.face == 0) continue;

        double t;
        if(!tr.Raytrace(rayPoint, rayDir, &t, NULL)) continue;
        if(t > faceT) {
            face  = tr.meta.face;
            faceT = t;
        }
    }

    return face;
}

Vector SMesh::GetCenterOfMass() const {
    Vector center = {};
    double vol = 0.0;
    for(int i = 0; i < l.n; i++) {
        const STriangle &tr = l[i];
        double tvol = tr.SignedVolume();
        center = center.Plus(tr.a.Plus(tr.b.Plus(tr.c)).ScaledBy(tvol / 4.0));
        vol += tvol;
    }
    return center.ScaledBy(1.0 / vol);
}


void SMesh::PrecomputeTransparency() {
    std::sort(l.begin(), l.end(),
              [&](const STriangle &sta, const STriangle &stb) {
        RgbaColor colora = sta.meta.color,
                  colorb = stb.meta.color;
        bool opaquea = colora.IsEmpty() || colora.alpha == 255,
             opaqueb = colorb.IsEmpty() || colorb.alpha == 255;

        if(!opaquea || !opaqueb) isTransparent = true;
        return (opaquea != opaqueb && opaquea == true);
    });
}

void SMesh::RemoveDegenerateTriangles() {
    for(auto &tr : l) {
        tr.tag = (int)tr.IsDegenerate();
    }
    l.RemoveTagged();
}

double SMesh::CalculateVolume() const {
    double vol = 0;
    for(STriangle tr : l) {
        // Translate to place vertex A at (x, y, 0)
        Vector trans = {tr.a.x, tr.a.y, 0};
        tr.a = (tr.a).Minus(trans);
        tr.b = (tr.b).Minus(trans);
        tr.c = (tr.c).Minus(trans);

        // Rotate to place vertex B on the y-axis. Depending on
        // whether the triangle is CW or CCW, C is either to the
        // right or to the left of the y-axis. This handles the
        // sign of our normal.
        Vector u = {-tr.b.y, tr.b.x, 0};
        u = u.WithMagnitude(1);
        Vector v = {tr.b.x, tr.b.y, 0};
        v = v.WithMagnitude(1);
        Vector n = {0, 0, 1};

        tr.a = (tr.a).DotInToCsys(u, v, n);
        tr.b = (tr.b).DotInToCsys(u, v, n);
        tr.c = (tr.c).DotInToCsys(u, v, n);

        n = tr.Normal().WithMagnitude(1);

        // Triangles on edge don't contribute
        if(fabs(n.z) < LENGTH_EPS) continue;

        // The plane has equation p dot n = a dot n
        double d = (tr.a).Dot(n);
        // nx*x + ny*y + nz*z = d
        // nz*z = d - nx*x - ny*y
        double A = -n.x/n.z, B = -n.y/n.z, C = d/n.z;

        double mac = tr.c.y/tr.c.x, mbc = (tr.c.y - tr.b.y)/tr.c.x;
        double xc = tr.c.x, yb = tr.b.y;

        // I asked Maple for
        //    int(int(A*x + B*y +C, y=mac*x..(mbc*x + yb)), x=0..xc);
        double integral =
            (1.0/3)*(
                A*(mbc-mac)+
                (1.0/2)*B*(mbc*mbc-mac*mac)
            )*(xc*xc*xc)+
            (1.0/2)*(A*yb+B*yb*mbc+C*(mbc-mac))*xc*xc+
            C*yb*xc+
            (1.0/2)*B*yb*yb*xc;

        vol += integral;
    }
    return vol;
}

double SMesh::CalculateSurfaceArea(const std::vector<uint32_t> &faces) const {
    double area = 0.0;
    for(uint32_t f : faces) {
        for(const STriangle &t : l) {
            if(f != t.meta.face) continue;
            area += t.Area();
        }
    }
    return area;
}

} // namespace SolveSpace
