#include "ssKdTree.h"

namespace SolveSpace
{
    // A kd-tree element needs to go on a side of a node if it's when KDTREE_EPS
    // of the boundary. So increasing this number never breaks anything, but may
    // result in more duplicated elements. So it's conservative to be sloppy here.
    #define KDTREE_EPS (20*LENGTH_EPS)


    //-----------------------------------------------------------------------------
    // Make a kd-tree of edges. This is used for O(log(n)) implementations of stuff
    // that would naively be O(n).
    //-----------------------------------------------------------------------------
    SKdNodeEdges *SKdNodeEdges::Alloc() {
        SKdNodeEdges *ne = (SKdNodeEdges *)Platform::AllocTemporary(sizeof(SKdNodeEdges));
        *ne = {};
        return ne;
    }
    SEdgeLl *SEdgeLl::Alloc() {
        SEdgeLl *sell = (SEdgeLl *)Platform::AllocTemporary(sizeof(SEdgeLl));
        *sell = {};
        return sell;
    }
    SKdNodeEdges *SKdNodeEdges::From(SEdgeList *sel) {
        SEdgeLl *sell = NULL;
        SEdge *se;
        for(se = sel->l.First(); se; se = sel->l.NextAfter(se)) {
            SEdgeLl *n = SEdgeLl::Alloc();
            n->se = se;
            n->next = sell;
            sell = n;
        }
        return SKdNodeEdges::From(sell);
    }
    SKdNodeEdges *SKdNodeEdges::From(SEdgeLl *sell) {
        SKdNodeEdges *n = SKdNodeEdges::Alloc();

        // Compute the midpoints (just mean, though median would be better) of
        // each component.
        Vector ptAve = {};
        SEdgeLl *flip;
        int totaln = 0;
        for(flip = sell; flip; flip = flip->next) {
            ptAve = ptAve.Plus(flip->se->a);
            ptAve = ptAve.Plus(flip->se->b);
            totaln++;
        }
        ptAve = ptAve.ScaledBy(1.0 / (2*totaln));

        // For each component, see how it splits.
        int ltln[3] = { 0, 0, 0 }, gtln[3] = { 0, 0, 0 };
        double badness[3];
        for(flip = sell; flip; flip = flip->next) {
            for(int i = 0; i < 3; i++) {
                if(flip->se->a.Element(i) < ptAve.Element(i) + KDTREE_EPS ||
                flip->se->b.Element(i) < ptAve.Element(i) + KDTREE_EPS)
                {
                    ltln[i]++;
                }
                if(flip->se->a.Element(i) > ptAve.Element(i) - KDTREE_EPS ||
                flip->se->b.Element(i) > ptAve.Element(i) - KDTREE_EPS)
                {
                    gtln[i]++;
                }
            }
        }
        for(int i = 0; i < 3; i++) {
            badness[i] = pow((double)ltln[i], 4) + pow((double)gtln[i], 4);
        }

        // Choose the least bad coordinate to split along.
        if(badness[0] < badness[1] && badness[0] < badness[2]) {
            n->which = 0;
        } else if(badness[1] < badness[2]) {
            n->which = 1;
        } else {
            n->which = 2;
        }
        n->c = ptAve.Element(n->which);

        if(totaln < 3 || totaln == gtln[n->which] || totaln == ltln[n->which]) {
            n->edges = sell;
            // and we're a leaf node
            return n;
        }

        // Sort the edges according to which side(s) of the split plane they're on.
        SEdgeLl *gtl = NULL, *ltl = NULL;
        for(flip = sell; flip; flip = flip->next) {
            if(flip->se->a.Element(n->which) < n->c + KDTREE_EPS ||
            flip->se->b.Element(n->which) < n->c + KDTREE_EPS)
            {
                SEdgeLl *selln = SEdgeLl::Alloc();
                selln->se = flip->se;
                selln->next = ltl;
                ltl = selln;
            }
            if(flip->se->a.Element(n->which) > n->c - KDTREE_EPS ||
            flip->se->b.Element(n->which) > n->c - KDTREE_EPS)
            {
                SEdgeLl *selln = SEdgeLl::Alloc();
                selln->se = flip->se;
                selln->next = gtl;
                gtl = selln;
            }
        }

        n->lt = SKdNodeEdges::From(ltl);
        n->gt = SKdNodeEdges::From(gtl);
        return n;
    }

    int SKdNodeEdges::AnyEdgeCrossings(Vector a, Vector b, int cnt,
            Vector *pi, SPointList *spl) const
    {
        int inters = 0;
        if(gt && lt) {
            if(a.Element(which) < c + KDTREE_EPS ||
            b.Element(which) < c + KDTREE_EPS)
            {
                inters += lt->AnyEdgeCrossings(a, b, cnt, pi, spl);
            }
            if(a.Element(which) > c - KDTREE_EPS ||
            b.Element(which) > c - KDTREE_EPS)
            {
                inters += gt->AnyEdgeCrossings(a, b, cnt, pi, spl);
            }
        } else {
            SEdgeLl *sell;
            for(sell = edges; sell; sell = sell->next) {
                SEdge *se = sell->se;
                if(se->tag == cnt) continue;
                if(se->EdgeCrosses(a, b, pi, spl)) {
                    inters++;
                }
                se->tag = cnt;
            }
        }
        return inters;
    }

    
SKdNode *SKdNode::Alloc()
    { return (SKdNode *)Platform::AllocTemporary(sizeof(SKdNode)); }

SKdNode *SKdNode::From(SMesh *m) {
    int i;
    STriangle *tra = (STriangle *)Platform::AllocTemporary((m->l.n) * sizeof(*tra));

    for(i = 0; i < m->l.n; i++) {
        tra[i] = m->l[i];
    }

    srand(0);
    int n = m->l.n;
    while(n > 1) {
        int k = rand() % n;
        n--;
        swap(tra[k], tra[n]);
    }

    STriangleLl *tll = NULL;
    for(i = 0; i < m->l.n; i++) {
        STriangleLl *tn = STriangleLl::Alloc();
        tn->tri = &(tra[i]);
        tn->next = tll;
        tll = tn;
    }

    return SKdNode::From(tll);
}

SKdNode *SKdNode::From(STriangleLl *tll) {
    int which = 0;
    SKdNode *ret = Alloc();

    int i;
    int gtc[3] = { 0, 0, 0 }, ltc[3] = { 0, 0, 0 }, allc = 0;
    double badness[3] = { 0, 0, 0 };
    double split[3] = { 0, 0, 0 };

    if(!tll) {
        goto leaf;
    }

    for(i = 0; i < 3; i++) {
        int tcnt = 0;
        STriangleLl *ll;
        for(ll = tll; ll; ll = ll->next) {
            split[i] += (ll->tri->a).Element(i);
            split[i] += (ll->tri->b).Element(i);
            split[i] += (ll->tri->c).Element(i);
            tcnt++;
        }
        split[i] /= (tcnt*3);

        for(ll = tll; ll; ll = ll->next) {
            STriangle *tr = ll->tri;

            double a = (tr->a).Element(i),
                   b = (tr->b).Element(i),
                   c = (tr->c).Element(i);

            if(a < split[i] + KDTREE_EPS ||
               b < split[i] + KDTREE_EPS ||
               c < split[i] + KDTREE_EPS)
            {
                ltc[i]++;
            }
            if(a > split[i] - KDTREE_EPS ||
               b > split[i] - KDTREE_EPS ||
               c > split[i] - KDTREE_EPS)
            {
                gtc[i]++;
            }
            if(i == 0) allc++;
        }
        badness[i] = pow((double)ltc[i], 4) + pow((double)gtc[i], 4);
    }
    if(badness[0] < badness[1] && badness[0] < badness[2]) {
        which = 0;
    } else if(badness[1] < badness[2]) {
        which = 1;
    } else {
        which = 2;
    }

    if(allc < 3 || allc == gtc[which] || allc == ltc[which]) {
        goto leaf;
    }

    STriangleLl *ll;
    STriangleLl *lgt, *llt; lgt = llt = NULL;
    for(ll = tll; ll; ll = ll->next) {
        STriangle *tr = ll->tri;

        double a = (tr->a).Element(which),
               b = (tr->b).Element(which),
               c = (tr->c).Element(which);

        if(a < split[which] + KDTREE_EPS ||
           b < split[which] + KDTREE_EPS ||
           c < split[which] + KDTREE_EPS)
        {
            STriangleLl *n = STriangleLl::Alloc();
            *n = *ll;
            n->next = llt;
            llt = n;
        }
        if(a > split[which] - KDTREE_EPS ||
           b > split[which] - KDTREE_EPS ||
           c > split[which] - KDTREE_EPS)
        {
            STriangleLl *n = STriangleLl::Alloc();
            *n = *ll;
            n->next = lgt;
            lgt = n;
        }
    }

    ret->which = which;
    ret->c = split[which];
    ret->gt = SKdNode::From(lgt);
    ret->lt = SKdNode::From(llt);
    return ret;

leaf:
    ret->tris = tll;
    return ret;
}

void SKdNode::ClearTags() const {
    if(gt && lt) {
        gt->ClearTags();
        lt->ClearTags();
    } else {
        STriangleLl *ll;
        for(ll = tris; ll; ll = ll->next) {
            ll->tri->tag = 0;
        }
    }
}

void SKdNode::AddTriangle(STriangle *tr) {
    if(gt && lt) {
        double ta = (tr->a).Element(which),
               tb = (tr->b).Element(which),
               tc = (tr->c).Element(which);
        if(ta < c + KDTREE_EPS ||
           tb < c + KDTREE_EPS ||
           tc < c + KDTREE_EPS)
        {
            lt->AddTriangle(tr);
        }
        if(ta > c - KDTREE_EPS ||
           tb > c - KDTREE_EPS ||
           tc > c - KDTREE_EPS)
        {
            gt->AddTriangle(tr);
        }
    } else {
        STriangleLl *tn = STriangleLl::Alloc();
        tn->tri = tr;
        tn->next = tris;
        tris = tn;
    }
}

void SKdNode::MakeMeshInto(SMesh *m) const {
    if(gt) gt->MakeMeshInto(m);
    if(lt) lt->MakeMeshInto(m);

    STriangleLl *ll;
    for(ll = tris; ll; ll = ll->next) {
        if(ll->tri->tag) continue;

        m->AddTriangle(ll->tri);
        ll->tri->tag = 1;
    }
}

void SKdNode::ListTrianglesInto(std::vector<STriangle *> *tl) const {
    if(gt) gt->ListTrianglesInto(tl);
    if(lt) lt->ListTrianglesInto(tl);

    STriangleLl *ll;
    for(ll = tris; ll; ll = ll->next) {
        if(ll->tri->tag) continue;

        tl->push_back(ll->tri);
        ll->tri->tag = 1;
    }
}

//-----------------------------------------------------------------------------
// If any triangles in the mesh have an edge that goes through v (but not
// a vertex at v), then split those triangles so that they now have a vertex
// there. The existing triangle is modified, and the new triangle appears
// in extras.
//-----------------------------------------------------------------------------
void SKdNode::SnapToVertex(Vector v, SMesh *extras) {
    if(gt && lt) {
        double vc = v.Element(which);
        if(vc < c + KDTREE_EPS) {
            lt->SnapToVertex(v, extras);
        }
        if(vc > c - KDTREE_EPS) {
            gt->SnapToVertex(v, extras);
        }
        // Nothing bad happens if the triangle to be split appears in both
        // branches; the first call will split the triangle, so that the
        // second call will do nothing, because the modified triangle will
        // already contain v
    } else {
        STriangleLl *ll;
        for(ll = tris; ll; ll = ll->next) {
            STriangle *tr = ll->tri;

            // Do a cheap bbox test first
            int k;
            bool mightHit = true;

            for(k = 0; k < 3; k++) {
                double trA = (tr->a).Element(k);
                double trB = (tr->b).Element(k);
                double trC = (tr->c).Element(k);
                double vk = v.Element(k);
                if(trA < vk - KDTREE_EPS &&
                   trB < vk - KDTREE_EPS &&
                   trC < vk - KDTREE_EPS)
                {
                    mightHit = false;
                    break;
                }
                if(trA > vk + KDTREE_EPS &&
                   trB > vk + KDTREE_EPS &&
                   trC > vk + KDTREE_EPS)
                {
                    mightHit = false;
                    break;
                }
            }
            if(!mightHit) continue;

            if(tr->a.Equals(v)) { tr->a = v; continue; }
            if(tr->b.Equals(v)) { tr->b = v; continue; }
            if(tr->c.Equals(v)) { tr->c = v; continue; }

            if(tr->IsDegenerate()) {
                continue;
            }

            if(v.OnLineSegment(tr->a, tr->b)) {
                STriangle nt = STriangle::From(tr->meta, tr->a, v, tr->c);
                extras->AddTriangle(&nt);
                tr->a = v;
                continue;
            }
            if(v.OnLineSegment(tr->b, tr->c)) {
                STriangle nt = STriangle::From(tr->meta, tr->b, v, tr->a);
                extras->AddTriangle(&nt);
                tr->b = v;
                continue;
            }
            if(v.OnLineSegment(tr->c, tr->a)) {
                STriangle nt = STriangle::From(tr->meta, tr->c, v, tr->b);
                extras->AddTriangle(&nt);
                tr->c = v;
                continue;
            }
        }
    }
}

//-----------------------------------------------------------------------------
// Snap to each vertex of each triangle of the given mesh. If the given mesh
// is identical to the mesh used to make this kd tree, then the result should
// be a vertex-to-vertex mesh.
//-----------------------------------------------------------------------------
void SKdNode::SnapToMesh(SMesh *m) {
    int i, j, k;
    for(i = 0; i < m->l.n; i++) {
        STriangle *tr = &(m->l[i]);
        if(tr->IsDegenerate()) {
            continue;
        }
        for(j = 0; j < 3; j++) {
            Vector v = tr->vertices[j];

            SMesh extra = {};
            SnapToVertex(v, &extra);

            for(k = 0; k < extra.l.n; k++) {
                STriangle *tra = (STriangle *)Platform::AllocTemporary(sizeof(*tra));
                *tra = extra.l[k];
                AddTriangle(tra);
            }
            extra.Clear();
        }
    }
}

//-----------------------------------------------------------------------------
// For all the edges in sel, split them against the given triangle, and test
// them for occlusion. sel is both our input and our output. tag indicates
// whether an edge is occluded.
//-----------------------------------------------------------------------------
void SKdNode::SplitLinesAgainstTriangle(SEdgeList *sel, STriangle *tr) const {
    SEdgeList seln = {};

    Vector tn = tr->Normal().WithMagnitude(1);
    double td = tn.Dot(tr->a);

    // Consider front-facing triangles only.
    if(tn.z > LENGTH_EPS) {
        // If the edge crosses our triangle's plane, then split into above
        // and below parts. Note that we must preserve auxA, which contains
        // the style associated with this line, as well as the tag, which
        // contains the occlusion status.
        SEdge *se;
        for(se = sel->l.First(); se; se = sel->l.NextAfter(se)) {
            double da = (se->a).Dot(tn) - td,
                   db = (se->b).Dot(tn) - td;
            if((da < -LENGTH_EPS && db > LENGTH_EPS) ||
               (db < -LENGTH_EPS && da > LENGTH_EPS))
            {
                Vector m = Vector::AtIntersectionOfPlaneAndLine(
                                        tn, td,
                                        se->a, se->b, NULL);
                seln.AddEdge(m, se->b, se->auxA, 0, se->tag);
                se->b = m;
            }
        }
        for(se = seln.l.First(); se; se = seln.l.NextAfter(se)) {
            sel->AddEdge(se->a, se->b, se->auxA, 0, se->tag);
        }
        seln.Clear();

        for(se = sel->l.First(); se; se = sel->l.NextAfter(se)) {
            Vector pt = ((se->a).Plus(se->b)).ScaledBy(0.5);
            if(pt.Dot(tn) - td > -LENGTH_EPS) {
                // Edge is in front of or on our plane (remember, tn.z > 0)
                // so it is exempt from further splitting
                se->auxB = 1;
            } else {
                // Edge is behind our plane, needs further splitting
                se->auxB = 0;
            }
        }

        // Considering only the (x, y) coordinates, split the edge against our
        // triangle.
        Point2d a = (tr->a).ProjectXy(),
                b = (tr->b).ProjectXy(),
                c = (tr->c).ProjectXy();

        Point2d n[3] = { (b.Minus(a)).Normal().WithMagnitude(1),
                         (c.Minus(b)).Normal().WithMagnitude(1),
                         (a.Minus(c)).Normal().WithMagnitude(1)  };

        double d[3] = { n[0].Dot(b),
                        n[1].Dot(c),
                        n[2].Dot(a)  };

        // Split all of the edges where they intersect the triangle edges
        int i;
        for(i = 0; i < 3; i++) {
            for(se = sel->l.First(); se; se = sel->l.NextAfter(se)) {
                if(se->auxB) continue;

                Point2d ap = (se->a).ProjectXy(),
                        bp = (se->b).ProjectXy();
                double da = n[i].Dot(ap) - d[i],
                       db = n[i].Dot(bp) - d[i];
                if((da < -LENGTH_EPS && db > LENGTH_EPS) ||
                   (db < -LENGTH_EPS && da > LENGTH_EPS))
                {
                    double dab = (db - da);
                    Vector spl = ((se->a).ScaledBy( db/dab)).Plus(
                                  (se->b).ScaledBy(-da/dab));
                    seln.AddEdge(spl, se->b, se->auxA, 0, se->tag);
                    se->b = spl;
                }
            }
            for(se = seln.l.First(); se; se = seln.l.NextAfter(se)) {
                // The split pieces are all behind the triangle, since only
                // edges behind the triangle got split. So their auxB is 0.
                sel->AddEdge(se->a, se->b, se->auxA, 0, se->tag);
            }
            seln.Clear();
        }

        for(se = sel->l.First(); se; se = sel->l.NextAfter(se)) {
            bool occluded;
            if(se->auxB) {
                // Lies above or on the triangle plane, so triangle doesn't
                // occlude it.
                occluded = false;
            } else {
                // Test the segment to see if it lies outside the triangle
                // (i.e., outside wrt at least one edge), and keep it only
                // then.
                Point2d pt = ((se->a).Plus(se->b).ScaledBy(0.5)).ProjectXy();
                occluded = true;
                for(i = 0; i < 3; i++) {
                    // If the test point lies on the boundary of our triangle,
                    // then we still discard the edge.
                    if(n[i].Dot(pt) - d[i] > LENGTH_EPS) occluded = false;
                }
            }

            if(occluded) {
                se->tag = 1;
            }
        }
    }
}

//-----------------------------------------------------------------------------
// Given an edge orig, occlusion test it against our mesh. We output an edge
// list in sel, where only invisible portions of the edge are tagged.
//-----------------------------------------------------------------------------
void SKdNode::OcclusionTestLine(SEdge orig, SEdgeList *sel, int cnt) const {
    if(gt && lt) {
        double ac = (orig.a).Element(which),
               bc = (orig.b).Element(which);
        // We can ignore triangles that are separated in x or y, but triangles
        // that are separated in z may still contribute
        if(ac < c + KDTREE_EPS ||
           bc < c + KDTREE_EPS ||
           which == 2)
        {
            lt->OcclusionTestLine(orig, sel, cnt);
        }
        if(ac > c - KDTREE_EPS ||
           bc > c - KDTREE_EPS ||
           which == 2)
        {
            gt->OcclusionTestLine(orig, sel, cnt);
        }
    } else {
        STriangleLl *ll;
        for(ll = tris; ll; ll = ll->next) {
            STriangle *tr = ll->tri;

            if(tr->tag == cnt) continue;

            SplitLinesAgainstTriangle(sel, tr);
            tr->tag = cnt;
        }
    }
}

//-----------------------------------------------------------------------------
// Search the mesh for a triangle with an edge from b to a (i.e., the mate
// for the edge from a to b), and increment info->count each time that we
// find one. If a triangle is found, then report whether it is front- or
// back-facing using info->frontFacing. And regardless of whether a mate is
// found, report whether the edge intersects the mesh with info->intersectsMesh;
// if coplanarIsInter then we count the edge as intersecting if it's coplanar
// with a triangle in the mesh, otherwise not.
//-----------------------------------------------------------------------------
void SKdNode::FindEdgeOn(Vector a, Vector b, int cnt, bool coplanarIsInter,
                         EdgeOnInfo *info) const
{
    if(gt && lt) {
        double ac = a.Element(which),
               bc = b.Element(which);
        if(ac < c + KDTREE_EPS ||
           bc < c + KDTREE_EPS)
        {
            lt->FindEdgeOn(a, b, cnt, coplanarIsInter, info);
        }
        if(ac > c - KDTREE_EPS ||
           bc > c - KDTREE_EPS)
        {
            gt->FindEdgeOn(a, b, cnt, coplanarIsInter, info);
        }
        return;
    }

    // We are a leaf node; so we iterate over all the triangles in our
    // linked list.
    STriangleLl *ll;
    for(ll = tris; ll; ll = ll->next) {
        STriangle *tr = ll->tri;

        if(tr->tag == cnt) continue;

        // Test if this triangle matches up with the given edge
        if((a.Equals(tr->b) && b.Equals(tr->a)) ||
           (a.Equals(tr->c) && b.Equals(tr->b)) ||
           (a.Equals(tr->a) && b.Equals(tr->c)))
        {
            info->count++;
            // Record whether this triangle is front- or back-facing.
            if(tr->Normal().z > LENGTH_EPS) {
                info->frontFacing = true;
            } else {
                info->frontFacing = false;
            }
            // Record the triangle
            info->tr = tr;
            // And record which vertices a and b correspond to
            info->ai = a.Equals(tr->a) ? 0 : (a.Equals(tr->b) ? 1 : 2);
            info->bi = b.Equals(tr->a) ? 0 : (b.Equals(tr->b) ? 1 : 2);
        } else if(((a.Equals(tr->a) && b.Equals(tr->b)) ||
                   (a.Equals(tr->b) && b.Equals(tr->c)) ||
                   (a.Equals(tr->c) && b.Equals(tr->a))))
        {
            // It's an edge of this triangle, okay.
        } else {
            // Check for self-intersection
            Vector n = (tr->Normal()).WithMagnitude(1);
            double d = (tr->a).Dot(n);
            double pa = a.Dot(n) - d, pb = b.Dot(n) - d;
            // It's an intersection if neither point lies in-plane,
            // and the edge crosses the plane (should handle in-plane
            // intersections separately but don't yet).
            if((pa < -LENGTH_EPS || pa > LENGTH_EPS) &&
               (pb < -LENGTH_EPS || pb > LENGTH_EPS) &&
               (pa*pb < 0))
            {
                // The edge crosses the plane of the triangle; now see if
                // it crosses inside the triangle.
                if(tr->ContainsPointProjd(b.Minus(a), a)) {
                    if(coplanarIsInter) {
                        info->intersectsMesh = true;
                    } else {
                        Vector p = Vector::AtIntersectionOfPlaneAndLine(
                                                n, d, a, b, NULL);
                        Vector ta = tr->a,
                               tb = tr->b,
                               tc = tr->c;
                        if((p.DistanceToLine(ta, tb.Minus(ta)) < LENGTH_EPS) ||
                           (p.DistanceToLine(tb, tc.Minus(tb)) < LENGTH_EPS) ||
                           (p.DistanceToLine(tc, ta.Minus(tc)) < LENGTH_EPS))
                        {
                            // Intersection lies on edge. This happens when
                            // our edge is from a triangle coplanar with
                            // another triangle in the mesh. We don't test
                            // the edge against triangles whose plane contains
                            // that edge, but we do end up testing against
                            // the coplanar triangle's neighbours, which we
                            // will intersect on their edges.
                        } else {
                            info->intersectsMesh = true;
                        }
                    }
                }
            }
        }

        // Ensure that we don't count this triangle twice if it appears
        // in two buckets of the kd tree.
        tr->tag = cnt;
    }
}

static bool CheckAndAddTrianglePair(std::set<std::pair<STriangle *, STriangle *>> *pairs,
                                    STriangle *a, STriangle *b)
{
    if(pairs->find(std::make_pair(a, b)) != pairs->end() ||
       pairs->find(std::make_pair(b, a)) != pairs->end())
        return true;

    pairs->emplace(a, b);
    return false;
}

//-----------------------------------------------------------------------------
// Pick certain classes of edges out from our mesh. These might be:
//    * naked edges (i.e., edges with no anti-parallel neighbor) and self-
//      intersecting edges (i.e., edges that cross another triangle)
//    * turning edges (i.e., edges where a front-facing triangle joins
//      a back-facing triangle)
//    * emphasized edges (i.e., edges where a triangle from one face joins
//      a triangle from a different face)
//-----------------------------------------------------------------------------
void SKdNode::MakeCertainEdgesInto(SEdgeList *sel, EdgeKind how, bool coplanarIsInter,
                                   bool *inter, bool *leaky, int auxA) const
{
    if(inter) *inter = false;
    if(leaky) *leaky = false;

    std::vector<STriangle *> tris;
    ClearTags();
    ListTrianglesInto(&tris);

    std::set<std::pair<STriangle *, STriangle *>> edgeTris;
    int cnt = 1234;
    for(STriangle *tr : tris) {
        for(int j = 0; j < 3; j++) {
            Vector a = tr->vertices[j];
            Vector b = tr->vertices[(j + 1) % 3];

            SKdNode::EdgeOnInfo info = {};
            FindEdgeOn(a, b, cnt, coplanarIsInter, &info);

            switch(how) {
                case EdgeKind::NAKED_OR_SELF_INTER:
                    // there should be one anti-parllel edge
                    if(info.count != 1) {
                        // but there may be multiple parallel coincident edges
                        SKdNode::EdgeOnInfo parallelInfo = {};
                        FindEdgeOn(b, a, -cnt, coplanarIsInter, &parallelInfo);
                        if (info.count != parallelInfo.count) {
                            sel->AddEdge(a, b, auxA);
                            if(leaky) *leaky = true;
                        }
                    }
                    if(info.intersectsMesh) {
                        sel->AddEdge(a, b, auxA);
                        if(inter) *inter = true;
                    }
                    break;

                case EdgeKind::SELF_INTER:
                    if(info.intersectsMesh) {
                        sel->AddEdge(a, b, auxA);
                        if(inter) *inter = true;
                    }
                    break;

                case EdgeKind::TURNING:
                    if((tr->Normal().z < LENGTH_EPS) &&
                       (info.count == 1) &&
                       info.frontFacing)
                    {
                        if(CheckAndAddTrianglePair(&edgeTris, tr, info.tr))
                            break;
                        // This triangle is back-facing (or on edge), and
                        // this edge has exactly one mate, and that mate is
                        // front-facing. So this is a turning edge.
                        sel->AddEdge(a, b, auxA);
                    }
                    break;

                case EdgeKind::EMPHASIZED:
                    if(info.count == 1 && tr->meta.face != info.tr->meta.face) {
                        if(CheckAndAddTrianglePair(&edgeTris, tr, info.tr))
                            break;
                        // The two triangles that join at this edge come from
                        // different faces; either really different faces,
                        // or one is from a face and the other is zero (i.e.,
                        // not from a face).
                        sel->AddEdge(a, b, auxA);
                    }
                    break;

                case EdgeKind::SHARP:
                    if(info.count == 1) {
                        Vector na0 = tr->normals[j].WithMagnitude(1.0);
                        Vector nb0 = tr->normals[(j + 1) % 3].WithMagnitude(1.0);
                        Vector na1 = info.tr->normals[info.ai].WithMagnitude(1.0);
                        Vector nb1 = info.tr->normals[info.bi].WithMagnitude(1.0);
                        if(!((na0.Equals(na1) && nb0.Equals(nb1)) ||
                             (na0.Equals(nb1) && nb0.Equals(na1)))) {
                            if(CheckAndAddTrianglePair(&edgeTris, tr, info.tr))
                                break;
                            // The two triangles that join at this edge meet at a sharp
                            // angle. This implies they come from different faces.
                            sel->AddEdge(a, b, auxA);
                        }
                    }
                    break;
            }

            cnt++;
        }
    }
}

void SKdNode::MakeOutlinesInto(SOutlineList *sol, EdgeKind edgeKind) const
{
    std::vector<STriangle *> tris;
    ClearTags();
    ListTrianglesInto(&tris);

    std::set<std::pair<STriangle *, STriangle *>> edgeTris;
    int cnt = 1234;
    for(STriangle *tr : tris) {
        for(int j = 0; j < 3; j++) {
            Vector a = tr->vertices[j];
            Vector b = tr->vertices[(j + 1) % 3];

            SKdNode::EdgeOnInfo info = {};
            FindEdgeOn(a, b, cnt, /*coplanarIsInter=*/false, &info);
            cnt++;
            if(info.count != 1) continue;
            if(CheckAndAddTrianglePair(&edgeTris, tr, info.tr))
                continue;

            int tag = 0;
            switch(edgeKind) {
                case EdgeKind::EMPHASIZED:
                    if(tr->meta.face != info.tr->meta.face) {
                        tag = 1;
                    }
                    break;

                case EdgeKind::SHARP: {
                        Vector na0 = tr->normals[j].WithMagnitude(1.0);
                        Vector nb0 = tr->normals[(j + 1) % 3].WithMagnitude(1.0);
                        Vector na1 = info.tr->normals[info.ai].WithMagnitude(1.0);
                        Vector nb1 = info.tr->normals[info.bi].WithMagnitude(1.0);
                        if(!((na0.Equals(na1) && nb0.Equals(nb1)) ||
                             (na0.Equals(nb1) && nb0.Equals(na1)))) {
                            tag = 1;
                        }
                    }
                    break;

                default:
                    ssassert(false, "Unexpected edge kind");
            }

            Vector nl = tr->Normal().WithMagnitude(1.0);
            Vector nr = info.tr->Normal().WithMagnitude(1.0);

            // We don't add edges with the same left and right
            // normals because they can't produce outlines.
            if(tag == 0 && nl.Equals(nr)) continue;
            sol->AddEdge(a, b, nl, nr, tag);
        }
    }
}

}//ss