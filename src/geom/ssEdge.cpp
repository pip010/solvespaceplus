#include "ssEdge.h"

namespace SolveSpace
{

    SEdge SEdge::From(Vector a, Vector b) {
        SEdge se = {};
        se.a = a;
        se.b = b;
        return se;
    }

    bool SEdge::EdgeCrosses(Vector ea, Vector eb, Vector *ppi, SPointList *spl) const {
        Vector d = eb.Minus(ea);
        double t_eps = LENGTH_EPS/d.Magnitude();

        double t, tthis;
        bool skew;
        Vector pi;
        bool inOrEdge0, inOrEdge1;

        Vector dthis = b.Minus(a);
        double tthis_eps = LENGTH_EPS/dthis.Magnitude();

        if((ea.Equals(a) && eb.Equals(b)) ||
        (eb.Equals(a) && ea.Equals(b)))
        {
            if(ppi) *ppi = a;
            if(spl) spl->Add(a);
            return true;
        }

        // Can't just test if distance between d and a equals distance between d and b;
        // they could be on opposite sides, since we don't have the sign.
        double m = sqrt(d.Magnitude()*dthis.Magnitude());
        if(sqrt(fabs(d.Dot(dthis))) > (m - LENGTH_EPS)) {
            // The edges are parallel.
            if(fabs(a.DistanceToLine(ea, d)) > LENGTH_EPS) {
                // and not coincident, so can't be intersecting
                return false;
            }
            // The edges are coincident. Make sure that neither endpoint lies
            // on the other
            bool inters = false;
            double t;
            t = a.Minus(ea).DivProjected(d);
            if(t > t_eps && t < (1 - t_eps)) inters = true;
            t = b.Minus(ea).DivProjected(d);
            if(t > t_eps && t < (1 - t_eps)) inters = true;
            t = ea.Minus(a).DivProjected(dthis);
            if(t > tthis_eps && t < (1 - tthis_eps)) inters = true;
            t = eb.Minus(a).DivProjected(dthis);
            if(t > tthis_eps && t < (1 - tthis_eps)) inters = true;

            if(inters) {
                if(ppi) *ppi = a;
                if(spl) spl->Add(a);
                return true;
            } else {
                // So coincident but disjoint, okay.
                return false;
            }
        }

        // Lines are not parallel, so look for an intersection.
        pi = Vector::AtIntersectionOfLines(ea, eb, a, b,
                                        &skew,
                                        &t, &tthis);
        if(skew) return false;

        inOrEdge0 = (t     > -t_eps)     && (t     < (1 + t_eps));
        inOrEdge1 = (tthis > -tthis_eps) && (tthis < (1 + tthis_eps));

        if(inOrEdge0 && inOrEdge1) {
            if(a.Equals(ea) || b.Equals(ea) ||
            a.Equals(eb) || b.Equals(eb))
            {
                // Not an intersection if we share an endpoint with an edge
                return false;
            }
            // But it's an intersection if a vertex of one edge lies on the
            // inside of the other (or if they cross away from either's
            // vertex).
            if(ppi) *ppi = pi;
            if(spl) spl->Add(pi);
            return true;
        }
        return false;
    }

    void SEdgeList::Clear() {
        l.Clear();
    }

    void SEdgeList::AddEdge(Vector a, Vector b, int auxA, int auxB, int tag) {
        SEdge e = {};
        e.tag = tag;
        e.a = a;
        e.b = b;
        e.auxA = auxA;
        e.auxB = auxB;
        l.Add(&e);
    }

    bool SEdgeList::AssembleContour(Vector first, Vector last, SContour *dest,
                                    SEdge *errorAt, bool keepDir, int start) const
    {
        int i;

        dest->AddPoint(first);
        dest->AddPoint(last);

        do {
            for(i = start; i < l.n; i++) {
                /// @todo fix const!
                SEdge *se = const_cast<SEdge*>(&(l[i]));
                if(se->tag) continue;

                if(se->a.Equals(last)) {
                    dest->AddPoint(se->b);
                    last = se->b;
                    se->tag = 1;
                    break;
                }
                // Don't allow backwards edges if keepDir is true.
                if(!keepDir && se->b.Equals(last)) {
                    dest->AddPoint(se->a);
                    last = se->a;
                    se->tag = 1;
                    break;
                }
            }
            if(i >= l.n) {
                // Couldn't assemble a closed contour; mark where.
                if(errorAt) {
                    errorAt->a = first;
                    errorAt->b = last;
                }
                return false;
            }

        } while(!last.Equals(first));

        return true;
    }

    bool SEdgeList::AssemblePolygon(SPolygon *dest, SEdge *errorAt, bool keepDir) const {
        dest->Clear();

        bool allClosed = true;
        Vector first = {};
        Vector last  = {};
        int i;
        for(i = 0; i < l.n; i++) {
            if(!l[i].tag) {
                first = l[i].a;
                last = l[i].b;
                /// @todo fix const!
                const_cast<SEdge*>(&(l[i]))->tag = 1;
                // Create a new empty contour in our polygon, and finish assembling
                // into that contour.
                dest->AddEmptyContour();
                if(!AssembleContour(first, last, dest->l.Last(), errorAt, keepDir, i+1)) {
                    allClosed = false;
                }
                // But continue assembling, even if some of the contours are open
            }
        }
        return allClosed;
    }

    //-----------------------------------------------------------------------------
    // Test if the specified edge crosses any of the edges in our list. Two edges
    // are not considered to cross if they share an endpoint (within LENGTH_EPS),
    // but they are considered to cross if they are coincident and overlapping.
    // If pi is not NULL, then a crossing is returned in that.
    //-----------------------------------------------------------------------------
    int SEdgeList::AnyEdgeCrossings(Vector a, Vector b, Vector *ppi, SPointList *spl) const {
        auto cnt = std::count_if(l.begin(), l.end(),
                                [&](SEdge const &se) { return se.EdgeCrosses(a, b, ppi, spl); });
        return static_cast<int>(cnt);
    }

    //-----------------------------------------------------------------------------
    // Returns true if the intersecting edge list contains an edge that shares
    // an endpoint with one of our edges.
    //-----------------------------------------------------------------------------
    bool SEdgeList::ContainsEdgeFrom(const SEdgeList *sel) const {
        for(const SEdge *se = l.First(); se; se = l.NextAfter(se)) {
            if(sel->ContainsEdge(se)) return true;
        }
        return false;
    }
    bool SEdgeList::ContainsEdge(const SEdge *set) const {
        for(const SEdge *se = l.First(); se; se = l.NextAfter(se)) {
            if((se->a).Equals(set->a) && (se->b).Equals(set->b)) return true;
            if((se->b).Equals(set->a) && (se->a).Equals(set->b)) return true;
        }
        return false;
    }

    //-----------------------------------------------------------------------------
    // Remove unnecessary edges:
    // - if two are anti-parallel then
    //     if both=true, remove both
    //     else remove only one.
    // - if two are parallel then remove one.
    //-----------------------------------------------------------------------------
    void SEdgeList::CullExtraneousEdges(bool both) {
        l.ClearTags();
        for(int i = 0; i < l.n; i++) {
            SEdge *se = &(l[i]);
            for(int j = i + 1; j < l.n; j++) {
                SEdge *set = &(l[j]);
                if((set->a).Equals(se->a) && (set->b).Equals(se->b)) {
                    // Two parallel edges exist; so keep only the first one.
                    set->tag = 1;
                }
                if((set->a).Equals(se->b) && (set->b).Equals(se->a)) {
                    // Two anti-parallel edges exist; if both=true, keep neither,
                    // otherwise keep only one.
                    if (both) se->tag = 1;
                    set->tag = 1;
                }
            }
        }
        l.RemoveTagged();
    }


    //-----------------------------------------------------------------------------
    // We have an edge list that contains only collinear edges, maybe with more
    // splits than necessary. Merge any collinear segments that join.
    //-----------------------------------------------------------------------------
    void SEdgeList::MergeCollinearSegments(Vector a, Vector b) {
        const Vector lineStart = a;
        const Vector lineDirection = b.Minus(a);
        std::sort(l.begin(), l.end(), [&](const SEdge &a, const SEdge &b) {
            double ta = (a.a.Minus(lineStart)).DivProjected(lineDirection);
            double tb = (b.a.Minus(lineStart)).DivProjected(lineDirection);

            return (ta < tb);
        });

        l.ClearTags();
        SEdge *prev = nullptr;
        for(auto &now : l) {
            if(prev != nullptr) {
                if((prev->b).Equals(now.a) && prev->auxA == now.auxA) {
                    // The previous segment joins up to us; so merge it into us.
                    prev->tag = 1;
                    now.a     = prev->a;
                }
            }
            prev = &now;
        }
        l.RemoveTagged();
    }
}//ss