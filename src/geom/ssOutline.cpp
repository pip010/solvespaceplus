//-----------------------------------------------------------------------------
// Operations on triangle meshes, like our mesh Booleans using the BSP, and
// the stuff to check for watertightness.
//
// Copyright 2008-2013 Jonathan Westhues.
//-----------------------------------------------------------------------------
#include "ssOutline.h"

#include <set>

namespace SolveSpace {


bool SOutline::IsVisible(Vector projDir) const {
    double ldot = nl.Dot(projDir);
    double rdot = nr.Dot(projDir);
    return (ldot > -LENGTH_EPS) == (rdot < LENGTH_EPS) ||
           (rdot > -LENGTH_EPS) == (ldot < LENGTH_EPS);
}

void SOutlineList::Clear() {
    l.Clear();
}

void SOutlineList::AddEdge(Vector a, Vector b, Vector nl, Vector nr, int tag) {
    SOutline so = {};
    so.a   = a;
    so.b   = b;
    so.nl  = nl;
    so.nr  = nr;
    so.tag = tag;
    l.Add(&so);
}

void SOutlineList::ListTaggedInto(SEdgeList *el, int auxA, int auxB) {
    for(const SOutline &so : l) {
        if(so.tag == 0) continue;
        el->AddEdge(so.a, so.b, auxA, auxB);
    }
}

void SOutlineList::MakeFromCopyOf(SOutlineList *sol) {
    for(SOutline *so = sol->l.First(); so; so = sol->l.NextAfter(so)) {
        l.Add(so);
    }
}


} // namespace SolveSpace
