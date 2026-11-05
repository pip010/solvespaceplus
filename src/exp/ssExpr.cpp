#include "ssExp.h"

namespace SolveSpace
{
    //-----------------------------------------------------------------------------
// The symbolic algebra system used to write our constraint equations;
// routines to build expressions in software or from a user-provided string,
// and to compute the partial derivatives that we'll use when write our
// Jacobian matrix.
//
// Copyright 2008-2013 Jonathan Westhues.
//-----------------------------------------------------------------------------


static inline Expr *AllocExpr() {
    return (Expr *)Platform::AllocTemporary(sizeof(Expr));
}

Expr *Expr::From(hParam p) {
    Expr *r = AllocExpr();
    r->op = Op::PARAM;
    r->parh = p;
    return r;
}

Expr *Expr::From(double v) {
    // Statically allocate common constants.
    // Note: this is only valid because AllocExpr() uses AllocTemporary(),
    // and Expr* is never explicitly freed.

    if(v == 0.0) {
        static Expr zero(0.0);
        return &zero;
    }

    if(v == 1.0) {
        static Expr one(1.0);
        return &one;
    }

    if(v == -1.0) {
        static Expr mone(-1.0);
        return &mone;
    }

    if(v == 0.5) {
        static Expr half(0.5);
        return &half;
    }

    if(v == -0.5) {
        static Expr mhalf(-0.5);
        return &mhalf;
    }

    Expr *r = AllocExpr();
    r->op = Op::CONSTANT;
    r->v = v;
    return r;
}

Expr *Expr::AnyOp(Op newOp, Expr *b) {
    Expr *r = AllocExpr();
    r->op = newOp;
    r->a = this;
    r->b = b;
    return r;
}

int Expr::Children() const {
    switch(op) {
        case Op::PARAM:
        case Op::PARAM_PTR:
        case Op::CONSTANT:
        case Op::VARIABLE:
            return 0;

        case Op::PLUS:
        case Op::MINUS:
        case Op::TIMES:
        case Op::DIV:
            return 2;

        case Op::NEGATE:
        case Op::SQRT:
        case Op::SQUARE:
        case Op::SIN:
        case Op::COS:
        case Op::ASIN:
        case Op::ACOS:
            return 1;
    }
    ssassert(false, "Unexpected operation");
}

int Expr::Nodes() const {
    switch(Children()) {
        case 0: return 1;
        case 1: return 1 + a->Nodes();
        case 2: return 1 + a->Nodes() + b->Nodes();
        default: ssassert(false, "Unexpected children count");
    }
}

Expr *Expr::DeepCopy() const {
    Expr *n = AllocExpr();
    *n = *this;
    int c = n->Children();
    if(c > 0) n->a = a->DeepCopy();
    if(c > 1) n->b = b->DeepCopy();
    return n;
}

Expr *Expr::DeepCopyWithParamsAsPointers(ParamList *firstTry, ParamList *thenTry,
                                         bool foldConstants) const {
    Expr *n = AllocExpr();
    if(op == Op::PARAM) {
        // A param that is referenced by its hParam gets rewritten to go
        // straight in to the parameter table with a pointer, or simply
        // into a constant if it's already known.
        Param *p = firstTry->FindByIdNoOops(parh);
        if(!p) p = thenTry->FindById(parh);
        if(p->known) {
            n->op = Op::CONSTANT;
            n->v = p->val;
        } else {
            n->op = Op::PARAM_PTR;
            n->parp = p;
        }
        return n;
    }

    *n = *this;
    int c = n->Children();
    if(c > 0) {
        n->a = a->DeepCopyWithParamsAsPointers(firstTry, thenTry, foldConstants);
        bool hasConstants = n->a->op == Op::CONSTANT;
        if(c > 1) {
            n->b = b->DeepCopyWithParamsAsPointers(firstTry, thenTry, foldConstants);
            hasConstants |= n->b->op == Op::CONSTANT;
        }
        if(hasConstants && foldConstants) {
            n = n->FoldConstants(false, 0);
        }
    }
    return n;
}

double Expr::Eval() const {
    switch(op) {
        case Op::PARAM:         return SK.GetParam(parh)->val;
        case Op::PARAM_PTR:     return parp->val;

        case Op::CONSTANT:      return v;
        case Op::VARIABLE:      ssassert(false, "Not supported yet");

        case Op::PLUS:          return a->Eval() + b->Eval();
        case Op::MINUS:         return a->Eval() - b->Eval();
        case Op::TIMES:         return a->Eval() * b->Eval();
        case Op::DIV:           return a->Eval() / b->Eval();

        case Op::NEGATE:        return -(a->Eval());
        case Op::SQRT:          return sqrt(a->Eval());
        case Op::SQUARE:        { double r = a->Eval(); return r*r; }
        case Op::SIN:           return sin(a->Eval());
        case Op::COS:           return cos(a->Eval());
        case Op::ACOS:          return acos(a->Eval());
        case Op::ASIN:          return asin(a->Eval());
    }
    ssassert(false, "Unexpected operation");
}

Expr *Expr::PartialWrt(hParam p) const {
    Expr *da, *db;

    switch(op) {
        case Op::PARAM_PTR: return From(p == parp->h ? 1 : 0);
        case Op::PARAM:     return From(p == parh ? 1 : 0);

        case Op::CONSTANT:  return From(0.0);
        case Op::VARIABLE:  ssassert(false, "Not supported yet");

        case Op::PLUS:      return (a->PartialWrt(p))->Plus(b->PartialWrt(p));
        case Op::MINUS:     return (a->PartialWrt(p))->Minus(b->PartialWrt(p));

        case Op::TIMES:
            da = a->PartialWrt(p);
            db = b->PartialWrt(p);
            return (a->Times(db))->Plus(b->Times(da));

        case Op::DIV:
            da = a->PartialWrt(p);
            db = b->PartialWrt(p);
            return ((da->Times(b))->Minus(a->Times(db)))->Div(b->Square());

        case Op::SQRT:
            return (From(0.5)->Div(a->Sqrt()))->Times(a->PartialWrt(p));

        case Op::SQUARE:
            return (From(2.0)->Times(a))->Times(a->PartialWrt(p));

        case Op::NEGATE:    return (a->PartialWrt(p))->Negate();
        case Op::SIN:       return (a->Cos())->Times(a->PartialWrt(p));
        case Op::COS:       return ((a->Sin())->Times(a->PartialWrt(p)))->Negate();

        case Op::ASIN:
            return (From(1)->Div((From(1)->Minus(a->Square()))->Sqrt()))
                        ->Times(a->PartialWrt(p));
        case Op::ACOS:
            return (From(-1)->Div((From(1)->Minus(a->Square()))->Sqrt()))
                        ->Times(a->PartialWrt(p));
    }
    ssassert(false, "Unexpected operation");
}

void Expr::ParamsUsedList(ParamSet *list) const {
    if(op == Op::PARAM || op == Op::PARAM_PTR) {
        // leaf: just add ourselves if we aren't already there
        hParam param = (op == Op::PARAM) ? parh : parp->h;
        list->insert(param);
        return;
    }

    int c = Children();
    if(c >= 1) {
        a->ParamsUsedList(list);
        if(c >= 2) b->ParamsUsedList(list);
    }
}

bool Expr::DependsOn(hParam p) const {
    if(op == Op::PARAM)     return (parh    == p);
    if(op == Op::PARAM_PTR) return (parp->h == p);

    int c = Children();
    if(c == 1)          return a->DependsOn(p);
    if(c == 2)          return a->DependsOn(p) || b->DependsOn(p);
    return false;
}

bool Expr::Tol(double a, double b) {
    return fabs(a - b) < 0.001;
}

bool Expr::IsZeroConst() const {
    return op == Op::CONSTANT && EXACT(v == 0.0);
}

Expr *Expr::FoldConstants(bool allocCopy, size_t depth) {
    Expr *n;
    
    if(!allocCopy) {
        n = this;
    } else {
        n = AllocExpr();
        *n = *this;
    }

    switch(op) {
        case Op::PARAM_PTR:
        case Op::PARAM:
        case Op::CONSTANT:
        case Op::VARIABLE:
            break;

        case Op::MINUS:
        case Op::TIMES:
        case Op::DIV:
        case Op::PLUS:
            if(depth > 0) {
                n->a = a->FoldConstants(allocCopy, depth - 1);
                n->b = b->FoldConstants(allocCopy, depth - 1);
            }
    
            // If both ops are known, then we can evaluate immediately
            if(n->a->op == Op::CONSTANT && n->b->op == Op::CONSTANT) {
                double nv = n->Eval();
                n->op = Op::CONSTANT;
                n->v = nv;
                break;
            }
            // x + 0 = 0 + x = x
            if(op == Op::PLUS && n->b->op == Op::CONSTANT && Tol(n->b->v, 0)) {
                *n = *(n->a); break;
            }
            if(op == Op::PLUS && n->a->op == Op::CONSTANT && Tol(n->a->v, 0)) {
                *n = *(n->b); break;
            }
            // 1*x = x*1 = x
            if(op == Op::TIMES && n->b->op == Op::CONSTANT && Tol(n->b->v, 1)) {
                *n = *(n->a); break;
            }
            if(op == Op::TIMES && n->a->op == Op::CONSTANT && Tol(n->a->v, 1)) {
                *n = *(n->b); break;
            }
            // 0*x = x*0 = 0
            if(op == Op::TIMES && n->b->op == Op::CONSTANT && Tol(n->b->v, 0)) {
                n->op = Op::CONSTANT; n->v = 0; break;
            }
            if(op == Op::TIMES && n->a->op == Op::CONSTANT && Tol(n->a->v, 0)) {
                n->op = Op::CONSTANT; n->v = 0; break;
            }

            break;

        case Op::SQRT:
        case Op::SQUARE:
        case Op::NEGATE:
        case Op::SIN:
        case Op::COS:
        case Op::ASIN:
        case Op::ACOS:
            if(depth > 0) {
                n->a = a->FoldConstants(allocCopy, depth - 1);
            }

            if(n->a->op == Op::CONSTANT) {
                double nv = n->Eval();
                n->op = Op::CONSTANT;
                n->v = nv;
            }
            break;
    }
    return n;
}

void Expr::Substitute(const SubstitutionMap &subMap) {
    ssassert(op != Op::PARAM_PTR, "Expected an expression that refer to params via handles");

    if(op == Op::PARAM) {
        auto it = subMap.find(parh);
        if(it != subMap.end()) {
            parh = it->second->h;
        }
    } else {
        int c = Children();
        if(c >= 1) {
            a->Substitute(subMap);
            if(c >= 2) b->Substitute(subMap);
        }
    }
}

//-----------------------------------------------------------------------------
// If the expression references only one parameter that appears in pl, then
// return that parameter. If no param is referenced, then return NO_PARAMS.
// If multiple params are referenced, then return MULTIPLE_PARAMS.
//-----------------------------------------------------------------------------
const hParam Expr::NO_PARAMS       = { 0 };
const hParam Expr::MULTIPLE_PARAMS = { std::numeric_limits<decltype(hParam::v)>::max() };
hParam Expr::ReferencedParams(ParamList *pl) const {
    if(op == Op::PARAM) {
        if(pl->FindByIdNoOops(parh)) {
            return parh;
        } else {
            return NO_PARAMS;
        }
    }
    ssassert(op != Op::PARAM_PTR, "Expected an expression that refer to params via handles");

    int c = Children();
    if(c == 0) {
        return NO_PARAMS;
    } else if(c == 1) {
        return a->ReferencedParams(pl);
    } else if(c == 2) {
        hParam pa, pb;
        pa = a->ReferencedParams(pl);
        pb = b->ReferencedParams(pl);
        if(pa == NO_PARAMS) {
            return pb;
        } else if(pb == NO_PARAMS) {
            return pa;
        } else if(pa == pb) {
            return pa; // either, doesn't matter
        } else {
            return MULTIPLE_PARAMS;
        }
    } else ssassert(false, "Unexpected children count");
}


//-----------------------------------------------------------------------------
// Routines to pretty-print an expression. Mostly for debugging.
//-----------------------------------------------------------------------------

std::string Expr::Print() const {
    char c;
    switch(op) {
        case Op::PARAM:     return ssprintf("param(%08x)", parh.v);
        case Op::PARAM_PTR: return ssprintf("param(p%08x)", parp->h.v);

        case Op::CONSTANT:  return ssprintf("%.3f", v);
        case Op::VARIABLE:  return "(var)";

        case Op::PLUS:      c = '+'; goto p;
        case Op::MINUS:     c = '-'; goto p;
        case Op::TIMES:     c = '*'; goto p;
        case Op::DIV:       c = '/'; goto p;
p:
            return "(" + a->Print() + " " + c + " " + b->Print() + ")";
            break;

        case Op::NEGATE:    return "(- " + a->Print() + ")";
        case Op::SQRT:      return "(sqrt " + a->Print() + ")";
        case Op::SQUARE:    return "(square " + a->Print() + ")";
        case Op::SIN:       return "(sin " + a->Print() + ")";
        case Op::COS:       return "(cos " + a->Print() + ")";
        case Op::ASIN:      return "(asin " + a->Print() + ")";
        case Op::ACOS:      return "(acos " + a->Print() + ")";
    }
    ssassert(false, "Unexpected operation");
}

} // namespace SolveSpace
