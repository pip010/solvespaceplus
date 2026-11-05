#include "ssExpVector.h"

namespace SolveSpace
{

ExprVector ExprVector::From(Expr *x, Expr *y, Expr *z) {
    ExprVector r = { x, y, z};
    return r;
}

ExprVector ExprVector::From(Vector vn) {
    ExprVector ve;
    ve.x = Expr::From(vn.x);
    ve.y = Expr::From(vn.y);
    ve.z = Expr::From(vn.z);
    return ve;
}

ExprVector ExprVector::From(hParam x, hParam y, hParam z) {
    ExprVector ve;
    ve.x = Expr::From(x);
    ve.y = Expr::From(y);
    ve.z = Expr::From(z);
    return ve;
}

ExprVector ExprVector::From(double x, double y, double z) {
    ExprVector ve;
    ve.x = Expr::From(x);
    ve.y = Expr::From(y);
    ve.z = Expr::From(z);
    return ve;
}

ExprVector ExprVector::Minus(ExprVector b) const {
    ExprVector r;
    r.x = x->Minus(b.x);
    r.y = y->Minus(b.y);
    r.z = z->Minus(b.z);
    return r;
}

ExprVector ExprVector::Plus(ExprVector b) const {
    ExprVector r;
    r.x = x->Plus(b.x);
    r.y = y->Plus(b.y);
    r.z = z->Plus(b.z);
    return r;
}

Expr *ExprVector::Dot(ExprVector b) const {
    Expr *r;
    r =         x->Times(b.x);
    r = r->Plus(y->Times(b.y));
    r = r->Plus(z->Times(b.z));
    return r;
}

ExprVector ExprVector::Cross(ExprVector b) const {
    ExprVector r;
    r.x = (y->Times(b.z))->Minus(z->Times(b.y));
    r.y = (z->Times(b.x))->Minus(x->Times(b.z));
    r.z = (x->Times(b.y))->Minus(y->Times(b.x));
    return r;
}

ExprVector ExprVector::ScaledBy(Expr *s) const {
    ExprVector r;
    r.x = x->Times(s);
    r.y = y->Times(s);
    r.z = z->Times(s);
    return r;
}

ExprVector ExprVector::WithMagnitude(Expr *s) const {
    Expr *m = Magnitude();
    return ScaledBy(s->Div(m));
}

Expr *ExprVector::Magnitude() const {
    Expr *r;
    r =         x->Square();
    r = r->Plus(y->Square());
    r = r->Plus(z->Square());
    return r->Sqrt();
}

Vector ExprVector::Eval() const {
    Vector r;
    r.x = x->Eval();
    r.y = y->Eval();
    r.z = z->Eval();
    return r;
}
    
} // namespace SolveSpace
