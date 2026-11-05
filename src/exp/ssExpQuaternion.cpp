#include "ssExpQuaternion.h"

namespace SolveSpace
{
    
ExprQuaternion ExprQuaternion::From(hParam w, hParam vx, hParam vy, hParam vz) {
    ExprQuaternion q;
    q.w  = Expr::From(w);
    q.vx = Expr::From(vx);
    q.vy = Expr::From(vy);
    q.vz = Expr::From(vz);
    return q;
}

ExprQuaternion ExprQuaternion::From(Expr *w, Expr *vx, Expr *vy, Expr *vz)
{
    ExprQuaternion q;
    q.w = w;
    q.vx = vx;
    q.vy = vy;
    q.vz = vz;
    return q;
}

ExprQuaternion ExprQuaternion::From(Quaternion qn) {
    ExprQuaternion qe;
    qe.w = Expr::From(qn.w);
    qe.vx = Expr::From(qn.vx);
    qe.vy = Expr::From(qn.vy);
    qe.vz = Expr::From(qn.vz);
    return qe;
}

ExprVector ExprQuaternion::RotationU() const {
    ExprVector u;
    Expr *two = Expr::From(2);

    u.x = w->Square();
    u.x = (u.x)->Plus(vx->Square());
    u.x = (u.x)->Minus(vy->Square());
    u.x = (u.x)->Minus(vz->Square());

    u.y = two->Times(w->Times(vz));
    u.y = (u.y)->Plus(two->Times(vx->Times(vy)));

    u.z = two->Times(vx->Times(vz));
    u.z = (u.z)->Minus(two->Times(w->Times(vy)));

    return u;
}

ExprVector ExprQuaternion::RotationV() const {
    ExprVector v;
    Expr *two = Expr::From(2);

    v.x = two->Times(vx->Times(vy));
    v.x = (v.x)->Minus(two->Times(w->Times(vz)));

    v.y = w->Square();
    v.y = (v.y)->Minus(vx->Square());
    v.y = (v.y)->Plus(vy->Square());
    v.y = (v.y)->Minus(vz->Square());

    v.z = two->Times(w->Times(vx));
    v.z = (v.z)->Plus(two->Times(vy->Times(vz)));

    return v;
}

ExprVector ExprQuaternion::RotationN() const {
    ExprVector n;
    Expr *two = Expr::From(2);

    n.x =              two->Times( w->Times(vy));
    n.x = (n.x)->Plus (two->Times(vx->Times(vz)));

    n.y =              two->Times(vy->Times(vz));
    n.y = (n.y)->Minus(two->Times( w->Times(vx)));

    n.z =               w->Square();
    n.z = (n.z)->Minus(vx->Square());
    n.z = (n.z)->Minus(vy->Square());
    n.z = (n.z)->Plus (vz->Square());

    return n;
}

ExprVector ExprQuaternion::Rotate(ExprVector p) const {
    // Express the point in the new basis
    return (RotationU().ScaledBy(p.x)).Plus(
            RotationV().ScaledBy(p.y)).Plus(
            RotationN().ScaledBy(p.z));
}

ExprQuaternion ExprQuaternion::Times(ExprQuaternion b) const {
    Expr *sa = w, *sb = b.w;
    ExprVector va = { vx, vy, vz };
    ExprVector vb = { b.vx, b.vy, b.vz };

    ExprQuaternion r;
    r.w = (sa->Times(sb))->Minus(va.Dot(vb));
    ExprVector vr = vb.ScaledBy(sa).Plus(
                    va.ScaledBy(sb).Plus(
                    va.Cross(vb)));
    r.vx = vr.x;
    r.vy = vr.y;
    r.vz = vr.z;
    return r;
}

Expr *ExprQuaternion::Magnitude() const {
    return ((w ->Square())->Plus(
            (vx->Square())->Plus(
            (vy->Square())->Plus(
            (vz->Square())))))->Sqrt();
}

} // namespace SolveSpace
