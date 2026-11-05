#pragma once

#include "ssExpr.h"
#include <ssVector.h>

namespace SolveSpace
{
    
    class ExprVector {
    public:
        Expr *x, *y, *z;

        static ExprVector From(Expr *x, Expr *y, Expr *z);
        static ExprVector From(Vector vn);
        static ExprVector From(hParam x, hParam y, hParam z);
        static ExprVector From(double x, double y, double z);

        ExprVector Plus(ExprVector b) const;
        ExprVector Minus(ExprVector b) const;
        Expr *Dot(ExprVector b) const;
        ExprVector Cross(ExprVector b) const;
        ExprVector ScaledBy(Expr *s) const;
        ExprVector WithMagnitude(Expr *s) const;
        Expr *Magnitude() const;

        Vector Eval() const;
    };

} // namespace SolveSpace
