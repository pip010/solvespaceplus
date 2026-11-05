#pragma once

#include "ssExp.h"
#include <ssQuaternion.h>

namespace SolveSpace
{
    class ExprQuaternion {
    public:
        Expr *w, *vx, *vy, *vz;

        static ExprQuaternion From(Expr *w, Expr *vx, Expr *vy, Expr *vz);
        static ExprQuaternion From(Quaternion qn);
        static ExprQuaternion From(hParam w, hParam vx, hParam vy, hParam vz);

        ExprVector RotationU() const;
        ExprVector RotationV() const;
        ExprVector RotationN() const;

        ExprVector Rotate(ExprVector p) const;
        ExprQuaternion Times(ExprQuaternion b) const;

        Expr *Magnitude() const;
    };

} // namespace SolveSpace