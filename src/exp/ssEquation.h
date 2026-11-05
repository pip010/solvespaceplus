#pragma once

#include "ssExpr.h"

#include <ssHandles.h>

namespace SolveSpace
{
    class Equation {
    public:
        int         tag;
        hEquation   h;

        Expr        *e;

        void Clear() {}
};

} // namespace SolveSpace
