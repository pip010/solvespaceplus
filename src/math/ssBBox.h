#pragma once

#include "ssVector.h"
#include "ssPoint.h"


namespace SolveSpace
{

 class BBox {
public:
    Vector minp;
    Vector maxp;

    static BBox From(const Vector &p0, const Vector &p1);

    Vector GetOrigin() const;
    Vector GetExtents() const;

    void Include(const Vector &v, double r = 0.0);
    bool Overlaps(const BBox &b1) const;
    bool Contains(const Point2d &p, double r = 0.0) const;
};

}//ss