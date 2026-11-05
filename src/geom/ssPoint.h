#pragma once

#include <ssVector.h>
#include <ssList.hpp>

namespace SovlveSpace
{
    enum class EarType : uint32_t {
        UNKNOWN = 0,
        NOT_EAR = 1,
        EAR     = 2
    };

    class SPoint {
    public:
        int     tag;

        EarType ear;

        Vector  p;
        Vector  auxv;
    };

    class SPointList {
    public:
        List<SPoint>    l;

        void Clear();
        bool ContainsPoint(Vector pt) const;
        int IndexForPoint(Vector pt) const;
        void IncrementTagFor(Vector pt);
        void Add(Vector pt);
    };

}//ss
