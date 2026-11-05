#pragma once

#include <functional>
#include <type_traits>
#include <cstdint>

// All of the hWhatever handles are a 32-bit ID, that is used to represent
// some data structure in the sketch.

namespace SolveSpace {

/// Trait indicating which types are handle types and should get the associated operators.
/// Specialize for each handle type and inherit from std::true_type.
template<class T>
struct IsHandleOracle : std::false_type {};

// Equality-compare any two instances of a handle type.
template<class T>
static inline typename std::enable_if<IsHandleOracle<T>::value, bool>::type
operator==(T const &lhs, T const &rhs) {
    return lhs.v == rhs.v;
}

// Inequality-compare any two instances of a handle type.
template<class T>
static inline typename std::enable_if<IsHandleOracle<T>::value, bool>::type
operator!=(T const &lhs, T const &rhs) {
    return !(lhs == rhs);
}

// Less-than-compare any two instances of a handle type.
template<class T>
static inline typename std::enable_if<IsHandleOracle<T>::value, bool>::type
operator<(T const &lhs, T const &rhs) {
    return lhs.v < rhs.v;
}

template<class T>
struct HandleHasher {
    static_assert(IsHandleOracle<T>::value, "Not a valid handle type");

    inline size_t operator()(const T &h) const {
        using Hasher = std::hash<decltype(T::v)>;
        return Hasher{}(h.v);
    }
};


enum class SolveResult : uint32_t;

enum class PolyError : uint32_t {
    GOOD              = 0,
    NOT_CLOSED        = 1,
    NOT_COPLANAR      = 2,
    SELF_INTERSECTING = 3,
    ZERO_LEN_EDGE     = 4
};

enum class StipplePattern : uint32_t {
    CONTINUOUS     = 0,
    SHORT_DASH     = 1,
    DASH           = 2,
    LONG_DASH      = 3,
    DASH_DOT       = 4,
    DASH_DOT_DOT   = 5,
    DOT            = 6,
    FREEHAND       = 7,
    ZIGZAG         = 8,

    LAST           = ZIGZAG
};



enum class Command : uint32_t;

// All of the hWhatever handles are a 32-bit ID, that is used to represent
// some data structure in the sketch.
class hEquation {
public:
    uint32_t v;

    inline bool isFromConstraint() const;
    inline hConstraint constraint() const;
};

template<>
struct IsHandleOracle<hEquation> : std::true_type {};

class hGroup {
public:
    // bits 15: 0   -- group index
    uint32_t v;

    inline hEntity entity(int i) const;
    inline hParam param(int i) const;
    inline hEquation equation(int i) const;
};

template<>
struct IsHandleOracle<hGroup> : std::true_type {};

class hParam;
class hRequest {
public:
    // bits 15: 0   -- request index
    uint32_t v;

    inline hEntity entity(int i) const;
    inline hParam param(int i) const;

    inline bool IsFromReferences() const;
};

template<>
struct IsHandleOracle<hRequest> : std::true_type {};

class hParam {
public:
    // bits 15: 0   -- param index
    //      31:16   -- request index
    uint32_t v;

    inline hRequest request() const;
};

template<>
struct IsHandleOracle<hParam> : std::true_type {};

class hEntity {
public:
    // bits 15: 0   -- entity index
    //      31:16   -- request index
    uint32_t v;

    inline bool isFromRequest() const;
    inline hRequest request() const;
    inline hGroup group() const;
    inline hEquation equation(int i) const;
};

template<>
struct IsHandleOracle<hEntity> : std::true_type {};



class hStyle {
public:
    uint32_t v;
};

template<>
struct IsHandleOracle<hStyle> : std::true_type {};

template<class T, class H> class IdList;

class Entity;
using EntityList = IdList<Entity,hEntity>;

struct EntityId {
    uint32_t v;     // entity ID, starting from 0
};

template<>
struct IsHandleOracle<EntityId> : std::true_type {};

struct EntityKey {
    hEntity     input;
    int         copyNumber;
    // (input, copyNumber) gets mapped to hGroup::entity(i)
};
struct EntityKeyHash {
    size_t operator()(const EntityKey &k) const {
        size_t h1 = std::hash<uint32_t>{}(k.input.v),
               h2 = std::hash<uint32_t>{}(k.copyNumber);
        return h1 ^ (h2 << 1);
    }
};
struct EntityKeyEqual {
    bool operator()(const EntityKey &a, const EntityKey &b) const {
        return std::tie(a.input, a.copyNumber) == std::tie(b.input, b.copyNumber);
    }
};
typedef std::unordered_map<EntityKey, EntityId, EntityKeyHash, EntityKeyEqual> EntityMap;



class hConstraint {
public:
    uint32_t v;

    inline hEquation equation(int i) const;
    inline hParam param(int i) const;
};

template<>
struct IsHandleOracle<hConstraint> : std::true_type {};

class hEquation {
public:
    uint32_t v;

    inline bool isFromConstraint() const;
    inline hConstraint constraint() const;
};


inline hEntity hGroup::entity(int i) const
    { hEntity r; r.v = 0x80000000 | (v << 16) | (uint32_t)i; return r; }
inline hParam hGroup::param(int i) const
    { hParam r; r.v = 0x80000000 | (v << 16) | (uint32_t)i; return r; }
inline hEquation hGroup::equation(int i) const
    { hEquation r; r.v = (v << 16) | 0x80000000 | (uint32_t)i; return r; }

inline bool hRequest::IsFromReferences() const {
    if(*this == Request::HREQUEST_REFERENCE_XY) return true;
    if(*this == Request::HREQUEST_REFERENCE_YZ) return true;
    if(*this == Request::HREQUEST_REFERENCE_ZX) return true;
    return false;
}
inline hEntity hRequest::entity(int i) const
    { hEntity r; r.v = (v << 16) | (uint32_t)i; return r; }
inline hParam hRequest::param(int i) const
    { hParam r; r.v = (v << 16) | (uint32_t)i; return r; }

inline bool hEntity::isFromRequest() const
    { if(v & 0x80000000) return false; else return true; }
inline hRequest hEntity::request() const
    { hRequest r; r.v = (v >> 16); return r; }
inline hGroup hEntity::group() const
    { hGroup r; r.v = (v >> 16) & 0x3fff; return r; }
inline hEquation hEntity::equation(int i) const
    { hEquation r; r.v = v | 0x40000000 | (uint32_t)i; return r; }

inline hRequest hParam::request() const
    { hRequest r; r.v = (v >> 16); return r; }


inline hEquation hConstraint::equation(int i) const
    { hEquation r; r.v = (v << 16) | (uint32_t)i; return r; }
inline hParam hConstraint::param(int i) const
    { hParam r; r.v = v | 0x40000000 | (uint32_t)i; return r; }

inline bool hEquation::isFromConstraint() const
    { if(v & 0xc0000000) return false; else return true; }
inline hConstraint hEquation::constraint() const
    { hConstraint r; r.v = (v >> 16); return r; }









class hSSurface {
public:
    uint32_t v;
};

template<>
struct IsHandleOracle<hSSurface> : std::true_type {};

class hSCurve {
public:
    uint32_t v;
};

template<>
struct IsHandleOracle<hSCurve> : std::true_type {};


} // namespace SolveSpace
