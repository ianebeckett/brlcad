#ifndef BVH_V2_TRI_H
#define BVH_V2_TRI_H

#include "vec.h"
#include "ray.h"
#include "bbox.h"

#include <limits>
#include <utility>
#include <optional>

namespace bvh::v2 {

template <typename T, size_t N, bool DO_STORE = true>
struct Tri {
    Vec<T, N, DO_STORE> p0, p1, p2;

    Tri() = default;

    BVH_ALWAYS_INLINE explicit Tri( T * p0, T * p1, T * p2 )
      : p0( p0 ), p1( p1 ), p2( p2) {
      static_assert( not DO_STORE );
    }

    template< bool DS = DO_STORE, class = std::enable_if_t< DS > >
    BVH_ALWAYS_INLINE Tri( Tri<T,N,false> const & other )
      : p0( other.p0 ), p1( other.p1 ), p2( other.p2 ) {}

    BVH_ALWAYS_INLINE Tri(const Vec<T, N, DO_STORE>& p0, const Vec<T, N, DO_STORE>& p1, const Vec<T, N, DO_STORE>& p2)
        : p0(p0), p1(p1), p2(p2)
    {}

    BVH_ALWAYS_INLINE BBox<T, N> get_bbox() const { return BBox<T, N>(p0).extend(p1).extend(p2); }
    BVH_ALWAYS_INLINE Vec<T, N, true> get_center() const { return (p0 + p1 + p2) * static_cast<T>(1. / 3.); }
};

template< typename T, size_t N >
using PtrTri = Tri<T,N,false>;

/// A 3d triangle, represented as two edges and a point, with an (unnormalized, left-handed) normal.
template <typename T>
struct PrecomputedTri {
    Vec<T, 3> p0, e1, e2, n;
    size_t prim_id = -1;

    PrecomputedTri() = default;

    BVH_ALWAYS_INLINE PrecomputedTri(const Vec<T, 3>& p0, const Vec<T, 3>& p1, const Vec<T, 3>& p2)
        : p0(p0), e1(p0 - p1), e2(p2 - p0), n(cross(e1, e2))
    {}

    template< bool B, bool B2>
    BVH_ALWAYS_INLINE PrecomputedTri(const Vec<T, 3, B>& p0, const Vec<T, 3, B>& p1, const Vec<T, 3, B>& p2, const Vec<T,3, B2>& nn)
        : p0(p0), e1(p0 - p1), e2(p2 - p0), n(nn)
    {}

    template< bool B >
    BVH_ALWAYS_INLINE PrecomputedTri(const Tri<T, 3, B>& triangle)
        : PrecomputedTri(triangle.p0, triangle.p1, triangle.p2)
    {}

    BVH_ALWAYS_INLINE Tri<T, 3> convert_to_tri() const { return Tri<T, 3>(p0, p0 - e1, e2 + p0); }
    BVH_ALWAYS_INLINE BBox<T, 3> get_bbox() const { return convert_to_tri().get_bbox(); }
    BVH_ALWAYS_INLINE Vec<T, 3> get_center() const { return convert_to_tri().get_center(); }

    /// Returns a pair containing the barycentric coordinates of the hit point if the given ray
    /// intersects the triangle, otherwise returns nothing. The distance at which the ray intersects
    /// the triangle is set in `ray.tmax`. The tolerance can be adjusted to account for numerical
    /// precision issues.
    template< bool B >
    BVH_ALWAYS_INLINE std::optional<std::pair<T, T>> intersect(
        Ray<T, 3, B>& ray,
        T & isect,
        T tolerance ) const; // = -std::numeric_limits<T>::epsilon()) const;
    // BVH_ALWAYS_INLINE std::optional<std::pair<T, T>> intersect(
    //     Ray<T, 3>& ray,
    //     T tolerance = -std::numeric_limits<T>::epsilon()) const {
    //   return intersect( ray, ray.tmax, tolerance );
    // }
};

template <typename T>
template< bool B >
std::optional<std::pair<T, T>> PrecomputedTri<T>::intersect(Ray<T, 3, B>& ray, T & isect, T tolerance) const {
    auto c = p0 - ray.org;
    auto r = cross(ray.dir, c);
    auto inv_det = static_cast<T>(1.) / dot(n, ray.dir);

    auto u = dot(r, e2) * inv_det;
    auto v = dot(r, e1) * inv_det;
    auto w = static_cast<T>(1.) - u - v;

    // These comparisons are designed to return false
    // when one of t, u, or v is a NaN
    if (u >= tolerance && v >= tolerance && w >= tolerance) {
        auto t = dot(n, c) * inv_det;
        if (t >= ray.tmin && t <= ray.tmax) {
            isect = t;
            return std::make_optional(std::pair<T, T> { u, v });
        }
    }

    return std::nullopt;
}

} // namespace bvh::v2

#endif
