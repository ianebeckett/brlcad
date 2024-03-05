#ifndef BVH_V2_VEC_H
#define BVH_V2_VEC_H

#include "utils.h"

#include <cstddef>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <type_traits>

namespace bvh::v2 {

template <typename T, size_t N, bool DO_STORE = true>
struct Vec {
    using value_type = std::conditional_t< DO_STORE, T[N], T* >;

    // T values[N];
    value_type values;

    Vec() = default;
    BVH_ALWAYS_INLINE explicit Vec( T * ptr ) {
      static_assert( not DO_STORE, "Can only set vec to pointer if not storing internally!" );
      values = ptr;
    }
    template <typename... Args>
    BVH_ALWAYS_INLINE Vec(T x, T y, Args&&... args) : values { x, y, static_cast<T>(std::forward<Args>(args))... } {
      static_assert( DO_STORE, "Cannot fill null pointer" );
    }
    BVH_ALWAYS_INLINE explicit Vec(T x) {
      static_assert( DO_STORE, "Cannot fill null pointer" );
      std::fill(values, values + N, x);
    }

    template <typename Compare>
    BVH_ALWAYS_INLINE size_t get_best_axis(Compare&& compare) const {
        size_t axis = 0;
        static_for<1, N>([&] (size_t i) {
            if (compare(values[i], values[axis]))
                axis = i;
        });
        return axis;
    }

    // Note: These functions are designed to be robust to NaNs
    BVH_ALWAYS_INLINE size_t get_largest_axis() const { return get_best_axis(std::greater<T>()); }
    BVH_ALWAYS_INLINE size_t get_smallest_axis() const { return get_best_axis(std::less<T>()); }

    BVH_ALWAYS_INLINE T& operator [] (size_t i) { return values[i]; }
    BVH_ALWAYS_INLINE T operator [] (size_t i) const { return values[i]; }

    template <typename F>
    BVH_ALWAYS_INLINE static Vec<T, N> generate(F&& f) {
        Vec<T, N> v;
        static_for<0, N>([&] (size_t i) { v[i] = f(i); });
        return v;
    }

    template< bool DNS = not DO_STORE >
    operator typename std::enable_if_t< DNS, Vec<T, N> >() const {
      Vec<T,N> vec;
      std::memcpy( vec.values, values, sizeof(T)*N );
      return vec;
    }

};

template <typename T, size_t N, bool B1, bool B2>
BVH_ALWAYS_INLINE Vec<T, N> operator + (const Vec<T, N, B1>& a, const Vec<T, N, B2>& b) {
    return Vec<T, N>::generate([&] (size_t i) { return a[i] + b[i]; });
}

template <typename T, size_t N, bool B1, bool B2>
BVH_ALWAYS_INLINE Vec<T, N> operator - (const Vec<T, N, B1>& a, const Vec<T, N, B2>& b) {
    return Vec<T, N>::generate([&] (size_t i) { return a[i] - b[i]; });
}

template <typename T, size_t N, bool B1>
BVH_ALWAYS_INLINE Vec<T, N> operator - (const Vec<T, N, B1>& a) {
    return Vec<T, N>::generate([&] (size_t i) { return -a[i]; });
}

template <typename T, size_t N, bool B1, bool B2>
BVH_ALWAYS_INLINE Vec<T, N> operator * (const Vec<T, N, B1>& a, const Vec<T, N, B2>& b) {
    return Vec<T, N>::generate([&] (size_t i) { return a[i] * b[i]; });
}

template <typename T, size_t N, bool B1, bool B2>
BVH_ALWAYS_INLINE Vec<T, N> operator / (const Vec<T, N, B1>& a, const Vec<T, N, B2>& b) {
    return Vec<T, N>::generate([&] (size_t i) { return a[i] / b[i]; });
}

template <typename T, size_t N, bool B1>
BVH_ALWAYS_INLINE Vec<T, N> operator * (const Vec<T, N, B1>& a, T b) {
    return Vec<T, N>::generate([&] (size_t i) { return a[i] * b; });
}

template <typename T, size_t N, bool B1>
BVH_ALWAYS_INLINE Vec<T, N> operator * (T a, const Vec<T, N, B1>& b) {
    return b * a;
}

template <typename T, size_t N, bool B1>
BVH_ALWAYS_INLINE Vec<T, N> operator / (T a, const Vec<T, N, B1>& b) {
    return Vec<T, N>::generate([&] (size_t i) { return a / b[i]; });
}

template <typename T, size_t N, bool B1, bool B2>
BVH_ALWAYS_INLINE Vec<T, N> robust_min(const Vec<T, N, B1>& a, const Vec<T, N, B2>& b) {
    return Vec<T, N>::generate([&] (size_t i) { return robust_min(a[i], b[i]); });
}

template <typename T, size_t N, bool B1, bool B2>
BVH_ALWAYS_INLINE Vec<T, N> robust_max(const Vec<T, N, B1>& a, const Vec<T, N, B2>& b) {
    return Vec<T, N>::generate([&] (size_t i) { return robust_max(a[i], b[i]); });
}

template <typename T, size_t N, bool B1, bool B2>
BVH_ALWAYS_INLINE T dot(const Vec<T, N, B1>& a, const Vec<T, N, B2>& b) {
    return std::transform_reduce(a.values, a.values + N, b.values, T(0));
}

template <typename T, bool b1, bool b2>
BVH_ALWAYS_INLINE Vec<T, 3> cross(const Vec<T, 3, b1>& a, const Vec<T, 3, b2>& b) {
    return Vec<T, 3>(
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]);
}

template <typename T, size_t N, bool B1, bool B2, bool B3>
BVH_ALWAYS_INLINE Vec<T, N> fast_mul_add(const Vec<T, N, B1>& a, const Vec<T, N, B2>& b, const Vec<T, N, B3>& c) {
    return Vec<T, N>::generate([&] (size_t i) { return fast_mul_add(a[i], b[i], c[i]); });
}

template <typename T, size_t N, bool B1>
BVH_ALWAYS_INLINE Vec<T, N> safe_inverse(const Vec<T, N, B1>& v) {
    return Vec<T, N>::generate([&] (size_t i) { return safe_inverse(v[i]); });
}

template <typename T, size_t N, bool B1, bool B2>
BVH_ALWAYS_INLINE T length(const Vec<T, N, B1>& v) {
    return std::sqrt(dot(v, v));
}

template <typename T, size_t N, bool B1>
BVH_ALWAYS_INLINE Vec<T, N> normalize(const Vec<T, N, B1>& v) {
    return v * (static_cast<T>(1.) / length(v));
}

} // namespace bvh::v2

#endif
