#ifndef QUATERNION_HPP
#define QUATERNION_HPP
#include <cstdint>
#include <functional>
#include <math.h>
#include <ostream>
#include <stdio.h>

#include "opflags.h"
#include "vecn.hpp"

namespace math3d {
namespace quaternion {
// holds quaternion related operations

/**Quaternion bases*/
enum QuaternionBase : std::uint8_t {
  kSCALAR = 1, // null value for quaternion base it has
               // no effect on computation
  kI = 2,      // i base for the second quaternion component
  kJ = 3,      // j base for the third quaternion component
  kK = 4       // k base for the fourth quaternion component
};

/**
  \brief Quaternion component
 */
template <class T> struct QuaternionComponent {
  QuaternionBase base;
  T r;

  QuaternionComponent() : r(0), base(kSCALAR) {}
  QuaternionComponent(QuaternionBase b, T a)
      : base(b), r(a) {}
};

template <class T> class Quaternion {
public:
  Quaternion()
      : coeffs{static_cast<T>(0), static_cast<T>(1),
               static_cast<T>(1), static_cast<T>(1)} {}
  //
  Quaternion(T c1, const QuaternionComponent<T> qs[3]) {
    coeffs[0] = c1;
    coeffs[static_cast<std::uint8_t>(qs[0].base) - 1] =
        qs[0].r;
    coeffs[static_cast<std::uint8_t>(qs[1].base) - 1] =
        qs[1].r;
    coeffs[static_cast<std::uint8_t>(qs[2].base) - 1] =
        qs[2].r;
  }
  Quaternion(T c1, const vecn::VecN<T, 3> &vs) {
    coeffs[0] = c1;
    T v_0;
    vs[0, v_0];
    T v_1;
    vs[1, v_1];
    T v_2;
    vs[2, v_2];
    coeffs[1] = v_0;
    coeffs[2] = v_1;
    coeffs[3] = v_2;
  }

  Quaternion(const QuaternionComponent<T> qs[4]) {
    coeffs[static_cast<std::uint8_t>(qs[0].base) - 1] =
        qs[0].r;
    coeffs[static_cast<std::uint8_t>(qs[1].base) - 1] =
        qs[1].r;
    coeffs[static_cast<std::uint8_t>(qs[2].base) - 1] =
        qs[2].r;
    coeffs[static_cast<std::uint8_t>(qs[3].base) - 1] =
        qs[3].r;
  }

  OpResult operator()(T &out) const { return scalar(out); }
  OpResult scalar(T &out) const {

    out = coeffs[0];
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "scalar", SUCCESS);
  }
  OpResult vector(vecn::VecN<T, 3> &out) const {
    T v[3];
    v[0] = coeffs[1];
    v[1] = coeffs[2];
    v[2] = coeffs[3];
    out = vecn::VecN<T, 3>(v);

    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "vector", SUCCESS);
  }
  OpResult operator()(vecn::VecN<T, 3> &out) const {
    return vector(out);
  }
  OpResult operator()(const QuaternionBase &b,
                      QuaternionComponent<T> &c) const {
    switch (b) {
    case kSCALAR: {
      c = QuaternionComponent(kSCALAR, coeffs[0]);
      break;
    }
    case kI: {
      c = QuaternionComponent(kI, coeffs[1]);
      break;
    }
    case kJ: {
      c = QuaternionComponent(kJ, coeffs[2]);
      break;
    }
    case kK: {
      c = QuaternionComponent(kK, coeffs[3]);
      break;
    }
    }

    return OpResult(
        __LINE__, __FILE__, __FUNCTION__,
        "(const QuaternionBase&, QuaternionComponent<T>&)",
        SUCCESS);
  }

  OpResult
  operator()(const QuaternionComponent<T> &c) const {
    QuaternionBase b = c.base;
    coeffs[static_cast<std::uint8_t>(c.base) - 1] = c.r;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "(const QuaternionComponent<T>&)",
                    SUCCESS);
  }
  OpResult multiply(T t, vecn::VecN<T, 3> &out) const {
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.multiply(t, out);
  }
  OpResult add(T t, vecn::VecN<T, 3> &out) const {
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.add(t, out);
  }
  OpResult subtract(T t, vecn::VecN<T, 3> &out) const {
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.subtract(t, out);
  }
  OpResult divide(T t, vecn::VecN<T, 3> &out) const {
    if (t == 0) {

      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "divide", ARG_ERROR);
    }
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.divide(t, out);
  }
  /** arithmetic operations with a vector on vector part*/
  OpResult multiply(vecn::VecN<T, 3> t,
                    vecn::VecN<T, 3> &out) const {
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.multiply(t, out);
  }
  OpResult add(vecn::VecN<T, 3> t,
               vecn::VecN<T, 3> &out) const {

    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.add(t, out);
  }
  OpResult subtract(vecn::VecN<T, 3> t,
                    vecn::VecN<T, 3> &out) const {

    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.subtract(t, out);
  }
  OpResult divide(vecn::VecN<T, 3> t,
                  vecn::VecN<T, 3> &out) const {
    for (std::size_t i = 0; i < 3; i++) {
      T v;
      t(i, v);
      if (v == 0) {
        return OpResult(__LINE__, __FILE__, __FUNCTION__,
                        "divide", ARG_ERROR);
      }
    }
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.divide(t, out);
  }
  /** dot product and cross product for two vec3*/
  OpResult dot(vecn::VecN<T, 3> t, T &out) const {

    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.dot(t, out);
  }
  OpResult cross(vecn::VecN<T, 3> t,
                 vecn::VecN<T, 3> &out) const {
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.cross(t, out);
  }

  /** Quaternion product as it is shown by Vince 2011 p.
   * 63
   Given a quaternion \f[q_a = [s_a, a]\f] and
   another quaternion \f[q_b = [s_b, b]\f]
   Their product is equal to:
   \f[s_a s_b - a \cdot b, s_a b + s_b a + a \times b \f]
   */
  OpResult hamilton_product(const Quaternion &q_b,
                            Quaternion<T> &out) const {
    // s_a, s_b, a, b
    T s_a = static_cast<T>(0);
    auto res = scalar(s_a);
    if (res.status != SUCCESS)
      return res;

    T s_b = static_cast<T>(0);
    res = q_b.scalar(s_b);
    if (res.status != SUCCESS)
      return res;

    vecn::VecN<T, 3> a;
    res = vector(a);
    if (res.status != SUCCESS)
      return res;

    vecn::VecN<T, 3> b;
    res = q_b.vector(b);
    if (res.status != SUCCESS)
      return res;

    // s_a * s_b
    T s_ab = s_a * s_b;

    // a \cdot b
    T a_dot_b = static_cast<T>(0);
    res = a.dot(b, a_dot_b);
    if (res.status != SUCCESS)
      return res;

    // a \times b
    vecn::VecN<T, 3> cross_ab;
    res = a.cross(b, cross_ab);
    if (res.status != SUCCESS)
      return res;

    // s_a * b + s_b * a + a \times b
    vecn::VecN<T, 3> tout;
    for (std::size_t i = 0; i < 3; i++) {
      T b_i;
      b[i, b_i];
      T a_i;
      a[i, a_i];
      T c_i;
      cross_ab[i, c_i];
      tout(i, s_a * b_i + s_b * a_i + c_i);
    }
    out = Quaternion(s_ab - a_dot_b, tout);
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "hamilton_product", SUCCESS);
  }
  OpResult conjugate(Quaternion<T> &out) const {
    T s = static_cast<T>(0);
    auto res = scalar(s);
    if (res.status != SUCCESS)
      return res;

    vecn::VecN<T, 3> vec;
    res = vector(vec);

    if (res.status != SUCCESS)
      return res;

    res = multiply(static_cast<T>(-1), vec);
    if (res.status != SUCCESS)
      return res;

    out = Quaternion(s, vec);
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "conjugate", SUCCESS);
  }
  /**
    \brief from Vince 2011 - Quaternions for Computer
    Graphics p. 69
   */
  OpResult normalized(Quaternion<T> &out) const {
    T nval = static_cast<T>(0);
    auto res = norm(nval);
    if (res.status != SUCCESS)
      return res;
    T inv_mag = static_cast<T>(1.0) / nval;

    res = scalar(nval);

    if (res.status != SUCCESS)
      return res;

    T scalar_part = nval * inv_mag;
    vecn::VecN<T, 3> vs;
    res = vector(vs);
    if (res.status != SUCCESS)
      return res;
    res = multiply(inv_mag, vs);

    if (res.status != SUCCESS)
      return res;
    out = Quaternion(scalar_part, vs);

    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "normalized", SUCCESS);
  }
  /**
    \brief from Vince 2011 - Quaternions for Computer
    Graphics p. 69
   */
  OpResult inversed(Quaternion<T> &out) const {
    //
    T det_out = static_cast<T>(0);
    auto res = det(det_out);
    if (res.status != SUCCESS)
      return res;

    T inv_mag2 = static_cast<T>(1.0) / det_out;
    Quaternion conj;

    res = conjugate(conj);
    if (res.status != SUCCESS)
      return res;

    T conj_scalar = static_cast<T>(0);
    res = conj.scalar(conj_scalar);
    if (res.status != SUCCESS)
      return res;

    T spart = conj_scalar * inv_mag2;

    vecn::VecN<T, 3> vs;
    res = multiply(inv_mag2, vs);
    if (res.status != SUCCESS)
      return res;

    out = Quaternion(spart, vs);
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "inversed", SUCCESS);
  }
  /**
   \brief from Vince 2011 - Quaternions for Computer
   Graphics p. 69
   */
  OpResult add(const Quaternion &q,
               Quaternion<T> &out) const {
    auto fn = [](T thisval, T tval) {
      return thisval + tval;
    };
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "add",
                    apply_el(q, fn, out));
  }
  /**
   \brief from Vince 2011 - Quaternions for Computer
   Graphics p. 69
  */
  OpResult subtract(const Quaternion &q,
                    Quaternion<T> &out) const {
    auto fn = [](T thisval, T tval) {
      return thisval - tval;
    };
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "subtract", apply_el(q, fn, out));
  }
  /**
    \brief from Vince 2011 - Quaternions for Computer
    Graphics p. 69
   */
  OpResult product(const Quaternion &q,
                   Quaternion<T> &out) const {
    return hamilton_product(q, out);
  }
  /**
    \brief from Vince 2011 - Quaternions for Computer
    Graphics p. 69
   */
  OpResult product(T r, Quaternion<T> &out) const {
    auto fn = [](T thisval, T tval) {
      return thisval * tval;
    };
    Quaternion q(r, vecn::VecN<T, 3>(r));
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "product", apply_el(q, fn, out));
  }
  OpResult power(std::size_t i, Quaternion<T> &out) const {
    Quaternion accumulant = *this;
    Quaternion result2 = *this;
    for (std::size_t j = 1; j < i; j++) {
      accumulant.hamilton_product(result2, accumulant);
    }
    out = accumulant;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "power", SUCCESS);
  }
  OpResult squared(Quaternion<T> &out) const {
    Quaternion r1 = *this;
    Quaternion r2 = *this;
    auto res = product(r1, out);
    if (res.status != SUCCESS)
      return res;

    res = product(r2, out);
    if (res.status != SUCCESS)
      return res;

    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "squared", SUCCESS);
  }
  /**
    \brief from Vince 2011 - Quaternions for Computer
    Graphics p. 69
   */
  OpResult norm(T &out) const {
    auto res = det(out);

    if (res.status != SUCCESS)
      return res;
    //
    out = sqrt(out);

    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "norm", SUCCESS);
  }

  /**
    \brief from Vince 2011 - Quaternions for Computer
    Graphics p. 25
   */
  OpResult determinant(T &out) const {
    T s = static_cast<T>(0);
    auto res = scalar(s);
    if (res.status != SUCCESS)
      return res;
    vecn::VecN<T, 3> vec;
    res = vector(vec);
    if (res.status != SUCCESS)
      return res;

    T a = s * s;
    T b;
    res = vec.dot(vec, b);
    if (res.status != SUCCESS)
      return res;
    out = a + b;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "determinant", SUCCESS);
  }
  OpResult det(T &out) const { return determinant(out); }
  OpResult magnitude(T &out) const { return norm(out); }

private:
  T coeffs[4];

  template <typename Func>
  opstatus_t apply_el(const Quaternion &q, const Func &fn,
                      Quaternion<T> &out) const {
    T s;
    scalar(s);
    vecn::VecN<T, 3> v;
    vector(v);
    T q_s;
    q.scalar(q_s);
    vecn::VecN<T, 3> q_v;
    q.vector(q_v);
    out(fn(s, q_s));
    vecn::VecN<T, 3> ovec;
    for (std::size_t i = 0; i < 3; ++i) {
      T v_i;
      v[i, v_i];
      T q_i;
      q_v[i, q_i];
      ovec(i, fn(v_i, q_i));
    }
    out(ovec);
    return SUCCESS;
  }
};

} // namespace quaternion

}; // namespace math3d

#endif
