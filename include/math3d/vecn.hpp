#ifndef VECN_HPP
#define VECN_HPP
#include <cstring>
#include <iostream>
#include <math.h>
#include <ostream>
#include <stdio.h>

#include "opflags.h"

namespace math3d {
namespace vecn {

template <typename T> struct VecNCell {
  T content;
  std::size_t index = 0;
};

template <typename T>
VecNCell<T> make_vecn_cell(const T &c, std::size_t i) {
  VecNCell<T> cell{.content = c, .index = i};
  return cell;
}

template <class T, std::size_t N> class VecN {

public:
  /*! Tested */
  VecN() {
    for (std::size_t i = 0; i < N; ++i) {
      data[i] = 0;
    }
  }
  /*! Tested */
  VecN(const T (&arr)[N]) {
    memcpy(data, arr, N * sizeof(T));
  }
  VecN(T s) {
    for (std::size_t i = 0; i < N; ++i) {
      data[i] = s;
    }
  }
  /*! Tested */
  OpResult operator()(std::size_t index, T &out) const {
    if (index >= N) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     "(size_t, T&)", INDEX_ERROR);
      return vflag;
    }
    out = data[index];
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "(size_t, T&)", SUCCESS);
    return vflag;
  }
  OpResult operator()(T (&out)[N]) const {
    memcpy(out, data, N * sizeof(T));
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "(T (&)[N])", SUCCESS);
    return vflag;
  }
  constexpr OpResult operator()(std::size_t &out) const {
    out = N;
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "(size_t)", SUCCESS);
    return vflag;
  }

  /*! Tested */
  OpResult operator()(const VecNCell<T> &cell) {
    if (cell.index >= N) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     "(const VecNCell<T>&)", INDEX_ERROR);
      return vflag;
    }
    data[cell.index] = cell.content;

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "(const VecNCell<T>&)", SUCCESS);
    return vflag;
  }

  /*! Tested */
  OpResult add(T v, T (&out)[N]) const {
    VecN<T, N> v_1(v);
    T v_in[N];
    v_1(v_in);
    //
    auto res = add(v_in, out);
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "add(T, T (&)[N])", res.status);
    return vflag;
  }

  /*! Tested */
  OpResult add(T v, VecN<T, N> &vout) const {
    T out[N];
    auto res = add(v, out);
    vout = VecN<T, N>(out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "add(T, VecN<T, N>&", res.status);
    return vflag;
  }
  /*! Tested */
  OpResult add(const T (&v)[N], T (&out)[N]) const {
    auto fn = [](T thisel, T argel) {
      return thisel + argel;
    };
    auto res = apply_el<decltype(fn)>(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "add(const T(&)[N], T (&)[N])",
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult add(const VecN<T, N> &v, VecN<T, N> &out) const {
    T v_in[N];
    v(v_in);
    T vout[N];
    auto res = add(v_in, vout);
    out = VecN<T, N>(vout);
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "add(const VecN<T, N>&, VecN<T, N>&)",
                   res.status);
    return vflag;
  }
  //
  OpResult subtract(T v, T (&out)[N]) const {
    VecN<T, N> v_1(v);
    T v_in[N];
    v_1(v_in);
    auto res = subtract(v_in, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "subtract(T, T (&)[N])", res.status);
    return vflag;
  }
  /*! Tested */
  OpResult subtract(T v, VecN<T, N> &vout) const {
    VecN<T, N> v_1(v);
    T v_in[N];
    v_1(v_in);
    T out[N];
    auto res = subtract(v_in, out);
    vout = VecN<T, N>(out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "subtract(T, VecN<T, N>&", res.status);
    return vflag;
  }
  /*! Tested */
  OpResult subtract(const T (&v)[N], T (&out)[N]) const {
    auto fn = [](T thisel, T argel) {
      return thisel - argel;
    };
    auto res = apply_el<decltype(fn)>(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "subtract(const T(&)[N], T (&)[N])",
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult subtract(const VecN<T, N> &v,
                    VecN<T, N> &out) const {
    T v_in[N];
    v(v_in);
    T vout[N];
    auto res = subtract(v_in, vout);
    out = VecN<T, N>(vout);

    OpResult vflag(
        __LINE__, __FILE__, __FUNCTION__,
        "subtract(const VecN<T, N>&, VecN<T, N>&)",
        res.status);
    return vflag;
  }
  //
  OpResult multiply(T v, T (&out)[N]) const {
    T v_in[N];
    VecN<T, N> vv(v);
    vv(v_in);

    auto res = multiply(v_in, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "multiply(T, T (&)[N])", res.status);
    return vflag;
  }
  /*! Tested */
  OpResult multiply(T v, VecN<T, N> &vout) const {
    T v_in[N];
    VecN<T, N> vv(v);
    vv(v_in);

    T out[N];
    auto res = multiply(v_in, out);
    vout = VecN<T, N>(out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "multiply(T, VecN<T, N>&", res.status);
    return vflag;
  }
  /*! Tested */
  OpResult multiply(const T (&v)[N], T (&out)[N]) const {
    auto fn = [](T thisel, T argel) {
      return thisel * argel;
    };
    auto res = apply_el<decltype(fn)>(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "multiply(const T(&)[N], T (&)[N])",
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult multiply(const VecN<T, N> &v,
                    VecN<T, N> &out) const {
    T v_in[N];
    v(v_in);
    T vout[N];
    auto res = multiply(v_in, vout);
    out = VecN<T, N>(vout);

    OpResult vflag(
        __LINE__, __FILE__, __FUNCTION__,
        "multiply(const VecN<T, N>&, VecN<T, N>&)",
        res.status);
    return vflag;
  }

  //
  OpResult divide(T v, T (&out)[N]) const {
    if (v == static_cast<T>(0)) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     "divide(T, T (&)[N])", ARG_ERROR);
      return vflag;
    }

    T v_in[N];
    VecN<T, N> v_(v);
    v_(v_in);

    auto res = divide(v_in, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "divide(T, T (&)[N])", res.status);
    return vflag;
  }
  /*! Tested */
  OpResult divide(T v, VecN<T, N> &vout) const {
    // check for zero division
    if (v == static_cast<T>(0)) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     "divide(T, VecN<T, N>&", ARG_ERROR);
      return vflag;
    }
    T v_in[N];
    VecN<T, N> v_(v);
    v_(v_in);
    T out[N];

    auto res = divide(v_in, out);
    vout = VecN<T, N>(out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "divide(T, VecN<T, N>&", res.status);
    return vflag;
  }
  /*! Tested */
  OpResult divide(const T (&v)[N], T (&out)[N]) const {
    for (std::size_t j = 0; j < N; j++) {
      if (v[j] == static_cast<T>(0)) {
        OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                       "divide(const T(&)[N], T (&)[N])",
                       ARG_ERROR);
        return vflag;
      }
    }
    auto fn = [](T thisel, T argel) {
      return thisel / argel;
    };
    auto res = apply_el<decltype(fn)>(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "divide(const T(&)[N], T (&)[N])",
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult divide(const VecN<T, N> &v,
                  VecN<T, N> &out) const {
    // check zero division
    T v_in[N];
    v(v_in);
    T vout[N];
    auto res = divide(v_in, vout);
    out = VecN<T, N>(vout);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "divide(const VecN<T, N>&, VecN<T, N>&)",
                   res.status);
    return vflag;
  }
  OpResult dot(const T &v, T &out) const {
    T v_in[N];
    VecN<T, N> v_(v);
    v_(v_in);
    dot(v_in, out);
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "dot(const T&, T &)", SUCCESS);
    return vflag;
  }
  OpResult dot(const T (&v)[N], T &out) const {
    out = static_cast<T>(0);
    for (std::size_t i = 0; i < N; i++) {
      out += data[i] * v[i];
    }

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "dot(const T(&)[N], T &)", SUCCESS);
    return vflag;
  }
  OpResult dot(const VecN<T, N> &v, T &out) const {
    T v_in[N];
    v(v_in);
    dot(v_in, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "dot(const VecN<T, N>&, T &)", SUCCESS);
    return vflag;
  }

  OpResult cross(const VecN<T, N> &v,
                 VecN<T, N> &out) const {
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "cross(const VecN<T, N>&, VecN<T, N> &)",
                   NOT_IMPLEMENTED);
    return vflag;
  }

private:
  /** holds the vector data*/
  T data[N];

  template <typename Func>
  OpResult apply_el(const T (&v)[N], const Func &fn,
                    T (&out)[N]) const {

    for (std::size_t i = 0; i < N; i++) {
      out[i] = fn(data[i], v[i]);
    }
    OpResult vflag(
        __LINE__, __FILE__, __FUNCTION__,
        "apply_el(const T(&)[N], const Func&, T(&)[N])",
        SUCCESS);
    return vflag;
  }
};

template <typename T, std::size_t BaseOrder, std::size_t N>
OpResult base(VecN<T, N> &vout) {
  if (BaseOrder >= N) {
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "base(VecN<T, N>&)", ARG_ERROR);
    return vflag;
  }
  VecN<T, N> out;
  out(make_vecn_cell<T>(1, BaseOrder));
  vout = out;
  OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                 "base(VecN<T, N>&)", SUCCESS);
  return vflag;
}

template <typename T>
OpResult cross(const VecN<T, 2> &v1, const VecN<T, 2> &v2,
               T &out) {
  T t1;
  v1(0, t1);
  T t2;
  v1(1, t2);
  T v_1;
  v2(0, v_1);
  T v_2;
  v2(1, v_2);
  T t1v2 = t1 * v_2;
  T t2v1 = t2 * v_1;
  out = t1v2 - t2v1;
  OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                 "cross(const VecN<T, 2>&, T&)", SUCCESS);
  return vflag;
}

template <typename T>
OpResult cross(const VecN<T, 3> &a, const VecN<T, 3> &b,
               VecN<T, 3> &out) {
  T tx;
  a(0, tx);
  T ty;
  a(1, ty);
  T tz;
  a(2, tz);
  //
  T vx;
  b(0, vx);
  T vy;
  b(1, vy);
  T vz;
  b(2, vz);
  //
  T x = ty * vz - tz * vy;
  T y = tz * vx - tx * vz;
  T z = tx * vy - ty * vx;
  out(make_vecn_cell(x, 0));
  out(make_vecn_cell(y, 1));
  out(make_vecn_cell(z, 2));
  OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                 "cross(const VecN<T, 3>&, VecN<T, 3>&)",
                 SUCCESS);
  return vflag;
}

} // namespace vecn
} // namespace math3d

#endif
