#ifndef VECN_HPP
#define VECN_HPP
#include <cstdint>
#include <cstring>
#include <functional>
#include <iostream>
#include <math.h>
#include <ostream>
#include <stdio.h>
#include <vector>

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
  /** holds the vector data*/
  T data[N];

public:
  /*! Tested */
  VecN() {}
  /*! Tested */
  VecN(const std::vector<T> &vd) {
    int nb_s = vd.size() - N;
    if (nb_s > 0) {
      // vector size is bigger than current vector
      for (std::size_t i = 0; i < N; i++) {
        data[i] = vd[i];
      }
    } else {
      // vector size is smaller than current vector
      for (std::size_t i = 0; i < N; i++) {
        if (i < vd.size()) {
          (*this)(make_vecn_cell<T>(vd[i], i));
        } else {
          (*this)(make_vecn_cell<T>(0, i));
        }
      }
    }
  } /*! Tested */
  VecN(const T (&arr)[N]) {
    memcpy(data, arr, N * sizeof(T));
  }
  VecN(T s) {
    for (std::size_t i = 0; i < N; ++i) {
      data[i] = s;
    }
  }
  /*! Tested */
  OpResult operator()(std::size_t &out) const {
    out = N;
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "(size_t)", SUCCESS);
    return vflag;
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
  /*! Tested */
  OpResult operator()(const VecNCell<T> &cell) {
    if (cell.index >= N) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     "(const VecNCell<T>&)", INDEX_ERROR);
      return vflag;
    }
    data[cell.index] = cell.el;

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "(const VecNCell<T>&)", SUCCESS);
    return vflag;
  }

  /*! Tested */
  OpResult add(T v, std::vector<T> &out) const {
    auto fn = [](T thisel, T argel) {
      return thisel + argel;
    };
    auto res = apply_el(v, fn, out);
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }

  /*! Tested */
  OpResult add(T v, VecN<T, N> &vout) const {
    auto fn = [](T thisel, T argel) {
      return thisel + argel;
    };
    auto res = apply_el(v, fn, vout);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult add(const std::vector<T> &v,
               std::vector<T> &out) const {
    auto fn = [](T thisel, T argel) {
      return thisel + argel;
    };
    auto res = apply_el(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult add(const VecN<T, N> &v, VecN<T, N> &out) const {
    auto fn = [](T thisel, T argel) {
      return thisel + argel;
    };
    auto res = apply_el(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  //
  OpResult subtract(T v, std::vector<T> &out) const {
    auto fn = [](T thisel, T argel) {
      return thisel - argel;
    };
    auto res = apply_el(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult subtract(T v, VecN<T, N> &vout) const {
    auto fn = [](T thisel, T argel) {
      return thisel - argel;
    };
    auto res = apply_el(v, fn, vout);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult subtract(const std::vector<T> &v,
                    std::vector<T> &out) const {
    auto fn = [](T thisel, T argel) {
      return thisel - argel;
    };
    auto res = apply_el(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult subtract(const VecN<T, N> &v,
                    VecN<T, N> &out) const {
    auto fn = [](T thisel, T argel) {
      return thisel - argel;
    };
    auto res = apply_el(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  //
  OpResult multiply(T v, std::vector<T> &out) const {
    auto fn = [](T thisel, T argel) {
      return thisel * argel;
    };
    auto res = apply_el(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult multiply(T v, VecN<T, N> &vout) const {
    auto fn = [](T thisel, T argel) {
      return thisel * argel;
    };
    auto res = apply_el(v, fn, vout);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult multiply(const std::vector<T> &v,
                    std::vector<T> &out) const {
    auto fn = [](T thisel, T argel) {
      return thisel * argel;
    };
    auto res = apply_el(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult multiply(const VecN<T, N> &v,
                    VecN<T, N> &out) const {
    auto fn = [](T thisel, T argel) {
      return thisel * argel;
    };
    auto res = apply_el(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }

  //
  OpResult divide(T v, std::vector<T> &out) const {
    if (v == static_cast<T>(0)) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     ARG_ERROR);
      return vflag;
    }

    auto fn = [](T thisel, T argel) {
      return thisel / argel;
    };
    auto res = apply_el(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult divide(T v, VecN<T, N> &vout) const {
    // check for zero division
    if (v == static_cast<T>(0)) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     ARG_ERROR);
      return vflag;
    }
    auto fn = [](T thisel, T argel) {
      return thisel / argel;
    };
    auto res = apply_el(v, fn, vout);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult divide(const std::vector<T> &v,
                  std::vector<T> &out) const {
    for (std::size_t j = 0; j < v.size(); j++) {
      if (v[j] == static_cast<T>(0)) {
        OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                       ARG_ERROR);
        return vflag;
      }
    }
    auto fn = [](T thisel, T argel) {
      return thisel / argel;
    };
    auto res = apply_el(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  /*! Tested */
  OpResult divide(const VecN<T, N> &v,
                  VecN<T, N> &out) const {
    std::size_t vsize = 0;
    v.size(vsize);
    // check zero division
    for (std::size_t j = 0; j < vsize; j++) {
      T vout = static_cast<T>(0);
      v(j, vout); // access
      if (vout == static_cast<T>(0)) {

        OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                       ARG_ERROR);
        return vflag;
      }
    }
    auto fn = [](T thisel, T argel) {
      return thisel / argel;
    };
    auto res = apply_el(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   res.status);
    return vflag;
  }
  OpResult dot(const T &v, T &out) const {
    out = static_cast<T>(0);
    for (std::size_t i = 0; i < N; i++) {
      out += data[i] * v;
    }

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
  OpResult dot(const std::vector<T> &v, T &out) const {
    if (v.size() != N) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     SIZE_ERROR);
      return vflag;
    }
    out = static_cast<T>(0);
    for (std::size_t i = 0; i < N; i++) {
      out += data[i] * v[i];
    }

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
  OpResult dot(const VecN<T, N> &v, T &out) const {

    std::size_t vsize = 0;
    v.size(vsize);
    // check size
    if (vsize != N) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     SIZE_ERROR);
      return vflag;
    }

    out = static_cast<T>(0);
    for (std::size_t i = 0; i < N; i++) {
      T tout;
      v(i, tout);
      out += data[i] * tout;
    }

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }

  OpResult cross(const VecN<T, N> &v,
                 VecN<T, N> &out) const {
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   NOT_IMPLEMENTED);
    return vflag;
  }

private:
  template<typename Func>
  OpResult apply_el(T v, const Func &fn,
                    T (&out)[N]) const {
    for (std::size_t i = 0; i < N; i++) {
      out[i] = fn(data[i], v);
    }
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "apply_el(T, <T(T,T)>, T(&)[N])",
                   SUCCESS);
    return vflag;
  }
  template<typename Func>
  OpResult apply_el(const T (&v)[N],
                    const Func &fn
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
  template <typename Func>
  OpResult apply_el(T v, const Func &fn,
                    VecN<T, N> &vout) const {
    std::vector<T> out;
    OpResult result = apply_el(v, fn, out);
    if (result.status != SUCCESS) {
      return result;
    }
    vout = VecN<T, N>(out);
    OpResult vflag(
        __LINE__, __FILE__, __FUNCTION__,
        "apply_el(T,const <T(T,T)>&, VecN<T, N>&)",
        SUCCESS);
    return vflag;
  }

  template <typename Func>
  OpResult apply_el(const VecN<T, N> &v,
                    const Func &fn,
                    VecN<T, N> &vout) const {

    for (std::size_t i = 0; i < N; i++) {
      T tout = static_cast<T>(0);
      v(i, tout);                                 // access
      vout(make_vecn_cell(fn(data[i], tout), i)); // setter
    }

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
};

template <typename T, std::size_t BaseOrder, std::size_t N>
OpResult base(VecN<T, N> &vout) {
  if (BaseOrder >= N) {
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   ARG_ERROR);
    return vflag;
  }
  VecN<T, N> out;
  out(BaseOrder, 1);
  vout = out;
  OpResult vflag(__LINE__, __FILE__, __FUNCTION__, SUCCESS);
  return vflag;
}

template <typename T> class VecN<T, 2> {

  OpResult cross(const VecN<T, 2> &v, T &out) const {
    T t1;
    this(0, t1);
    T t2;
    this(1, t2);
    T v1;
    v(0, v1);
    T v2;
    v(1, v2);
    T t1v2 = t1 * v2;
    T t2v1 = t2 * v1;
    out = t1v2 - t2v1;
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
};

template <typename T> class VecN<T, 3> {

  OpResult cross(const VecN<T, 3> &v,
                 VecN<T, 3> &out) const {
    T tx;
    (*this)(0, tx);
    T ty;
    (*this)(1, ty);
    T tz;
    (*this)(2, tz);
    //
    T vx;
    v(0, vx);
    T vy;
    v(1, vy);
    T vz;
    v(2, vz);
    //
    T x = ty * vz - tz * vy;
    T y = tz * vx - tx * vz;
    T z = tx * vy - ty * vx;
    out(make_vecn_cell(x, 0));
    out(make_vecn_cell(y, 1));
    out(make_vecn_cell(z, 2));
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
};

} // namespace vecn
} // namespace math3d

#endif
