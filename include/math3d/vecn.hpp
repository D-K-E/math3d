#ifndef VECN_HPP
#define VECN_HPP
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <math.h>
#include <ostream>
#include <stdio.h>
#include <vector>

#include "opflags.h"

namespace math3d {
namespace vecn {

template <class T, unsigned int N> class VecN {
  /** holds the vector data*/
  std::array<T, N> data;

public:
  /*! Tested */
  VecN() {}
  /*! Tested */
  VecN(const std::vector<T> &vd) {
    int nb_s = vd.size() - N;
    if (nb_s > 0) {
      // vector size is bigger than current vector
      for (unsigned int i = 0; i < N; i++) {
        data[i] = vd[i];
      }
    } else {
      // vector size is smaller than current vector
      for (unsigned int i = 0; i < N; i++) {
        if (i < vd.size()) {
          this(i, vd[i]);
        } else {
          this(i, static_cast<T>(0));
        }
      }
    }
  } /*! Tested */
  VecN(const std::array<T, N> &arr) : data(arr) {}
  VecN(T s): data.fill(s) {
    
  }
  /*! Tested */
  OpResult operator[](unsigned int &out) const {
    out = static_cast<unsigned int>(data.size());
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
  /*! Tested */
  OpResult operator[](unsigned int index, T &out) const {
    if (index >= data.size()) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     INDEX_ERROR);
      return vflag;
    }
    out = data[index];
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
  OpResult operator[](std::array<T, N> &out) const {
    out = data;
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
  /*! Tested */
  OpResult operator()(unsigned int index, const T &el) {
    if (index >= data.size()) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     INDEX_ERROR);
      return vflag;
    }
    data[index] = el;

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
  /*! Tested */
  static OpResult base(unsigned int nb_dimensions,
                       unsigned int base_order,
                       std::vector<T> &out) {
    if (base_order >= nb_dimensions) {

      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     ARG_ERROR);
      return vflag;
    }
    //
    if (out.size() != nb_dimensions) {
      out.clear();
      out.resize(static_cast<std::size_t>(nb_dimensions));
    }
    for (unsigned int i = 0; i < nb_dimensions; i++) {
      out[i] = static_cast<T>(0);
    }
    out[base_order] = static_cast<T>(1);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
  OpResult apply_el(T v, const std::function<T(T, T)> &fn,
                    std::vector<T> &out) const {
    if (out.size() != data.size()) {
      out.resize(data.size());
    }
    for (unsigned int i = 0; i < data.size(); i++) {
      out[i] = fn(data[i], v);
    }
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
  OpResult apply_el(const std::vector<T> &v,
                    const std::function<T(T, T)> &fn,
                    std::vector<T> &out) const {
    if (v.size() != data.size()) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     SIZE_ERROR);
      return vflag;
    }
    if (data.size() != out.size()) {
      out.resize(data.size());
    }
    for (unsigned int i = 0; i < data.size(); i++) {
      out[i] = fn(data[i], v[i]);
    }
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
  OpResult apply_el(T v, const std::function<T(T, T)> &fn,
                    VecN<T, N> &vout) const {
    std::vector<T> out;
    OpResult result = apply_el(v, fn, out);
    if (result.status != SUCCESS) {
      return result;
    }
    vout = VecN<T, N>(out);
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
  OpResult apply_el(const VecN<T, N> &v,
                    const std::function<T(T, T)> &fn,
                    VecN<T, N> &vout) const {

    for (unsigned int i = 0; i < data.size(); i++) {
      T tout = static_cast<T>(0);
      v[i, tout];                 // access
      vout(i, fn(data[i], tout)); // setter
    }

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
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
    for (unsigned int j = 0; j < v.size(); j++) {
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
    unsigned int vsize = 0;
    v.size(vsize);
    // check zero division
    for (unsigned int j = 0; j < vsize; j++) {
      T vout = static_cast<T>(0);
      v[j, vout]; // access
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
    for (unsigned int i = 0; i < data.size(); i++) {
      out += data[i] * v;
    }

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
  OpResult dot(const std::vector<T> &v, T &out) const {
    if (v.size() != data.size()) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     SIZE_ERROR);
      return vflag;
    }
    out = static_cast<T>(0);
    for (unsigned int i = 0; i < data.size(); i++) {
      out += data[i] * v[i];
    }

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
  OpResult dot(const VecN<T, N> &v, T &out) const {

    unsigned int vsize = 0;
    v.size(vsize);
    // check size
    if (vsize != data.size()) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                     SIZE_ERROR);
      return vflag;
    }

    out = static_cast<T>(0);
    for (unsigned int i = 0; i < data.size(); i++) {
      T tout = static_cast<T>(0);
      v[i, tout];
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
};

template <typename T, unsigned int BaseOrder,
          unsigned int N>
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
    this[0, t1];
    T t2;
    this[1, t2];
    T v1;
    v[0, v1];
    T v2;
    v[1, v2];
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
    this[0, t1];
    T ty;
    this[1, t2];
    T tz;
    this[2, t2];
    //
    T vx;
    v[0, v1];
    T vy;
    v[1, v2];
    T vz;
    v[2, v2];
    //
    T x = ty * vz - tz * vy;
    T y = tz * vx - tx * vz;
    T z = tx * vy - ty * vx;
    out(0, x);
    out(1, y);
    out(2, z);
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   SUCCESS);
    return vflag;
  }
};

} // namespace vecn
} // namespace math3d

#endif
