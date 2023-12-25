#ifndef MATH3D_HPP
#define MATH3D_HPP

#include "lu.hpp"
#include "matn.hpp"
#include "opflags.hpp"
#include "quaternion.hpp"
#include "vecn.hpp"

namespace math3d {

std::ostream &operator<<(std::ostream &out,
                         opstatus_t flag) {
  switch (flag) {
  case SUCCESS: {
    out << "SUCCESS" << std::endl;
    break;
  }
  case SIZE_ERROR: {
    out << "SIZE_ERROR "
        << " :: " << __FILE__ << " :: " << __LINE__;
    break;
  }
  case INDEX_ERROR: {
    out << "INDEX_ERROR "
        << " :: " << __FILE__ << " :: " << __LINE__;
    break;
  }
  case ARG_ERROR: {
    out << "ARG_ERROR "
        << " :: " << __FILE__ << " :: " << __LINE__;
    break;
  }
  case NOT_IMPLEMENTED: {
    out << "NOT_IMPLEMENTED "
        << " :: " << __FILE__ << " :: " << __LINE__;
    break;
  }
  case LU_ERROR: {
    out << "LU_DECOMPOSITION_ERROR "
        << " :: " << __FILE__ << " :: " << __LINE__;
    break;
  }
  }
  return out;
}

template <typename T, unsigned int R, unsigned int C = R>
std::ostream &operator<<(std::ostream &out,
                         const matn::MatN<T, R, C> &m) {
  std::array<T, R * C> arr;
  m[arr];
  for (std::size_t i = 0; i < arr.size(); i++) {
    if (i % C == 0) {
      out << std::endl;
    }
    if (arr[i] >= 0) {
      out << " " << arr[i] << " ";
    } else {
      out << arr[i] << " ";
    }
  }
  return out;
}

template <typename T, unsigned int R>
std::ostream &operator<<(std::ostream &out,
                         const vecn::VecN<T, R> &m) {
  std::array<T, R> arr;
  m[arr];
  for (std::size_t i = 0; i < arr.size(); i++) {
    if (i % R == 0) {
      out << std::endl;
    }
    if (arr[i] >= 0) {
      out << " " << arr[i] << " ";
    } else {
      out << arr[i] << " ";
    }
  }
  return out;
}

template <typename T>
std::ostream &
operator<<(std::ostream &out,
           const quaternion::QuaternionBase<T> &c) {
  switch (c.base) {
  case kSCALAR_BASE: {
    out << "SCALAR_BASE::" << c.r << std::endl;
    break;
  }
  case kI: {
    out << "I_BASE::" << c.r << std::endl;
    break;
  }
  case kJ: {
    out << "J_BASE::" << c.r << std::endl;
    break;
  }
  case kK: {
    out << "K_BASE::" << c.r << std::endl;
    break;
  }
  }
  return out;
}
template <typename T>
std::ostream &
operator<<(std::ostream &out,
           const quaternion::QuaternionComponent<T> &c) {
  switch (c.base) {
  case kSCALAR_BASE: {
    out << c.r << std::endl;
    break;
  }
  case kI: {
    out << c.r << "i";
    break;
  }
  case kJ: {
    out << c.r << "j";
    break;
  }
  case kK: {
    out << c.r << "k";
    break;
  }
  }
  return out;
}

template <typename T>
std::ostream &
operator<<(std::ostream &out,
           const quaternion::Quaternion<T> &q) {
  quaternion::QuaternionComponent<T> c(kSCALAR_BASE, 0);
  q[c];
  out << c << " + ";
  c = quaternion::QuaternionComponent<T>(kI, 0);
  q[c];
  out << c << " + ";
  c = quaternion::QuaternionComponent<T>(kJ, 0);
  q[c];
  out << c << " + ";
  c = quaternion::QuaternionComponent<T>(kK, 0);
  q[c];
  out << c << std::endl;
  return out;
}

template <typename T>
matn::MatN<T, 3, 3>
to_skew_mat(const vecn::VecN<T, 3> &vec) {
  //
  T a3;
  vec[2, a3];
  T a2;
  vec[1, a2];
  T a1;
  vec[0, a1];

  matn::MatN<T, 3, 3> m;

  std::array<T, 3> c1{0, a3, -a2};
  std::array<T, 3> c2{-a3, 0, a1};
  std::array<T, 3> c3{a2, -a1, 0};
  m.set_column(0, c1);
  m.set_column(1, c2);
  m.set_column(2, c3);
  return m;
}

/**convert quaternion to rotation matrix
  from Vince, 2011, Quaternion ..., p. 123
 */
template <typename T>
matn::MatN<T, 3, 3>
toRotMat3x3(const quaternion::Quaternion<T> &q_in) {
  quaternion::Quaternion<T> q;
  q_in.normalized(q);
  //
  T s;
  q.scalar(s);

  //
  vecn::VecN<T, 3> xyz;
  q.vector(xyz);
  T x;
  xyz[0, x];
  T y;
  xyz[1, y];
  T z;
  xyz[2, z];

  // to rotation matrix
  T c1_1 = 1 - (2 * (std::pow(y, 2) + std::pow(z, 2)));
  T c1_2 = 2 * (x * y + s * z);
  T c1_3 = 2 * (x * z - s * y);
  //
  T c2_1 = 2 * (x * y - s * z);
  T c2_2 = 1 - (2 * (std::pow(x, 2) + std::pow(z, 2)));
  T c2_3 = 2 * (y * z + s * x);
  //
  T c3_1 = 2 * (x * z + s * y);
  T c3_2 = 2 * (y * z - s * x);
  T c3_3 = 1 - 2 * (std::pow(x, 2) + std::pow(y, 2));

  std::array<T, 3> c1{c1_1, c1_2, c1_3};
  std::array<T, 3> c2{c2_1, c2_2, c2_3};
  std::array<T, 3> c3{c3_1, c3_2, c3_3};
  matn::MatN<T, 3, 3> m;
  m.set_column(0, c1);
  m.set_column(1, c2);
  m.set_column(2, c3);
  return m;
}

template <typename T, unsigned int N>
matn::MatN<T, N, N>
invert_mat(const matn::MatN<T, N, N> &m) {
  //
  lu::LUdecomp lu(m);
  matn::MatN<T, N, N> B;
  matn::MatN<T, N, N>::identity(N, B);
  matn::MatN<T, N, N> inv_m;
  lu.solve_mat(B, inv_m);
  return inv_m;
}

} // namespace math3d

#endif
