#ifndef MATH3D_HPP
#define MATH3D_HPP

#include "core.h"
#include "utils.h"

namespace math3d {

template <typename T>
OpResult to_skew_mat(const vecn::VecN<T, 3> &vec,
                     matn::MatN<T, 3, 3> &out) {
  //
  T a3;
  vec(2, a3);
  T a2;
  vec(1, a2);
  T a1;
  vec(0, a1);

  matn::MatN<T, 3, 3> m;

  T c1[] = {0, a3, -a2};
  T c2[] = {-a3, 0, a1};
  T c3[] = {a2, -a1, 0};
  m.set_column(0, c1);
  m.set_column(1, c2);
  m.set_column(2, c3);
  out = m;
  return OpResult(__LINE__, __FILE__, __FUNCTION__,
                  "to_skew_mat", SUCCESS);
}

/**convert quaternion to rotation matrix
  from Vince, 2011, Quaternion ..., p. 123
 */
template <typename T>
OpResult toRotMat3x3(const quaternion::Quaternion<T> &q_in,
                     matn::MatN<T, 3, 3> &out) {
  quaternion::Quaternion<T> q;
  q_in.normalized(q);
  //
  T s;
  q.scalar(s);

  //
  vecn::VecN<T, 3> xyz;
  q.vector(xyz);
  T x;
  xyz(0, x);
  T y;
  xyz(1, y);
  T z;
  xyz(2, z);

  // to rotation matrix
  T c1_1 = 1 - (2 * ((y * y) + (z * z)));
  T c1_2 = 2 * (x * y + s * z);
  T c1_3 = 2 * (x * z - s * y);
  //
  T c2_1 = 2 * (x * y - s * z);
  T c2_2 = 1 - 2 * ((x * x) + (z * z));
  T c2_3 = 2 * (y * z + s * x);
  //
  T c3_1 = 2 * (x * z + s * y);
  T c3_2 = 2 * (y * z - s * x);
  T c3_3 = 1 - 2 * ((x * x) + (y * y));

  T c1[] = {c1_1, c1_2, c1_3};
  T c2[] = {c2_1, c2_2, c2_3};
  T c3[] = {c3_1, c3_2, c3_3};
  matn::MatN<T, 3, 3> m;
  m.set_column(0, c1);
  m.set_column(1, c2);
  m.set_column(2, c3);
  out = m;
  return OpResult(__LINE__, __FILE__, __FUNCTION__,
                  "toRotMat3x3", SUCCESS);
}

template <typename T, std::size_t N>
OpResult invert_mat(const matn::MatN<T, N, N> &m,
                    matn::MatN<T, N, N> &out) {
  //
  lu::LUdecomp<T, N> lu_d(m);
  matn::MatN<T, N, N> B;
  matn::identity<T, N>(B);
  matn::MatN<T, N, N> inv_m;
  lu_d.solve_mat(B, inv_m);
  out = inv_m;
  return OpResult(__LINE__, __FILE__, __FUNCTION__,
                  "invert_mat", SUCCESS);
}

} // namespace math3d

#endif
