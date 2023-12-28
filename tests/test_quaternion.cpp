#include <gtest/gtest.h>
#include <math3d/quaternion.hpp>

typedef float real;
using namespace math3d::quaternion;
using namespace math3d;

/*! @{
  Test accessors of Quaternion object. These functions help
  us to access
  real parts, that is the scalar and the vector part of the
  Quaternion object.
 */

TEST(QuaternionTest, test_scalar1) {
  Quaternion<real> q;
  real t = static_cast<real>(4561);
  auto res = q.scalar(t);

  EXPECT_EQ(t, 0);
  EXPECT_EQ(res.status, SUCCESS);
}
TEST(QuaternionTest, test_vector1) {
  Quaternion<real> q;
  vecn::VecN<real, 3> v_s;
  auto res = q.vector(v_s);
  real vs[3];
  v_s(vs);
  real comp[3] = {1, 1, 1};
  EXPECT_EQ(res.status, SUCCESS);
  EXPECT_EQ(vs[0], comp[0]);
  EXPECT_EQ(vs[1], comp[1]);
  EXPECT_EQ(vs[2], comp[2]);
};

/*! @{ test component accessors */
TEST(QuaternionTest, test_get_component_0) {
  Quaternion<real> q(2, vecn::VecN<real, 3>(3));
  QuaternionComponent<real> comp;
  auto result = q(kSCALAR, comp);
  EXPECT_EQ(result.status, SUCCESS);
  EXPECT_EQ(comp.r, 2);
  EXPECT_EQ(comp.base, kSCALAR);
}

TEST(QuaternionTest, test_get_component_1) {
  Quaternion<real> q(2, vecn::VecN<real, 3>(3));
  QuaternionComponent<real> comp;
  auto result = q(kI, comp);
  EXPECT_EQ(result.status, SUCCESS);
  EXPECT_EQ(comp.r, 3);
  EXPECT_EQ(comp.base, kI);
}

TEST(QuaternionTest, test_get_component_2) {
  Quaternion<real> q(2, vecn::VecN<real, 3>(3));
  QuaternionComponent<real> comp;
  auto result = q(kJ, comp);
  EXPECT_EQ(result.status, SUCCESS);
  EXPECT_EQ(comp.r, 3);
  EXPECT_EQ(comp.base, kJ);
}

TEST(QuaternionTest, test_get_component_3) {
  Quaternion<real> q(2, vecn::VecN<real, 3>(3));
  QuaternionComponent<real> comp;
  auto result = q(kK, comp);
  EXPECT_EQ(result.status, SUCCESS);
  EXPECT_EQ(comp.r, 3);
  EXPECT_EQ(comp.base, kK);
}
/*! @} */
/*! @} */

/*! @{ Test constructors of the Quaternion object*/

TEST(QuaternionTest, test_constructor_0) {
  Quaternion<real> q;
  real s = static_cast<real>(1);
  auto res = q.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec2;
  res = q.vector(vec2);
  EXPECT_EQ(res.status, SUCCESS);
  real vec[3];
  vec2(vec);
  //
  EXPECT_EQ(s, 0);
  EXPECT_EQ(vec[0], 1);
  EXPECT_EQ(vec[1], 1);
  EXPECT_EQ(vec[2], 1);
}

TEST(QuaternionTest, test_constructor_1) {
  real v[3] = {3, 4, 5};
  Quaternion<real> q(2, vecn::VecN<real, 3>(v));

  real s = static_cast<real>(1);
  auto res = q.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec2;
  res = q.vector(vec2);
  EXPECT_EQ(res.status, SUCCESS);
  //
  real vec[3];
  vec2(vec);
  EXPECT_EQ(s, 2);
  EXPECT_EQ(vec[0], 3);
  EXPECT_EQ(vec[1], 4);
  EXPECT_EQ(vec[2], 5);
}

TEST(QuaternionTest, test_constructor_2) {
  real c_1 = 2;
  real cs[3];
  cs[0] = 3;
  cs[1] = 4;
  cs[2] = 5;
  Quaternion<real> q(c_1, vecn::VecN<real, 3>(cs));

  real s = static_cast<real>(1);
  auto res = q.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec2;
  res = q.vector(vec2);
  EXPECT_EQ(res.status, SUCCESS);
  real vec[3];
  vec2(vec);
  //
  EXPECT_EQ(s, 2);
  EXPECT_EQ(vec[0], 3);
  EXPECT_EQ(vec[1], 4);
  EXPECT_EQ(vec[2], 5);
}

TEST(QuaternionTest, test_constructor_3) {
  QuaternionComponent<real> c1(kSCALAR, 2);
  QuaternionComponent<real> c2(kI, 3);
  QuaternionComponent<real> c3(kJ, 4);
  QuaternionComponent<real> c4(kK, 5);
  QuaternionComponent<real> cs[] = {c1, c2, c3, c4};
  Quaternion<real> q(cs);

  real s = static_cast<real>(1);
  auto res = q.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec1;
  res = q.vector(vec1);
  EXPECT_EQ(res.status, SUCCESS);
  real vec[3];
  vec1(vec);
  //
  EXPECT_EQ(s, 2);
  EXPECT_EQ(vec[0], 3);
  EXPECT_EQ(vec[1], 4);
  EXPECT_EQ(vec[2], 5);
}

TEST(QuaternionTest, test_constructor_4) {
  real c1 = 2;
  QuaternionComponent<real> c2(kI, 3);
  QuaternionComponent<real> c3(kJ, 4);
  QuaternionComponent<real> c4(kK, 5);
  QuaternionComponent<real> cs[] = {c2, c3, c4};
  Quaternion<real> q(c1, cs);

  real s = static_cast<real>(1);
  auto res = q.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec1;
  res = q.vector(vec1);
  EXPECT_EQ(res.status, SUCCESS);
  real vec[3];
  vec1(vec);
  //
  EXPECT_EQ(s, 2);
  EXPECT_EQ(vec[0], 3);
  EXPECT_EQ(vec[1], 4);
  EXPECT_EQ(vec[2], 5);
}

TEST(QuaternionTest, test_constructor_5) {
  real c1 = 2;
  vecn::VecN<real, 3> cs({3, 4, 5});
  Quaternion<real> q(c1, cs);
  real s = static_cast<real>(1);
  auto res = q.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec1;
  res = q.vector(vec1);
  EXPECT_EQ(res.status, SUCCESS);
  real vec[3];
  vec1(vec);
  //
  EXPECT_EQ(s, 2);
  EXPECT_EQ(vec[0], 3);
  EXPECT_EQ(vec[1], 4);
  EXPECT_EQ(vec[2], 5);
}

/*! @} */

/*! @{ Test vector operations */

TEST(QuaternionTest, test_vector_multiplication_0) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({3, 4, 5}));
  vecn::VecN<real, 3> v1;

  auto res = q.multiply(2, v1);
  EXPECT_EQ(res.status, SUCCESS);

  real s = static_cast<real>(1);
  res = q.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);
  real v[3];
  v1(v);

  //
  EXPECT_EQ(s, 2);
  EXPECT_EQ(v[0], 6);
  EXPECT_EQ(v[1], 8);
  EXPECT_EQ(v[2], 10);
}
TEST(QuaternionTest, test_vector_multiplication_1) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({3, 4, 5}));
  vecn::VecN<real, 3> v1;
  auto res = q.multiply(2, v1);
  EXPECT_EQ(res.status, SUCCESS);
  real v[3];
  v1(v);
  //
  EXPECT_EQ(v[0], 6);
  EXPECT_EQ(v[1], 8);
  EXPECT_EQ(v[2], 10);
}
TEST(QuaternionTest, test_vector_multiplication_2) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({3, 4, 5}));
  vecn::VecN<real, 3> v1;
  real t[3];
  t[0] = 2;
  t[1] = 3;
  t[2] = 1;
  auto result = q.multiply(t, v1);
  real v[3];
  v1(v);
  EXPECT_EQ(v[0], 6);
  EXPECT_EQ(v[1], 12);
  EXPECT_EQ(v[2], 5);
  EXPECT_EQ(result.status, SUCCESS);
}
TEST(QuaternionTest, test_vector_addition_0) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({3, 4, 5}));
  vecn::VecN<real, 3> v_out;
  auto result = q.add(2, v_out);
  real v[3];
  v_out(v);
  real s = static_cast<real>(10);
  q.scalar(s);
  EXPECT_EQ(s, 2);
  EXPECT_EQ(v[0], 5);
  EXPECT_EQ(v[1], 6);
  EXPECT_EQ(v[2], 7);
  EXPECT_EQ(result.status, SUCCESS);
}
TEST(QuaternionTest, test_vector_addition_1) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({3, 4, 5}));
  vecn::VecN<real, 3> v_out;
  auto result = q.add(2, v_out);
  real v[3];
  v_out(v);
  EXPECT_EQ(v[0], 5);
  EXPECT_EQ(v[1], 6);
  EXPECT_EQ(v[2], 7);
  EXPECT_EQ(result.status, SUCCESS);
}
TEST(QuaternionTest, test_vector_addition_2) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({3, 4, 5}));
  vecn::VecN<real, 3> v_in;
  vecn::VecN<real, 3> t({2, 3, 1});
  auto result = q.add(t, v_in);
  real v[3];
  v_in(v);
  EXPECT_EQ(v[0], 5);
  EXPECT_EQ(v[1], 7);
  EXPECT_EQ(v[2], 6);
  EXPECT_EQ(result.status, SUCCESS);
}
TEST(QuaternionTest, test_vector_subtraction_0) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({3, 4, 5}));
  vecn::VecN<real, 3> v_in;
  q.vector(v_in);
  auto result = q.subtract(2, v_in);
  real s = static_cast<real>(22);
  q.scalar(s);
  real v[3];
  v_in(v);
  EXPECT_EQ(s, 2);
  EXPECT_EQ(v[0], 1);
  EXPECT_EQ(v[1], 2);
  EXPECT_EQ(v[2], 3);
  EXPECT_EQ(result.status, SUCCESS);
}
TEST(QuaternionTest, test_vector_subtraction_1) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({3, 4, 5}));
  vecn::VecN<real, 3> v_in;
  auto result = q.subtract(2, v_in);
  real v[3];
  v_in(v);
  EXPECT_EQ(v[0], 1);
  EXPECT_EQ(v[1], 2);
  EXPECT_EQ(v[2], 3);
  EXPECT_EQ(result.status, SUCCESS);
}
TEST(QuaternionTest, test_vector_subtraction_2) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({3, 4, 5}));
  vecn::VecN<real, 3> v_in;
  vecn::VecN<real, 3> t({2, 3, 1});
  auto result = q.subtract(t, v_in);
  real v[3];
  v_in(v);
  EXPECT_EQ(v[0], 1);
  EXPECT_EQ(v[1], 1);
  EXPECT_EQ(v[2], 4);
  EXPECT_EQ(result.status, SUCCESS);
}
TEST(QuaternionTest, test_vector_division_0) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({4, 4, 6}));
  vecn::VecN<real, 3> v_in;
  q.vector(v_in);
  auto result = q.divide(2, v_in);

  real s = static_cast<real>(1);
  result = q.scalar(s);
  EXPECT_EQ(result.status, SUCCESS);
  real v[3];
  v_in(v);

  EXPECT_EQ(s, 2);
  EXPECT_EQ(v[0], 2);
  EXPECT_EQ(v[1], 2);
  EXPECT_EQ(v[2], 3);
  EXPECT_EQ(result.status, SUCCESS);
}
TEST(QuaternionTest, test_vector_division_1) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({2, 4, 6}));
  vecn::VecN<real, 3> v_in;
  auto result = q.divide(2, v_in);
  real v[3];
  v_in(v);
  EXPECT_EQ(v[0], 1);
  EXPECT_EQ(v[1], 2);
  EXPECT_EQ(v[2], 3);
  EXPECT_EQ(result.status, SUCCESS);
}
TEST(QuaternionTest, test_vector_division_2) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({4, 4, 5}));
  vecn::VecN<real, 3> v_in;
  vecn::VecN<real, 3> t({2, 4, 1});
  auto result = q.divide(t, v_in);
  real v[3];
  v_in(v);
  EXPECT_EQ(v[0], 2);
  EXPECT_EQ(v[1], 1);
  EXPECT_EQ(v[2], 5);
  EXPECT_EQ(result.status, SUCCESS);
}
TEST(QuaternionTest, test_vector_division_3) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({3, 4, 5}));
  vecn::VecN<real, 3> v;
  q.vector(v);
  auto result = q.divide(static_cast<real>(0), v);
  EXPECT_EQ(result.status, ARG_ERROR);
}
TEST(QuaternionTest, test_vector_division_5) {
  Quaternion<real> q(2, vecn::VecN<real, 3>({3, 4, 5}));
  vecn::VecN<real, 3> v;
  vecn::VecN<real, 3> t({2, 0, 1});
  auto result = q.divide(t, v);
  EXPECT_EQ(result.status, ARG_ERROR);
}
TEST(QuaternionTest, test_vector_dot_0) {
  Quaternion<real> q(1, vecn::VecN<real, 3>({9, 2, 7}));
  vecn::VecN<real, 3> t({4, 8, 10});
  real out;
  auto result = q.dot(t, out);
  EXPECT_EQ(out, 122);
  EXPECT_EQ(result.status, SUCCESS);
}
TEST(QuaternionTest, test_vector_dot_1) {
  Quaternion<real> q(0, vecn::VecN<real, 3>({9, 2, 7}));
  vecn::VecN<real, 3> t({4, 8, 10});
  real s = static_cast<real>(55156);
  auto result = q.dot(t, s);
  EXPECT_EQ(s, 122);
  EXPECT_EQ(result.status, SUCCESS);
}
TEST(QuaternionTest, test_vector_cross_0) {
  Quaternion<real> q(1, vecn::VecN<real, 3>({2, 3, 4}));
  vecn::VecN<real, 3> vout;
  vecn::VecN<real, 3> t({5, 6, 7});
  auto result = q.cross(t, vout);
  real out[3];
  vout(out);
  EXPECT_EQ(result.status, SUCCESS);
  EXPECT_EQ(out[0], static_cast<real>(-3));
  EXPECT_EQ(out[1], static_cast<real>(6));
  EXPECT_EQ(out[2], static_cast<real>(-3));
}
/*! @} */

/*! @{ Test hamilton product operation from Vince 2011 -
 * Quaternions for
 * Computer Graphics, p. 70 */
TEST(QuaternionTest, test_hamilton_product) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  Quaternion<real> q_b(1, vecn::VecN<real, 3>({-2, 5, -6}));
  Quaternion<real> q_out;
  q_a.hamilton_product(q_b, q_out);
  // Quaternion q_ab = q_a * q_b;

  real s = static_cast<real>(0);
  auto res = q_out.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> v_in;
  res = q_out.vector(v_in);
  EXPECT_EQ(res.status, SUCCESS);
  real v[3];
  v_in(v);

  EXPECT_EQ(s, static_cast<real>(-41));
  EXPECT_EQ(v[0], static_cast<real>(-4));
  EXPECT_EQ(v[1], static_cast<real>(9));
  EXPECT_EQ(v[2], static_cast<real>(-20));
}
/*! @} */
/*! @{ Test conjugate operation from Vince 2011, p. 70 */
TEST(QuaternionTest, test_conjugate_0) {
  Quaternion<real> q_b(2, vecn::VecN<real, 3>({-2, 3, -4}));
  Quaternion<real> q_out;
  q_b.conjugate(q_out);

  real s = static_cast<real>(0);
  auto res = q_out.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> v_in;
  res = q_out.vector(v_in);
  EXPECT_EQ(res.status, SUCCESS);
  real v[3];
  v_in(v);

  EXPECT_EQ(s, static_cast<real>(2));
  EXPECT_EQ(v[0], static_cast<real>(2));
  EXPECT_EQ(v[1], static_cast<real>(-3));
  EXPECT_EQ(v[2], static_cast<real>(4));
}
TEST(QuaternionTest, test_conjugate_1) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  Quaternion<real> q_b(1, vecn::VecN<real, 3>({-2, 5, -6}));
  Quaternion<real> q_ab_prod;
  q_a.hamilton_product(q_b, q_ab_prod);

  Quaternion<real> q_out1;
  q_ab_prod.conjugate(q_out1);

  //
  Quaternion<real> q_b_conj;
  q_b.conjugate(q_b_conj);

  Quaternion<real> q_a_conj;
  q_a.conjugate(q_a_conj);

  Quaternion<real> q_out2;
  q_b_conj.hamilton_product(q_a_conj, q_out2);

  //
  real s = static_cast<real>(0);
  real sq = static_cast<real>(0);

  q_out1.scalar(s);
  q_out2.scalar(sq);

  //
  vecn::VecN<real, 3> v_in;
  vecn::VecN<real, 3> vq_in;
  q_out1.vector(v_in);
  q_out2.vector(vq_in);
  real v[3];
  v_in(v);
  real vq[3];
  vq_in(vq);

  EXPECT_EQ(s, sq);
  EXPECT_EQ(v[0], vq[0]);
  EXPECT_EQ(v[1], vq[1]);
  EXPECT_EQ(v[2], vq[2]);
}

/*! @{ arithmetic ops for Quaternions */
TEST(QuaternionTest, test_add_operator) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  Quaternion<real> q_b(1, vecn::VecN<real, 3>({-2, 5, -6}));
  Quaternion<real> q_ab;
  q_a.add(q_b, q_ab);

  real s = static_cast<real>(1);
  auto res = q_ab.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec_in;
  res = q_ab.vector(vec_in);
  EXPECT_EQ(res.status, SUCCESS);
  real vec[3];
  vec_in(vec);

  EXPECT_EQ(s, 3);
  EXPECT_EQ(vec[0], -4);
  EXPECT_EQ(vec[1], 8);
  EXPECT_EQ(vec[2], -10);
}
TEST(QuaternionTest, test_subtract_operator) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  Quaternion<real> q_b(1, vecn::VecN<real, 3>({-2, 5, -6}));
  Quaternion<real> q_ab;
  q_a.subtract(q_b, q_ab);

  real s = static_cast<real>(1);
  auto res = q_ab.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec_in;
  res = q_ab.vector(vec_in);
  EXPECT_EQ(res.status, SUCCESS);
  real vec[3];
  vec_in(vec);

  EXPECT_EQ(s, 1);
  EXPECT_EQ(vec[0], 0);
  EXPECT_EQ(vec[1], -2);
  EXPECT_EQ(vec[2], 2);
}
TEST(QuaternionTest, test_product_operator_0) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  Quaternion<real> q_b(1, vecn::VecN<real, 3>({-2, 5, -6}));
  Quaternion<real> q_ab;
  q_a.product(q_b, q_ab);

  real s = static_cast<real>(1);
  auto res = q_ab.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec_i;
  res = q_ab.vector(vec_i);
  EXPECT_EQ(res.status, SUCCESS);
  real vec[3];
  vec_i(vec);

  EXPECT_EQ(s, static_cast<real>(-41));
  EXPECT_EQ(vec[0], static_cast<real>(-4));
  EXPECT_EQ(vec[1], static_cast<real>(9));
  EXPECT_EQ(vec[2], static_cast<real>(-20));
}
TEST(QuaternionTest, test_product_operator_1) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  Quaternion<real> q_ab;
  q_a.product(2, q_ab);

  real s = static_cast<real>(1);
  auto res = q_ab.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec_in;
  res = q_ab.vector(vec_in);
  EXPECT_EQ(res.status, SUCCESS);
  real vec[3];
  vec_in(vec);

  EXPECT_EQ(s, 4);
  EXPECT_EQ(vec[0], -4);
  EXPECT_EQ(vec[1], 6);
  EXPECT_EQ(vec[2], -8);
}
/*! @} */

/*! @{ Test determinant of Quaternion */
TEST(QuaternionTest, test_norm_squared_v1) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  real dv = static_cast<real>(156);
  q_a.norm_squared(dv);
  // 4 + 4 + 9 + 16
  EXPECT_EQ(dv, static_cast<real>(33));
}

/*! Test determinant of Quaternion */
TEST(QuaternionTest, test_norm_squared_v2) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  real dv = static_cast<real>(5);
  q_a.norm_squared(dv);
  // 4 + 4 + 9 + 16
  EXPECT_EQ(dv, static_cast<real>(33));
}
/*! @} */

/*! @{ Test magnitude of Quaternion */
TEST(QuaternionTest, test_magnitude) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  real dv = static_cast<real>(5);
  q_a.magnitude(dv);
  // 4 + 4 + 9 + 16
  EXPECT_EQ(dv, static_cast<real>(sqrt(33)));
}

/*! Test norm of Quaternion */
TEST(QuaternionTest, test_norm) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  real n = static_cast<real>(0);
  auto res = q_a.norm(n);
  EXPECT_EQ(res.status, SUCCESS);
  // 4 + 4 + 9 + 16
  EXPECT_EQ(n, static_cast<real>(sqrt(33)));
}

/*! @} */

/*! @{ Test exponent of Quaternion */
TEST(QuaternionTest, test_squared) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  Quaternion<real> q_a2;
  q_a.squared(q_a2);

  //
  real s = static_cast<real>(1);
  auto res = q_a2.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec_in;
  res = q_a2.vector(vec_in);
  EXPECT_EQ(res.status, SUCCESS);

  real vec[3];
  vec_in(vec);
  EXPECT_EQ(s, -25);
  EXPECT_EQ(vec[0], -8);
  EXPECT_EQ(vec[1], 12);
  EXPECT_EQ(vec[2], -16);
}

TEST(QuaternionTest, test_power) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  Quaternion<real> q_a2;
  q_a.power(2, q_a2);

  real s = static_cast<real>(1);
  auto res = q_a2.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec_in;
  res = q_a2.vector(vec_in);
  EXPECT_EQ(res.status, SUCCESS);
  real vec[3];
  vec_in(vec);

  EXPECT_EQ(s, -25); // might be -25
  EXPECT_EQ(vec[0], -8);
  EXPECT_EQ(vec[1], 12);
  EXPECT_EQ(vec[2], -16);
}

/*! @} */

/*! @{ Test inverse of Quaternion Vince 2011 -
 * Quaternions for Computer Graphics, p. 71
 */
TEST(QuaternionTest, test_inversed) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  Quaternion<real> inv;
  q_a.inversed(inv);

  real s = static_cast<real>(1);
  auto res = inv.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec_in;
  res = inv.vector(vec_in);
  EXPECT_EQ(res.status, SUCCESS);
  real vec[3];
  vec_in(vec);

  EXPECT_EQ(s, static_cast<real>(
                   static_cast<real>(1.0 / 33) * 2));
  EXPECT_EQ(vec[0], static_cast<real>(
                        static_cast<real>(1.0 / 33) * 2));
  EXPECT_EQ(vec[1], static_cast<real>(
                        static_cast<real>(1.0 / 33) * -3));
  EXPECT_EQ(vec[2], static_cast<real>(
                        static_cast<real>(1.0 / 33) * 4));
}
/*! @} */

/*! @{ Test normalization method of Quaternion */
TEST(QuaternionTest, test_normalized) {
  Quaternion<real> q_a(2, vecn::VecN<real, 3>({-2, 3, -4}));
  Quaternion<real> n_q;
  q_a.normalized(n_q);

  real s = static_cast<real>(1);
  auto res = n_q.scalar(s);
  EXPECT_EQ(res.status, SUCCESS);

  vecn::VecN<real, 3> vec_in;
  res = n_q.vector(vec_in);
  EXPECT_EQ(res.status, SUCCESS);

  EXPECT_EQ(s, static_cast<real>(
                   static_cast<real>(1.0 / sqrt(33)) * 2));
  real vec[3];
  vec_in(vec);
  EXPECT_EQ(vec[0],
            static_cast<real>(
                static_cast<real>(1.0 / sqrt(33)) * -2));
  EXPECT_EQ(vec[1],
            static_cast<real>(
                static_cast<real>(1.0 / sqrt(33)) * 3));
  EXPECT_EQ(vec[2],
            static_cast<real>(
                static_cast<real>(1.0 / sqrt(33)) * -4));
}
/*! @} */
