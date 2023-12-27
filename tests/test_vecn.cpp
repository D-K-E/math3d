// test file for quaternion
#include <gtest/gtest.h>
#include <math3d/vecn.hpp>

/*! @{
 */

typedef float real;
using namespace math3d::vecn;
using namespace math3d;

/*! @{ testing constructors for the vector
 */

TEST(VecNTest, test_empty_constructor) {
  VecN<real, 2> vec;
  std::size_t vsize;
  vec(vsize);
  EXPECT_EQ(vsize, static_cast<std::size_t>(2));
}

TEST(VecNTest, test_vector_constructor) {
  real tvec[5];
  real v1 = static_cast<real>(-21);
  real v2 = static_cast<real>(1);
  real v3 = static_cast<real>(2);
  real v4 = static_cast<real>(1.53);
  real v5 = static_cast<real>(3);
  tvec[0] = v1;
  tvec[1] = v2;
  tvec[2] = v3;
  tvec[3] = v4;
  tvec[4] = v5;
  VecN<real, 5> v(tvec);
  real t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0;
  // test size
  std::size_t vsize = 6156;
  v(vsize);
  EXPECT_EQ(vsize, 5);

  // test values
  // first value
  EXPECT_EQ(v(0, t1).status, SUCCESS);
  EXPECT_EQ(t1, v1);

  // second value
  EXPECT_EQ(v(1, t2).status, SUCCESS);
  EXPECT_EQ(t2, v2);
  // third value

  EXPECT_EQ(v(2, t3).status, SUCCESS);
  EXPECT_EQ(t3, v3);

  // fourth value
  EXPECT_EQ(v(3, t4).status, SUCCESS);
  EXPECT_EQ(t4, v4);

  // fifth value
  EXPECT_EQ(v(4, t5).status, SUCCESS);
  EXPECT_EQ(t5, v5);

  // sixth value
  EXPECT_EQ(v(5, t6).status, INDEX_ERROR);
}

TEST(VecNTest, test_dimension_constructor) {
  VecN<real, 5> v(5);
  real r1 = static_cast<real>(5);
  // test values
  real t1 = 0;
  // first value
  EXPECT_EQ(v(0, t1).status, SUCCESS);
  EXPECT_EQ(t1, r1);

  // second value
  EXPECT_EQ(v(1, t1).status, SUCCESS);
  EXPECT_EQ(t1, r1);
  // third value

  EXPECT_EQ(v(2, t1).status, SUCCESS);
  EXPECT_EQ(t1, r1);

  // fourth value
  EXPECT_EQ(v(3, t1).status, SUCCESS);
  EXPECT_EQ(t1, r1);

  // fifth value

  EXPECT_EQ(v(make_vecn_cell(t1, 4)).status, SUCCESS);
  EXPECT_EQ(t1, r1);

  // sixth value
  EXPECT_EQ(v(5, t1).status, INDEX_ERROR);
  EXPECT_EQ(v(make_vecn_cell(t1, 5)).status, INDEX_ERROR);
}

TEST(VecNTest, test_vecn_flag_fn_name) {
  VecN<real, 2> v;
  std::size_t vsize = 5;
  auto vflag = v(vsize);
  const std::string fname = vflag.fn_name;
  const std::string compval = "operator()";
  EXPECT_EQ(fname, compval);
}
TEST(VecNTest, test_debug_check_true) {
  VecN<real, 2> v;
  std::size_t vsize = 5;
  bool vflag = CHECK(v(vsize));
  EXPECT_EQ(vflag, true);
}

TEST(VecNTest, test_debug_check_m_true) {
  VecN<real, 2> v;
  std::size_t vsize = 5;
  OpResult res = v(vsize);
  CHECK_MATH3D(v(vsize), res);
  EXPECT_EQ(res.success, true);
  std::string cname = "v(vsize)";
  std::string rcname = res.call_name;
  EXPECT_EQ(rcname, cname);
}

/*! @} */

/*! @{ testing base dimensions for the vector
 */
TEST(VecNTest, test_to_base_dimensions_vector) {
  //
  VecN<real, 6> vout{};
  base<real, 1, 6>(vout);
  real out[6];
  vout(out);
  EXPECT_EQ(out[1], static_cast<real>(1));
}
TEST(VecNTest, test_to_base_dimensions_vecn) {
  //
  VecN<real, 6> vout;
  base<real, 1, 6>(vout);
  real tout = 0;
  EXPECT_EQ(vout(1, tout).status, SUCCESS);
  EXPECT_EQ(tout, static_cast<real>(1));
}

/*! @} */

/*! @{ testing summation for the vector
 */
TEST(VecNTest, test_add_scalar_vector) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 3;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);
  real out[5];
  auto result = v.add(1, out).status;
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(out[0], 2);
  EXPECT_EQ(out[1], 3);
  EXPECT_EQ(out[2], 4);
  EXPECT_EQ(out[3], 3);
  EXPECT_EQ(out[4], 2);
}
TEST(VecNTest, test_add_scalar_vecn) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 3;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);
  VecN<real, 5> out;
  auto result = v.add(1, out).status;
  EXPECT_EQ(result, SUCCESS);
  real t = static_cast<real>(0);
  out(0, t);
  EXPECT_EQ(t, 2);
  out(1, t);
  EXPECT_EQ(t, 3);
  out(2, t);
  EXPECT_EQ(t, 4);
  out(3, t);
  EXPECT_EQ(t, 3);
  out(4, t);
  EXPECT_EQ(t, 2);
}

TEST(VecNTest, test_add_vector_to_vector_t) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 3;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);
  real out[5];
  real avec[5];
  avec[0] = 1;
  avec[1] = 2;
  avec[2] = 2;
  avec[3] = 2;
  avec[4] = 1;
  auto result = v.add(avec, out).status;
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(out[0], static_cast<real>(2));
  EXPECT_EQ(out[1], static_cast<real>(4));
  EXPECT_EQ(out[2], static_cast<real>(5));
  EXPECT_EQ(out[3], static_cast<real>(4));
  EXPECT_EQ(out[4], static_cast<real>(2));
}
TEST(VecNTest, test_add_vecn_to_vecn_t) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 3;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);

  // output vector
  VecN<real, 5> out;

  real avec[5];
  // prepare vector to add
  avec[0] = 1;
  avec[1] = 2;
  avec[2] = 2;
  avec[3] = 2;
  avec[4] = 1;

  VecN<real, 5> avecn(avec);
  auto result = v.add(avecn, out).status;

  real t = static_cast<real>(0);
  out(0, t);
  //
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(t, static_cast<real>(2));

  out(1, t);
  EXPECT_EQ(t, static_cast<real>(4));

  out(2, t);
  EXPECT_EQ(t, static_cast<real>(5));

  out(3, t);
  EXPECT_EQ(t, static_cast<real>(4));

  out(4, t);
  EXPECT_EQ(t, static_cast<real>(2));
}
/*! @} */

/*! @{ testing subtraction for the vector
 */
TEST(VecNTest, test_subtract_scalar_vector) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 3;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);
  real out[5];
  auto result = v.subtract(1, out).status;
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(out[0], 0);
  EXPECT_EQ(out[1], 1);
  EXPECT_EQ(out[2], 2);
  EXPECT_EQ(out[3], 1);
  EXPECT_EQ(out[4], 0);
}
TEST(VecNTest, test_subtract_scalar_vecn) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 3;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);
  VecN<real, 5> out;
  auto result = v.subtract(1, out).status;
  EXPECT_EQ(result, SUCCESS);
  real t = static_cast<real>(0);
  out(0, t);
  EXPECT_EQ(t, 0);
  out(1, t);
  EXPECT_EQ(t, 1);
  out(2, t);
  EXPECT_EQ(t, 2);
  out(3, t);
  EXPECT_EQ(t, 1);
  out(4, t);
  EXPECT_EQ(t, 0);
}
TEST(VecNTest, test_subtract_vector_to_vector_t) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 3;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);
  real out[5];
  real avec[5];
  avec[0] = 1;
  avec[1] = 2;
  avec[2] = 2;
  avec[3] = 2;
  avec[4] = 1;
  auto result = v.subtract(avec, out).status;
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(out[0], static_cast<real>(0));
  EXPECT_EQ(out[1], static_cast<real>(0));
  EXPECT_EQ(out[2], static_cast<real>(1));
  EXPECT_EQ(out[3], static_cast<real>(0));
  EXPECT_EQ(out[4], static_cast<real>(0));
}
TEST(VecNTest, test_subtract_vecn_to_vecn_t) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 3;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);

  // output vector
  VecN<real, 5> out;

  real avec[5];
  // prepare vector to subtract
  avec[0] = 1;
  avec[1] = 2;
  avec[2] = 2;
  avec[3] = 2;
  avec[4] = 1;

  VecN<real, 5> avecn(avec);
  auto result = v.subtract(avecn, out).status;

  real t = static_cast<real>(0);
  out(0, t);
  //
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(t, static_cast<real>(0));

  out(1, t);
  EXPECT_EQ(t, static_cast<real>(0));

  out(2, t);
  EXPECT_EQ(t, static_cast<real>(1));

  out(3, t);
  EXPECT_EQ(t, static_cast<real>(0));

  out(4, t);
  EXPECT_EQ(t, static_cast<real>(0));
}

/*! @} */

/*! @{ testing element-wise multiplication for the vector
 */
TEST(VecNTest, test_multiply_scalar_vector) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 3;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);
  real out[5];
  auto result = v.multiply(2, out).status;
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(out[0], 2);
  EXPECT_EQ(out[1], 4);
  EXPECT_EQ(out[2], 6);
  EXPECT_EQ(out[3], 4);
  EXPECT_EQ(out[4], 2);
}
TEST(VecNTest, test_multiply_scalar_vecn) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 3;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);
  VecN<real, 5> out;
  auto result = v.multiply(2, out).status;
  EXPECT_EQ(result, SUCCESS);
  real t = static_cast<real>(0);
  out(0, t);
  EXPECT_EQ(t, 2);
  out(1, t);
  EXPECT_EQ(t, 4);
  out(2, t);
  EXPECT_EQ(t, 6);
  out(3, t);
  EXPECT_EQ(t, 4);
  out(4, t);
  EXPECT_EQ(t, 2);
}
TEST(VecNTest, test_multiply_vector_to_vector_t) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 3;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);
  real out[5];
  real avec[5];
  avec[0] = 1;
  avec[1] = 2;
  avec[2] = 2;
  avec[3] = 2;
  avec[4] = 1;
  auto result = v.multiply(avec, out).status;
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(out[0], static_cast<real>(1));
  EXPECT_EQ(out[1], static_cast<real>(4));
  EXPECT_EQ(out[2], static_cast<real>(6));
  EXPECT_EQ(out[3], static_cast<real>(4));
  EXPECT_EQ(out[4], static_cast<real>(1));
}
TEST(VecNTest, test_multiply_vecn_to_vecn_t) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 3;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);

  // output vector
  VecN<real, 5> out;

  real avec[5];
  // prepare vector to multiply
  avec[0] = 1;
  avec[1] = 2;
  avec[2] = 2;
  avec[3] = 2;
  avec[4] = 1;

  VecN<real, 5> avecn(avec);
  auto result = v.multiply(avecn, out).status;

  real t = static_cast<real>(0);
  out(0, t);
  //
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(t, static_cast<real>(1));

  out(1, t);
  EXPECT_EQ(t, static_cast<real>(4));

  out(2, t);
  EXPECT_EQ(t, static_cast<real>(6));

  out(3, t);
  EXPECT_EQ(t, static_cast<real>(4));

  out(4, t);
  EXPECT_EQ(t, static_cast<real>(1));
}

/*! @} */
/*! @{ testing element-wise division for the vector
 */
TEST(VecNTest, test_divide_scalar_vector) {
  real inv[5];
  inv[0] = 2;
  inv[1] = 4;
  inv[2] = 6;
  inv[3] = 2;
  inv[4] = 0;
  VecN<real, 5> v(inv);
  real out[5];
  auto result = v.divide(2, out).status;
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(out[0], 1);
  EXPECT_EQ(out[1], 2);
  EXPECT_EQ(out[2], 3);
  EXPECT_EQ(out[3], 1);
  EXPECT_EQ(out[4], 0);
}
TEST(VecNTest, test_divide_scalar_vector_zero) {
  real inv[5];
  inv[0] = 2;
  inv[1] = 4;
  inv[2] = 6;
  inv[3] = 2;
  inv[4] = 0;
  VecN<real, 5> v(inv);
  real out[5];
  auto result = v.divide(0, out).status;
  EXPECT_EQ(result, ARG_ERROR);
}
TEST(VecNTest, test_divide_scalar_vecn) {
  real inv[5];
  inv[0] = 2;
  inv[1] = 4;
  inv[2] = 6;
  inv[3] = 2;
  inv[4] = 0;
  VecN<real, 5> v(inv);
  VecN<real, 5> out;
  auto result = v.divide(2, out).status;
  EXPECT_EQ(result, SUCCESS);
  real t = static_cast<real>(0);
  out(0, t);
  EXPECT_EQ(t, 1);
  out(1, t);
  EXPECT_EQ(t, 2);
  out(2, t);
  EXPECT_EQ(t, 3);
  out(3, t);
  EXPECT_EQ(t, 1);
  out(4, t);
  EXPECT_EQ(t, 0);
}
TEST(VecNTest, test_divide_scalar_vecn_zero) {
  real inv[5];
  inv[0] = 2;
  inv[1] = 4;
  inv[2] = 6;
  inv[3] = 2;
  inv[4] = 0;
  VecN<real, 5> v(inv);
  VecN<real, 5> out;
  auto result = v.divide(0, out).status;
  EXPECT_EQ(result, ARG_ERROR);
}
TEST(VecNTest, test_divide_vector_to_vector_t) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 4;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);
  real out[5];
  real avec[5];
  avec[0] = 1;
  avec[1] = 2;
  avec[2] = 2;
  avec[3] = 2;
  avec[4] = 1;
  auto result = v.divide(avec, out).status;
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(out[0], static_cast<real>(1));
  EXPECT_EQ(out[1], static_cast<real>(1));
  EXPECT_EQ(out[2], static_cast<real>(2));
  EXPECT_EQ(out[3], static_cast<real>(1));
  EXPECT_EQ(out[4], static_cast<real>(1));
}
TEST(VecNTest, test_divide_vector_to_vector_zero) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 4;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);
  real out[5];
  real avec[5];
  avec[0] = 1;
  avec[1] = 2;
  avec[2] = 0;
  avec[3] = 2;
  avec[4] = 1;
  auto result = v.divide(avec, out).status;
  EXPECT_EQ(result, ARG_ERROR);
}
TEST(VecNTest, test_divide_vecn_to_vecn_t) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 4;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);

  // output vector
  VecN<real, 5> out;

  real avec[5];
  // prepare vector to divide
  avec[0] = 1;
  avec[1] = 2;
  avec[2] = 2;
  avec[3] = 2;
  avec[4] = 1;

  VecN<real, 5> avecn(avec);
  auto result = v.divide(avecn, out).status;

  real t = static_cast<real>(0);
  out(0, t);
  //
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(t, static_cast<real>(1));

  out(1, t);
  EXPECT_EQ(t, static_cast<real>(1));

  out(2, t);
  EXPECT_EQ(t, static_cast<real>(2));

  out(3, t);
  EXPECT_EQ(t, static_cast<real>(1));

  out(4, t);
  EXPECT_EQ(t, static_cast<real>(1));
}
TEST(VecNTest, test_divide_vecn_to_vecn_zero) {
  real inv[5];
  inv[0] = 1;
  inv[1] = 2;
  inv[2] = 4;
  inv[3] = 2;
  inv[4] = 1;
  VecN<real, 5> v(inv);

  // output vector
  VecN<real, 5> out;

  real avec[5];
  // prepare vector to divide
  avec[0] = 1;
  avec[1] = 2;
  avec[2] = -10;
  avec[3] = 2;
  avec[4] = 0;

  VecN<real, 5> avecn(avec);
  auto result = v.divide(avecn, out).status;

  real t = static_cast<real>(0);
  out(0, t);
  //
  EXPECT_EQ(result, ARG_ERROR);
}

/*! @} */

/*! @{ testing element-wise division for the vector
 */

TEST(VecNTest, test_dot_scalar_t) {
  real inv[5];
  inv[0] = 1; // 2
  inv[1] = 2; // 4
  inv[2] = 3; // 6
  inv[3] = 2; // 4
  inv[4] = 1; // 2
  // sum == 18
  VecN<real, 5> v(inv);
  real out;
  auto result = v.dot(static_cast<real>(2), out).status;
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(out, static_cast<real>(18));
}
TEST(VecNTest, test_dot_scalar_zero) {
  real inv[5];
  inv[0] = 1; // 2
  inv[1] = 2; // 4
  inv[2] = 3; // 6
  inv[3] = 2; // 4
  inv[4] = 1; // 2
  // sum == 18
  VecN<real, 5> v(inv);
  real out;
  auto result = v.dot(static_cast<real>(0), out).status;
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(out, static_cast<real>(0));
}
TEST(VecNTest, test_dot_vector_t) {
  real inv[5];
  inv[0] = 1; // 1
  inv[1] = 2; // 4
  inv[2] = 3; // 3
  inv[3] = 2; // 0
  inv[4] = 1; // 0
  // sum == 8
  real avec[5];
  VecN<real, 5> v(inv);
  avec[0] = 1;
  avec[1] = 2;
  avec[2] = 1;
  avec[3] = 0;
  avec[4] = 0;
  real out = static_cast<real>(0);
  auto result = v.dot(avec, out).status;
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(out, static_cast<real>(8));
}
TEST(VecNTest, test_dot_vecn_t) {
  real inv[5];
  inv[0] = 1; // 1
  inv[1] = 2; // 4
  inv[2] = 3; // 3
  inv[3] = 2; // 0
  inv[4] = 1; // 0
  // sum == 8
  real avec[5];
  VecN<real, 5> v(inv);
  avec[0] = 1;
  avec[1] = 2;
  avec[2] = 1;
  avec[3] = 0;
  avec[4] = 0;
  VecN<real, 5> av(avec);
  real out = static_cast<real>(0);
  auto result = v.dot(av, out).status;
  EXPECT_EQ(result, SUCCESS);
  EXPECT_EQ(out, static_cast<real>(8));
}

/*! @} */
