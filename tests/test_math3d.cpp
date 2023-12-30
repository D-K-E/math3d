#include <gtest/gtest.h>
#include <math3d/math3d.hpp>

typedef float real;
using namespace math3d::lu;
using namespace math3d::matn;
using namespace math3d::vecn;
using namespace math3d::quaternion;
using namespace math3d;

TEST(Math3DTest, test_invert_mat) {

  MatN<real, 3, 3> m;
  m.set_column(0, {25, 64, 144});
  m.set_column(1, {5, 8, 12});
  m.set_column(2, {1, 1, 1});
  MatN<real, 3, 3> X;
  invert_mat(m, X);
  real c1[3], c2[3], c3[3];
  X.get_column(0, c1);
  X.get_column(1, c2);
  X.get_column(2, c3);
  //
  EXPECT_EQ(std::round(c1[0] * 100.0f) / 100.0f, 0.05f);
  EXPECT_EQ(std::round(c1[1] * 100.0f) / 100.0f, -0.95f);
  EXPECT_EQ(std::round(c1[2] * 100.0f) / 100.0f, 4.57f);
  //
  EXPECT_EQ(std::round(c2[0] * 100.0f) / 100.0f, -0.08f);
  EXPECT_EQ(std::round(c2[1] * 100.0f) / 100.0f, 1.42f);
  EXPECT_EQ(std::round(c2[2] * 100.0f) / 100.0f, -5.0f);
  //
  EXPECT_EQ(std::round(c3[0] * 100.0f) / 100.0f, 0.04f);
  EXPECT_EQ(std::round(c3[1] * 100.0f) / 100.0f, -0.46f);
  EXPECT_EQ(std::round(c3[2] * 100.0f) / 100.0f, 1.43f);
}

TEST(Math3DTest, test_toRotMat3x3) {
  /*
0.35 + 0.2i + 0.3j + 0.1k

0.238095  0.190476  0.952381
0.72381  0.619048 -0.304762
-0.647619  0.761905  0.0095236
   */

  MatN<real, 3, 3> m;
  Quaternion<real> q(0.35f,
                     VecN<real, 3>({0.2f, 0.3f, 0.1f}));
  toRotMat3x3(q, m);
  real c1[3], c2[3], c3[3];
  m.get_column(0, c1);
  m.get_column(1, c2);
  m.get_column(2, c3);
  //
  EXPECT_EQ(std::round(c1[0] * 100.0f) / 100.0f, 0.24f);
  EXPECT_EQ(std::round(c1[1] * 100.0f) / 100.0f, 0.72f);
  EXPECT_EQ(std::round(c1[2] * 100.0f) / 100.0f, -0.65f);
  //
  EXPECT_EQ(std::round(c2[0] * 100.0f) / 100.0f, 0.19f);
  EXPECT_EQ(std::round(c2[1] * 100.0f) / 100.0f, 0.62f);
  EXPECT_EQ(std::round(c2[2] * 100.0f) / 100.0f, 0.76f);
  //
  EXPECT_EQ(std::round(c3[0] * 100.0f) / 100.0f, 0.95f);
  EXPECT_EQ(std::round(c3[1] * 100.0f) / 100.0f, -0.3f);
  EXPECT_EQ(std::round(c3[2] * 100.0f) / 100.0f, 0.01f);
}

TEST(Math3DTest, test_debug_check_m_true) {
  VecN<real, 2> v;
  std::size_t vsize = 5;
  OpResult res;
  CHECK_MATH3D(v(vsize), res);
  EXPECT_EQ(res.success, true);
  std::string cname = "v(vsize)";
  std::string rcname = res.call_name;
  EXPECT_EQ(rcname, cname);
}

TEST(MatnTest, test_divide_scalar_value_false_check_macro) {

  real mv[] = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  auto r = CHECK(m.divide(0, out));
  EXPECT_EQ(r, false);
  //
}

TEST(Math3DTest, test_debug_check_true) {
  VecN<real, 2> v;
  std::size_t vsize = 5;
  bool vflag = CHECK(v(vsize));
  EXPECT_EQ(vflag, true);
}
