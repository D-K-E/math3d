#include <gtest/gtest.h>
#include <math3d/lu.hpp>

typedef float real;
using namespace math3d::lu;
using namespace math3d::matn;
using namespace math3d::vecn;
using namespace math3d;

/**
  all test values are
  from:https://math.libretexts.org/Bookshelves/Linear_Algebra/Introduction_to_Matrix_Algebra_(Kaw)/01%3A_Chapters/1.07%3A_LU_Decomposition_Method_for_Solving_Simultaneous_Linear_Equations
 */

TEST(LUTest, test_L) {
  //

  MatN<real, 3, 3> m;
  m.set_column(0, {25, 64, 144});
  m.set_column(1, {5, 8, 12});
  m.set_column(2, {1, 1, 1});
  LUdecomp lu_d(m);
  MatN<real, 3, 3> L;
  lu_d.lower(L);
  real row1[3], row2[3], row3[3];
  L.get_row(0, row1);
  L.get_row(1, row2);
  L.get_row(2, row3);
  //
  EXPECT_EQ(row1[0], 1);
  EXPECT_EQ(row1[1], 0);
  EXPECT_EQ(row1[2], 0);
  //
  EXPECT_EQ(std::round(row2[0] * 100.0f) / 100.0f, 2.56f);
  EXPECT_EQ(std::round(row2[1] * 100.0f) / 100.0f, 1.0f);
  EXPECT_EQ(std::round(row2[2] * 100.0f) / 100.0f, 0);
  //
  EXPECT_EQ(std::round(row3[0] * 100.0f) / 100.0f, 5.76f);
  EXPECT_EQ(std::round(row3[1] * 100.0f) / 100.0f, 3.5f);
  EXPECT_EQ(std::round(row3[2] * 100.0f) / 100.0f, 1.0);
}
TEST(LUTest, test_U) {

  MatN<real, 3, 3> m;
  m.set_column(0, {25, 64, 144});
  m.set_column(1, {5, 8, 12});
  m.set_column(2, {1, 1, 1});
  LUdecomp lu_d(m);
  MatN<real, 3, 3> U;
  lu_d.upper(U);
  real row1[3], row2[3], row3[3];
  U.get_row(0, row1);
  U.get_row(1, row2);
  U.get_row(2, row3);
  //
  EXPECT_EQ(row1[0], 25);
  EXPECT_EQ(row1[1], 5);
  EXPECT_EQ(row1[2], 1);
  //
  EXPECT_EQ(std::round(row2[0] * 100.0f) / 100.0f, 0);
  EXPECT_EQ(std::round(row2[1] * 100.0f) / 100.0f, -4.8f);
  EXPECT_EQ(std::round(row2[2] * 100.0f) / 100.0f, -1.56f);
  //
  EXPECT_EQ(std::round(row3[0] * 100.0f) / 100.0f, 0);
  EXPECT_EQ(std::round(row3[1] * 100.0f) / 100.0f, 0);
  EXPECT_EQ(std::round(row3[2] * 100.0f) / 100.0f, 0.7f);
}
TEST(LUTest, test_solve_forward) {

  MatN<real, 3, 3> m;
  m.set_column(0, {25, 64, 144});
  m.set_column(1, {5, 8, 12});
  m.set_column(2, {1, 1, 1});
  VecN<real, 3> b({106.8f, 177.2f, 279.2f});
  LUdecomp lu_d(m);
  VecN<real, 3> x;
  lu_d.solve_forward(b, x);
  real row[3];
  x(row);
  //
  EXPECT_EQ(std::round(row[0] * 100.0f) / 100.0f, 106.8f);
  EXPECT_EQ(std::round(row[1] * 100.0f) / 100.0f, -96.21f);
  EXPECT_EQ(std::round(row[2] * 100.0f) / 100.0f, 0.76f);
  //
}
TEST(LUTest, test_solve_backward) {

  MatN<real, 3, 3> m;
  m.set_column(0, {25, 64, 144});
  m.set_column(1, {5, 8, 12});
  m.set_column(2, {1, 1, 1});

  VecN<real, 3> b({106.8f, -96.21f, 0.76f});
  LUdecomp<real, 3> lu_d(m);
  VecN<real, 3> x;
  lu_d.solve_backward(b, x);
  real row[3];
  x(row);
  //
  EXPECT_EQ(std::round(row[0] * 100.0f) / 100.0f, 0.29f);
  EXPECT_EQ(std::round(row[1] * 100.0f) / 100.0f, 19.69f);
  EXPECT_EQ(std::round(row[2] * 100.0f) / 100.0f, 1.09f);
  //
}

TEST(LUTest, test_solve) {

  MatN<real, 3, 3> m;
  m.set_column(0, {25, 64, 144});
  m.set_column(1, {5, 8, 12});
  m.set_column(2, {1, 1, 1});
  VecN<real, 3> b({106.8f, 177.2f, 279.2f});
  LUdecomp lu_d(m);
  VecN<real, 3> x;
  lu_d.solve(b, x);
  real row[3];
  x(row);
  //
  EXPECT_EQ(std::round(row[0] * 100.0f) / 100.0f, 0.29f);
  EXPECT_EQ(std::round(row[1] * 100.0f) / 100.0f, 19.69f);
  EXPECT_EQ(std::round(row[2] * 100.0f) / 100.0f, 1.09f);
  //
}

TEST(LUTest, test_solve_mat) {

  MatN<real, 3, 3> m;
  m.set_column(0, {25, 64, 144});
  m.set_column(1, {5, 8, 12});
  m.set_column(2, {1, 1, 1});
  MatN<real, 3, 3> B;
  MatN<real, 3, 3>::identity<3>(B);
  LUdecomp lu_d(m);
  MatN<real, 3, 3> X;
  lu_d.solve_mat<3>(B, X);
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
