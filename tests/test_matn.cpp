#include <math3d/matn.hpp>

//
#include <chrono>
#include <gtest/gtest.h>
#include <ratio>

/*! @{
 */

typedef float real;
using namespace math3d::matn;
using namespace math3d;

/*! @{ Test constructors of the matrix
 */
TEST(MatnTest, test_empty_constructor) {
  MatN<real, 0, 0> m;
  std::size_t vsize = 0;
  auto r = m(vsize);
  EXPECT_EQ(vsize, static_cast<std::size_t>(0));
  EXPECT_EQ(r.status, SUCCESS);
}
TEST(MatnTest, test_double_vec_constructor) {

  std::vector<std::vector<real>> mv;
  std::vector<real> r1 = {0, 1, 2};
  std::vector<real> r2 = {0, 2, 4};
  mv.push_back(r1);
  mv.push_back(r2);

  MatN<real, 2, 3> m(mv);
  std::size_t vsize = 0;

  auto r = m(vsize);
  EXPECT_EQ(r.status, SUCCESS);
  EXPECT_EQ(vsize, static_cast<std::size_t>(6));

  std::size_t rows, cols;
  r = m(rows, cols);
  EXPECT_EQ(cols, static_cast<std::size_t>(3));
  EXPECT_EQ(r.status, SUCCESS);

  EXPECT_EQ(rows, static_cast<std::size_t>(2));
  EXPECT_EQ(r.status, SUCCESS);

  real arr[6];
  r = m(arr);
  EXPECT_EQ(r.status, SUCCESS);
  EXPECT_DOUBLE_EQ(arr[0], static_cast<real>(0));
  EXPECT_DOUBLE_EQ(arr[1], static_cast<real>(1));
  EXPECT_DOUBLE_EQ(arr[2], static_cast<real>(2));
  EXPECT_DOUBLE_EQ(arr[3], static_cast<real>(0));
  EXPECT_DOUBLE_EQ(arr[4], static_cast<real>(2));
  EXPECT_DOUBLE_EQ(arr[5], static_cast<real>(4));
}
TEST(MatnTest, test_double_vec_constructor_bigger_row) {

  std::vector<std::vector<real>> mv;
  std::vector<real> r1 = {0, 1, 2};
  std::vector<real> r2 = {0, 2, 4};
  std::vector<real> r3 = {1, 2, 4};
  mv.push_back(r1);
  mv.push_back(r2);
  mv.push_back(r3);

  MatN<real, 2, 3> m(mv);
  std::size_t vsize = 0;

  auto r = m(vsize);
  EXPECT_EQ(vsize, static_cast<std::size_t>(6));
  EXPECT_EQ(r.status, SUCCESS);

  std::size_t rows, cols;
  r = m(rows, cols);
  EXPECT_EQ(cols, static_cast<std::size_t>(3));
  EXPECT_EQ(r.status, SUCCESS);

  EXPECT_EQ(rows, static_cast<std::size_t>(2));

  real arr[6];
  r = m(arr);
  EXPECT_EQ(r.status, SUCCESS);
  EXPECT_DOUBLE_EQ(arr[0], static_cast<real>(0));
  EXPECT_DOUBLE_EQ(arr[1], static_cast<real>(1));
  EXPECT_DOUBLE_EQ(arr[2], static_cast<real>(2));
  EXPECT_DOUBLE_EQ(arr[3], static_cast<real>(0));
  EXPECT_DOUBLE_EQ(arr[4], static_cast<real>(2));
  EXPECT_DOUBLE_EQ(arr[5], static_cast<real>(4));
}
TEST(MatnTest, test_single_vec_constructor) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  std::size_t vsize = 0;

  auto r = m(vsize);
  EXPECT_EQ(vsize, static_cast<std::size_t>(6));
  EXPECT_EQ(r.status, SUCCESS);

  std::size_t rows, cols;
  r = m(rows, cols);
  EXPECT_EQ(cols, static_cast<std::size_t>(3));
  EXPECT_EQ(r.status, SUCCESS);

  EXPECT_EQ(rows, static_cast<std::size_t>(2));
  EXPECT_EQ(r.status, SUCCESS);

  real arr[6];
  r = m(arr);
  EXPECT_EQ(r.status, SUCCESS);
  EXPECT_DOUBLE_EQ(arr[0], static_cast<real>(0));
  EXPECT_DOUBLE_EQ(arr[1], static_cast<real>(1));
  EXPECT_DOUBLE_EQ(arr[2], static_cast<real>(2));
  EXPECT_DOUBLE_EQ(arr[3], static_cast<real>(0));
  EXPECT_DOUBLE_EQ(arr[4], static_cast<real>(2));
  EXPECT_DOUBLE_EQ(arr[5], static_cast<real>(4));
}
TEST(MatnTest, test_single_vec_constructor_bigger) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4, 56, 35};

  MatN<real, 2, 3> m(mv);
  std::size_t vsize = 0;

  auto r = m(vsize);
  EXPECT_EQ(vsize, static_cast<std::size_t>(6));
  EXPECT_EQ(r.status, SUCCESS);

  std::size_t rows, cols;
  r = m(rows, cols);
  EXPECT_EQ(cols, static_cast<std::size_t>(3));
  EXPECT_EQ(r.status, SUCCESS);

  EXPECT_EQ(rows, static_cast<std::size_t>(2));

  real arr[6];
  r = m(arr);
  EXPECT_EQ(r.status, SUCCESS);
  EXPECT_DOUBLE_EQ(arr[0], static_cast<real>(0));
  EXPECT_DOUBLE_EQ(arr[1], static_cast<real>(1));
  EXPECT_DOUBLE_EQ(arr[2], static_cast<real>(2));
  EXPECT_DOUBLE_EQ(arr[3], static_cast<real>(0));
  EXPECT_DOUBLE_EQ(arr[4], static_cast<real>(2));
  EXPECT_DOUBLE_EQ(arr[5], static_cast<real>(4));
}
TEST(MatnTest, test_single_vec_constructor_smaller) {

  std::vector<real> mv = {0, 1, 2};

  MatN<real, 2, 3> m(mv);
  std::size_t vsize = 0;

  auto r = m(vsize);
  EXPECT_EQ(vsize, static_cast<std::size_t>(6));
  EXPECT_EQ(r.status, SUCCESS);

  std::size_t rows, cols;
  r = m(rows, cols);
  EXPECT_EQ(cols, static_cast<std::size_t>(3));
  EXPECT_EQ(r.status, SUCCESS);

  EXPECT_EQ(rows, static_cast<std::size_t>(2));

  real arr[6];
  r = m(arr);
  EXPECT_EQ(r.status, SUCCESS);
  EXPECT_DOUBLE_EQ(arr[0], static_cast<real>(0));
  EXPECT_DOUBLE_EQ(arr[1], static_cast<real>(1));
  EXPECT_DOUBLE_EQ(arr[2], static_cast<real>(2));
  EXPECT_DOUBLE_EQ(arr[3], static_cast<real>(0));
  EXPECT_DOUBLE_EQ(arr[4], static_cast<real>(0));
  EXPECT_DOUBLE_EQ(arr[5], static_cast<real>(0));
}
/*! @}
 */

TEST(MatnTest, test_static_from_row_cols) {

  MatN<real, 2, 3> m;
  auto r = MatN<real>::from_row_cols(m);
  EXPECT_EQ(r.status, SUCCESS);

  std::size_t vsize = 0;

  r = m(vsize);
  EXPECT_EQ(vsize, static_cast<std::size_t>(6));
  EXPECT_EQ(r.status, SUCCESS);

  std::size_t rows, cols;
  r = m(rows, cols);
  EXPECT_EQ(cols, static_cast<std::size_t>(3));
  EXPECT_EQ(r.status, SUCCESS);

  EXPECT_EQ(rows, static_cast<std::size_t>(2));

  real arr[6];
  r = m(arr);
  EXPECT_EQ(r.status, SUCCESS);
  EXPECT_EQ(arr[0], static_cast<real>(0));
  EXPECT_EQ(arr[1], static_cast<real>(0));
  EXPECT_EQ(arr[2], static_cast<real>(0));
  EXPECT_EQ(arr[3], static_cast<real>(0));
  EXPECT_EQ(arr[4], static_cast<real>(0));
  EXPECT_EQ(arr[5], static_cast<real>(0));
}
TEST(MatnTest, test_static_from_row_cols_fill) {

  MatN<real, 2, 3> m;
  MatN<real, 10, 10>::from_row_cols(20, m);
  std::size_t vsize = 0;

  m(vsize);
  EXPECT_EQ(vsize, static_cast<std::size_t>(6));

  std::size_t rows, cols;
  m(rows, cols);
  EXPECT_EQ(cols, static_cast<std::size_t>(3));

  EXPECT_EQ(rows, static_cast<std::size_t>(2));

  real out = 0;
  auto r = m(0, out);
  EXPECT_EQ(out, static_cast<real>(20));
  EXPECT_EQ(r.status, SUCCESS);
  real arr[6];
  r = m(arr);
  EXPECT_EQ(r.status, SUCCESS);
  EXPECT_EQ(arr[0], static_cast<real>(20));
  EXPECT_EQ(arr[1], static_cast<real>(20));
  EXPECT_EQ(arr[2], static_cast<real>(20));
  EXPECT_EQ(arr[3], static_cast<real>(20));
  EXPECT_EQ(arr[4], static_cast<real>(20));
  EXPECT_EQ(arr[5], static_cast<real>(20));
}
TEST(MatnTest, test_fill) {

  real vsize = static_cast<real>(2);
  MatN<real, 2, 3> m(vsize);

  real out = 0;
  auto res = m(0, out);
  EXPECT_EQ(res.status, SUCCESS);
  EXPECT_EQ(out, static_cast<real>(2));
  res = m(3, out);
  EXPECT_EQ(out, static_cast<real>(2));
  EXPECT_EQ(res.status, SUCCESS);
}
TEST(MatnTest, test_transpose) {

  MatN<real, 2, 3> m;
  MatN<real, 10>::from_row_cols(m);
  MatN<real, 3, 2> out;
  MatN<real, 10>::from_row_cols(out);

  auto res = m.transpose(out);
  EXPECT_EQ(res.status, SUCCESS);
  std::size_t rows, cols;
  //
  res = out(rows, cols);
  EXPECT_EQ(res.status, SUCCESS);
  EXPECT_EQ(cols, 2);
  //
  EXPECT_EQ(res.status, SUCCESS);
  EXPECT_EQ(rows, 3);
  //
}
TEST(MatnTest, test_col_nb) {

  MatN<real, 2, 3> m;
  MatN<real, 1>::from_row_cols(m);
  std::size_t rows, cnb;
  auto r = m(rows, cnb);
  EXPECT_EQ(r.status, SUCCESS);
  EXPECT_EQ(cnb, static_cast<std::size_t>(3));
  //
}
TEST(MatnTest, test_row_nb) {
  MatN<real, 2, 3> m;
  MatN<real, 1>::from_row_cols(m);
  std::size_t rnb, cols;
  auto r = m(rnb, cols);
  EXPECT_EQ(r.status, SUCCESS);
  EXPECT_EQ(rnb, static_cast<std::size_t>(2));
  //
}
TEST(MatnTest, test_get_size) {
  MatN<real, 2, 3> m;
  MatN<real, 1>::from_row_cols(m);
  std::size_t snb = 0;
  auto r = m(snb);
  EXPECT_EQ(r.status, SUCCESS);
  EXPECT_EQ(snb, static_cast<std::size_t>(6));
  //
}
TEST(MatnTest, test_add_scalar_value) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  auto r = m.add(2, out);
  EXPECT_EQ(r.status, SUCCESS);
  real tval = static_cast<real>(0);
  EXPECT_EQ(SUCCESS, (out(0, tval)).status);
  EXPECT_DOUBLE_EQ(tval, 2);
  EXPECT_EQ(SUCCESS, (out(1, tval)).status);
  EXPECT_DOUBLE_EQ(tval, 3);
  EXPECT_EQ(SUCCESS, out(2, tval).status);
  EXPECT_DOUBLE_EQ(tval, 4);
  EXPECT_EQ(SUCCESS, out(3, tval).status);
  EXPECT_DOUBLE_EQ(tval, 2);
  EXPECT_EQ(SUCCESS, out(4, tval).status);
  EXPECT_DOUBLE_EQ(tval, 4);
  EXPECT_EQ(SUCCESS, out(5, tval).status);
  EXPECT_DOUBLE_EQ(tval, 6);
  //
}
TEST(MatnTest, test_add_mat_value) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  MatN<real, 2, 3> tout(static_cast<real>(2));
  // filled matrix
  auto r = m.add(tout, out);
  EXPECT_EQ(r.status, SUCCESS);
  real tval = static_cast<real>(0);
  EXPECT_EQ(SUCCESS, out(0, tval).status);
  EXPECT_DOUBLE_EQ(tval, 2);
  EXPECT_EQ(SUCCESS, out(1, tval).status);
  EXPECT_DOUBLE_EQ(tval, 3);
  EXPECT_EQ(SUCCESS, out(2, tval).status);
  EXPECT_DOUBLE_EQ(tval, 4);
  EXPECT_EQ(SUCCESS, out(3, tval).status);
  EXPECT_DOUBLE_EQ(tval, 2);
  EXPECT_EQ(SUCCESS, out(4, tval).status);
  EXPECT_DOUBLE_EQ(tval, 4);
  EXPECT_EQ(SUCCESS, out(5, tval).status);
  EXPECT_DOUBLE_EQ(tval, 6);
  //
}
TEST(MatnTest, test_subtract_scalar_value) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  auto r = m.subtract(2, out);
  EXPECT_EQ(r.status, SUCCESS);
  real tval = static_cast<real>(0);
  EXPECT_EQ(SUCCESS, out(0, tval).status);
  EXPECT_DOUBLE_EQ(tval, -2);
  EXPECT_EQ(SUCCESS, out(1, tval).status);
  EXPECT_DOUBLE_EQ(tval, -1);
  EXPECT_EQ(SUCCESS, out(2, tval).status);
  EXPECT_DOUBLE_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(3, tval).status);
  EXPECT_DOUBLE_EQ(tval, -2);
  EXPECT_EQ(SUCCESS, out(4, tval).status);
  EXPECT_DOUBLE_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(5, tval).status);
  EXPECT_DOUBLE_EQ(tval, 2);
  //
}
TEST(MatnTest, test_subtract_mat_value) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;

  std::vector<real> mv2 = {0, 1, 2, 0, 2, 4};
  MatN<real, 2, 3> tout(mv2);
  // filled matrix
  auto r = m.subtract(tout, out);
  EXPECT_EQ(r.status, SUCCESS);
  real tval = static_cast<real>(0);
  EXPECT_EQ(SUCCESS, out(0, tval).status);
  EXPECT_DOUBLE_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(1, tval).status);
  EXPECT_DOUBLE_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(2, tval).status);
  EXPECT_DOUBLE_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(3, tval).status);
  EXPECT_DOUBLE_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(4, tval).status);
  EXPECT_DOUBLE_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(5, tval).status);
  EXPECT_DOUBLE_EQ(tval, 0);
  //
}
TEST(MatnTest, test_hadamard_product_scalar_value) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  auto r = m.hadamard_product(2, out);
  EXPECT_EQ(r.status, SUCCESS);
  real tval = static_cast<real>(0);
  EXPECT_EQ(SUCCESS, out(0, tval).status);
  EXPECT_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(1, tval).status);
  EXPECT_EQ(tval, 2);
  EXPECT_EQ(SUCCESS, out(2, tval).status);
  EXPECT_EQ(tval, 4);
  EXPECT_EQ(SUCCESS, out(3, tval).status);
  EXPECT_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(4, tval).status);
  EXPECT_EQ(tval, 4);
  EXPECT_EQ(SUCCESS, out(5, tval).status);
  EXPECT_EQ(tval, 8);
  //
}
TEST(MatnTest, test_hadamard_product_mat_value) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;

  std::vector<real> mv2 = {0, 1, 2, 0, 2, 4};
  MatN<real, 2, 3> tout(mv2);
  // filled matrix
  auto r = m.hadamard_product(tout, out);
  EXPECT_EQ(r.status, SUCCESS);
  real tval = static_cast<real>(0);
  EXPECT_EQ(SUCCESS, out(0, tval).status);
  EXPECT_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(1, tval).status);
  EXPECT_EQ(tval, 1);
  EXPECT_EQ(SUCCESS, out(2, tval).status);
  EXPECT_EQ(tval, 4);
  EXPECT_EQ(SUCCESS, out(3, tval).status);
  EXPECT_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(4, tval).status);
  EXPECT_EQ(tval, 4);
  EXPECT_EQ(SUCCESS, out(5, tval).status);
  EXPECT_EQ(tval, 16);
  //
}
TEST(MatnTest, test_divide_scalar_value_false) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  auto r = m.divide(0, out);
  EXPECT_EQ(r.status, ARG_ERROR);
  //
}
TEST(MatnTest, test_divide_scalar_value_false_check_macro) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  auto r = CHECK(m.divide(0, out));
  EXPECT_EQ(r, false);
  //
}
TEST(MatnTest, test_divide_mat_false) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;

  std::vector<real> mv2 = {0, 1, 2, 0, 2, 4};
  MatN<real, 2, 3> tout(mv2);
  // filled matrix
  auto r = m.divide(tout, out);
  EXPECT_EQ(r.status, ARG_ERROR);
  //
}
TEST(MatnTest, test_divide_scalar_value_true) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  auto r = m.divide(2, out);
  EXPECT_EQ(r.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  EXPECT_EQ(SUCCESS, out(0, tval).status);
  EXPECT_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(1, tval).status);
  EXPECT_EQ(tval, static_cast<real>(0.5));
  EXPECT_EQ(SUCCESS, out(2, tval).status);
  EXPECT_EQ(tval, 1);
  EXPECT_EQ(SUCCESS, out(3, tval).status);
  EXPECT_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(4, tval).status);
  EXPECT_EQ(tval, 1);
  EXPECT_EQ(SUCCESS, out(5, tval).status);
  EXPECT_EQ(tval, 2);
}
TEST(MatnTest, test_divide_mat_true) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;

  std::vector<real> mv2 = {1, 1, 2, 1, 2, 4};
  MatN<real, 2, 3> tout(mv2);
  // filled matrix
  auto r = m.divide(tout, out);
  EXPECT_EQ(r.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  EXPECT_EQ(SUCCESS, out(0, tval).status);
  EXPECT_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(1, tval).status);
  EXPECT_EQ(tval, 1);
  EXPECT_EQ(SUCCESS, out(2, tval).status);
  EXPECT_EQ(tval, 1);
  EXPECT_EQ(SUCCESS, out(3, tval).status);
  EXPECT_EQ(tval, 0);
  EXPECT_EQ(SUCCESS, out(4, tval).status);
  EXPECT_EQ(tval, 1);
  EXPECT_EQ(SUCCESS, out(5, tval).status);
  EXPECT_EQ(tval, 1);
}
TEST(MatnTest, test_multiply_mat) {

  // values from S. Lang, Introduction to Linear Algebra,
  // 1986, p. 48

  std::vector<real> Amat_values = {2, 1, 5, 1, 3, 2};
  std::vector<real> Bmat_values = {3, 4, -1, 2, 2, 1};

  MatN<real, 2, 3> A(Amat_values);
  MatN<real, 3, 2> B(Bmat_values);
  MatN<real, 2, 2> out;
  //

  auto start = std::chrono::steady_clock::now();
  auto res = A.multiply(B, out);
  auto stop = std::chrono::steady_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(
          stop - start);

  RecordProperty("duration in microseconds %i",
                 (int)duration.count());

  EXPECT_EQ(res.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  EXPECT_EQ(SUCCESS, out(0, tval).status);
  EXPECT_EQ(tval, 15);
  EXPECT_EQ(SUCCESS, out(1, tval).status);
  EXPECT_EQ(tval, 15);
  EXPECT_EQ(SUCCESS, out(2, tval).status);
  EXPECT_EQ(tval, 4);
  EXPECT_EQ(SUCCESS, out(3, tval).status);
  EXPECT_EQ(tval, 12);
}
TEST(MatnTest, test_gaxpy) {

  // values from Golub, Van Loan 2013, p.5

  real Amat_values[] = {1, 2, 3, 4, 5, 6};
  real x[] = {7, 8};
  real y[] = {0.0, 0.0, 0.0};

  MatN<real, 3, 2> A(Amat_values);
  auto start = std::chrono::high_resolution_clock::now();
  //
  auto res = A.gaxpy(x, y);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(
          stop - start);

  RecordProperty("duration in microseconds %i",
                 (int)duration.count());

  EXPECT_EQ(res.status, SUCCESS);

  //
  EXPECT_EQ(y[0], 23);
  EXPECT_EQ(y[1], 53);
  EXPECT_EQ(y[2], 83);
}
TEST(MatnTest, test_outer_product) {

  // values from Golub, Van Loan 2013, p.7

  real x[] = {1, 2, 3};
  real y[] = {4, 5};
  MatN<real, 3, 2> out(0);
  MatN<real, 3, 2> m(0);
  //
  auto start = std::chrono::steady_clock::now();
  auto res = m.outer_product<3, 2>(x, y, out);
  auto stop = std::chrono::steady_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(
          stop - start);

  RecordProperty("duration in microseconds %i",
                 (int)duration.count());

  EXPECT_EQ(res.status, SUCCESS);

  real d[6];
  out(d);

  //
  EXPECT_EQ(d[0], 4);
  EXPECT_EQ(d[1], 5);
  EXPECT_EQ(d[2], 8);
  EXPECT_EQ(d[3], 10);
  EXPECT_EQ(d[4], 12);
  EXPECT_EQ(d[5], 15);
}
TEST(MatnTest, test_multiply_scalar) {

  // values from S. Lang, Introduction to Linear Algebra,
  // 1986, p. 48

  std::vector<real> Amat_values = {2, 1, 5, 1, 3, 2};
  real B = 2;

  MatN<real, 2, 3> A(Amat_values);
  MatN<real, 2, 2> out;
  //
  auto res = A.multiply(B, out);
  //
  EXPECT_EQ(res.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  EXPECT_EQ(SUCCESS, out(0, tval).status);
  EXPECT_EQ(tval, 16);
  EXPECT_EQ(SUCCESS, out(1, tval).status);
  EXPECT_EQ(tval, 16);
  EXPECT_EQ(SUCCESS, out(2, tval).status);
  EXPECT_EQ(tval, 12);
  EXPECT_EQ(SUCCESS, out(3, tval).status);
  EXPECT_EQ(tval, 12);
}
TEST(MatnTest, test_dot_mat) {

  // values from S. Lang, Introduction to Linear Algebra,
  // 1986, p. 48

  std::vector<real> Amat_values = {2, 1, 5, 1, 3, 2};
  std::vector<real> Bmat_values = {3, 4, -1, 2, 2, 1};

  MatN<real, 2, 3> A(Amat_values);
  MatN<real, 3, 2> B(Bmat_values);
  MatN<real, 2, 2> out;
  //
  auto res = A.dot(B, out);
  EXPECT_EQ(res.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  EXPECT_EQ(SUCCESS, out(0, tval).status);
  EXPECT_EQ(tval, 15);
  EXPECT_EQ(SUCCESS, out(1, tval).status);
  EXPECT_EQ(tval, 15);
  EXPECT_EQ(SUCCESS, out(2, tval).status);
  EXPECT_EQ(tval, 4);
  EXPECT_EQ(SUCCESS, out(3, tval).status);
  EXPECT_EQ(tval, 12);
}
TEST(MatnTest, test_dot_scalar) {

  // values from S. Lang, Introduction to Linear Algebra,
  // 1986, p. 48

  std::vector<real> Amat_values = {2, 1, 5, 1, 3, 2};
  real B = 2;

  MatN<real, 2, 3> A(Amat_values);
  MatN<real, 2, 2> out;
  //
  auto res = A.dot(B, out);
  //
  EXPECT_EQ(res.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  EXPECT_EQ(SUCCESS, out(0, tval).status);
  EXPECT_EQ(tval, 16);
  EXPECT_EQ(SUCCESS, out(1, tval).status);
  EXPECT_EQ(tval, 16);
  EXPECT_EQ(SUCCESS, out(2, tval).status);
  EXPECT_EQ(tval, 12);
  EXPECT_EQ(SUCCESS, out(3, tval).status);
  EXPECT_EQ(tval, 12);
}
TEST(MatnTest, test_dot_vector) {
  // values from
  // https://people.math.umass.edu/~havens/m235Lectures/Lecture05.pdf

  real Amat_values[] = {1,  -4, 7,  -2, 5,
                                   -8, 3,  -6, 9};
  real Bv[] = {2, 1, -1};

  MatN<real, 3, 3> A(Amat_values);
  MatN<real, 3, 1> B(Bv);
  MatN<real, 3, 1> out;
  //
  auto res = A.dot(B, out);
  //
  EXPECT_EQ(res.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  EXPECT_EQ(SUCCESS, out(0, tval).status);
  EXPECT_EQ(tval, -9);
  EXPECT_EQ(SUCCESS, out(1, tval).status);
  EXPECT_EQ(tval, 9);
  EXPECT_EQ(SUCCESS, out(2, tval).status);
  EXPECT_EQ(tval, -9);
}
TEST(MatnTest, test_add_rows) {

  real new_rows[] = {1, -4, 7,
                    -2, 5, -8,
                    3, -6, 9};
  real Avals[] = {2, 1, -1};

  MatN<real, 1, 3> A(Avals);
  MatN<real, 4, 3> out;
  //
  auto res = A.add_rows<9>(new_rows, out);
  //
  EXPECT_EQ(res.status, SUCCESS);
  real check[12];
  out(check);
  //
  EXPECT_EQ(check[0], 2);
  EXPECT_EQ(check[1], 1);
  EXPECT_EQ(check[2], -1);
  EXPECT_EQ(check[3], 1);
  EXPECT_EQ(check[4], -4);
  EXPECT_EQ(check[5], 7);
  EXPECT_EQ(check[6], -2);
  EXPECT_EQ(check[7], 5);
  EXPECT_EQ(check[8], -8);
  EXPECT_EQ(check[9], 3);
  EXPECT_EQ(check[10], -6);
  EXPECT_EQ(check[11], 9);
}
TEST(MatnTest, test_add_columns) {

  real new_columns[] = {1,  -4, 7,  -2, 5,
                                     -8, 3,  -6, 9};
  real Avals[] = {2, 1, -1};

  MatN<real, 3, 1> A(Avals);
  MatN<real, 3, 4> out;
  //
  auto res =
      A.add_columns<9>(new_columns, out);
  //
  EXPECT_EQ(res.status, SUCCESS);
  real check[12];
  out(check);
  //

  EXPECT_EQ(check[0], 2);
  EXPECT_EQ(check[1], 1);
  EXPECT_EQ(check[2], -2);
  EXPECT_EQ(check[3], 3);
  EXPECT_EQ(check[4], 1);
  EXPECT_EQ(check[5], -4);
  EXPECT_EQ(check[6], 5);
  EXPECT_EQ(check[7], -6);
  EXPECT_EQ(check[8], -1);
  EXPECT_EQ(check[9], 7);
  EXPECT_EQ(check[10], -8);
  EXPECT_EQ(check[11], 9);
}
