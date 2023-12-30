#include <math3d/math3d.hpp>

typedef float real;

int main(void) {
  // empty constructor
  math3d::matn::MatN<real, 3, 2> m1;

  // single array constructor
  math3d::matn::MatN<real, 2, 3> m2({0, 1, 2, 0, 2, 4});
  std::cout << m2 << std::endl;
  /** prints:
0  1  2
0  2  4
   */

  // single value constructor
  math3d::matn::MatN<real, 2, 3> m3(2.0f);
  std::cout << m3 << std::endl;
  /**
2  2  2
2  2  2
   */

  // access to certain properties
  std::size_t snb;
  m2(snb); // get size
  std::cout << snb << std::endl;
  // 6

  std::size_t rnb, cols;
  auto r = m2(rnb, cols); // get row and col number
  std::cout << rnb << "x" << cols << std::endl;
  // 2x3

  // get column
  real c1[2];
  m2.get_column(1, c1);
  for (std::size_t i = 0; i < 2; ++i) {
    std::cout << c1[i] << " ";
  }
  std::cout << std::endl;
  // prints: 1 2

  // get row
  real r1[3];
  m2.get_row(1, r1);
  for (std::size_t i = 0; i < 3; ++i) {
    std::cout << r1[i] << " ";
  }
  std::cout << std::endl;
  // prints: 0 2 4

  // set value
  m2({.content = 1.1, .row = 0, .column = 1});
  std::cout << m2 << std::endl;
  /*
0  1.1  2
0  2  4
   */

  // set column
  c1[0] = 4.1f;
  c1[1] = 0.2f;
  m2.set_column(1, c1);
  std::cout << m2 << std::endl;
  /*
0  4.1  2
0  0.2  4
   */

  // set row
  r1[0] = -1.1;
  r1[1] = 2.01;
  r1[2] = 3.1f;
  m2.set_row(1, r1);
  std::cout << m2 << std::endl;
  /*
0  4.1  2
-1.1  2.01  3.1
   */

  // set 2 rows
  real r2[3 * 2] = {1, 9, 6, 3, 0, 4};
  math3d::matn::MatN<real, 3 + 3, 2> padded;
  m1.add_rows<6>(r2, padded);
  std::cout << padded << std::endl;
  /*
0  0
0  0
0  0
1  9
6  3
0  4
   */

  // transpose
  math3d::matn::MatN<real, 3, 2> m2_T;
  m2.transpose(m2_T);
  std::cout << m2_T << std::endl;
  /**
 0  4.1
-1.1  2.01
 3.1  0
   */

  // some basic arithmetic
  math3d::matn::MatN<real, 2, 3> mout;
  m2.add(m3, mout);
  std::cout << mout << std::endl;
  /*
 2  6.1  4
 0.9  4.01  5.1
   */

  m2.subtract(m3, mout);
  std::cout << mout << std::endl;
  /*
 -2  2.1  0
 -3.1  0.009  1.1
   */

  m2.hadamard_product(m3, mout);
  std::cout << mout << std::endl;
  /*
  0  8.2  4
 -2.2  4.02  6.2
   */

  m2.divide(m3, mout);
  std::cout << mout << std::endl;
  /*
 0  2.05  1
 -0.55  1.005  1.55
   */

  math3d::matn::MatN<real, 2, 2> mout2;
  m2.multiply(m2_T, mout2);
  std::cout << mout2 << std::endl;
  /*
 1.69  8.241
 7.399  -0.4699
   */

  // make the identity matrix
  math3d::matn::identity(mout2);
  std::cout << mout2 << std::endl;
  /*
 1  0
 0  1
   */
  return 0;
}
