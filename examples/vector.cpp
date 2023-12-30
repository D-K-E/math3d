#include <math3d/math3d.hpp>

typedef float real;

int main(void) {
  //
  math3d::vecn::VecN<real, 5> v1({-21, 1, 2, 1.53, 3});
  math3d::vecn::VecN<real, 5> v2(3);

  // access to different properties of vector

  // size
  std::size_t size;
  v1(size);
  std::cout << (size == 5) << std::endl;
  // prints out: 1

  // all content
  real a[5];
  v1(a);
  for (std::size_t i = 0; i < size; ++i) {
    std::cout << a[i] << " ";
  }
  std::cout << std::endl;
  // prints out: -21 1 2 1.53 3

  // single member
  real b;
  v1(3, b);
  std::cout << (b == 1.53f) << std::endl;
  // prints out: 1

  // set a new value
  v1(math3d::vecn::make_vecn_cell(0.1f, 1));
  // or
  v1({.content = 0.1f, .index = 1});

  // some basic math
  math3d::vecn::VecN<real, 5> out;
  v1.add(v2, out);
  std::cout << out << std::endl;
  // prints out: { -18, 3.1, 5, 4.53, 6 }

  //
  v1.subtract(v2, out);
  std::cout << out << std::endl;
  // prints out: { -24, -2.9, -1, -1.47, 0 }

  // elementwise multiplication
  v1.multiply(v2, out);
  std::cout << out << std::endl;
  // prints out: { -63, 0.3, 6, 4.59, 9 }

  v1.divide(v2, out);
  std::cout << out << std::endl;
  // prints out: { -7, 0.0333333, 0.666667, 0.51, 1 }

  // what if there is a 0 division?
  math3d::OpResult result = v1.divide(0, out);
  std::cout << (result.status == math3d::ARG_ERROR)
            << std::endl;
  // prints out: 1

  // notice the `out` doesn't get modified
  std::cout << out << std::endl;
  // prints out: { -7, 0.0333333, 0.666667, 0.51, 1 }

  // dot product
  real out2;
  v1.dot(v2, out2);
  std::cout << out2 << std::endl;
  // prints out: -43.11

  // some information about the call
  math3d::OpResult res;
  INFO_VERBOSE_MATH3D(v1.dot(v2, out2), res);
  // prints to std::cerr some useful information in xml
  // format
  /* <OpResult status='SUCCESS' function='dot'
   *   file='/usr/share/math3d/examples/vector.cpp'
   *   callsite='v1.dot(v2, out2)' line='70' duration='0'
   *   durationUnit='microsecond'/>
   */

  // create a base vector
  math3d::vecn::base<real, 1, 5>(out);
  std::cout << out << std::endl;
  // prints out: { 0, 1, 0, 0, 0 }

  math3d::vecn::base<real, 3, 5>(out);
  std::cout << out << std::endl;
  // prints out: { 0, 0, 0, 1, 0 }

  return 0;
}
