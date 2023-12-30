# math3d
yet another 3d math library

This one doesn't throw any exceptions. It doesn't allocate memory on the heap.
It doesn't contain any virtual function, nor a recursive function. It also
doesn't use any function pointers. No concurrency and threading libraries are
used either. All methods are `const` qualified except for setters.

Hence it is highly compatible with OpenCL-C++ and CUDA at the same time.
For CUDA, one would need to add `__host__`, `__device__` attributes in front
of the methods. 
For OpenCL-C++ you probably don't need anything, just copy-paste the stuff
into your kernel and it should be good to go.

We also provide 2 optional headers for facilitating integration of this
library to cuda/opencl projects:

- `include/math3d/core.cuh`: contains all 4 objects (`VecN`, `MatN`, `Quaternion`,
  `LUdecomp`) in a single header. It is destined to be used in CUDA projects.

- `include/math3d/core_cl.h`: contains all 4 objects (`VecN`, `MatN`, `Quaternion`,
  `LUdecomp`) in a single header. It is destined to be used in OpenCL-C++
  projects.

# Usage

Here are some usage examples.

The following snippet outlines several operations that can be done with
vectors.

```c++
// vector example
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
```

The following snippet outlines several operations that can be done with
matrices:
```c++
// matrix example
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
```
