#include "core.h"
#include <chrono>
#include <iostream>
#include <string>

std::ostream &operator<<(std::ostream &out,
                         math3d::opstatus_t flag);

template <typename T, size_t R, size_t C = R>
std::ostream &
operator<<(std::ostream &out,
           const math3d::matn::MatN<T, R, C> &m) {
  T arr[R * C];
  m(arr);
  for (std::size_t i = 0; i < (R * C); i++) {
    if (i % C == 0) {
      out << std::endl;
    }
    out << " " << arr[i] << " ";
  }
  out << std::endl;
  return out;
}

template <typename T, size_t R>
std::ostream &
operator<<(std::ostream &out,
           const math3d::vecn::VecN<T, R> &m) {
  T arr[R];
  m(arr);
  out << "{ " << arr[0];
  for (std::size_t i = 1; i < R; i++) {
    out << ", " << arr[i];
  }
  out << " }";
  return out;
}

template <typename T>
std::ostream &
operator<<(std::ostream &out,
           const math3d::quaternion::QuaternionBase &b) {
  switch (b) {
  case math3d::quaternion::kSCALAR: {
    out << "scalar_base" << std::endl;
    break;
  }
  case math3d::quaternion::kI: {
    out << "i_base" << std::endl;
    break;
  }
  case math3d::quaternion::kJ: {
    out << "j_base::" << std::endl;
    break;
  }
  case math3d::quaternion::kK: {
    out << "k_base::" << std::endl;
    break;
  }
  }
  return out;
}
template <typename T>
std::ostream &operator<<(
    std::ostream &out,
    const math3d::quaternion::QuaternionComponent<T> &c) {
  switch (c.base) {
  case math3d::quaternion::kSCALAR: {
    out << c.r << std::endl;
    break;
  }
  case math3d::quaternion::kI: {
    out << c.r << "i";
    break;
  }
  case math3d::quaternion::kJ: {
    out << c.r << "j";
    break;
  }
  case math3d::quaternion::kK: {
    out << c.r << "k";
    break;
  }
  }
  return out;
}

template <typename T>
std::ostream &
operator<<(std::ostream &out,
           const math3d::quaternion::Quaternion<T> &q) {
  math3d::quaternion::QuaternionComponent<T> c;
  q(math3d::quaternion::kSCALAR, c);
  out << c << " + ";
  q(math3d::quaternion::kI, c);
  out << c << " + ";
  q(math3d::quaternion::kJ, c);
  out << c << " + ";
  q(math3d::quaternion::kK, c);
  out << c << std::endl;
  return out;
}

namespace math3d {
bool CHECK(OpResult r);
OpResult INFO(const OpResult &res);
OpResult INFO_VERBOSE(const OpResult &res);
} // namespace math3d

/**Checks if operation was successful*/
#define CHECK_MATH3D(call, res)                            \
  do {                                                     \
    res = call;                                            \
    res.call_name = #call;                                 \
    res.line_info = __LINE__;                              \
    res.file_name = __FILE__;                              \
  } while (0)

/**Prints information about the operation*/
#define INFO_MATH3D(call, res)                             \
  do {                                                     \
    res = call;                                            \
    res.call_name = #call;                                 \
    res.line_info = __LINE__;                              \
    res.file_name = __FILE__;                              \
    if (res.status != SUCCESS) {                           \
      res = math3d::INFO(res);                             \
    }                                                      \
  } while (0)

/**Prints everything about operation result*/
#define INFO_VERBOSE_MATH3D(call, res)                     \
  do {                                                     \
    auto start = std::chrono::steady_clock::now();         \
    res = call;                                            \
    auto stop = std::chrono::steady_clock::now();          \
    auto duration = std::chrono::duration_cast<            \
        std::chrono::microseconds>(stop - start);          \
    auto str_inf = std::to_string(                         \
        static_cast<int>(duration.count()));               \
    res.duration_info = str_inf.c_str();                   \
    res.call_name = #call;                                 \
    res.line_info = __LINE__;                              \
    res.file_name = __FILE__;                              \
    INFO_VERBOSE(res);                                     \
  } while (0)
