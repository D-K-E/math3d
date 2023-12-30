#ifndef CORE_CUH
#define CORE_CUH

#include <cstdint>

namespace math3d {

enum opstatus_t : uint8_t {
  SUCCESS = 1,
  INDEX_ERROR = 2,
  SIZE_ERROR = 3,
  ARG_ERROR = 4,
  LU_ERROR = 5,
  NOT_IMPLEMENTED = 6
};

struct OpResult {
  //
  size_t line_info = 0;
  const char *file_name{nullptr};
  const char *fn_name{nullptr};
  const char *call_name{nullptr};
  const char *duration_info{nullptr};

  opstatus_t status = NOT_IMPLEMENTED;
  bool success = false;

  __host__ __device__ OpResult() = default;
  __host__ __device__ ~OpResult() = default;

  __host__ __device__ OpResult(size_t line, const char *fname,
                               const char *funcname, const char *cname,
                               opstatus_t op)
      : line_info(line), file_name(fname), fn_name(funcname), call_name(cname),
        status(op), success(op == SUCCESS) {}
};

namespace vecn {

template <typename T> struct VecNCell {
  T content;
  std::size_t index = 0;
};

template <typename T>
__host__ __device__ VecNCell<T> make_vecn_cell(const T &c, std::size_t i) {
  VecNCell<T> cell{.content = c, .index = i};
  return cell;
}

template <class T, std::size_t N> class VecN {

public:
  /*! Tested */
  __host__ __device__ VecN() {
    for (std::size_t i = 0; i < N; ++i) {
      data[i] = 0;
    }
  }
  /*! Tested */
  __host__ __device__ VecN(const T (&arr)[N]) {
    memcpy(data, arr, N * sizeof(T));
  }
  __host__ __device__ VecN(T s) {
    for (std::size_t i = 0; i < N; ++i) {
      data[i] = s;
    }
  }
  /*! Tested */
  __host__ __device__ OpResult operator()(std::size_t index, T &out) const {
    if (index >= N) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "(size_t, T&)",
                     INDEX_ERROR);
      return vflag;
    }
    out = data[index];
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "(size_t, T&)", SUCCESS);
    return vflag;
  }
  __host__ __device__ OpResult operator()(T (&out)[N]) const {
    memcpy(out, data, N * sizeof(T));
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "(T (&)[N])", SUCCESS);
    return vflag;
  }
  __host__ __device__ constexpr OpResult operator()(std::size_t &out) const {
    out = N;
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "(size_t)", SUCCESS);
    return vflag;
  }

  /*! Tested */
  __host__ __device__ OpResult operator()(const VecNCell<T> &cell) {
    if (cell.index >= N) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "(const VecNCell<T>&)",
                      INDEX_ERROR);
    }
    data[cell.index] = cell.content;

    return OpResult(__LINE__, __FILE__, __FUNCTION__, "(const VecNCell<T>&)",
                    SUCCESS);
  }

  /*! Tested */
  __host__ __device__ OpResult add(T v, T (&out)[N]) const {
    VecN<T, N> v_1(v);
    T v_in[N];
    v_1(v_in);
    //
    auto res = add(v_in, out);
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "add(T, T (&)[N])",
                   res.status);
    return vflag;
  }

  /*! Tested */
  __host__ __device__ OpResult add(T v, VecN<T, N> &vout) const {
    T out[N];
    auto res = add(v, out);
    vout = VecN<T, N>(out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "add(T, VecN<T, N>&",
                   res.status);
    return vflag;
  }
  /*! Tested */
  __host__ __device__ OpResult add(const T (&v)[N], T (&out)[N]) const {
    auto fn = [](T thisel, T argel) { return thisel + argel; };
    auto res = apply_el<decltype(fn)>(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "add(const T(&)[N], T (&)[N])", res.status);
    return vflag;
  }
  /*! Tested */
  __host__ __device__ OpResult add(const VecN<T, N> &v, VecN<T, N> &out) const {
    T v_in[N];
    v(v_in);
    T vout[N];
    auto res = add(v_in, vout);
    out = VecN<T, N>(vout);
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "add(const VecN<T, N>&, VecN<T, N>&)", res.status);
    return vflag;
  }
  //
  __host__ __device__ OpResult subtract(T v, T (&out)[N]) const {
    VecN<T, N> v_1(v);
    T v_in[N];
    v_1(v_in);
    auto res = subtract(v_in, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "subtract(T, T (&)[N])",
                   res.status);
    return vflag;
  }
  /*! Tested */
  __host__ __device__ OpResult subtract(T v, VecN<T, N> &vout) const {
    VecN<T, N> v_1(v);
    T v_in[N];
    v_1(v_in);
    T out[N];
    auto res = subtract(v_in, out);
    vout = VecN<T, N>(out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "subtract(T, VecN<T, N>&",
                   res.status);
    return vflag;
  }
  /*! Tested */
  __host__ __device__ OpResult subtract(const T (&v)[N], T (&out)[N]) const {
    auto fn = [](T thisel, T argel) { return thisel - argel; };
    auto res = apply_el<decltype(fn)>(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "subtract(const T(&)[N], T (&)[N])", res.status);
    return vflag;
  }
  /*! Tested */
  __host__ __device__ OpResult subtract(const VecN<T, N> &v,
                                        VecN<T, N> &out) const {
    T v_in[N];
    v(v_in);
    T vout[N];
    auto res = subtract(v_in, vout);
    out = VecN<T, N>(vout);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "subtract(const VecN<T, N>&, VecN<T, N>&)", res.status);
    return vflag;
  }
  //
  __host__ __device__ OpResult multiply(T v, T (&out)[N]) const {
    T v_in[N];
    VecN<T, N> vv(v);
    vv(v_in);

    auto res = multiply(v_in, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "multiply(T, T (&)[N])",
                   res.status);
    return vflag;
  }
  /*! Tested */
  __host__ __device__ OpResult multiply(T v, VecN<T, N> &vout) const {
    T v_in[N];
    VecN<T, N> vv(v);
    vv(v_in);

    T out[N];
    auto res = multiply(v_in, out);
    vout = VecN<T, N>(out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "multiply(T, VecN<T, N>&",
                   res.status);
    return vflag;
  }
  /*! Tested */
  __host__ __device__ OpResult multiply(const T (&v)[N], T (&out)[N]) const {
    auto fn = [](T thisel, T argel) { return thisel * argel; };
    auto res = apply_el<decltype(fn)>(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "multiply(const T(&)[N], T (&)[N])", res.status);
    return vflag;
  }
  /*! Tested */
  __host__ __device__ OpResult multiply(const VecN<T, N> &v,
                                        VecN<T, N> &out) const {
    T v_in[N];
    v(v_in);
    T vout[N];
    auto res = multiply(v_in, vout);
    out = VecN<T, N>(vout);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "multiply(const VecN<T, N>&, VecN<T, N>&)", res.status);
    return vflag;
  }

  //
  __host__ __device__ OpResult divide(T v, T (&out)[N]) const {
    if (v == static_cast<T>(0)) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "divide(T, T (&)[N])",
                     ARG_ERROR);
      return vflag;
    }

    T v_in[N];
    VecN<T, N> v_(v);
    v_(v_in);

    auto res = divide(v_in, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "divide(T, T (&)[N])",
                   res.status);
    return vflag;
  }
  /*! Tested */
  __host__ __device__ OpResult divide(T v, VecN<T, N> &vout) const {
    // check for zero division
    if (v == static_cast<T>(0)) {
      OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "divide(T, VecN<T, N>&",
                     ARG_ERROR);
      return vflag;
    }
    T v_in[N];
    VecN<T, N> v_(v);
    v_(v_in);
    T out[N];

    auto res = divide(v_in, out);
    vout = VecN<T, N>(out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "divide(T, VecN<T, N>&",
                   res.status);
    return vflag;
  }
  /*! Tested */
  __host__ __device__ OpResult divide(const T (&v)[N], T (&out)[N]) const {
    for (std::size_t j = 0; j < N; j++) {
      if (v[j] == static_cast<T>(0)) {
        OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                       "divide(const T(&)[N], T (&)[N])", ARG_ERROR);
        return vflag;
      }
    }
    auto fn = [](T thisel, T argel) { return thisel / argel; };
    auto res = apply_el<decltype(fn)>(v, fn, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "divide(const T(&)[N], T (&)[N])", res.status);
    return vflag;
  }
  /*! Tested */
  __host__ __device__ OpResult divide(const VecN<T, N> &v,
                                      VecN<T, N> &out) const {
    // check zero division
    T v_in[N];
    v(v_in);
    T vout[N];
    auto res = divide(v_in, vout);
    out = VecN<T, N>(vout);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "divide(const VecN<T, N>&, VecN<T, N>&)", res.status);
    return vflag;
  }
  __host__ __device__ OpResult dot(const T &v, T &out) const {
    T v_in[N];
    VecN<T, N> v_(v);
    v_(v_in);
    dot(v_in, out);
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "dot(const T&, T &)",
                   SUCCESS);
    return vflag;
  }
  __host__ __device__ OpResult dot(const T (&v)[N], T &out) const {
    out = static_cast<T>(0);
    for (std::size_t i = 0; i < N; i++) {
      out += data[i] * v[i];
    }

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "dot(const T(&)[N], T &)",
                   SUCCESS);
    return vflag;
  }
  __host__ __device__ OpResult dot(const VecN<T, N> &v, T &out) const {
    T v_in[N];
    v(v_in);
    dot(v_in, out);

    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "dot(const VecN<T, N>&, T &)", SUCCESS);
    return vflag;
  }

  __host__ __device__ OpResult cross(const VecN<T, N> &v,
                                     VecN<T, N> &out) const {
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "cross(const VecN<T, N>&, VecN<T, N> &)", NOT_IMPLEMENTED);
    return vflag;
  }

private:
  /** holds the vector data*/
  T data[N];

  template <typename Func>
  __host__ __device__ OpResult apply_el(const T (&v)[N], const Func &fn,
                                        T (&out)[N]) const {

    for (std::size_t i = 0; i < N; i++) {
      out[i] = fn(data[i], v[i]);
    }
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                   "apply_el(const T(&)[N], const Func&, T(&)[N])", SUCCESS);
    return vflag;
  }
};

template <typename T, std::size_t BaseOrder, std::size_t N>
__host__ __device__ OpResult base(VecN<T, N> &vout) {
  if (BaseOrder >= N) {
    OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "base(VecN<T, N>&)",
                   ARG_ERROR);
    return vflag;
  }
  VecN<T, N> out;
  out(make_vecn_cell<T>(1, BaseOrder));
  vout = out;
  OpResult vflag(__LINE__, __FILE__, __FUNCTION__, "base(VecN<T, N>&)",
                 SUCCESS);
  return vflag;
}

template <typename T>
__host__ __device__ OpResult cross(const VecN<T, 2> &v1, const VecN<T, 2> &v2,
                                   T &out) {
  T t1;
  v1(0, t1);
  T t2;
  v1(1, t2);
  T v_1;
  v2(0, v_1);
  T v_2;
  v2(1, v_2);
  T t1v2 = t1 * v_2;
  T t2v1 = t2 * v_1;
  out = t1v2 - t2v1;
  OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                 "cross(const VecN<T, 2>&, T&)", SUCCESS);
  return vflag;
}

template <typename T>
__host__ __device__ OpResult cross(const VecN<T, 3> &a, const VecN<T, 3> &b,
                                   VecN<T, 3> &out) {
  T tx;
  a(0, tx);
  T ty;
  a(1, ty);
  T tz;
  a(2, tz);
  //
  T vx;
  b(0, vx);
  T vy;
  b(1, vy);
  T vz;
  b(2, vz);
  //
  T x = ty * vz - tz * vy;
  T y = tz * vx - tx * vz;
  T z = tx * vy - ty * vx;
  out(make_vecn_cell(x, 0));
  out(make_vecn_cell(y, 1));
  out(make_vecn_cell(z, 2));
  OpResult vflag(__LINE__, __FILE__, __FUNCTION__,
                 "cross(const VecN<T, 3>&, VecN<T, 3>&)", SUCCESS);
  return vflag;
}

} // namespace vecn

namespace quaternion {
// holds quaternion related operations

/**Quaternion bases*/
enum QuaternionBase : std::uint8_t {
  kSCALAR = 1, // null value for quaternion base it has
               // no effect on computation
  kI = 2,      // i base for the second quaternion component
  kJ = 3,      // j base for the third quaternion component
  kK = 4       // k base for the fourth quaternion component
};

/**
  \brief Quaternion component
 */
template <class T> struct QuaternionComponent {
  QuaternionBase base;
  T r;

  __host__ __device__ QuaternionComponent() : r(0), base(kSCALAR) {}
  __host__ __device__ QuaternionComponent(QuaternionBase b, T a)
      : base(b), r(a) {}
};

template <class T> class Quaternion {
public:
  __host__ __device__ Quaternion()
      : coeffs{static_cast<T>(0), static_cast<T>(1), static_cast<T>(1),
               static_cast<T>(1)} {}
  //
  __host__ __device__ Quaternion(T c1, const QuaternionComponent<T> qs[3]) {
    coeffs[0] = c1;
    coeffs[static_cast<std::uint8_t>(qs[0].base) - 1] = qs[0].r;
    coeffs[static_cast<std::uint8_t>(qs[1].base) - 1] = qs[1].r;
    coeffs[static_cast<std::uint8_t>(qs[2].base) - 1] = qs[2].r;
  }
  __host__ __device__ Quaternion(T c1, const vecn::VecN<T, 3> &vs) {
    coeffs[0] = c1;
    T v_0;
    vs(0, v_0);
    T v_1;
    vs(1, v_1);
    T v_2;
    vs(2, v_2);
    coeffs[1] = v_0;
    coeffs[2] = v_1;
    coeffs[3] = v_2;
  }

  __host__ __device__ Quaternion(const QuaternionComponent<T> qs[4]) {
    coeffs[static_cast<std::uint8_t>(qs[0].base) - 1] = qs[0].r;
    coeffs[static_cast<std::uint8_t>(qs[1].base) - 1] = qs[1].r;
    coeffs[static_cast<std::uint8_t>(qs[2].base) - 1] = qs[2].r;
    coeffs[static_cast<std::uint8_t>(qs[3].base) - 1] = qs[3].r;
  }

  __host__ __device__ OpResult operator()(T &out) const { return scalar(out); }
  __host__ __device__ OpResult scalar(T &out) const {

    out = coeffs[0];
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "scalar", SUCCESS);
  }
  __host__ __device__ OpResult vector(vecn::VecN<T, 3> &out) const {
    T v[3];
    v[0] = coeffs[1];
    v[1] = coeffs[2];
    v[2] = coeffs[3];
    out = vecn::VecN<T, 3>(v);

    return OpResult(__LINE__, __FILE__, __FUNCTION__, "vector", SUCCESS);
  }
  __host__ __device__ OpResult operator()(vecn::VecN<T, 3> &out) const {
    return vector(out);
  }
  __host__ __device__ OpResult operator()(const QuaternionBase &b,
                                          QuaternionComponent<T> &c) const {
    switch (b) {
    case kSCALAR: {
      c = QuaternionComponent(kSCALAR, coeffs[0]);
      break;
    }
    case kI: {
      c = QuaternionComponent(kI, coeffs[1]);
      break;
    }
    case kJ: {
      c = QuaternionComponent(kJ, coeffs[2]);
      break;
    }
    case kK: {
      c = QuaternionComponent(kK, coeffs[3]);
      break;
    }
    }

    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "(const QuaternionBase&, QuaternionComponent<T>&)",
                    SUCCESS);
  }

  __host__ __device__ OpResult
  operator()(const QuaternionComponent<T> &c) const {
    QuaternionBase b = c.base;
    coeffs[static_cast<std::uint8_t>(c.base) - 1] = c.r;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "(const QuaternionComponent<T>&)", SUCCESS);
  }
  __host__ __device__ OpResult multiply(T t, vecn::VecN<T, 3> &out) const {
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.multiply(t, out);
  }
  __host__ __device__ OpResult add(T t, vecn::VecN<T, 3> &out) const {
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.add(t, out);
  }
  __host__ __device__ OpResult subtract(T t, vecn::VecN<T, 3> &out) const {
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.subtract(t, out);
  }
  __host__ __device__ OpResult divide(T t, vecn::VecN<T, 3> &out) const {
    if (t == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "divide", ARG_ERROR);
    }
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.divide(t, out);
  }
  /** arithmetic operations with a vector on vector part*/
  __host__ __device__ OpResult multiply(vecn::VecN<T, 3> t,
                                        vecn::VecN<T, 3> &out) const {
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.multiply(t, out);
  }
  __host__ __device__ OpResult add(vecn::VecN<T, 3> t,
                                   vecn::VecN<T, 3> &out) const {

    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.add(t, out);
  }
  __host__ __device__ OpResult subtract(vecn::VecN<T, 3> t,
                                        vecn::VecN<T, 3> &out) const {

    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.subtract(t, out);
  }
  __host__ __device__ OpResult divide(vecn::VecN<T, 3> t,
                                      vecn::VecN<T, 3> &out) const {
    for (std::size_t i = 0; i < 3; i++) {
      T v;
      t(i, v);
      if (v == 0) {
        return OpResult(__LINE__, __FILE__, __FUNCTION__, "divide", ARG_ERROR);
      }
    }
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.divide(t, out);
  }
  /** dot product and cross product for two vec3*/
  __host__ __device__ OpResult dot(const vecn::VecN<T, 3> &t, T &out) const {

    vecn::VecN<T, 3> vec;
    vector(vec);
    return vec.dot(t, out);
  }
  __host__ __device__ OpResult cross(const vecn::VecN<T, 3> &t,
                                     vecn::VecN<T, 3> &out) const {
    vecn::VecN<T, 3> vec;
    vector(vec);
    return vecn::cross(vec, t, out);
  }

  /** Quaternion product as it is shown by Vince 2011 p.
   * 63
   Given a quaternion \f[q_a = [s_a, a]\f] and
   another quaternion \f[q_b = [s_b, b]\f]
   Their product is equal to:
   \f[s_a s_b - a \cdot b, s_a b + s_b a + a \times b \f]
   */
  __host__ __device__ OpResult hamilton_product(const Quaternion &q_b,
                                                Quaternion<T> &out) const {
    // s_a, s_b, a, b
    T s_a = static_cast<T>(0);
    auto res = scalar(s_a);
    if (res.status != SUCCESS) {
      return res;
    }

    T s_b = static_cast<T>(0);
    res = q_b.scalar(s_b);
    if (res.status != SUCCESS) {
      return res;
    }

    vecn::VecN<T, 3> a;
    res = vector(a);
    if (res.status != SUCCESS) {
      return res;
    }

    vecn::VecN<T, 3> b;
    res = q_b.vector(b);
    if (res.status != SUCCESS)
      return res;

    // s_a * s_b
    T s_ab = s_a * s_b;

    // a \cdot b
    T a_dot_b = static_cast<T>(0);
    res = a.dot(b, a_dot_b);
    if (res.status != SUCCESS)
      return res;

    // a \times b
    vecn::VecN<T, 3> cross_ab;
    res = vecn::cross(a, b, cross_ab);
    if (res.status != SUCCESS) {
      return res;
    }

    // s_a * b + s_b * a + a \times b
    T out_v[3];
    T b_s[3];
    T a_s[3];
    T ab_s[3];
    b(b_s);
    a(a_s);
    cross_ab(ab_s);
    for (std::size_t i = 0; i < 3; i++) {
      out_v[i] = s_a * b_s[i] + s_b * a_s[i] + ab_s[i];
    }
    vecn::VecN<T, 3> tout(out_v);
    out = Quaternion(s_ab - a_dot_b, tout);
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "hamilton_product",
                    SUCCESS);
  }
  __host__ __device__ OpResult conjugate(Quaternion<T> &out) const {
    T s = static_cast<T>(0);
    auto res = scalar(s);
    if (res.status != SUCCESS)
      return res;

    vecn::VecN<T, 3> vec;
    res = vector(vec);

    if (res.status != SUCCESS)
      return res;

    res = multiply(static_cast<T>(-1), vec);
    if (res.status != SUCCESS)
      return res;

    out = Quaternion(s, vec);
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "conjugate", SUCCESS);
  }
  /**
    \brief from Vince 2011 - Quaternions for Computer
    Graphics p. 69
   */
  __host__ __device__ OpResult normalized(Quaternion<T> &out) const {
    T nval = static_cast<T>(0);
    auto res = norm(nval);
    if (res.status != SUCCESS)
      return res;
    T inv_mag = static_cast<T>(1.0) / nval;

    res = scalar(nval);

    if (res.status != SUCCESS)
      return res;

    T scalar_part = nval * inv_mag;
    vecn::VecN<T, 3> vs;
    res = vector(vs);
    if (res.status != SUCCESS)
      return res;
    res = multiply(inv_mag, vs);

    if (res.status != SUCCESS)
      return res;
    out = Quaternion(scalar_part, vs);

    return OpResult(__LINE__, __FILE__, __FUNCTION__, "normalized", SUCCESS);
  }
  /**
    \brief from Vince 2011 - Quaternions for Computer
    Graphics p. 69
   */
  __host__ __device__ OpResult inversed(Quaternion<T> &out) const {
    //
    T q_norm = static_cast<T>(0);
    auto res = norm_squared(q_norm);
    if (res.status != SUCCESS) {
      return res;
    }

    Quaternion conj;
    res = conjugate(conj);
    if (res.status != SUCCESS) {
      return res;
    }
    T c_s;
    vecn::VecN<T, 3> c_v;
    conj(c_s);
    conj(c_v);
    vecn::VecN<T, 3> v;
    c_v.divide(q_norm, v);
    out = Quaternion(c_s / q_norm, v);
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "inversed", SUCCESS);
  }
  /**
   \brief from Vince 2011 - Quaternions for Computer
   Graphics p. 69
   */
  __host__ __device__ OpResult add(const Quaternion &q,
                                   Quaternion<T> &out) const {
    auto fn = [](T thisval, T tval) { return thisval + tval; };
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "add",
                    apply_el(q, fn, out));
  }
  /**
   \brief from Vince 2011 - Quaternions for Computer
   Graphics p. 69
  */
  __host__ __device__ OpResult subtract(const Quaternion &q,
                                        Quaternion<T> &out) const {
    auto fn = [](T thisval, T tval) { return thisval - tval; };
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "subtract",
                    apply_el(q, fn, out));
  }
  /**
    \brief from Vince 2011 - Quaternions for Computer
    Graphics p. 69
   */
  __host__ __device__ OpResult product(const Quaternion &q,
                                       Quaternion<T> &out) const {
    return hamilton_product(q, out);
  }
  /**
    \brief from Vince 2011 - Quaternions for Computer
    Graphics p. 69
   */
  __host__ __device__ OpResult product(T r, Quaternion<T> &out) const {
    auto fn = [](T thisval, T tval) { return thisval * tval; };
    Quaternion q(r, vecn::VecN<T, 3>(r));
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "product",
                    apply_el(q, fn, out));
  }
  __host__ __device__ OpResult power(std::size_t i, Quaternion<T> &out) const {
    Quaternion accumulant = *this;
    Quaternion result2 = *this;
    for (std::size_t j = 1; j < i; j++) {
      accumulant.hamilton_product(result2, accumulant);
    }
    out = accumulant;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "power", SUCCESS);
  }
  __host__ __device__ OpResult squared(Quaternion<T> &out) const {
    Quaternion r1 = *this;
    Quaternion r2 = *this;
    auto res = product(r1, out);
    if (res.status != SUCCESS)
      return res;

    res = product(r2, out);
    if (res.status != SUCCESS)
      return res;

    return OpResult(__LINE__, __FILE__, __FUNCTION__, "squared", SUCCESS);
  }
  /**
    \brief from Vince 2011 - Quaternions for Computer
    Graphics p. 69
   */
  __host__ __device__ OpResult norm(T &out) const {
    auto res = norm_squared(out);

    if (res.status != SUCCESS)
      return res;
    //
    out = sqrt(out);

    return OpResult(__LINE__, __FILE__, __FUNCTION__, "norm", SUCCESS);
  }

  /**
    \brief from Vince 2011 - Quaternions for Computer
    Graphics p. 25
   */
  __host__ __device__ OpResult norm_squared(T &out) const {
    T s = static_cast<T>(0);
    auto res = scalar(s);
    if (res.status != SUCCESS)
      return res;
    vecn::VecN<T, 3> vec;
    res = vector(vec);
    if (res.status != SUCCESS)
      return res;

    T a = s * s;
    T b;
    T v[3];
    vec(v);
    res = vec.dot(v, b);
    if (res.status != SUCCESS)
      return res;
    out = a + b;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "norm_squared", SUCCESS);
  }
  __host__ __device__ OpResult magnitude(T &out) const { return norm(out); }

private:
  template <typename Func>
  __host__ __device__ opstatus_t apply_el(const Quaternion &q, const Func &fn,
                                          Quaternion<T> &out) const {
    T s;
    scalar(s);
    vecn::VecN<T, 3> v;
    vector(v);
    T vs[3];
    v(vs);

    //
    T q_s;
    q.scalar(q_s);
    vecn::VecN<T, 3> q_v;
    q.vector(q_v);
    T qs[3];
    q_v(qs);
    //
    T rs[3];
    for (std::size_t i = 0; i < 3; ++i) {
      rs[i] = fn(vs[i], qs[i]);
    }
    T r = fn(s, q_s);
    out = Quaternion(r, vecn::VecN<T, 3>(rs));
    return SUCCESS;
  }

  //
  T coeffs[4];
};

} // namespace quaternion

namespace matn {

template <typename T> struct MatNCell {
  T content;
  std::size_t row = 0;
  std::size_t column = 0;
  int index = -1;
};

template <typename T>
__host__ __device__ MatNCell<T>
make_matn_cell(const T &c, std::size_t row_param, std::size_t col_param) {
  MatNCell<T> cell{.content = c, .row = row_param, .column = col_param};
  return cell;
}

template <typename T>
__host__ __device__ MatNCell<T> make_matn_cell(const T &c, std::size_t i) {
  MatNCell<T> cell{
      .content = c, .row = 0, .column = 0, .index = static_cast<int>(i)};
  return cell;
}

template <class T = float, std::size_t RowNb = 1, std::size_t ColNb = RowNb>
class MatN {
  /** holds the vector data*/
  T data[ColNb * RowNb];
  static const std::size_t Size = ColNb * RowNb;

public:
  __host__ __device__ MatN() {
    for (std::size_t i = 0; i < Size; ++i) {
      data[i] = 0;
    }
  }
  template <size_t N> __host__ __device__ MatN(const T (&vd)[N]) {
    if (N == Size) {
      memcpy(data, vd, Size * sizeof(T));
    } else if (N < Size) {
      for (size_t i = 0; i < N; ++i) {
        data[i] = vd[i];
      }
      for (size_t i = N; i < Size; ++i) {
        data[i] = 0;
      }
    } else {
      for (size_t i = 0; i < Size; ++i) {
        data[i] = vd[i];
      }
    }
  }

  __host__ __device__ MatN(T fill_value) {
    for (std::size_t i = 0; i < Size; ++i) {
      data[i] = fill_value;
    }
  }

  static __host__ __device__ OpResult
  from_row_cols(T v, MatN<T, RowNb, ColNb> &out) {
    MatN<T, RowNb, ColNb> mat(v);
    out = mat;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "from_row_cols", SUCCESS);
  }

  /**\brief Create matrix based on argument matrix*/
  static __host__ __device__ OpResult zeros(MatN<T, RowNb, ColNb> &out) {
    MatN<T, RowNb, ColNb>::from_row_cols(0, out);
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "zeros", SUCCESS);
  }
  // tested
  template <std::size_t OutRowNb = RowNb, std::size_t OutColNb = ColNb>
  __host__ __device__ OpResult fill(T v,
                                    MatN<T, OutRowNb, OutColNb> &out) const {
    std::size_t s = 0;
    out(s);
    for (std::size_t i = 0; i < s; i++) {
      out(i, v);
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "fill", SUCCESS);
  }
  // tested
  __host__ __device__ OpResult transpose(MatN<T, ColNb, RowNb> &out) const {

    for (std::size_t i = 0; i < RowNb; i++) {
      for (std::size_t j = 0; j < ColNb; j++) {
        T tout = static_cast<T>(0);
        (*this)(i, j, tout);                // getter
        out(make_matn_cell<T>(tout, i, j)); // setter
      }
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "transpose", SUCCESS);
  }

  __host__ __device__ OpResult operator()(std::size_t &rows,
                                          std::size_t &cols) const {
    rows = RowNb;
    cols = ColNb;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "get_rows_cols", SUCCESS);
  }

  // tested
  __host__ __device__ OpResult operator()(std::size_t &out) const {
    out = Size;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "get_size", SUCCESS);
  }
  __host__ __device__ OpResult operator()(std::size_t row, std::size_t col,
                                          T &out) const {
    std::size_t index = row * ColNb + col;
    OpResult r = (*this)(index, out);
    if (r.status != SUCCESS)
      return r;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "get", SUCCESS);
  }
  __host__ __device__ OpResult operator()(std::size_t index, T &out) const {
    if (index >= Size)
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "get", INDEX_ERROR);
    out = data[index];
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "get", SUCCESS);
  }
  __host__ __device__ OpResult operator()(T (&out)[RowNb * ColNb]) const {
    memcpy(out, data, Size * sizeof(T));
    // for (std::size_t i = 0; i < size; ++i){
    //     out[i] = data[i];
    // }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "get", SUCCESS);
  }
  __host__ __device__ OpResult operator()(const MatNCell<T> &cell) {
    std::size_t row = cell.row;
    std::size_t col = cell.column;
    T el = cell.content;
    int cell_index = cell.index;
    std::size_t index;
    if (cell_index < 0) {
      index = row * ColNb + col;
    } else {
      index = static_cast<std::size_t>(cell_index);
    }
    //
    if (index >= Size) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "set", INDEX_ERROR);
    }
    data[index] = el;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "set", SUCCESS);
  }
  __host__ __device__ OpResult get_column(std::size_t index,
                                          T (&out)[RowNb]) const {
    if (index >= ColNb) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "get_column",
                      INDEX_ERROR);
    }
    for (std::size_t i = 0; i < RowNb; i++) {
      out[i] = data[i * ColNb + index];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "get_column", SUCCESS);
  }
  __host__ __device__ OpResult set_column(std::size_t index,
                                          const T (&idata)[RowNb]) {
    if (index >= ColNb) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "set_column",
                      INDEX_ERROR);
    }
    for (std::size_t i = 0; i < RowNb; i++) {
      data[i * ColNb + index] = idata[i];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "set_column", SUCCESS);
  }
  __host__ __device__ OpResult get_row(std::size_t index,
                                       T (&out)[ColNb]) const {
    if (index >= RowNb) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "get_row", INDEX_ERROR);
    }
    for (std::size_t i = 0; i < ColNb; i++) {
      out[i] = data[index * ColNb + i];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "get_row", SUCCESS);
  }
  __host__ __device__ OpResult set_row(std::size_t index,
                                       const T (&idata)[ColNb]) {
    if (index >= RowNb) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "set_row", INDEX_ERROR);
    }
    for (std::size_t i = 0; i < ColNb; i++) {
      data[index * ColNb + i] = idata[i];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "set_row", SUCCESS);
  }

  /**Obtain submatrix TODO*/
  __host__ __device__ OpResult submat(std::size_t row_start,
                                      std::size_t col_start,
                                      MatN<T, RowNb, ColNb> &out) const {
    std::size_t row_size = RowNb - row_start;
    std::size_t col_size = ColNb - col_start;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "submat",
                    NOT_IMPLEMENTED);
  }
  __host__ __device__ OpResult add(const MatN<T, RowNb, ColNb> &v,
                                   MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv + val; };
    return apply(v, fn, out);
  }
  __host__ __device__ OpResult add(T v, MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv + val; };
    return apply(v, fn, out);
  }
  __host__ __device__ OpResult subtract(const MatN<T, RowNb, ColNb> &v,
                                        MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv - val; };
    return apply(v, fn, out);
  }
  __host__ __device__ OpResult subtract(T v, MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv - val; };
    return apply(v, fn, out);
  }
  __host__ __device__ OpResult hadamard_product(
      const MatN<T, RowNb, ColNb> &v, MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv * val; };
    return apply(v, fn, out);
  }
  __host__ __device__ OpResult
  hadamard_product(T v, MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv * val; };
    return apply(v, fn, out);
  }
  __host__ __device__ OpResult divide(const MatN<T, RowNb, ColNb> &v,
                                      MatN<T, RowNb, ColNb> &out) const {
    std::size_t osize = 0;
    v(osize);
    for (std::size_t i = 0; i < osize; i++) {
      T tout = static_cast<T>(0);
      v(i, tout); // getter
      if (tout == static_cast<T>(0)) {
        // zero division risk
        return OpResult(__LINE__, __FILE__, __FUNCTION__, "divide", ARG_ERROR);
      }
    }
    auto fn = [](T matv, T val) { return matv / val; };
    return apply(v, fn, out);
  }
  __host__ __device__ OpResult divide(T v, MatN<T, RowNb, ColNb> &out) const {
    if (v == static_cast<T>(0)) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "divide", ARG_ERROR);
    }
    auto fn = [](T matv, T val) { return matv / val; };
    return apply(v, fn, out);
  }
  /**Declares inner vector product*/
  template <std::size_t N = RowNb>
  __host__ __device__ OpResult vdot(const T (&x)[N], const T (&y)[N], T &out) {
    if (N == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "vdot", SIZE_ERROR);
    }

    out = static_cast<T>(0);
    for (std::size_t i = 0; i < N; i++) {
      out += x[i] * y[i];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "vdot", SUCCESS);
  }

  /**Declares inner vector product with scalars*/
  template <std::size_t N = RowNb>
  __host__ __device__ OpResult vdot_s(const T (&x)[N], const T &a,
                                      T (&out)[N]) {

    if (N == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "vdot_s", SIZE_ERROR);
    }
    for (std::size_t i = 0; i < N; i++) {
      out[i] = x[i] * a;
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "vdot_s", SUCCESS);
  }
  /**Implements saxpy algorithm from Golub, Van Loan 2013,
   * p. 4 alg.1.1.2*/
  template <std::size_t N = RowNb>
  __host__ __device__ OpResult saxpy(const T &a, const T (&x)[N],
                                     T (&y)[N]) const {
    if (N == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "saxpy", SIZE_ERROR);
    }
    for (std::size_t i = 0; i < N; i++) {
      y[i] += x[i] * a; //
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "saxpy", SUCCESS);
  }
  /**
    Implements gaxpy algorithm from Golub, Van Loan 2013, p.
    4 alg.1.1.3

    as specified in p. 6-7
   */
  __host__ __device__ OpResult gaxpy(const T (&x)[ColNb], T (&y)[RowNb]) const {
    for (std::size_t j = 0; j < ColNb; j++) {
      T c_j[RowNb];
      get_column(j, c_j);
      saxpy<RowNb>(x[j], c_j, y);
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "gaxpy", SUCCESS);
  }
  /**
     Implements outer product update from Golub, Van Loan
     2013, p. 7 as a series of saxpy operations
    */
  template <std::size_t Rn, std::size_t Cn>
  __host__ __device__ OpResult outer_product(const T (&x)[Rn], const T (&y)[Cn],
                                             MatN<T, Rn, Cn> &out) const {
    for (std::size_t i = 0; i < Rn; i++) {
      T A_i[Cn];
      out.get_row(i, A_i);
      saxpy<Cn>(x[i], y, A_i);
      out.set_row(i, A_i);
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "outer_product", SUCCESS);
  }
  template <std::size_t OutColNb = RowNb>
  __host__ __device__ OpResult multiply(T v,
                                        MatN<T, RowNb, OutColNb> &out) const {
    // m x n \cdot  vmat (n x l) = out (m x l)
    // RowNb x ColNb \codt (n x l) = out (OutRowNb x
    // OutColNb)
    MatN<T, ColNb, OutColNb> vmat(v);

    auto r = multiply(vmat, out);
    if (r.status != SUCCESS)
      return r;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "multiply", SUCCESS);
  }
  /*matrix to matrix multiplication*/
  template <std::size_t OutColNb = RowNb>
  __host__ __device__ OpResult dot(const MatN<T, ColNb, OutColNb> &v,
                                   MatN<T, RowNb, OutColNb> &out) const {
    return multiply(v, out);
  }
  /*matrix to scalar multiplication*/
  template <std::size_t OutColNb = RowNb>
  __host__ __device__ OpResult dot(T v, MatN<T, RowNb, OutColNb> &out) const {
    return multiply(v, out);
  }
  /*matrix to vector multiplication*/
  __host__ __device__ OpResult dot(const T (&v)[ColNb],
                                   MatN<T, RowNb, 1> &out) const {
    MatN<T, ColNb, 1> vmat(v);
    auto r = multiply<1>(vmat, out);
    if (r.status != SUCCESS)
      return r;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "dot", SUCCESS);
  }

  /**
    m x n \cdot  vmat (n x l) = out (m x l)
    RowNb x ColNb \codt (OutRowNb x OutColNb) = out (RowNb x
    OutColNb)

    We are using the kij (row outer product) variant from
    Golub, van Loan 2013, p. 11 alg. 1.1.8 due to
    implementing this algorithm in C++. For fortran etc one
    should use jki since it access matrices by column.  For
    a comparison of algorithms see table 1.1.1 in p. 9

    tested
   */
  template <std::size_t OutColNb = RowNb>
  __host__ __device__ OpResult multiply(const MatN<T, ColNb, OutColNb> &B,
                                        MatN<T, RowNb, OutColNb> &out) const {

    // fill out matrix with zero
    out = MatN<T, RowNb, OutColNb>(static_cast<T>(0));
    for (std::size_t k = 0; k < ColNb; k++) {
      // x vector
      T A_k[RowNb];
      get_column(k, A_k);

      // y vector
      T B_k[OutColNb];
      B.get_row(k, B_k);

      // compute their outer product
      outer_product<RowNb, OutColNb>(A_k, B_k, out);
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "multiply", SUCCESS);
  }
  /**
    add row
   */
  __host__ __device__ OpResult add_row(const T (&r_data)[ColNb],
                                       MatN<T, RowNb + 1, ColNb> &out) const {
    return add_rows<ColNb>(r_data, out);
  }
  /**
    add rows if the incoming data has a size of multiple of
    number of columns
    of this array
  */
  template <std::size_t InRow = ColNb>
  __host__ __device__ OpResult
  add_rows(const T (&r_data)[InRow],
           MatN<T, RowNb + (InRow / ColNb), ColNb> &out) const {
    if ((InRow % ColNb) != 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "add_rows", SIZE_ERROR);
    }
    // fill output matrix with zeros
    out = MatN<T, RowNb + (InRow / ColNb), ColNb>(0);

    // fill with the output matrix with current matrix
    // elements
    std::size_t i = 0;
    std::size_t j = 0;
    for (i = 0; i < RowNb; i++) {
      for (j = 0; j < ColNb; j++) {
        T value = static_cast<T>(0);
        (*this)(i, j, value);                // getter
        out(make_matn_cell<T>(value, i, j)); // setter
      }
    }

    // fill from r_data the remaining values
    std::size_t nb_of_rows_to_add = static_cast<std::size_t>(InRow / ColNb);
    for (i = 0; i <= nb_of_rows_to_add; i++) {
      std::size_t row = RowNb + i;
      for (std::size_t j = 0; j < ColNb; j++) {
        T row_val = r_data[i * ColNb + j];
        out(make_matn_cell(row_val, row, j)); // setter
      }
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "add_rows", SUCCESS);
  }
  /**
    add column
  */
  __host__ __device__ OpResult
  add_column(const T (&r_data)[RowNb], MatN<T, RowNb, ColNb + 1> &out) const {
    return add_columns<RowNb>(r_data, out);
  }
  template <std::size_t InCol = RowNb>
  __host__ __device__ OpResult
  add_columns(const T (&c_data)[InCol],
              MatN<T, RowNb, ColNb + (InCol / RowNb)> &out) const {
    if ((InCol % RowNb) != 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "add_columns",
                      SIZE_ERROR);
    }
    // fill output matrix with zeros
    MatN<T, RowNb, ColNb + (InCol / RowNb)>::zeros(out);

    // fill with the output matrix with current matrix
    // elements
    std::size_t i = 0;
    std::size_t j = 0;
    for (i = 0; i < RowNb; i++) {
      for (j = 0; j < ColNb; j++) {
        T value = static_cast<T>(0);
        (*this)(i, j, value);             // getter
        out(make_matn_cell(value, i, j)); // setter
      }
    }
    // fill from c_data the remaining values
    std::size_t nb_of_cols_to_add = static_cast<std::size_t>(InCol / RowNb);

    // even if there are zero columns to add the output
    // should be one
    for (i = 0; i < nb_of_cols_to_add; i++) {
      std::size_t col = ColNb + i;
      for (j = 0; j < RowNb; j++) {
        T col_val = c_data[i * RowNb + j];
        out(make_matn_cell(col_val, j, col));
      }
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "add_columns", SUCCESS);
  }

private:
  template <typename Func>
  __host__ __device__ OpResult apply(const MatN<T, RowNb, ColNb> &vmat,
                                     const Func &fn,
                                     MatN<T, RowNb, ColNb> &out) const {
    for (std::size_t i = 0; i < Size; i++) {
      T tout = static_cast<T>(0);
      vmat(i, tout); // getter
      T val = fn(data[i], tout);
      auto r = out(make_matn_cell<T>(val, i));
      if (r.status != SUCCESS)
        return r;
    }

    return OpResult(__LINE__, __FILE__, __FUNCTION__, "apply", SUCCESS);
  }
  template <typename Func>
  __host__ __device__ OpResult apply(const Func &fn,
                                     MatN<T, RowNb, ColNb> &out) const {
    for (std::size_t i = 0; i < Size; i++) {
      T val = fn(data[i]);
      auto r = out(make_matn_cell<T>(val, i));
      if (r.status != SUCCESS)
        return r;
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "apply", SUCCESS);
  }
  template <typename Func>
  __host__ __device__ OpResult apply(const T &v, const Func &fn,
                                     MatN<T, RowNb, ColNb> &out) const {
    for (std::size_t i = 0; i < Size; i++) {
      T val = fn(data[i], v);
      auto r = out(make_matn_cell<T>(val, i));
      if (r.status != SUCCESS)
        return r;
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "apply", SUCCESS);
  }
};

template <typename T, std::size_t N>
__host__ __device__ OpResult identity(MatN<T, N, N> &out) {
  MatN<T, N, N> mat;
  for (std::size_t i = 0; i < N; i++) {
    mat(make_matn_cell(static_cast<T>(1), i, i));
  }
  out = mat;
  return OpResult(__LINE__, __FILE__, __FUNCTION__, "identity", SUCCESS);
}

} // namespace matn

namespace lu {

template <typename T, std::size_t N> struct LUdecomp {
  matn::MatN<T, N, N> LU;
  /**
    decompose the given matrix into two matrices, L, U, so
    that A=LU
    implemented from Golub and Van Loan 2013, Matrix
    computations, p. 116, algorithm 3.2.1
   */
  __host__ __device__ LUdecomp(const matn::MatN<T, N, N> &in_m) {
    T mdata[N * N];
    in_m(mdata);

    T ludata[N * N];
    memcpy(ludata, mdata, (N * N) * sizeof(T));

    //
    LU = matn::MatN<T, N, N>(ludata);
    for (std::size_t k = 0; k < (N - 1); ++k) {
      //

      for (std::size_t r = k + 1; r < N; ++r) {
        //
        T r_k;
        LU(r, k, r_k);
        T k_k;
        LU(k, k, k_k);

        // set L
        LU(matn::make_matn_cell(r_k / k_k, r, k));
      }
      for (std::size_t i = k + 1; i < N; ++i) {
        for (std::size_t j = k + 1; j < N; ++j) {
          //
          T i_k;
          LU(i, k, i_k);
          T k_j;
          LU(k, j, k_j);

          //
          T i_j;
          LU(i, j, i_j);
          LU(matn::make_matn_cell(i_j - i_k * k_j, i, j));
        }
      }
    }
  }

  // tested
  __host__ __device__ OpResult upper(matn::MatN<T, N, N> &U) const {
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = i; j < N; ++j) {
        T udata;
        LU(i, j, udata);
        U({.content = udata, .row = i, .column = j});
      }
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "upper", SUCCESS);
  }
  // tested
  __host__ __device__ OpResult lower(matn::MatN<T, N, N> &L) const {
    matn::identity(L);
    for (std::size_t i = 0; i < (N - 1); ++i) {
      for (std::size_t j = i + 1; j < N; ++j) {
        T udata;
        LU(j, i, udata);
        L({.content = udata, .row = j, .column = i});
      }
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "lower", SUCCESS);
  }

  // from Golub and Van Loan 2013, Matrix computations, p.
  // 108
  // tested
  __host__ __device__ OpResult solve_forward(const vecn::VecN<T, N> &b,
                                             vecn::VecN<T, N> &x) const {
    matn::MatN<T, N, N> L;
    OpResult res = lower(L);
    T b_0;
    b(0, b_0);
    T L_0;
    L(0, 0, L_0);
    if (L_0 == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "solve_forward",
                      ARG_ERROR);
    }
    x({.content = b_0 / L_0, .index = 0});
    for (std::size_t i = 1; i < N; ++i) {
      T x_i;
      b(i, x_i);
      x({.content = x_i, .index = i});

      // compute vdot
      T out = 0;
      for (std::size_t j = 0; j <= (i - 1); ++j) {
        T ij;
        L(i, j, ij);
        T x_j;
        x(j, x_j);
        out += ij * x_j;
      }

      x(i, x_i);
      T diff = x_i - out;

      T i_i;
      L(i, i, i_i);

      if (i_i == 0) {
        return OpResult(__LINE__, __FILE__, __FUNCTION__, "solve_forward",
                        ARG_ERROR);
      }
      x({.content = diff / i_i, .index = i});
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "solve_forward", SUCCESS);
  }

  // from Golub and Van Loan 2013, Matrix computations, p.
  // 108
  // tested
  __host__ __device__ OpResult solve_backward(const vecn::VecN<T, N> &b,
                                              vecn::VecN<T, N> &x) const {
    matn::MatN<T, N, N> U;
    OpResult res = upper(U);
    T u_n;
    U(N - 1, N - 1, u_n);
    if (u_n == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__, "solve_backward",
                      ARG_ERROR);
    }
    T b_n;
    b(N - 1, b_n);
    x({.content = b_n / u_n, .index = N - 1});

    std::size_t i = N - 2;

    while (i >= 0) {
      T b_i;
      b(i, b_i);
      x({.content = b_i, .index = i});

      // compute vdot
      T out = 0;
      for (std::size_t j = i + 1; j < N; ++j) {
        T ij;
        U(i, j, ij);
        T x_j;
        x(j, x_j);
        out += ij * x_j;
      }
      //
      T x_i;
      x(i, x_i);
      T i_i;
      U(i, i, i_i);

      if (i_i == 0) {
        return OpResult(__LINE__, __FILE__, __FUNCTION__, "solve_backward",
                        ARG_ERROR);
      }

      x({.content = (x_i - out) / i_i, .index = i});

      // size_t acts weird if we don't decrement it this way
      if (i == 0) {
        break;
      } else {
        --i;
      }
    }

    return OpResult(__LINE__, __FILE__, __FUNCTION__, "solve_backward",
                    SUCCESS);
  }
  __host__ __device__ OpResult solve(const vecn::VecN<T, N> &b,
                                     vecn::VecN<T, N> &x) {
    vecn::VecN<T, N> in_x;
    OpResult res = solve_forward(b, in_x);
    if (res.status != SUCCESS) {
      return res;
    }
    res = solve_backward(in_x, x);
    if (res.status != SUCCESS) {
      return res;
    }

    return OpResult(__LINE__, __FILE__, __FUNCTION__, "solve", SUCCESS);
  }

  // from Golub and Van Loan 2013, Matrix computations, p.
  // 108
  template <std::size_t Q>
  __host__ __device__ OpResult solve_mat(const matn::MatN<T, N, Q> &B,
                                         matn::MatN<T, N, Q> &X) {
    for (std::size_t j = 0; j < Q; ++j) {
      //
      T b_j[N];
      B.get_column(j, b_j);
      vecn::VecN<T, N> b(b_j);

      T x_j[N];
      X.get_column(j, x_j);
      //
      vecn::VecN<T, N> x(x_j);
      OpResult res = solve(b, x);
      if (res.status != SUCCESS) {
        return res;
      }
      // pass data to x_j
      x(x_j);
      X.set_column(j, x_j);
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "solve_mat", SUCCESS);
  }
};

} // namespace lu
template <typename T>
__host__ __device__ OpResult to_skew_mat(const vecn::VecN<T, 3> &vec,
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
  return OpResult(__LINE__, __FILE__, __FUNCTION__, "to_skew_mat", SUCCESS);
}

/**convert quaternion to rotation matrix
  from Vince, 2011, Quaternion ..., p. 123
 */
template <typename T>
__host__ __device__ OpResult toRotMat3x3(const quaternion::Quaternion<T> &q_in,
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
  return OpResult(__LINE__, __FILE__, __FUNCTION__, "toRotMat3x3", SUCCESS);
}

template <typename T, std::size_t N>
__host__ __device__ OpResult invert_mat(const matn::MatN<T, N, N> &m,
                                        matn::MatN<T, N, N> &out) {
  //
  lu::LUdecomp<T, N> lu_d(m);
  matn::MatN<T, N, N> B;
  matn::identity<T, N>(B);
  matn::MatN<T, N, N> inv_m;
  lu_d.solve_mat(B, inv_m);
  out = inv_m;
  return OpResult(__LINE__, __FILE__, __FUNCTION__, "invert_mat", SUCCESS);
}

} // namespace math3d

#endif
