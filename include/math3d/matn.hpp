#ifndef MATN_HPP
#define MATN_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <math.h>
#include <ostream>
#include <stdio.h>
#include <vector>

#include "opflags.h"

namespace math3d {

namespace matn {

template <class T = float, unsigned int RowNb = 1,
          unsigned int ColNb = RowNb>
class MatN {
  /** holds the vector data*/
  std::array<T, ColNb * RowNb> data;
  static const unsigned int nb_rows = RowNb;
  static const unsigned int nb_cols = ColNb;
  static const unsigned int size = ColNb * RowNb;

public:
  MatN() : data.fill(0) {}
  MatN(const std::array<T, RowNb * ColNb> &vd) : data(vd) {}
  MatN(const std::vector<std::vector<T>> &vdata) {

    int rdiff = vdata.size() - RowNb;
    int cdiff = vdata[0].size() - ColNb;
    unsigned int i = 0;
    unsigned int j = 0;

    if (rdiff >= 0 && cdiff >= 0) {
      // row and col size is bigger than current matrix
      for (i = 0; i < RowNb; i++) {
        for (j = 0; j < ColNb; j++) {
          this(i, j, vdata[i][j]);
        }
      }
    } else if (rdiff < 0 && cdiff >= 0) {
      // row size smaller col size is bigger than current
      // matrix
      for (i = 0; i < RowNb; i++) {
        for (j = 0; j < ColNb; j++) {
          T val = static_cast<T>(0);
          if (i < vdata.size()) {
            val = vdata[i][j];
          }
          this(i, j, val);
        }
      }
    } else if (rdiff < 0 && cdiff < 0) {
      // row and col sizes are smaller than current matrix
      for (i = 0; i < RowNb; i++) {
        for (j = 0; j < ColNb; j++) {
          T val = static_cast<T>(0);
          if (i < vdata.size() && j < vdata[0].size()) {
            val = vdata[i][j];
          }
          this(i, j, val);
        }
      }
    } else {
      // row size is bigger and col size is smaller than
      // current matrix
      for (i = 0; i < RowNb; i++) {
        for (j = 0; j < ColNb; j++) {
          T val = static_cast<T>(0);
          if (j < vdata[0].size()) {
            val = vdata[i][j];
          }
          this(i, j, val);
        }
      }
    }
    //
  }
  MatN(const std::vector<T> &vdata) {
    int size_nb = static_cast<unsigned int>(vdata.size()) -
                  RowNb * ColNb;
    unsigned int i = 0;
    if (size_nb < 0) {
      // vector size is smaller than matrix size
      for (i = 0; i < size; i++) {
        if (i < vdata.size()) {
          this(i, vdata[i]);
        } else {
          this(i, static_cast<T>(0));
        }
      }
    } else {
      // vector size is bigger than matrix size
      for (i = 0; i < size; i++) {
        this(i, vdata[i]);
      }
    }

    // lu_decomposition = LUdcmp<T, RowNb, ColNb>(data);
  }
  MatN(T fill_value) {
    for (unsigned int i = 0; i < size; i++) {
      this(i, fill_value);
    }
  }
  /**\brief Create matrix based on argument matrix*/
  template <unsigned int OutRowNb = RowNb,
            unsigned int OutColNb = ColNb>
  static OpResult
  from_row_cols(MatN<T, OutRowNb, OutColNb> &out) {
    out = MatN<T, OutRowNb, OutColNb>(static_cast<T>(0));
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "from_row_cols", SUCCESS);
  }
  template <unsigned int OutRowNb = RowNb,
            unsigned int OutColNb = ColNb>
  static OpResult
  from_row_cols(T v, MatN<T, OutRowNb, OutColNb> &out) {
    MatN<T, OutRowNb, OutColNb> mat(v);
    out = mat;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "from_row_cols", SUCCESS);
  }
  template <unsigned int OutRowNb = RowNb,
            unsigned int OutColNb = ColNb>
  static OpResult
  identity(unsigned int nb,
           MatN<T, OutRowNb, OutColNb> &out) {
    MatN<T, OutRowNb, OutColNb> mat;
    auto r = from_row_cols<OutRowNb, OutColNb>(mat);
    if (r.status != SUCCESS)
      return r;
    for (unsigned int i = 0; i < nb; i++) {
      mat(i, i, static_cast<T>(1));
    }
    out = mat;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "identity", SUCCESS);
  }
  OpResult apply(const MatN<T, RowNb, ColNb> &vmat,
                 const std::function<T(T, T)> &fn,
                 MatN<T, RowNb, ColNb> &out) const {
    for (unsigned int i = 0; i < data.size(); i++) {
      T tout = static_cast<T>(0);
      vmat[i, tout]; // getter
      T val = fn(data[i], tout);
      auto r = out(i, val);
      if (r.status != SUCCESS)
        return r;
    }

    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "apply", SUCCESS);
  }
  OpResult apply(const std::function<T(T)> &fn,
                 MatN<T, RowNb, ColNb> &out) const {
    for (unsigned int i = 0; i < data.size(); i++) {
      T val = fn(data[i]);
      auto r = out(i, val);
      if (r.status != SUCCESS)
        return r;
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "apply", SUCCESS);
  }
  OpResult apply(const T &v,
                 const std::function<T(T, T)> &fn,
                 MatN<T, RowNb, ColNb> &out) const {
    for (unsigned int i = 0; i < data.size(); i++) {
      T val = fn(data[i], v);
      auto r = out(i, val);
      if (r.status != SUCCESS)
        return r;
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "apply", SUCCESS);
  }
  // tested
  template <unsigned int OutRowNb = RowNb,
            unsigned int OutColNb = ColNb>
  OpResult fill(T v,
                MatN<T, OutRowNb, OutColNb> &out) const {
    unsigned int s = 0;
    out.get_size(s);
    for (unsigned int i = 0; i < s; i++) {
      out(i, v);
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "fill", SUCCESS);
  }
  // tested
  OpResult transpose(MatN<T, ColNb, RowNb> &out) const {

    for (unsigned int i = 0; i < nb_rows; i++) {
      for (unsigned int j = 0; j < nb_cols; j++) {
        T tout = static_cast<T>(0);
        this[i, j, tout]; // getter
        out(j, i, tout);  // setter
      }
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "transpose", SUCCESS);
  }

  OpResult operator[](unsigned int &rows,
                      unsigned int &cols) const {
    rows = nb_rows;
    cols = nb_cols;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "get_rows_cols", SUCCESS);
  }

  // tested
  OpResult operator[](unsigned int &out) const {
    out = size;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "get_size", SUCCESS);
  }
  OpResult operator[](unsigned int row, unsigned int col,
                      T &out) const {
    unsigned int index = row * nb_cols + col;
    OpResult r = this[index, out];
    if (r.status != SUCCESS)
      return r;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "get",
                    SUCCESS);
  }
  OpResult operator[](unsigned int index, T &out) const {
    if (index >= data.size())
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "get", INDEX_ERROR);
    out = data[index];
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "get",
                    SUCCESS);
  }
  OpResult
  operator[](std::array<T, RowNb * ColNb> &out) const {
    out = data;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "get",
                    SUCCESS);
  }
  OpResult operator()(unsigned int row, unsigned int col,
                      const T &el) {
    unsigned int index = row * nb_cols + col;
    auto r = this(index, el);
    if (r.status != SUCCESS)
      return r;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "set",
                    SUCCESS);
  }
  OpResult operator()(unsigned int index, const T &el) {
    if (index >= data.size())
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "set", INDEX_ERROR);

    data[index] = el;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "set",
                    SUCCESS);
  }
  OpResult get_column(unsigned int index,
                      std::array<T, RowNb> &out) const {
    if (index >= ColNb) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "get_column", INDEX_ERROR);
    }
    for (unsigned int i = 0; i < RowNb; i++) {
      out[i] = data[i * ColNb + index];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "get_column", SUCCESS);
  }
  OpResult set_column(unsigned int index,
                      const std::array<T, RowNb> &idata) {
    if (index >= ColNb) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "set_column", INDEX_ERROR);
    }
    for (unsigned int i = 0; i < RowNb; i++) {
      data[i * ColNb + index] = idata[i];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "set_column", SUCCESS);
  }
  OpResult get_row(unsigned int index,
                   std::array<T, ColNb> &out) const {
    if (index >= RowNb) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "get_row", INDEX_ERROR);
    }
    for (unsigned int i = 0; i < ColNb; i++) {
      out[i] = data[index * ColNb + i];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "get_row", SUCCESS);
  }
  OpResult set_row(unsigned int index,
                   const std::array<T, ColNb> &idata) {
    if (index >= RowNb) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "set_row", INDEX_ERROR);
    }
    for (unsigned int i = 0; i < ColNb; i++) {
      data[index * ColNb + i] = idata[i];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "set_row", SUCCESS);
  }

  /**Obtain submatrix TODO*/
  OpResult submat(unsigned int row_start,
                  unsigned int col_start,
                  MatN<T, RowNb, ColNb> &out) const {
    unsigned int row_size = nb_rows - row_start;
    unsigned int col_size = nb_cols - col_start;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "submat", NOT_IMPLEMENTED);
  }
  OpResult add(const MatN<T, RowNb, ColNb> &v,
               MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv + val; };
    return apply(v, fn, out);
  }
  OpResult add(T v, MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv + val; };
    return apply(v, fn, out);
  }
  OpResult subtract(const MatN<T, RowNb, ColNb> &v,
                    MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv - val; };
    return apply(v, fn, out);
  }
  OpResult subtract(T v, MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv - val; };
    return apply(v, fn, out);
  }
  OpResult
  hadamard_product(const MatN<T, RowNb, ColNb> &v,
                   MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv * val; };
    return apply(v, fn, out);
  }
  OpResult
  hadamard_product(T v, MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv * val; };
    return apply(v, fn, out);
  }
  OpResult divide(const MatN<T, RowNb, ColNb> &v,
                  MatN<T, RowNb, ColNb> &out) const {
    unsigned int osize = 0;
    v[osize];
    for (unsigned int i = 0; i < osize; i++) {
      T tout = static_cast<T>(0);
      v[i, tout]; // getter
      if (tout == static_cast<T>(0)) {
        // zero division risk
        return OpResult(__LINE__, __FILE__, __FUNCTION__,
                        "divide", ARG_ERROR);
      }
    }
    auto fn = [](T matv, T val) { return matv / val; };
    return apply(v, fn, out);
  }
  OpResult divide(T v, MatN<T, RowNb, ColNb> &out) const {
    if (v == static_cast<T>(0)) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "divide", ARG_ERROR);
    }
    auto fn = [](T matv, T val) { return matv / val; };
    return apply(v, fn, out);
  }
  /**Declares inner vector product*/
  template <unsigned int N = RowNb>
  OpResult vdot(const std::array<T, N> &x,
                const std::array<T, N> &y, T &out) {
    if (N == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "vdot", SIZE_ERROR);
    }

    out = static_cast<T>(0);
    for (unsigned int i = 0; i < N; i++) {
      out += x[i] * y[i];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "vdot", SUCCESS);
  }

  /**Declares inner vector product with scalars*/
  template <unsigned int N = RowNb>
  OpResult vdot_s(const std::array<T, N> &x, const T &a,
                  const std::array<T, N> &out) {

    if (N == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "vdot_s", SIZE_ERROR);
    }
    for (unsigned int i = 0; i < N; i++) {
      out[i] = x[i] * a;
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "vdot_s", SUCCESS);
  }
  /**Implements saxpy algorithm from Golub, Van Loan 2013,
   * p. 4 alg.1.1.2*/
  template <unsigned int N = RowNb>
  OpResult saxpy(const T &a, const std::array<T, N> &x,
                 std::array<T, N> &y) const {
    if (N == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "saxpy", SIZE_ERROR);
    }
    for (unsigned int i = 0; i < N; i++) {
      y[i] += x[i] * a; //
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "saxpy", SUCCESS);
  }
  /**
    Implements gaxpy algorithm from Golub, Van Loan 2013, p.
    4 alg.1.1.3

    as specified in p. 6-7
   */
  OpResult gaxpy(const std::array<T, ColNb> &x,
                 std::array<T, RowNb> &y) const {
    for (unsigned int j = 0; j < ColNb; j++) {
      std::array<T, RowNb> c_j;
      get_column(j, c_j);
      saxpy<RowNb>(x[j], c_j, y);
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "gaxpy", SUCCESS);
  }
  /**
     Implements outer product update from Golub, Van Loan
     2013, p. 7 as a series of saxpy operations
    */
  template <unsigned int Rn, unsigned int Cn>
  OpResult outer_product(const std::array<T, Rn> &x,
                         const std::array<T, Cn> &y,
                         MatN<T, Rn, Cn> &out) const {
    for (unsigned int i = 0; i < Rn; i++) {
      std::array<T, Cn> A_i;
      out.get_row(i, A_i);
      saxpy<Cn>(x[i], y, A_i);
      out.set_row(i, A_i);
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "outer_product", SUCCESS);
  }
  template <unsigned int OutColNb = RowNb>
  OpResult multiply(T v,
                    MatN<T, RowNb, OutColNb> &out) const {
    // m x n \cdot  vmat (n x l) = out (m x l)
    // RowNb x ColNb \codt (n x l) = out (OutRowNb x
    // OutColNb)
    MatN<T, ColNb, OutColNb> vmat(v);

    auto r = multiply(vmat, out);
    if (r.status != SUCCESS)
      return r;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "multiply", SUCCESS);
  }
  /*matrix to matrix multiplication*/
  template <unsigned int OutColNb = RowNb>
  OpResult dot(const MatN<T, ColNb, OutColNb> &v,
               MatN<T, RowNb, OutColNb> &out) const {
    return multiply(v, out);
  }
  /*matrix to scalar multiplication*/
  template <unsigned int OutColNb = RowNb>
  OpResult dot(T v, MatN<T, RowNb, OutColNb> &out) const {
    return multiply(v, out);
  }
  /*matrix to vector multiplication*/
  OpResult dot(const std::array<T, ColNb> &v,
               MatN<T, RowNb, 1> &out) const {
    MatN<T, ColNb, 1> vmat(v);
    auto r = multiply<1>(vmat, out);
    if (r.status != SUCCESS)
      return r;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "dot",
                    SUCCESS);
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
  template <unsigned int OutColNb = RowNb>
  OpResult multiply(const MatN<T, ColNb, OutColNb> &B,
                    MatN<T, RowNb, OutColNb> &out) const {

    // fill out matrix with zero
    out = MatN<T, RowNb, OutColNb>(static_cast<T>(0));
    for (unsigned int k = 0; k < ColNb; k++) {
      // x vector
      std::array<T, RowNb> A_k;
      get_column(k, A_k);

      // y vector
      std::array<T, OutColNb> B_k;
      B.get_row(k, B_k);

      // compute their outer product
      outer_product<RowNb, OutColNb>(A_k, B_k, out);
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "multiply", SUCCESS);
  }
  /**
    add row
   */
  OpResult add_row(const std::array<T, ColNb> &r_data,
                   MatN<T, RowNb + 1, ColNb> &out) const {
    return add_rows<ColNb>(r_data, out);
  }
  /**
    add rows if the incoming data has a size of multiple of
    number of columns
    of this array
  */
  template <unsigned int InRow = ColNb>
  OpResult add_rows(
      const std::array<T, InRow> &r_data,
      MatN<T, RowNb + (InRow / ColNb), ColNb> &out) const {
    if ((InRow % ColNb) != 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "add_rows", SIZE_ERROR);
    }
    // fill output matrix with zeros
    from_row_cols(out);

    // fill with the output matrix with current matrix
    // elements
    unsigned int i = 0;
    unsigned int j = 0;
    for (i = 0; i < RowNb; i++) {
      for (j = 0; j < ColNb; j++) {
        T value = static_cast<T>(0);
        this[i, j, value]; // getter
        out(i, j, value);  // setter
      }
    }

    // fill from r_data the remaining values
    unsigned int nb_of_rows_to_add =
        static_cast<unsigned int>(InRow / ColNb);
    for (i = 0; i <= nb_of_rows_to_add; i++) {
      unsigned int row = RowNb + i;
      for (unsigned int j = 0; j < ColNb; j++) {
        T row_val = r_data[i * ColNb + j];
        out(row, j, row_val); // setter
      }
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "add_rows", SUCCESS);
  }
  /**
    add column
  */
  OpResult
  add_column(const std::array<T, RowNb> &r_data,
             MatN<T, RowNb, ColNb + 1> &out) const {
    return add_columns<RowNb>(r_data, out);
  }
  template <unsigned int InCol = RowNb>
  OpResult add_columns(
      const std::array<T, InCol> &c_data,
      MatN<T, RowNb, ColNb + (InCol / RowNb)> &out) const {
    if ((InCol % RowNb) != 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "add_columns", SIZE_ERROR);
    }
    // fill output matrix with zeros
    from_row_cols(out);

    // fill with the output matrix with current matrix
    // elements
    unsigned int i = 0;
    unsigned int j = 0;
    for (i = 0; i < RowNb; i++) {
      for (j = 0; j < ColNb; j++) {
        T value = static_cast<T>(0);
        this[i, j, value]; // getter
        out(i, j, value);  // setter
      }
    }
    // fill from c_data the remaining values
    unsigned int nb_of_cols_to_add =
        static_cast<unsigned int>(InCol / RowNb);

    // even if there are zero columns to add the output
    // should be one
    for (i = 0; i < nb_of_cols_to_add; i++) {
      unsigned int col = ColNb + i;
      for (j = 0; j < RowNb; j++) {
        T col_val = c_data[i * RowNb + j];
        out(j, col, col_val);
      }
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "add_columns", SUCCESS);
  }
  OpResult
  to_double_vec(std::vector<std::vector<T>> &ovec) const {
    //
    std::vector<std::vector<T>> out(
        nb_rows, std::vector<T>(nb_cols));
    for (unsigned int i = 0; i < nb_rows; i++) {
      for (unsigned int j = 0; j < nb_cols; j++) {
        this[i, j, out[i][j]];
      }
    }
    ovec = out;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "to_double_vec", SUCCESS);
  }
};

} // namespace matn

bool CHECK(OpResult r) { return r.status == SUCCESS; }

/**Checks if operation was successful*/
#define CHECK_MATN(call, res)                              \
  do {                                                     \
    res = call;                                            \
    res.call_name = #call;                                 \
    res.line_info = __LINE__;                              \
    res.file_name = __FILE__;                              \
  } while (0)

OpResult INFO(const MResult res) {
  std::cerr << res.status << " :: " << res.file_name
            << " :: " << res.line_info
            << " :: " << res.fn_name
            << " :: " << res.call_name << std::endl;
  return res;
}

/**Prints information about the operation*/
#define INFO_MATN(call, res)                               \
  do {                                                     \
    res = call;                                            \
    res.call_name = #call;                                 \
    res.line_info = __LINE__;                              \
    res.file_name = __FILE__;                              \
    if (res.status != SUCCESS) {                           \
      res = INFO(res);                                     \
    }                                                      \
  } while (0)

OpResult INFO_VERBOSE(MResult res) {
  if (res.status == SUCCESS) {
    std::cerr << "SUCCESS "
              << " :: " << res.file_name
              << " :: " << res.fn_name
              << " :: " << res.call_name
              << " :: " << res.duration_info
              << " microseconds" << std::endl;
    return res;
  }
  return INFO(res);
}

/**Prints everything about operation result*/
#define INFO_VERBOSE_MATN(call, res)                       \
  do {                                                     \
    auto start = std::chrono::steady_clock::now();         \
    res = call;                                            \
    auto stop = std::chrono::steady_clock::now();          \
    auto duration = std::chrono::duration_cast<            \
        std::chrono::microseconds>(stop - start);          \
    res.duration_info =                                    \
        to_string(static_cast<int>(duration.count()));     \
    res.call_name = #call;                                 \
    res.line_info = __LINE__;                              \
    res.file_name = __FILE__;                              \
    INFO_VERBOSE(res);                                     \
  } while (0)

} // namespace math3d

#endif
