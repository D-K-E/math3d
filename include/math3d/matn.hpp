#ifndef MATN_HPP
#define MATN_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iostream>
#include <math.h>
#include <ostream>
#include <stdio.h>
#include <vector>

#include "opflags.h"

namespace math3d {

namespace matn {

template <typename T> struct MatNCell {
  T content;
  std::size_t row = 0;
  std::size_t column = 0;
  int index = -1;
};

template <typename T>
MatNCell<T> make_matn_cell(const T &c,
                           std::size_t row_param,
                           std::size_t col_param) {
  MatNCell<T> cell{
      .content = c, .row = row_param, .column = col_param};
  return cell;
}

template <typename T>
MatNCell<T> make_matn_cell(const T &c, std::size_t i) {
  MatNCell<T> cell{.content = c,
                   .row = 0,
                   .column = 0,
                   .index = static_cast<int>(i)};
  return cell;
}

template <typename T>
MatNCell<T> make_matn_cell(
    const std::tuple<T, std::size_t, std::size_t> &tpl) {
  auto [c, row_param, col_param] = tpl;
  MatNCell<T> cell{
      .content = c, .row = row_param, .column = col_param};
  return cell;
}

template <class T = float, std::size_t RowNb = 1,
          std::size_t ColNb = RowNb>
class MatN {
  /** holds the vector data*/
  T data[ColNb * RowNb];
  static const std::size_t Size = ColNb * RowNb;

public:
  MatN() {
    for (std::size_t i = 0; i < Size; ++i) {
      data[i] = 0;
    }
  }
  MatN(const T (&vd)[RowNb * ColNb]) {
    memcpy(data, vd, Size * sizeof(T));
  }
  MatN(const std::vector<std::vector<T>> &vdata) {

    int rdiff = vdata.size() - RowNb;
    int cdiff = vdata[0].size() - ColNb;
    std::size_t i = 0;
    std::size_t j = 0;

    if (rdiff >= 0 && cdiff >= 0) {
      // row and col size is bigger than current matrix
      for (i = 0; i < RowNb; i++) {
        for (j = 0; j < ColNb; j++) {
          (*this)(make_matn_cell<T>(vdata[i][j], i, j));
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
          (*this)(make_matn_cell<T>(val, i, j));
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
          (*this)(make_matn_cell<T>(val, i, j));
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
          (*this)(make_matn_cell<T>(val, i, j));
        }
      }
    }
    //
  }
  MatN(const std::vector<T> &vdata) {
    int size_nb = vdata.size() - Size;
    std::size_t i = 0;
    if (size_nb < 0) {
      // vector size is smaller than matrix size
      for (i = 0; i < Size; i++) {
        if (i < vdata.size()) {
          (*this)(make_matn_cell<T>(vdata[i], i));
        } else {
          (*this)(make_matn_cell<T>(0, i));
        }
      }
    } else {
      // vector size is bigger than matrix size
      for (i = 0; i < Size; i++) {
        (*this)(make_matn_cell<T>(vdata[i], i));
      }
    }
  }

  MatN(T fill_value) {
    for (std::size_t i = 0; i < Size; ++i) {
      data[i] = fill_value;
    }
  }
  /**\brief Create matrix based on argument matrix*/
  template <std::size_t OutRowNb = RowNb,
            std::size_t OutColNb = ColNb>
  static OpResult
  from_row_cols(MatN<T, OutRowNb, OutColNb> &out) {
    out = MatN<T, OutRowNb, OutColNb>(static_cast<T>(0));
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "from_row_cols", SUCCESS);
  }
  template <std::size_t OutRowNb = RowNb,
            std::size_t OutColNb = ColNb>
  static OpResult
  from_row_cols(T v, MatN<T, OutRowNb, OutColNb> &out) {
    MatN<T, OutRowNb, OutColNb> mat(v);
    out = mat;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "from_row_cols", SUCCESS);
  }
  template <std::size_t OutRowNb = RowNb,
            std::size_t OutColNb = ColNb>
  static OpResult
  identity(std::size_t nb,
           MatN<T, OutRowNb, OutColNb> &out) {
    MatN<T, OutRowNb, OutColNb> mat;
    auto r = from_row_cols<OutRowNb, OutColNb>(mat);
    if (r.status != SUCCESS)
      return r;
    for (std::size_t i = 0; i < nb; i++) {
      mat(i, i, static_cast<T>(1));
    }
    out = mat;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "identity", SUCCESS);
  }
  OpResult apply(const MatN<T, RowNb, ColNb> &vmat,
                 const std::function<T(T, T)> &fn,
                 MatN<T, RowNb, ColNb> &out) const {
    for (std::size_t i = 0; i < Size; i++) {
      T tout = static_cast<T>(0);
      vmat(i, tout); // getter
      T val = fn(data[i], tout);
      auto r = out(make_matn_cell<T>(val, i));
      if (r.status != SUCCESS)
        return r;
    }

    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "apply", SUCCESS);
  }
  OpResult apply(const std::function<T(T)> &fn,
                 MatN<T, RowNb, ColNb> &out) const {
    for (std::size_t i = 0; i < Size; i++) {
      T val = fn(data[i]);
      auto r = out(make_matn_cell<T>(val, i));
      if (r.status != SUCCESS)
        return r;
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "apply", SUCCESS);
  }
  OpResult apply(const T &v,
                 const std::function<T(T, T)> &fn,
                 MatN<T, RowNb, ColNb> &out) const {
    for (std::size_t i = 0; i < Size; i++) {
      T val = fn(data[i], v);
      auto r = out(make_matn_cell<T>(val, i));
      if (r.status != SUCCESS)
        return r;
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "apply", SUCCESS);
  }
  // tested
  template <std::size_t OutRowNb = RowNb,
            std::size_t OutColNb = ColNb>
  OpResult fill(T v,
                MatN<T, OutRowNb, OutColNb> &out) const {
    std::size_t s = 0;
    out(s);
    for (std::size_t i = 0; i < s; i++) {
      out(i, v);
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "fill", SUCCESS);
  }
  // tested
  OpResult transpose(MatN<T, ColNb, RowNb> &out) const {

    for (std::size_t i = 0; i < RowNb; i++) {
      for (std::size_t j = 0; j < ColNb; j++) {
        T tout = static_cast<T>(0);
        (*this)(i, j, tout);                // getter
        out(make_matn_cell<T>(tout, i, j)); // setter
      }
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "transpose", SUCCESS);
  }

  OpResult operator()(std::size_t &rows,
                      std::size_t &cols) const {
    rows = RowNb;
    cols = ColNb;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "get_rows_cols", SUCCESS);
  }

  // tested
  OpResult operator()(std::size_t &out) const {
    out = Size;
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "get_size", SUCCESS);
  }
  OpResult operator()(std::size_t row, std::size_t col,
                      T &out) const {
    std::size_t index = row * ColNb + col;
    OpResult r = (*this)(index, out);
    if (r.status != SUCCESS)
      return r;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "get",
                    SUCCESS);
  }
  OpResult operator()(std::size_t index, T &out) const {
    if (index >= Size)
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "get", INDEX_ERROR);
    out = data[index];
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "get",
                    SUCCESS);
  }
  OpResult operator()(T (&out)[RowNb * ColNb]) const {
    memcpy(out, data, Size * sizeof(T));
    // for (std::size_t i = 0; i < size; ++i){
    //     out[i] = data[i];
    // }
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "get",
                    SUCCESS);
  }
  OpResult operator()(const MatNCell<T> &cell) {
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
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "set", INDEX_ERROR);
    }
    data[index] = el;
    return OpResult(__LINE__, __FILE__, __FUNCTION__, "set",
                    SUCCESS);
  }
  OpResult get_column(std::size_t index,
                      T (&out)[RowNb]) const {
    if (index >= ColNb) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "get_column", INDEX_ERROR);
    }
    for (std::size_t i = 0; i < RowNb; i++) {
      out[i] = data[i * ColNb + index];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "get_column", SUCCESS);
  }
  OpResult set_column(std::size_t index,
                      const T (&idata)[RowNb]) {
    if (index >= ColNb) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "set_column", INDEX_ERROR);
    }
    for (std::size_t i = 0; i < RowNb; i++) {
      data[i * ColNb + index] = idata[i];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "set_column", SUCCESS);
  }
  OpResult get_row(std::size_t index,
                   T (&out)[ColNb]) const {
    if (index >= RowNb) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "get_row", INDEX_ERROR);
    }
    for (std::size_t i = 0; i < ColNb; i++) {
      out[i] = data[index * ColNb + i];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "get_row", SUCCESS);
  }
  OpResult set_row(std::size_t index,
                   const T (&idata)[ColNb]) {
    if (index >= RowNb) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "set_row", INDEX_ERROR);
    }
    for (std::size_t i = 0; i < ColNb; i++) {
      data[index * ColNb + i] = idata[i];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "set_row", SUCCESS);
  }

  /**Obtain submatrix TODO*/
  OpResult submat(std::size_t row_start,
                  std::size_t col_start,
                  MatN<T, RowNb, ColNb> &out) const {
    std::size_t row_size = RowNb - row_start;
    std::size_t col_size = ColNb - col_start;
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
    std::size_t osize = 0;
    v(osize);
    for (std::size_t i = 0; i < osize; i++) {
      T tout = static_cast<T>(0);
      v(i, tout); // getter
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
  template <std::size_t N = RowNb>
  OpResult vdot(const T (&x)[N], const T (&y)[N], T &out) {
    if (N == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "vdot", SIZE_ERROR);
    }

    out = static_cast<T>(0);
    for (std::size_t i = 0; i < N; i++) {
      out += x[i] * y[i];
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "vdot", SUCCESS);
  }

  /**Declares inner vector product with scalars*/
  template <std::size_t N = RowNb>
  OpResult vdot_s(const T (&x)[N], const T &a,
                  T (&out)[N]) {

    if (N == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "vdot_s", SIZE_ERROR);
    }
    for (std::size_t i = 0; i < N; i++) {
      out[i] = x[i] * a;
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "vdot_s", SUCCESS);
  }
  /**Implements saxpy algorithm from Golub, Van Loan 2013,
   * p. 4 alg.1.1.2*/
  template <std::size_t N = RowNb>
  OpResult saxpy(const T &a, const T (&x)[N],
                 T (&y)[N]) const {
    if (N == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "saxpy", SIZE_ERROR);
    }
    for (std::size_t i = 0; i < N; i++) {
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
  OpResult gaxpy(const T (&x)[ColNb], T (&y)[RowNb]) const {
    for (std::size_t j = 0; j < ColNb; j++) {
      T c_j[RowNb];
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
  template <std::size_t Rn, std::size_t Cn>
  OpResult outer_product(const T (&x)[Rn], const T (&y)[Cn],
                         MatN<T, Rn, Cn> &out) const {
    for (std::size_t i = 0; i < Rn; i++) {
      T A_i[Cn];
      out.get_row(i, A_i);
      saxpy<Cn>(x[i], y, A_i);
      out.set_row(i, A_i);
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "outer_product", SUCCESS);
  }
  template <std::size_t OutColNb = RowNb>
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
  template <std::size_t OutColNb = RowNb>
  OpResult dot(const MatN<T, ColNb, OutColNb> &v,
               MatN<T, RowNb, OutColNb> &out) const {
    return multiply(v, out);
  }
  /*matrix to scalar multiplication*/
  template <std::size_t OutColNb = RowNb>
  OpResult dot(T v, MatN<T, RowNb, OutColNb> &out) const {
    return multiply(v, out);
  }
  /*matrix to vector multiplication*/
  OpResult dot(const T (&v)[ColNb],
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
  template <std::size_t OutColNb = RowNb>
  OpResult multiply(const MatN<T, ColNb, OutColNb> &B,
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
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "multiply", SUCCESS);
  }
  /**
    add row
   */
  OpResult add_row(const T (&r_data)[ColNb],
                   MatN<T, RowNb + 1, ColNb> &out) const {
    return add_rows<ColNb>(r_data, out);
  }
  /**
    add rows if the incoming data has a size of multiple of
    number of columns
    of this array
  */
  template <std::size_t InRow = ColNb>
  OpResult add_rows(
      const T (&r_data)[InRow],
      MatN<T, RowNb + (InRow / ColNb), ColNb> &out) const {
    if ((InRow % ColNb) != 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "add_rows", SIZE_ERROR);
    }
    // fill output matrix with zeros
    from_row_cols(out);

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
    std::size_t nb_of_rows_to_add =
        static_cast<std::size_t>(InRow / ColNb);
    for (i = 0; i <= nb_of_rows_to_add; i++) {
      std::size_t row = RowNb + i;
      for (std::size_t j = 0; j < ColNb; j++) {
        T row_val = r_data[i * ColNb + j];
        out(make_matn_cell(row_val, row, j)); // setter
      }
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "add_rows", SUCCESS);
  }
  /**
    add column
  */
  OpResult
  add_column(const T (&r_data)[RowNb],
             MatN<T, RowNb, ColNb + 1> &out) const {
    return add_columns<RowNb>(r_data, out);
  }
  template <std::size_t InCol = RowNb>
  OpResult add_columns(
      const T (&c_data)[InCol],
      MatN<T, RowNb, ColNb + (InCol / RowNb)> &out) const {
    if ((InCol % RowNb) != 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "add_columns", SIZE_ERROR);
    }
    // fill output matrix with zeros
    from_row_cols(out);

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
    std::size_t nb_of_cols_to_add =
        static_cast<std::size_t>(InCol / RowNb);

    // even if there are zero columns to add the output
    // should be one
    for (i = 0; i < nb_of_cols_to_add; i++) {
      std::size_t col = ColNb + i;
      for (j = 0; j < RowNb; j++) {
        T col_val = c_data[i * RowNb + j];
        out(make_matn_cell(col_val, j, col));
      }
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "add_columns", SUCCESS);
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

OpResult INFO(const OpResult res) {
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

OpResult INFO_VERBOSE(OpResult res) {
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
