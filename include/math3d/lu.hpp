#ifndef LU_HPP
#define LU_HPP
#include "matn.hpp"
#include "vecn.hpp"
namespace math3d {
namespace lu {

template <typename T, std::size_t N> struct LUdecomp {
  matn::MatN<T, N, N> LU;
  /**
    decompose the given matrix into two matrices, L, U, so
    that A=LU
    implemented from Golub and Van Loan 2013, Matrix
    computations, p. 116, algorithm 3.2.1
   */
  LUdecomp(const matn::MatN<T, N, N> &in_m) {
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
  OpResult upper(matn::MatN<T, N, N> &U) const {
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = i; j < N; ++j) {
        T udata;
        LU(i, j, udata);
        U({.content = udata, .row = i, .column = j});
      }
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "upper", SUCCESS);
  }
  // tested
  OpResult lower(matn::MatN<T, N, N> &L) const {
    matn::MatN<T, N, N>::identity(L);
    for (std::size_t i = 0; i < (N - 1); ++i) {
      for (std::size_t j = i + 1; j < N; ++j) {
        T udata;
        LU(j, i, udata);
        L({.content = udata, .row = j, .column = i});
      }
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "lower", SUCCESS);
  }

  // from Golub and Van Loan 2013, Matrix computations, p.
  // 108
  // tested
  OpResult solve_forward(const vecn::VecN<T, N> &b,
                         vecn::VecN<T, N> &x) const {
    matn::MatN<T, N, N> L;
    OpResult res = lower(L);
    T b_0;
    b(0, b_0);
    T L_0;
    L(0, 0, L_0);
    if (L_0 == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "solve_forward", ARG_ERROR);
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
        return OpResult(__LINE__, __FILE__, __FUNCTION__,
                        "solve_forward", ARG_ERROR);
      }
      x({.content = diff / i_i, .index = i});
    }
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "solve_forward", SUCCESS);
  }

  // from Golub and Van Loan 2013, Matrix computations, p.
  // 108
  // tested
  OpResult solve_backward(const vecn::VecN<T, N> &b,
                          vecn::VecN<T, N> &x) const {
    matn::MatN<T, N, N> U;
    OpResult res = upper(U);
    T u_n;
    U(N - 1, N - 1, u_n);
    if (u_n == 0) {
      return OpResult(__LINE__, __FILE__, __FUNCTION__,
                      "solve_backward", ARG_ERROR);
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
        return OpResult(__LINE__, __FILE__, __FUNCTION__,
                        "solve_backward", ARG_ERROR);
      }

      x({.content = (x_i - out) / i_i, .index = i});

      // size_t acts weird if we don't decrement it this way
      if (i == 0) {
        break;
      } else {
        --i;
      }
    }

    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "solve_backward", SUCCESS);
  }
  OpResult solve(const vecn::VecN<T, N> &b,
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

    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "solve", SUCCESS);
  }

  // from Golub and Van Loan 2013, Matrix computations, p.
  // 108
  template <unsigned int Q>
  OpResult solve_mat(const matn::MatN<T, N, Q> &B,
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
    return OpResult(__LINE__, __FILE__, __FUNCTION__,
                    "solve_mat", SUCCESS);
  }
};

} // namespace lu
} // namespace math3d
#endif
