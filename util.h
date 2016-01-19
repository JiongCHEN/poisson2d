#ifndef UTIL_H
#define UTIL_H

#include <eigen3/Eigen/Sparse>

namespace poisson {

template <typename T1, typename T2>
void rm_spmat_row_col(Eigen::SparseMatrix<T1> &A, const std::vector<T2> &g2l) {
  T2 new_size = 0;
  for (T2 i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != static_cast<T2>(-1) )
      ++new_size;
  }
  std::vector<Eigen::Triplet<T1>> trips;
  for (T2 j = 0; j < A.outerSize(); ++j) {
    for (typename Eigen::SparseMatrix<T1>::InnerIterator it(A, j); it; ++it) {
      if ( g2l[it.row()] != static_cast<T2>(-1) && g2l[it.col()] != static_cast<T2>(-1) )
        trips.push_back(Eigen::Triplet<T1>(g2l[it.row()], g2l[it.col()], it.value()));
    }
  }
  A.resize(new_size, new_size);
  A.reserve(trips.size());
  A.setFromTriplets(trips.begin(), trips.end());
}

template <typename T1, typename T2>
void remove_vector_row(Eigen::Matrix<T1, -1, 1> &b, const std::vector<T2> &g2l) {
  T2 new_size = 0;
  for (T2 i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != static_cast<T2>(-1) )
      ++new_size;
  }
  Eigen::Matrix<T1, -1, 1> sub;
  sub.resize(new_size);
#pragma omp parallel for
  for (T2 i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != static_cast<T2>(-1) )
      sub[g2l[i]] = b[i];
  }
  b = sub;
}

template <typename T1, typename T2>
void recover_vector_row(const Eigen::Matrix<T1, -1, 1> &l, const std::vector<T2> &g2l, Eigen::Matrix<T1, -1, 1> &g) {
#pragma omp parallel for
  for (T2 i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != static_cast<T2>(-1) )
      g[i] = l[g2l[i]];
  }
}

}
#endif
