#ifndef PROBLEM_H
#define PROBLEM_H

#include <unordered_set>
#include <eigen3/Eigen/Sparse>

#include "grid.h"
#include "func.h"
#include "util.h"

namespace poisson {

/// @brief problem description framework
template <typename T1, typename T2>
class problem
{
public:
  ~problem() {}
  virtual void set_initial_cond(grid2d<T1, T2> *g) = 0;
  virtual void config_operator(const grid2d<T1, T2> *g) = 0;
  virtual void config_rhs(const grid2d<T1, T2> *g) = 0;
  virtual void set_dirichlet_bc(grid2d<T1,T2> *g) = 0;
  virtual void set_neumann_bc(grid2d<T1, T2> *g) = 0;
  virtual int solve(grid2d<T1, T2> *g) = 0;
};

template <typename T1, typename T2>
class laplace_equation : public problem<T1, T2>
{
public:
  void set_initial_cond(grid2d<T1, T2> *g) {
      return;
  }
  void config_operator(const grid2d<T1, T2> *g) {
    laplacian(g, L_);
  }
  void set_dirichlet_bc(grid2d<T1, T2> *g) {
    fixDoF_.clear();
    for (T2 i = 0; i < g->nx(); ++i) {
      for (T2 j = 0; j < g->ny(); ++j) {
        if ( i == 0 ) {
          (*g)(i, j) = 0;
          fixDoF_.insert(g->idx(i, j));
        }
        if ( j == 0 ) {
          (*g)(i, j) = 0;
          fixDoF_.insert(g->idx(i, j));
        }
        if ( i == g->nx()-1 ) {
          (*g)(i, j) = 0;
          fixDoF_.insert(g->idx(i, j));
        }
        if ( j == g->ny()-1 ) {
          (*g)(i, j) = sin(M_PI*i*g->dx()/g->width());
          fixDoF_.insert(g->idx(i, j));
        }
      }
    }
    g2l_.resize(g->dim());
    T2 cnt = 0;
    for (T2 i = 0; i < g2l_.size(); ++i) {
      if ( fixDoF_.find(i) != fixDoF_.end() )
        g2l_[i] = static_cast<T2>(-1);
      else
        g2l_[i] = cnt++;
    }
  }
  void set_neumann_bc(grid2d<T1, T2> *g) {
    return;
  }
  void config_rhs(const grid2d<T1, T2> *g) {
    rhs_.setZero(g->dim());
  }
  int solve(grid2d<T1, T2> *g) {
    Eigen::Map<Eigen::Matrix<T1, -1, 1>> X(g->data(), g->dim());
    rhs_ -= L_*X;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<T1>> solver;
    if ( !fixDoF_.empty() ) {
      rm_spmat_row_col<T1, T2>(L_, g2l_);
      remove_vector_row<T1, T2>(rhs_, g2l_);
    }
    solver.compute(L_);
    assert(solver.info() == Eigen::Success);
    Eigen::Matrix<T1, -1, 1> dx = solver.solve(rhs_);
    assert(solver.info() == Eigen::Success);
    Eigen::Matrix<T1, -1, 1> DX;
    DX.setZero(g->dim());
    if ( !fixDoF_.empty() )
      recover_vector_row<T1, T2>(dx, g2l_, DX);
    else
      DX = dx;
    X += DX;
    return 0;
  }
private:
  void laplacian(const grid2d<T1, T2> *g, Eigen::SparseMatrix<T1> &L) const {
    std::vector<Eigen::Triplet<T1>> trips;
    const T1 dx2 = g->dx()*g->dx(), dy2 = g->dy()*g->dy();
    for (T2 i = 0; i < g->nx(); ++i) {
      for (T2 j = 0; j < g->ny(); ++j) {
        T2 mid = g->idx(i, j);
        T2 I, J, left, right, up, down; {
          trips.push_back(Eigen::Triplet<T1>(mid, mid, -2.0/dx2-2.0/dy2));
        } {
          I = i-1; J = j;
          g->clamp(I, J);
          left = g->idx(I, J);
          trips.push_back(Eigen::Triplet<T1>(mid, left, 1.0/dx2));
        } {
          I = i+1; J = j;
          g->clamp(I, J);
          right = g->idx(I, J);
          trips.push_back(Eigen::Triplet<T1>(mid, right, 1.0/dx2));
        } {
          I = i; J = j-1;
          g->clamp(I, J);
          down = g->idx(I, J);
          trips.push_back(Eigen::Triplet<T1>(mid, down, 1.0/dy2));
        } {
          I = i; J = j+1;
          g->clamp(I, J);
          up = g->idx(I, J);
          trips.push_back(Eigen::Triplet<T1>(mid, up, 1.0/dy2));
        }
      }
    }
    L.resize(g->dim(), g->dim());
    L.setFromTriplets(trips.begin(), trips.end());
  }
private:
  std::vector<T2> g2l_;
  Eigen::SparseMatrix<T1> L_;
  Eigen::Matrix<T1, -1, 1> rhs_;
  std::unordered_set<T2> fixDoF_;
};

template <typename T1, typename T2>
class poisson_equation : public problem<T1, T2>
{
};

/// solve for the equation $\frac{\partial u}{\partial t}=D\Delta u$
/// on domain [0, L]x[0, L]x[0, T]
/// with initial condition:
/// u(x, y, 0) = exp(-20(x-L/2)^2-20(x-L/2)^2)
/// and boundary conditions
/// u(0, y, t) = 0, u(L, y, t) = 0, u(x, 0, t) = 0, u(x, L, t) = 0.
template <typename T1, typename T2>
class heat_equation : public problem<T1, T2>
{
public:
    heat_equation(const T1 h, const T1 T) : h_(h), T_(T), clk_(0) {}
    void set_initial_cond(grid2d<T1, T2> *g) {
      const T1 dx = g->dx(), dy = g->dy();
      for (T2 i = 0; i < g->nx(); ++i)
        for (T2 j = 0; j < g->ny(); ++j)
          (*g)(i, j) = exp(-20*sq(i*dx-g->width()/2.0)-20*sq(j*dy-g->height()/2.0));
    }
    void config_operator(const grid2d<T1, T2> *g) {
      I_.resize(g->dim(), g->dim());
      I_.setIdentity();
      laplacian(g, L_);
    }
    void config_rhs(const grid2d<T1, T2> *g) {
      rhs_.resize(g->dim());
      std::copy(g->data(), g->data()+g->dim(), rhs_.data());
    }
    void set_dirichlet_bc(grid2d<T1, T2> *g) {
      fixDoF_.clear();
      for (T2 i = 0; i < g->nx(); ++i) {
        for (T2 j = 0; j < g->ny(); ++j) {
          if ( i == 0 ) {
            (*g)(i, j) = 0;
            fixDoF_.insert(g->idx(i, j));
          }
          if ( j == 0 ) {
            (*g)(i, j) = 0;
            fixDoF_.insert(g->idx(i, j));
          }
          if ( i == g->nx()-1 ) {
            (*g)(i, j) = 0;
            fixDoF_.insert(g->idx(i, j));
          }
          if ( j == g->ny()-1 ) {
            (*g)(i, j) = 0;
            fixDoF_.insert(g->idx(i, j));
          }
        }
      }
      g2l_.resize(g->dim());
      T2 cnt = 0;
      for (T2 i = 0; i < g2l_.size(); ++i) {
        if ( fixDoF_.find(i) != fixDoF_.end() )
          g2l_[i] = static_cast<T2>(-1);
        else
          g2l_[i] = cnt++;
      }
    }
    void set_neumann_bc(grid2d<T1, T2> *g) {
      return;
    }
    int solve(grid2d<T1, T2> *g) {
      Eigen::Map<Eigen::Matrix<T1, -1, 1>> X(g->data(), g->dim());
      Eigen::SparseMatrix<T1> LHS = I_-h_*L_;
      rhs_ -= LHS*X;
      Eigen::SimplicialCholesky<Eigen::SparseMatrix<T1>> solver;
      if ( !fixDoF_.empty() ) {
        rm_spmat_row_col<T1, T2>(LHS, g2l_);
        remove_vector_row<T1, T2>(rhs_, g2l_);
      }
      solver.compute(LHS);
      assert(solver.info() == Eigen::Success);
      Eigen::Matrix<T1, -1, 1> dx = solver.solve(rhs_);
      assert(solver.info() == Eigen::Success);
      Eigen::Matrix<T1, -1, 1> DX;
      DX.setZero(g->dim());
      if ( !fixDoF_.empty() )
        recover_vector_row<T1, T2>(dx, g2l_, DX);
      else
        DX = dx;
      X += DX;
      clk_ += h_;
      return 0;
    }
    bool stop() const {
      return clk_ >= T_;
    }
private:
    T1 sq(const T1 x) const {
      return x*x;
    }
    void laplacian(const grid2d<T1, T2> *g, Eigen::SparseMatrix<T1> &L) const {
      std::vector<Eigen::Triplet<T1>> trips;
      const T1 dx2 = g->dx()*g->dx(), dy2 = g->dy()*g->dy();
      for (T2 i = 0; i < g->nx(); ++i) {
        for (T2 j = 0; j < g->ny(); ++j) {
          T2 mid = g->idx(i, j);
          T2 I, J, left, right, up, down; {
            trips.push_back(Eigen::Triplet<T1>(mid, mid, -2.0/dx2-2.0/dy2));
          } {
            I = i-1; J = j;
            g->clamp(I, J);
            left = g->idx(I, J);
            trips.push_back(Eigen::Triplet<T1>(mid, left, 1.0/dx2));
          } {
            I = i+1; J = j;
            g->clamp(I, J);
            right = g->idx(I, J);
            trips.push_back(Eigen::Triplet<T1>(mid, right, 1.0/dx2));
          } {
            I = i; J = j-1;
            g->clamp(I, J);
            down = g->idx(I, J);
            trips.push_back(Eigen::Triplet<T1>(mid, down, 1.0/dy2));
          } {
            I = i; J = j+1;
            g->clamp(I, J);
            up = g->idx(I, J);
            trips.push_back(Eigen::Triplet<T1>(mid, up, 1.0/dy2));
          }
        }
      }
      L.resize(g->dim(), g->dim());
      L.setFromTriplets(trips.begin(), trips.end());
    }
private:
    const T1 h_, T_;
    T1 clk_;
    std::vector<T2> g2l_;
    Eigen::SparseMatrix<T1> L_, I_;
    Eigen::Matrix<T1, -1, 1> rhs_;
    std::unordered_set<T2> fixDoF_;
};

}
#endif
