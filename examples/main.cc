#include <iostream>
#include <fstream>
#include <memory>

#include "src/problem.h"
#include "src/func.h"

using namespace std;
using namespace poisson;
using namespace Eigen;

#define XDIM   16
#define YDIM   16
#define WIDTH  1.0
#define HEIGHT 1.0

template <typename T1, typename T2>
static int dump(const char *filename, const grid2d<T1, T2> *g) {
  ofstream ofs(filename);
  if ( ofs.fail() ) {
    cerr << "[INFO] cannot open " << filename << endl;
    return __LINE__;
  }
  ofs << "# X Y Z\n";
  T1 x, y, v;
  for (T2 i = 0; i < g->nx(); ++i) {
    for (T2 j = 0; j < g->ny(); ++j) {
      x = i*g->dx();
      y = j*g->dy();
      v = (*g)(i, j);
      ofs << x << " " << y << " " << v << "\n";
    }
    ofs << "\n";
  }
  ofs.close();
  return 0;
}

template <typename T>
class stress_distribution : public func2d<T> {
public:
  stress_distribution(const T w0, const T a, const T b)
    : w0_(w0), a_(a), b_(b) {}
  T operator ()(const T x, const T y) {
    return w0_/sinh(M_PI*b_/a_)*sin(M_PI*x/a_)*sinh(M_PI*y/a_);
  }
private:
  T w0_, a_, b_;
};

int main(int argc, char *argv[])
{
  if ( argc != 2 ) {
    cerr << "# Usage: main grid.dat\n";
    return __LINE__;
  }
  shared_ptr<grid2d<double, int>> grid = make_shared<grid2d<double, int>>(XDIM, YDIM, WIDTH, HEIGHT);
  shared_ptr<problem<double, int>> prb = make_shared<laplace_equation<double, int>>();
  prb->config_operator(grid.get());
  prb->set_dirichlet_bc(grid.get());
  prb->config_rhs(grid.get());
  prb->solve(grid.get());
  Map<const VectorXd> x(grid->data(), grid->dim());

  // sampling analytic solution
  shared_ptr<func2d<double>> fx = make_shared<stress_distribution<double>>(1.0, WIDTH, HEIGHT);
  VectorXd xstar(grid->dim());
  for (int i = 0; i < grid->nx(); ++i) {
    for (int j = 0; j < grid->ny(); ++j) {
      double x = i*grid->dx(), y = j*grid->dy();
      xstar[grid->idx(i, j)] = (*fx)(x, y);
    }
  }
  cout << "[INFO] error infinity norm: " << (xstar-x).lpNorm<Infinity>() << endl;

  dump(argv[1], grid.get());
  
  cout << "[INFO] done\n";
  return 0;
}
