#include <iostream>
#include <fstream>

#include "problem.h"
#include "func.h"

using namespace std;
using namespace poisson;
using namespace Eigen;

#define XDIM   128
#define YDIM   128
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
    cerr << "# Usage: grid.dat\n";
    return __LINE__;
  }
  grid2d<double, int> *grid = new grid2d<double, int>(XDIM, YDIM, WIDTH, HEIGHT);
  problem<double, int> *prb = new laplace_equation<double, int>();
  prb->config_operator(grid);
  prb->set_dirichlet_bc(grid);
  prb->config_rhs(grid);
  prb->solve(grid);

  // analytic solution
  func2d<double> *fx = new stress_distribution<double>(1.0, WIDTH, HEIGHT);
  VectorXd xstar(grid->dim());
  for (int i = 0; i < grid->nx(); ++i) {
    for (int j = 0; j < grid->ny(); ++j) {
      double x = i*grid->dx(), y = j*grid->dy();
      xstar[grid->idx(i, j)] = (*fx)(x, y);
    }
  }
  Map<VectorXd> x(grid->data(), grid->dim());
  cout << "[INFO] error infinity norm: " << (xstar-x).lpNorm<Infinity>() << endl;

  dump(argv[1], grid);

  delete prb;
  delete grid;
  delete fx;

  cout << "[INFO] done\n";
  return 0;
}
