#include <iostream>
#include <fstream>
#include <memory>
#include <eigen3/Eigen/Dense>

#include "src/grid.h"
#include "src/io.h"

using namespace std;
using namespace poisson;

static inline double gaussian_distr(const double x, const double y) {
  return 2.0*exp( -(std::pow(x-0.25, 2)+std::pow(y-0.5, 2))/0.01 );
}

static inline double hor_velocity(const double x, const double y) {
  return y-0.5;
}

static inline double ver_velocity(const double x, const double y) {
  return -(x-0.5);
}

const double WIDTH    = 1;
const double HEIGHT   = 1;
const size_t XDIM     = 100;
const size_t YDIM     = 100;
const double TIMESTEP = 0.1;
const size_t FRMS     = 500;

class advection_equation
{
public:
  typedef Eigen::Vector2d Vec2;

  advection_equation()
      : h_(TIMESTEP) {
    T_ = make_shared<grid2d<double, size_t>>(XDIM, YDIM, WIDTH, HEIGHT);
    u_ = make_shared<grid2d<double, size_t>>(XDIM, YDIM, WIDTH, HEIGHT);
    v_ = make_shared<grid2d<double, size_t>>(XDIM, YDIM, WIDTH, HEIGHT);
    T_->apply_func(gaussian_distr);
    u_->apply_func(hor_velocity);
    v_->apply_func(ver_velocity);
  }
  int solve(const double curr_time) {
    grid2d<double, size_t> next_u, next_v, next_T, half_u, half_v;

    next_T = *T_;
    
    // time invariant velocity field
    next_u = *u_;  
    next_v = *v_;

    // assumption
    half_u = *u_; // (u_+next_u)/2.;
    half_v = *v_; // (v_+next_v)/2.;
    
    // solve the internal points
    #pragma omp parallel for
    for (size_t i = 1; i < T_->nx()-1; ++i) {
      for (size_t j = 1; j < T_->ny()-1; ++j) {        
        Vec2 curr_p = T_->pos<Vec2>(i, j), prev_p;
        Vec2 curr_v = Vec2((*u_)(i, j), (*v_)(i, j));

        int iter = 5;
        while ( iter-- ) {
          prev_p = curr_p-0.5*h_*curr_v;
          curr_v[0] = half_u.interplinear(prev_p[0], prev_p[1]);
          curr_v[1] = half_v.interplinear(prev_p[0], prev_p[1]);
        }

        prev_p = curr_p-h_*curr_v;

        next_T(i, j) = T_->interpcubic(prev_p[0], prev_p[1]);
      }
    }

    *T_ = next_T;

    return 0;
  }
public:
  double h_;
  std::shared_ptr<grid2d<double, size_t>> T_, u_, v_;
};

int main(int argc, char *argv[])
{  
  advection_equation eqn;

  char outfile[256];
  for (size_t i = 0; i < FRMS; ++i) {
    cout << "# frame " << i << endl;
    sprintf(outfile, "./adv_frm_%zu.vtk", i);
    ofstream ofs(outfile);
    grid2vtk(ofs, XDIM, YDIM, WIDTH, HEIGHT);
    point_data(ofs, eqn.T_->data(), eqn.T_->dim(), "adv");
    ofs.close();

    eqn.solve(i*TIMESTEP);
  }

  cout << "# done" << endl;
  return 0;
}
