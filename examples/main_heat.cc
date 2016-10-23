#include <iostream>
#include <fstream>

#include "src/problem.h"
#include "src/func.h"
#include "src/io.h"

using namespace std;
using namespace poisson;
using namespace Eigen;

#define XDIM   100
#define YDIM   100
#define WIDTH  1.0
#define HEIGHT 1.0
#define TIMESTEP 0.001
#define T 1

int main(int argc, char *argv[])
{
  grid2d<double, int> *grid = new grid2d<double, int>(XDIM, YDIM, WIDTH, HEIGHT);
  heat_equation<double, int> *prb = new heat_equation<double, int>(TIMESTEP, T);

  prb->config_operator(grid);
  prb->set_initial_cond(grid);
//  prb->set_dirichlet_bc(grid); // time invariant here

  int count = 0;
  char outfile[128];
  while ( true ) {
    cout << "[Info] frame " << count << endl;
    sprintf(outfile, "./frame_%d.vtk", count);
    ++count;
    ofstream os(outfile);
    grid2vtk(os, XDIM, YDIM, WIDTH, HEIGHT);
    point_data(os, grid->data(), grid->dim(), "heat");
    os.close();
    prb->config_rhs(grid);
    prb->solve(grid);
    if ( prb->stop() )
      break;
  }

  delete prb;
  delete grid;
  cout << "[INFO] done\n";
  return 0;
}
