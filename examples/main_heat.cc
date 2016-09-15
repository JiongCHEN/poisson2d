#include <iostream>
#include <fstream>

#include "src/problem.h"
#include "src/func.h"

using namespace std;
using namespace poisson;
using namespace Eigen;

#define XDIM   100
#define YDIM   100
#define WIDTH  1.0
#define HEIGHT 1.0
#define TIMESTEP 0.001
#define T 1

template <typename OS, typename FLOAT, typename INT>
static void grid2vtk(OS &os, const INT xdim, const INT ydim, const FLOAT width, const FLOAT height)
{
  os << "# vtk DataFile Version 2.0\nSample rectilinear grid\nASCII\nDATASET RECTILINEAR_GRID\n";
  os << "DIMENSIONS" << " " << xdim << " " << ydim << " " << 1 << endl;
  FLOAT dx = width/(xdim-1), dy = height/(ydim-1);
  os << "X_COORDINATES " << xdim << " float\n";
  for (INT i = 0; i < xdim; ++i) {
    os << i*dx << endl;
  }
  os << "Y_COORDINATES " << ydim << " float\n";
  for (INT i = 0; i < ydim; ++i) {
    os << i*dy << endl;
  }
  os << "Z_COORDINATES " << 1 << " float\n";
  os << 0 << endl;
}

template <typename OS, typename Iterator, typename INT>
static void vtk_data(OS &os, Iterator first, INT size, const char *value_name, const char *table_name = "my_table")
{
  os << "SCALARS " << value_name << " float\nLOOKUP_TABLE " << table_name << "\n";
  for(size_t i = 0; i < size; ++i, ++first)
    os << *first << "\n";
}

template <typename OS, typename Iterator, typename INT>
static void point_data(OS &os, Iterator first, INT size, const char *value_name, const char *table_name = "my_table")
{
  os << "POINT_DATA " << size << "\n";
  vtk_data(os, first, size, value_name, table_name);
}

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
