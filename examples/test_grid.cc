#include <iostream>
#include <fstream>
#include <algorithm>

#include "src/grid.h"

using namespace std;
using namespace poisson;

int main(int argc, char *argv[])
{
  grid2d<double, size_t> g(8, 8, 1, 1);

  size_t cnt = 0;
  for (auto it = g.begin(); it != g.end(); ++it)
    *it = static_cast<grid2d<double, size_t>::value_type>(cnt++);
  
  std::sort(g.begin(), g.end(), std::greater<double>());

  for (auto it = g.begin(); it != g.end(); ++it)
    cout << *it << endl;
  
  cout << "[INFO] done" << endl;
  return 0;
}
