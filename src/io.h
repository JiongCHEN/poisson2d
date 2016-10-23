#ifndef IO_H
#define IO_H

#include <iostream>
#include <fstream>

namespace poisson {

template <typename OS, typename FLOAT, typename INT>
void grid2vtk(OS &os, const INT xdim, const INT ydim, const FLOAT width, const FLOAT height)
{
  os << "# vtk DataFile Version 2.0\nSample rectilinear grid\nASCII\nDATASET RECTILINEAR_GRID\n";
  os << "DIMENSIONS" << " " << xdim << " " << ydim << " " << 1 << std::endl;
  FLOAT dx = width/(xdim-1), dy = height/(ydim-1);
  os << "X_COORDINATES " << xdim << " float\n";
  for (INT i = 0; i < xdim; ++i) {
    os << i*dx << std::endl;
  }
  os << "Y_COORDINATES " << ydim << " float\n";
  for (INT i = 0; i < ydim; ++i) {
    os << i*dy << std::endl;
  }
  os << "Z_COORDINATES " << 1 << " float\n";
  os << 0 << std::endl;
}

template <typename OS, typename Iterator, typename INT>
void vtk_data(OS &os, Iterator first, INT size, const char *value_name, const char *table_name = "my_table")
{
  os << "SCALARS " << value_name << " float\nLOOKUP_TABLE " << table_name << "\n";
  for(size_t i = 0; i < size; ++i, ++first)
    os << *first << "\n";
}

template <typename OS, typename Iterator, typename INT>
void point_data(OS &os, Iterator first, INT size, const char *value_name, const char *table_name = "my_table")
{
  os << "POINT_DATA " << size << "\n";
  vtk_data(os, first, size, value_name, table_name);
}

}

#endif
