#ifndef GRID_H
#define GRID_H

#include <algorithm>

namespace poisson {

/// @brief 2D regular grid on [0, w]x[0, h]
template <typename T1, typename T2>
class grid2d
{
public:
  grid2d(const T2 nx, const T2 ny, const T1 w, const T1 h)
    : nx_(nx), ny_(ny), dim_(nx*ny), w_(w), h_(h), dx_(w/(nx-1)), dy_(h/(ny-1)) {
    data_ = (T1 *)malloc(sizeof(T1)*dim_);
    memset(data_, 0, sizeof(T1)*dim_);
  }
  ~grid2d() {
    free(data_);
  }
  T2 nx() const {
    return nx_;
  }
  T2 ny() const {
    return ny_;
  }
  T2 dim() const {
    return dim_;
  }
  T1 width() const {
    return w_;
  }
  T1 height() const {
    return h_;
  }
  T1 dx() const {
    return dx_;
  }
  T1 dy() const {
    return dy_;
  }
  T1* data() const {
    return data_;
  }
  void clamp(T2 &i, T2 &j) const {
    i = std::max(0, std::min(nx_-1, i));
    j = std::max(0, std::min(ny_-1, j));
  }
  T2 idx(const T2 i, const T2 j) const {
    return j*nx_+i;
  }
  T1& operator ()(const T2 i, const T2 j) {
    return data_[idx(i, j)];
  }
  T1 operator ()(const T2 i, const T2 j) const {
    return data_[idx(i, j)];
  }
  T1 sample(const T1 x, const T1 y) {
    T1 xx = x/dx_, yy = y/dy_;
    T1 ax = xx-floor(xx), ay = yy-floor(yy);
    T1 bx = 1.0-ax, by = 1.0-ay;
    T2 x0 = (T2)floor(xx), y0 = (T2)floor(yy);
    T2 x1 = x0+1, y1 = y0+1;
    clamp(x0, y0);
    clamp(x1, y1);
    return by*(bx*(*this)(x0, y0)+ax*(*this)(x1, y0))+
           ay*(bx*(*this)(x0, y1)+ax*(*this)(x1, y1));
  }
private:
  const T2 nx_, ny_, dim_;
  const T1 w_, h_;
  const T1 dx_, dy_;
  T1* data_;
};

}
#endif
