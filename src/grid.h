#ifndef GRID_H
#define GRID_H

#include <cmath>
#include <cstring>
#include <algorithm>

#include "util.h"

namespace poisson {

/// @brief 2D regular grid on [0, w]x[0, h]
template <typename T1, typename T2>
class grid2d
{
public:
  typedef T1           value_type;
  typedef T2           size_type;
  typedef T1*          iterator;
  typedef const T1*    const_iterator;
  typedef T1&          reference;
  typedef const T1&    const_reference;

  grid2d()
      : nx_(0), ny_(0), dim_(0), w_(0), h_(0), dx_(0), dy_(0), data_(nullptr) {}
  grid2d(const T2 nx, const T2 ny, const T1 w, const T1 h)
    : nx_(nx), ny_(ny), dim_(nx*ny), w_(w), h_(h), dx_(w/(nx-1)), dy_(h/(ny-1)) {
    data_ = new T1[dim_];
    std::fill(data_, data_+dim_, 0);
  }
  virtual ~grid2d() {
    delete[] data_;
  }
  iterator begin() const {
    return data_;
  }
  iterator end() const {
    return data_+dim_;
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
  iterator data() const {
    return data_;
  }
  void clamp(T2 &i, T2 &j) const {
    i = std::max(static_cast<T2>(0), std::min(nx_-1, i));
    j = std::max(static_cast<T2>(0), std::min(ny_-1, j));
  }
  T2 idx(const T2 i, const T2 j) const {
    return j*nx_+i;
  }
  reference operator ()(const T2 i, const T2 j) {
    return data_[idx(i, j)];
  }
  T1 operator ()(const T2 i, const T2 j) const {
    return data_[idx(i, j)];
  }
  template <class Vec2d>
  Vec2d pos(const T2 i, const T2 j) {
    return Vec2d(i*dx_, j*dy_);
  }
  T1 interplinear(const T1 x, const T1 y) {
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
  T1 interpcubic(const T1 x, const T1 y) {
    T1 xx = x/dx_, yy = y/dy_;

    const T2 x1 = (T2)floor(xx);
    const T2 x2 = x1+1;
    const T2 x3 = x1+2;
    const T2 x0 = x1-1;

    const T2 y1 = (T2)floor(yy);
    const T2 y2 = y1+1;
    const T2 y3 = y1+2;
    const T2 y0 = y1-1;

    if ( x0 < 0 || y0 < 0 || x3 >= nx_ || y3 >= ny_ )
      return interplinear(x, y);

    const T1 ax = xx-floor(xx);
    const T1 bx = yy-floor(yy);

    const T1 p0[] = {(*this)(x0, y0), (*this)(x1, y0), (*this)(x2, y0), (*this)(x3, y0)};
    const T1 p1[] = {(*this)(x0, y1), (*this)(x1, y1), (*this)(x2, y1), (*this)(x3, y1)};
    const T1 p2[] = {(*this)(x0, y2), (*this)(x1, y2), (*this)(x2, y2), (*this)(x3, y2)};
    const T1 p3[] = {(*this)(x0, y3), (*this)(x1, y3), (*this)(x2, y3), (*this)(x3, y3)};

    const T1 p[] = {cubic_interp(ax, p0), cubic_interp(ax, p1), cubic_interp(ax, p2), cubic_interp(ax, p3)};
    return cubic_interp(bx, p);
  }
  bool is_valid(const T2 i, const T2 j) const {
    return ( i >= 0 && i < nx_ && j >= 0 && j < ny_  );
  }
  void apply_func(T1 (*fx)(const T1, const T1)) {
    for (T2 i = 0; i < nx_; ++i)
      for (T2 j = 0; j < ny_; ++j)
        (*this)(i, j) = fx(i*dx_, j*dy_);
  }
  grid2d<T1, T2>& operator =(const grid2d<T1, T2> &other) {
    nx_  = other.nx();
    ny_  = other.ny();
    dim_ = other.dim();
    w_   = other.width();
    h_   = other.height();
    dx_  = other.dx();
    dy_  = other.dy();
    if ( data_ )
      delete[] data_;
    data_ = new T1[dim_];
    std::copy(other.data(), other.data()+other.dim(), data_);
    return *this;
  }
private:
  T2 nx_, ny_, dim_;
  T1 w_, h_;
  T1 dx_, dy_;
  T1* data_;
};

}
#endif
