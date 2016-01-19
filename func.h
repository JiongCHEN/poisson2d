#ifndef FUNC_H
#define FUNC_H

#include <cmath>

namespace poisson {

template <typename T>
class func2d {
public:
  virtual T operator ()(const T x, const T y) = 0;
};

template <typename T>
class your_func : public func2d<T> {
public:
  T operator ()(const T x, const T y) {
    return 0.0;
  }
};

template <typename T>
class distance_to_origin : public func2d<T> {
public:
  T operator ()(const T x, const T y) {
    return sqrt(x*x+y*y);
  }
};

}
#endif
