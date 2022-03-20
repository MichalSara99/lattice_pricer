/**
 * @file lattice_interpolators.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief  Interpolators
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_INTERPOLATORS_HPP_)
#define _LATTICE_INTERPOLATORS_HPP_

#include "lattice_logging.hpp"
#include <iostream>
#include <sstream>

namespace lattice {

/**
 * @brief Linear interpolation object
 *
 * @tparam container
 */
template <typename container> class linear_interpolator {
private:
  container x_;
  container y_;

public:
  typedef container container_type;
  typedef typename container::value_type container_value_type;

  linear_interpolator() : x_(), y_() {}
  linear_interpolator(linear_interpolator const &cpy)
      : x_(cpy.x_), y_(cpy.y_) {}
  linear_interpolator(linear_interpolator &&other) noexcept
      : x_(std::move(other.x_)), y_(std::move(other.y_)) {}
  ~linear_interpolator() {}

  linear_interpolator &operator=(linear_interpolator const &cpy) {
    if (this != &cpy) {
      x_ = cpy.x_;
      y_ = cpy.y_;
    }
    return *this;
  }
  linear_interpolator &operator=(linear_interpolator &&other) noexcept {
    if (this != &other) {
      x_ = std::move(other.x_);
      y_ = std::move(other.y_);
    }
    return *this;
  }

  void set_points(container const &xpoints, container const &ypoints,
                  bool sort_x_points = false) {
    x_ = xpoints;
    y_ = ypoints;

    if (sort_x_points)
      std::sort(x_.begin(), x_.end());

    for (std::size_t i = 0; i < x_.size(); ++i) {
      for (std::size_t j = 0; j < x_.size(); ++j) {
        if (x_[i] == xpoints[j]) {
          y_[i] = ypoints[j];
          break;
        }
      }
    }
  }

  container_value_type const get_value(container_value_type const &x) const {
    container_value_type x0{}, y0{}, x1{}, y1{};
    if (x < x_[0]) {
      std::stringstream ss{};
      ss << "Lower constant extrapolation for x ( " << x << " ) occured.\n";
      logger::get().warning(std::cout, ss.str());
      return y_[0];
    }
    if ((x > x_[x_.size() - 1])) {
      std::stringstream ss{};
      ss << "Upper constant extrapolation for x ( " << x << " ) occured.\n";
      logger::get().warning(std::cout, ss.str());
      return y_[x_.size() - 1];
    }

    for (std::size_t i = 0; i < x_.size(); ++i) {
      if (x_[i] < x) {
        x0 = x_[i];
        y0 = y_[i];
      } else if (x_[i] >= x) {
        x1 = x_[i];
        y1 = y_[i];
        break;
      }
    }
    return ((y0 * (x - x1) / (x0 - x1)) + (y1 * (x - x0) / (x1 - x0)));
  }
};

} // namespace lattice

#endif ///_LATTICE_INTERPOLATORS_HPP_