/**
 * @file lattice_risk_factors.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief Risk factors
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_RISK_FACTORS_HPP_)
#define _LATTICE_RISK_FACTORS_HPP_

#include "../utilities/lattice_enums.hpp"
#include <cmath>
#include <functional>

namespace lattice {

/**
 * @brief Discounting factor struct
 *
 * @tparam T
 */
template <typename T> struct discounting_factor {

  static std::function<T(T, T)> function(discounting_style style) {
    if (style == discounting_style::Continuous) {
      return [=](T rate, T delta) -> T {
        return std::exp(static_cast<T>(-1.0) * rate * delta);
      };
    }
    return [=](T rate, T delta) -> T {
      return (static_cast<T>(1.0) / (static_cast<T>(1.0) + rate * delta));
    };
  }
};

} // namespace lattice

#endif ///_LATTICE_RISK_FACTORS_HPP_