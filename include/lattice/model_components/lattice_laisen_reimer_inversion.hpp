#if !defined(_LATTICE_LAISEN_REIMER_INVERSION_HPP_)
#define _LATTICE_LAISEN_REIMER_INVERSION_HPP_

#include "../utilities/lattice_miscellaneous.hpp"

namespace lattice {

/**
 * @brief Peizer-Pratt 1 inversion object
 * for Laisen-Reimer model object
 *
 * @tparam T
 */
template <typename T = double> struct peizer_pratt_1 {
private:
  std::size_t n_;

public:
  explicit peizer_pratt_1(std::size_t number_time_points)
      : n_{number_time_points} {
    // numberTimePoints must be always odd number:
    assert(n_ % 2 != 0);
  }

  T operator()(T x) const {
    auto const first_bra = (x / (n_ + static_cast<T>(1.0 / 3.0)));
    auto const second_bra = (n_ + static_cast<T>(1.0 / 6.0));
    auto const quater = static_cast<T>(0.25);
    auto const sqrt =
        std::sqrt(quater - quater * std::exp(static_cast<T>(-1.0) * first_bra *
                                             first_bra * second_bra));
    return (static_cast<T>(0.5) + sign(x) * sqrt);
  }
};

/**
 * @brief Peizer-Pratt 2 inversion object
 * for Laisen-Reimer model object
 *
 * @tparam T
 */
template <typename T = double> struct peizer_pratt_2 {
private:
  std::size_t n_;

public:
  explicit peizer_pratt_2(std::size_t number_time_points)
      : n_{number_time_points} {
    // numberTimePoints must be always odd number:
    assert(n_ % 2 != 0);
  }

  T operator()(T x) const {
    auto const first_bra = (x / (n_ + static_cast<T>(1.0 / 3.0) +
                                 (static_cast<T>(0.1) / (n_ + 1.0))));
    auto const second_bra = (n_ + static_cast<T>(1.0 / 6.0));
    auto const quater = static_cast<T>(0.25);
    auto const sqrt =
        std::sqrt(quater - quater * std::exp(static_cast<T>(-1.0) * first_bra *
                                             first_bra * second_bra));
    return (static_cast<T>(0.5) + sign(x) * sqrt);
  }
};

} // namespace lattice

#endif ///_LATTICE_LAISEN_REIMER_INVERSION_HPP_