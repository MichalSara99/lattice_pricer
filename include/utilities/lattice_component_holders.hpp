/**
 * @file lattice_component_holders.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief  component holders
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_COMPONENT_HOLDERS_HPP_)
#define _LATTICE_COMPONENT_HOLDERS_HPP_

#include <type_traits>

namespace lattice {

/**
 * @brief Static delta time holder
 *
 * @tparam delta_time_t
 */
template <typename delta_time_t> struct delta_time_holder {
private:
  static auto const _delta_time_impl(std::size_t idx,
                                     delta_time_t const &delta_time,
                                     std::true_type) {
    return delta_time.at(idx);
  }

  static delta_time_t const _delta_time_impl(std::size_t idx,
                                             delta_time_t const &delta_time,
                                             std::false_type) {
    return delta_time;
  }

public:
  static auto const delta_time(std::size_t idx,
                               delta_time_t const &delta_time) {
    return _delta_time_impl(idx, delta_time, std::is_compound<delta_time_t>());
  }
};

/**
 * @brief Static risk-free rate holder
 *
 * @tparam risk_free_rate_t
 */
template <typename risk_free_rate_t> struct risk_free_rate_holder {
private:
  static auto const _rate_impl(std::size_t idx,
                               risk_free_rate_t const &risk_free_rate,
                               std::true_type) {
    return risk_free_rate.at(idx);
  }

  static auto const _rate_impl(std::size_t idx,
                               risk_free_rate_t const &risk_free_rate,
                               std::false_type) {
    return risk_free_rate;
  }

public:
  static auto const rate(std::size_t idx,
                         risk_free_rate_t const &risk_free_rate) {
    return _rate_impl(idx, risk_free_rate,
                      std::is_compound<risk_free_rate_t>());
  }
};

} // namespace lattice

#endif ///_LATTICE_COMPONENT_HOLDERS_HPP_