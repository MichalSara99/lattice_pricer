/**
 * @file lattice_leisen_reimer.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief Leisen-Reimer Model
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_LEISEN_REIMER_HPP_)
#define _LATTICE_LEISEN_REIMER_HPP_

#include "../utilities/lattice_enums.hpp"
#include "lattice_binomial_model.hpp"
#include "lattice_model_params.hpp"
#include <functional>

namespace lattice {

/**
 * @brief Leisen-Reimer model object
 *
 * @tparam T
 */
template <typename T = double>
class leisen_reimer : public binomial_model<1, T> {
private:
  model_params<1, asset_class::Equity, T> params_;
  std::size_t number_periods_;
  std::function<T(T)> inversion_;

public:
  explicit leisen_reimer(model_params<1, asset_class::Equity, T> const &params,
                         std::size_t number_periods,
                         std::function<T(T)> const &inversion_formula)
      : params_{params}, number_periods_{number_periods},
        inversion_{inversion_formula} {}

  // Forward generator:
  std::tuple<T, T> operator()(T value, T dt, std::size_t leaf_idx,
                              std::size_t time_idx,
                              bool is_mean_reverting) const override {
    T const sig = params_.volatility;
    T const r = params_.risk_free_rate;
    T const q = params_.dividend_rate;
    T const s = params_.spot;
    T const k = params_.strike;
    T const d1 =
        ((std::log(s / k) + (r - q + static_cast<T>(0.5) * sig * sig)) *
         (number_periods_ * dt)) /
        (sig * std::sqrt(number_periods_ * dt));
    T const d2 = d1 - sig * std::sqrt(number_periods_ * dt);
    T const e = std::exp((r - q) * dt);
    T const h1 = inversion_(d1);
    T const h2 = inversion_(d2);
    T const p = h2;
    T const up = e * (h1 / h2);
    T const down = (e - p * up) / (static_cast<T>(1.0) - p);
    return std::make_tuple(down * value, up * value);
  }

  // Backward generator:
  T operator()(T curr_value, T down_value, T up_value, T dt) const override {
    T const sig = params_.volatility;
    T const r = params_.risk_free_rate;
    T const q = params_.dividend_rate;
    T const s = params_.spot;
    T const k = params_.strike;
    T const d1 =
        ((std::log(s / k) + (r - q + static_cast<T>(0.5) * sig * sig)) *
         (number_periods_ * dt)) /
        (sig * std::sqrt(number_periods_ * dt));
    T const d2 = d1 - sig * std::sqrt(number_periods_ * dt);
    T const h2 = inversion_(d2);
    T const prob = h2;
    T const disc = std::exp(static_cast<T>(-1.0) * r * dt);
    return (disc *
            (prob * up_value + (static_cast<T>(1.0) - prob) * down_value));
  }

  static std::string const name() { return std::string{"Leisen-Reimer model"}; }

  static constexpr asset_class class_of_asset() { return asset_class::Equity; }
};

} // namespace lattice
#endif //_LATTICE_LEISEN_REIMER_HPP_