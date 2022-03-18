/**
 * @file lattice_tian.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief Tian model
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_TIAN_HPP_)
#define _LATTICE_TIAN_HPP_

#include "../utilities/lattice_enums.hpp"
#include "lattice_binomial_model.hpp"
#include "lattice_model_params.hpp"

namespace lattice {

/**
 * @brief Tian model object
 *
 * @tparam T
 */
template <typename T = double> class tian : public binomial_model<1, T> {
private:
  model_params<1, asset_class::Equity, T> params_;

public:
  explicit tian(model_params<1, asset_class::Equity, T> const &params)
      : params_{params} {}

  // Forward generator:
  std::tuple<T, T> operator()(T value, T dt, std::size_t leaf_idx,
                              std::size_t time_idx,
                              bool is_mean_reverting) const override {
    T const q = params_.dividend_rate;
    T const sig = params_.volatility;
    T const r = params_.risk_free_rate;
    T const half = static_cast<T>(0.5);
    T const one = static_cast<T>(1.0);
    T const v = std::exp(sig * sig * dt);
    T const x =
        std::sqrt(v * v + static_cast<T>(2.0) * v - static_cast<T>(3.0));
    T const up = (half * std::exp((r - q) * dt) * v * (v + one + x));
    T const down = (half * std::exp((r - q) * dt) * v * (v + one - x));
    return std::make_tuple(down * value, up * value);
  }

  // Backward generator:
  T operator()(T curr_value, T down_value, T up_value, T dt) const override {
    T const q = params_.dividend_rate;
    T const sig = params_.volatility;
    T const r = params_.risk_free_rate;
    T const half = static_cast<T>(0.5);
    T const one = static_cast<T>(1.0);
    T const v = std::exp(sig * sig * dt);
    T const x =
        std::sqrt(v * v + static_cast<T>(2.0) * v - static_cast<T>(3.0));
    T const up = (half * std::exp((r - q) * dt) * v * (v + one + x));
    T const down = (half * std::exp((r - q) * dt) * v * (v + one - x));
    T const disc = std::exp(-one * r * dt);
    T const prob = ((std::exp((r - q) * dt) - down) / (up - down));
    return (disc * (prob * upValue + (one - prob) * downValue));
  }

  static std::string const name() { return std::string{"Tian model"}; }

  static constexpr asset_class class_of_asset() { return asset_class::Equity; }
};

} // namespace lattice

#endif //_LATTICE_TIAN_HPP_