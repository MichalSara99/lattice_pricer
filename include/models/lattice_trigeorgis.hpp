/**
 * @file lattice_trigeorgis.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief Trigeorgis model
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_TRIGEORGIS_HPP_)
#define _LATTICE_TRIGEORGIS_HPP_

#include "../utilities/lattice_enums.hpp"
#include "lattice_binomial_model.hpp"
#include "lattice_model_params.hpp"

namespace lattice {

/**
 * @brief Trigeorgis model object
 *
 * @tparam T
 */
template <typename T = double> class trigeorgis : public binomial_model<1, T> {
private:
  T gamma_;
  model_params<1, asset_class::Equity, T> params_;

public:
  explicit trigeorgis(model_params<1, asset_class::Equity, T> const &params)
      : params_{params}, gamma_{params.risk_free_rate - params.dividend_rate -
                                static_cast<T>(0.5) * params.volatility *
                                    params.volatility} {}

  // Forward generator:
  std::tuple<T, T> operator()(T value, T dt, std::size_t leaf_idx,
                              std::size_t time_idx,
                              bool is_mean_reverting) const override {
    T const sig = params_.volatility;
    T const x = std::sqrt(sig * sig * dt + gamma_ * gamma_ * dt * dt);
    T const up = std::exp(x);
    T const down = (static_cast<T>(1.0) / up);
    return std::make_tuple(down * value, up * value);
  }

  // Backward generator:
  T operator()(T curr_value, T down_value, T up_value, T dt) override {
    T const sig = params_.volatility;
    T const r = params_.risk_free_rate;
    T const x = std::sqrt(sig * sig * dt + gamma_ * gamma_ * dt * dt);
    T const disc = std::exp(static_cast<T>(-1.0) * r * dt);
    T const prob =
        static_cast<T>(0.5) * (static_cast<T>(1.0) + (gamma_ * (dt / x)));
    return (disc * (prob * upValue + (static_cast<T>(1.0) - prob) * downValue));
  }

  static std::string const name() { return std::string{"Trigeorgis model"}; }

  static constexpr asset_class class_of_asset() { return asset_class::Equity; }
};

} // namespace lattice
#endif ///_LATTICE_TRIGEORGIS_HPP_