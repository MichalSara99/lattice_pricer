/**
 * @file lattice_boyle.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief  Boyle model
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_BOYLE_HPP_)
#define _LATTICE_BOYLE_HPP_

#include "../common/lattice_risk_factors.hpp"
#include "../utilities/lattice_enums.hpp"
#include "lattice_model_params.hpp"
#include "lattice_trinomial_model.hpp"
#include <functional>

namespace lattice {

/**
 * @brief Boyle model object
 *
 * @tparam T
 */
template <typename T = double> class boyle : public trinomial_model<1, T> {
private:
  model_params<1, asset_class::Equity, T> params_;
  // some constants:
  const T chalf_{0.5};
  const T cone_{1.0};
  const T ctwo_{2.0};

public:
  explicit boyle(model_params<1, asset_class::Equity, T> const &params)
      : params_{params} {}

  // Forward generator
  std::tuple<T, T, T>
  operator()(T value, T dt, std::size_t leaf_idx, std::size_t time_idx,
             bool is_mean_reverting = false) const override {
    T const sig = params_.volatility;
    T const expon = sig * std::sqrt(ctwo_ * dt);
    T const up = std::exp(expon);
    T const mid = cone_;
    T const down = cone_ / up;
    return std::make_tuple(down * value, mid * value, up * value);
  }

  // Backward generator
  T operator()(T curr_value, T down_value, T mid_value, T up_value, T dt,
               std::size_t revert_branches_size, std::size_t nodes_size,
               std::size_t leaf_idx) const override {
    T const q = params_.dividend_rate;
    T const r = params_.risk_free_rate;
    T const sig = params_.volatility;
    T const mu = r - q;
    T const e_sig = std::exp(sig * std::sqrt(chalf_ * dt));
    T const one_over_esig = cone_ / e_sig;
    T const e_mu = std::exp(mu * chalf_ * dt);
    T const p_u =
        std::pow(((e_mu - one_over_esig) / (e_sig - one_over_esig)), ctwo_);
    T const p_d = std::pow(((e_sig - e_mu) / (e_sig - one_over_esig)), ctwo_);
    T const p_m = cone_ - (p_u + p_d);
    T const disc = std::exp(-cone_ * r * dt);
    return (disc * (p_u * up_value + p_m * mid_value + p_d * down_value));
  }

  static std::string const name() { return std::string{"Boyle model"}; }

  static constexpr asset_class class_of_asset() { return asset_class::Equity; }
};

} // namespace lattice

#endif ///_LATTICE_BOYLE_HPP_