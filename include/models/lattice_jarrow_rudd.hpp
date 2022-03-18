/**
 * @file lattice_jarrow_rudd.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief Jarrow-Rudd model
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_JARROW_RUDD_HPP_)
#define _LATTICE_JARROW_RUDD_HPP_

#include "../utilities/lattice_enums.hpp"
#include "lattice_binomial_model.hpp"
#include "lattice_model_params.hpp"

namespace lattice {

/**
 * @brief Jarrow-Rudd model object
 *
 * @tparam T
 */
template <typename T = double> class jarrow_rudd : public binomial_model<1, T> {
private:
  T prob_;
  model_params<1, asset_class::Equity, T> params_;

public:
  explicit jarrow_rudd(model_params<1, asset_class::Equity, T> const &params)
      : params_{params}, prob_{static_cast<T>(0.5)} {}

  // Forward generator
  std::tuple<T, T> operator()(T value, T dt, std::size_t leaf_idx,
                              std::size_t time_idx,
                              bool is_mean_reverting) const override {
    T const sig = params_.volatility;
    T const d = (params_.risk_free_rate - params_.dividend_rate -
                 static_cast<T>(0.5) * sig * sig) *
                dt;
    T const x1 = d + sig * std::sqrt(dt);
    T const x2 = d - sig * std::sqrt(dt);
    T const up = std::exp(x1);
    T const down = std::exp(x2);
    return std::make_tuple(down * value, up * value);
  }

  // Backward generator
  T operator()(T curr_value, T down_value, T up_value, T dt) const override {
    T const r = params_.risk_free_rate;
    T const disc = std::exp(static_cast<T>(-1.0) * r * dt);
    return (disc *
            (prob_ * up_value + (static_cast<T>(1.0) - prob_) * down_value));
  }

  static std::string const name() { return std::string{"Jarrow-Rudd model"}; }
  static constexpr asset_class class_of_asset() { return asset_class::Equity; }
};

} // namespace lattice

#endif ///_LATTICE_JARROW_RUDD_HPP_