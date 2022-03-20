/**
 * @file lattice_black_derman_toy.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief  BDT model
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_BLACK_DERMAN_TOY_HPP_)
#define _LATTICE_BLACK_DERMAN_TOY_HPP_

#include "../common/lattice_risk_factors.hpp"
#include "../utilities/lattice_enums.hpp"
#include "lattice_binomial_model.hpp"
#include "lattice_model_params.hpp"
#include <functional>

namespace lattice {

/**
 * @brief Black-Derman-Toy model object
 *
 * | d ln(r(t)) = theta(t)*dt + sigma*dw(t) |
 *
 * @tparam T
 */
template <typename T = double>
class black_derman_toy : public binomial_model<1, T> {
private:
  // typedef discounting factor:
  typedef discounting_factor<T> DCF;

  std::function<T(T, T)> dcf_;
  discounting_style ds_;
  T prob_;
  std::vector<T> theta_;
  model_params<1, asset_class::InterestRate, T> params_;

  std::tuple<T, T> _branching(T theta, T sig, T sqrtdt,
                              std::size_t leaf_idx) const {
    T const up = std::exp(sig * (leaf_idx + 1) * sqrtdt);
    T const down = std::exp(sig * leaf_idx * sqrtdt);
    return std::make_tuple(down * theta, up * theta);
  }

public:
  explicit black_derman_toy(
      model_params<1, asset_class::InterestRate, T> const &params,
      std::vector<T> const &theta,
      discounting_style style = discounting_style::Continuous)
      : params_{params}, theta_{theta}, prob_{T{0.5}}, ds_{style},
        dcf_{DCF::function(style)} {}

  discounting_style style_of_discounting() const { return ds_; }

  // Returns risk-neutral probability:
  T node_probability() const { return this->prob_; }

  // Forward generator
  std::tuple<T, T> operator()(T value, T dt, std::size_t leaf_idx,
                              std::size_t time_idx,
                              bool is_mean_reverting) const override {
    T const sig = params_.volatility;
    T const theta = theta_.at(time_idx - 1);
    T const sqrtdt = std::sqrt(dt);
    return _branching(theta, sig, sqrtdt, leaf_idx);
  }

  // Backward generator
  T operator()(T curr_value, T down_value, T up_value, T dt) const override {
    T const disc = dcf_(curr_value, dt);
    return (disc *
            (prob_ * up_value + (static_cast<T>(1.0) - prob_) * down_value));
  }

  static std::string const name() {
    return std::string{"Black-Derman-Toy model"};
  }

  static constexpr asset_class class_of_asset() {
    return asset_class::InterestRate;
  }

  // TODO: this should be moved to general calibrationm engine:

  // Calibration minimizer:
  auto calibration_minimizer() const {
    LASSERT(false, "analytical theta via continuous/discrete discounting is "
                   "not available");
  }

  // Calibration objective function:
  auto calibration_objective() const {
    T const sig = params_.volatility;

    return [=](T theta, T dt, T market_discount,
               std::vector<T> const &prev_rate_states,
               std::vector<T> const &arrow_debreu_states,
               std::function<T(T, T)> const &dcf) -> T {
      T const sqrtdt = std::sqrt(dt);
      T sum{0.0};
      T rate{0.0};
      std::size_t const states_size = arrow_debreu_states.size();
      for (std::size_t i = 0; i < states_size; ++i) {
        rate = theta * std::exp(sig * i * sqrtdt);
        sum += arrow_debreu_states.at(i) * dcf(rate, dt);
      }
      return ((sum - market_discount) * (sum - market_discount));
    };
  }

  // Calibration forward function:
  auto calibrationForwardGenerator() const {
    T const sig = params_.volatility;

    return
        [=](T theta, T value, T dt, std::size_t leaf_idx) -> std::tuple<T, T> {
          T const sqrtdt = std::sqrt(dt);
          return _branching(theta, sig, sqrtdt, leaf_idx);
        };
  }
};
} // namespace lattice

#endif //_LATTICE_BLACK_DERMAN_TOY_HPP_