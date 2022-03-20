/**
 * @file lattice_ho_lee.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief  Ho-Lee model
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(LATTICE_HO_LEE_HPP_)
#define LATTICE_HO_LEE_HPP_

#include "../common/lattice_risk_factors.hpp"
#include "../utilities/lattice_enums.hpp"
#include "lattice_binomial_model.hpp"
#include "lattice_model_params.hpp"
#include <functional>

namespace lattice {

/**
 * @brief Ho-Lee model object
 *
 * |  dr(t) = theta(t)*dt + sigma*dw(t)  |
 *
 * @tparam T
 */
template <typename T = double> class ho_lee : public binomial_model<1, T> {
private:
  // typedef discounting factor:
  typedef discounting_factor<T> DCF;

  std::function<T(T, T)> dcf_;
  discounting_style ds_;
  T prob_;
  std::vector<T> theta_;
  model_params<1, asset_class::InterestRate, T> params_;

  std::tuple<T, T> _branching(T mean, T sigdt) const {
    T const down = mean - sigdt;
    T const up = mean + sigdt;
    return std::make_tuple(down, up);
  }

public:
  explicit ho_lee(model_params<1, asset_class::InterestRate, T> const &params,
                  std::vector<T> const &theta,
                  discounting_style style = discounting_style::Continuous)
      : params_{params}, theta_{theta}, prob_{0.5}, ds_{style},
        dcf_{DCF::function(style)} {}

  discounting_style style_of_discounting() const { return ds_; }

  // Returns risk-neutral probability:
  T node_probability() const { return this->prob_; }

  // Forward generator
  std::tuple<T, T> operator()(T value, T dt, std::size_t leaf_idx,
                              std::size_t time_idx,
                              bool is_mean_reverting) const override {
    T const sig = params_.volatility;
    T const sqrtdt = std::sqrt(dt);
    T const sigdt = sig * sqrtdt;
    T const mean = value + theta_.at(time_idx - 1) * dt;
    return _branching(mean, sigdt);
  }

  // Backward generator
  T operator()(T curr_value, T down_value, T up_value, T dt) const override {
    T const disc = dcf_(curr_value, dt);
    return (disc *
            (prob_ * up_value + (static_cast<T>(1.0) - prob_) * down_value));
  }

  static std::string const name() { return std::string{"Ho-Lee model"}; }

  static constexpr asset_class class_of_asset() {
    return asset_class::InterestRate;
  }

  // TODO: move this to the calibration engine:

  // Calibration minimizer:
  auto calibration_minimizer() const {
    T const sig = params_.volatility;

    // analytical theta via continuous discounting:
    return [=](T theta, T dt, T market_discount,
               std::vector<T> const &prev_rate_states,
               std::vector<T> const &arrow_debreu_states) -> T {
      const std::size_t prev_node_size = prev_rate_states.size();
      std::vector<T> rate_states(prev_node_size + 1);
      T const sqrtdt = std::sqrt(dt);
      T sum{0.0};

      for (std::size_t l = 0; l < prev_node_size; ++l) {
        rate_states[l] = prev_rate_states.at(l) - sig * sqrtdt;     // down l
        rate_states[l + 1] = prev_rate_states.at(l) + sig * sqrtdt; // up l + 1
      }

      std::size_t const states_size = arrow_debreu_states.size();
      for (std::size_t i = 0; i < states_size; ++i) {
        sum +=
            arrow_debreu_states.at(i) * std::exp(-1.0 * rate_states.at(i) * dt);
      }
      sum = (sum / market_discount);
      return std::log(sum) / (dt * dt);
    };
  }

  // Calibration objective function:
  auto calibration_objective() const {

    T const sig = params_.volatility;

    return [=](T theta, T dt, T market_discount,
               std::vector<T> const &prev_rate_states,
               std::vector<T> const &arrow_debreu_states,
               std::function<T(T, T)> const &dcf) -> T {
      const std::size_t prev_node_size = prev_rate_states.size();
      std::vector<T> rate_states(prev_node_size + 1);
      T const sqrtdt = std::sqrt(dt);
      T sum{0.0};

      for (std::size_t l = 0; l < prev_node_size; ++l) {
        rate_states[l] =
            prev_rate_states.at(l) + (theta * dt) - sig * sqrtdt; // down l
        rate_states[l + 1] =
            prev_rate_states.at(l) + (theta * dt) + sig * sqrtdt; // up l + 1
      }

      std::size_t const states_size = arrow_debreu_states.size();
      for (std::size_t i = 0; i < states_size; ++i) {
        sum += arrow_debreu_states.at(i) * dcf(rate_states.at(i), dt);
      }
      return ((sum - market_discount) * (sum - market_discount));
    };
  }

  // Calibration forward function:
  auto calibration_forward_generator() const {
    T const sig = params_.volatility;

    return
        [=](T theta, T value, T dt, std::size_t leaf_idx) -> std::tuple<T, T> {
          T const sqrtdt = std::sqrt(dt);
          T const sigdt = sig * sqrtdt;
          T const mean = value + theta * dt;
          return _branching(mean, sigdt);
        };
  }
};

} // namespace lattice

#endif /// LATTICE_HO_LEE_HPP_