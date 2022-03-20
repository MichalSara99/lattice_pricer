/**
 * @file lattice_black_karasinski.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief BK Model
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_BLACK_KARASINSKI_HPP_)
#define _LATTICE_BLACK_KARASINSKI_HPP_

#include "../common/lattice_risk_factors.hpp"
#include "../utilities/lattice_enums.hpp"
#include "lattice_model_params.hpp"
#include "lattice_trinomial_model.hpp"
#include <functional>

namespace lattice {

/**
 * @brief Black-Karasinski model object
 *
 *  |  d ln(r(t)) = (theta(t) - a*ln(r(t)))*dt + sigma * dw(t) |
 *
 * @tparam T
 */
template <typename T = double>
class black_karasinski : public trinomial_model<1, T> {
private:
  // typedef discounting factor:
  typedef discounting_factor<T> DCF;

  std::function<T(T, T)> dcf_;
  discounting_style ds_;
  std::vector<T> theta_;
  model_params<1, asset_class::InterestRate, T> params_;

  std::tuple<T, T, T> _risk_neutral_prob(std::size_t revert_branches_size,
                                         std::size_t nodes_size,
                                         std::size_t leaf_idx, T dt) const {
    T const sig = params_.volatility;
    T const a = params_.reversion_speed;
    T const sqrtdt = std::sqrt(static_cast<T>(3.0) * dt);
    T const dr = sig * sqrtdt;
    T const denom = dr * dr;
    T const one = static_cast<T>(1.0);
    T const half = static_cast<T>(0.5);
    std::size_t const decrement{(nodes_size - 1) / 2};
    long const i = (leaf_idx - decrement);
    T const eps = sign(i) * std::floor(std::abs(i) / revert_branches_size);
    T const nu = (-one * a * dt * i + eps) * dr;
    T const numer = (nu * nu + sig * sig * dt);
    return std::make_tuple(half * ((numer / denom) - (nu / dr)),
                           (one - (numer / denom)),
                           half * ((numer / denom) + (nu / dr)));
  }

  std::tuple<T, T, T> _branching(T mean, T dr, std::size_t leaf_idx,
                                 bool is_mean_reverting = false) const {
    T const one = static_cast<T>(1.0);
    T const two = static_cast<T>(2.0);
    if ((leaf_idx == 0) && (is_mean_reverting == true)) {
      // we are at the bottom of the tree -> going up with branching:
      return std::make_tuple(std::exp(mean) /*low*/,
                             std::exp(mean + one * dr) /*mid*/,
                             std::exp(mean + two * dr) /*high*/);
    }
    if (is_mean_reverting) {
      // we are at the top of the tree -> going down with branching:
      return std::make_tuple(std::exp(mean - two * dr) /*low*/,
                             std::exp(mean - one * dr) /*mid*/,
                             std::exp(mean) /*high*/);
    }
    // Normal branching here:
    return std::make_tuple(std::exp(mean - one * dr) /*low*/,
                           std::exp(mean) /*mid*/,
                           std::exp(mean + one * dr) /*high*/);
  }

public:
  explicit black_karasinski(
      model_params<1, asset_class::InterestRate, T> const &params,
      std::vector<T> const &theta,
      discounting_style style = discounting_style::Continuous)
      : params_{params}, theta_{theta}, ds_{style}, dcf_{DCF::function(style)} {
  }

  discounting_style style_of_discounting() const { return ds_; }

  // Returns tuple of probabilities:
  std::tuple<T, T, T> node_probability(std::size_t revert_branches_size,
                                       std::size_t nodes_size,
                                       std::size_t leaf_idx, T dt) const {
    return _risk_neutral_prob(revert_branches_size, nodes_size, leaf_idx, dt);
  }

  // Forward generator
  std::tuple<T, T, T>
  operator()(T value, T dt, std::size_t leaf_idx, std::size_t time_idx,
             bool is_mean_reverting = false) const override {
    T const sig = params_.volatility;
    T const a = params_.reversion_speed;
    T const sqrtdt = std::sqrt(static_cast<T>(3.0) * dt);
    T const dr = sig * sqrtdt;
    T const mean = std::log(value) + theta_.at(time_idx - 1);
    return _branching(mean, dr, leaf_idx, is_mean_reverting);
  }

  // Backward generator
  T operator()(T curr_value, T down_value, T mid_value, T up_value, T dt,
               std::size_t revert_branches_size, std::size_t nodes_size,
               std::size_t leaf_idx) const override {
    T const disc = dcf_(curr_value, dt);
    std::tuple<T, T, T> const probs =
        node_probability(revert_branches_size, nodes_size, leaf_idx, dt);
    return (disc *
            (std::get<0>(probs) * down_value + std::get<1>(probs) * mid_value +
             std::get<2>(probs) * up_value));
  }

  static std::string const name() {
    return std::string{"Black-Karasinski model"};
  }

  static constexpr asset_class class_of_asset() {
    return asset_class::InterestRate;
  }

  // TODO: Move this to some calibration engine:

  // Calibration objective function:
  std::function<T(T, T, T, std::vector<T> const &, std::vector<T> const &,
                  std::function<T(T, T)> const &)>
  calibration_objective(branching_style branch_style) const {
    T const sig = params_.volatility;
    T const a = params_.reversion_speed;

    if (branch_style == branching_style::Normal) {
      // Calibration normal:
      return [=](T theta, T dt, T market_discount,
                 std::vector<T> const &prev_rate_states,
                 std::vector<T> const &arrow_debreu_states,
                 std::function<T(T, T)> const &dcf) -> T {
        T const sqrtdt = std::sqrt(static_cast<T>(3.0) * dt);
        T const dr = sig * sqrtdt;
        T mean{};
        const std::size_t prev_node_size = prev_rate_states.size();
        std::vector<T> rate_states(prev_node_size + 2);

        for (std::size_t l = 0; l < prev_node_size; ++l) {
          mean = std::log(prev_rate_states.at(l)) + theta;
          rate_states[l] = std::exp(mean - dr);     // down
          rate_states[l + 1] = std::exp(mean);      // mid
          rate_states[l + 2] = std::exp(mean + dr); // high
        }

        T sum{};
        std::size_t const states_size = arrow_debreu_states.size();
        for (std::size_t i = 0; i < states_size; ++i) {
          sum += arrow_debreu_states.at(i) * dcf(rate_states.at(i), dt);
        }
        return ((sum - market_discount) * (sum - market_discount));
      };
    } else {

      // Calibration reverting:
      return [=](T theta, T dt, T market_discount,
                 std::vector<T> const &prev_rate_states,
                 std::vector<T> const &arrow_debreu_states,
                 std::function<T(T, T)> const &dcf) -> T {
        T const sqrtdt = std::sqrt(static_cast<T>(3.0) * dt);
        T const dr = sig * sqrtdt;
        T mean{};
        const std::size_t prev_node_size = prev_rate_states.size();
        std::vector<T> rate_states(prev_node_size);

        // Branching upward here:
        mean = std::log(prev_rate_states.at(0)) + theta;
        rate_states[0] = std::exp(mean);      // down
        rate_states[1] = std::exp(mean + dr); // mid
        rate_states[2] =
            std::exp(mean + static_cast<T>(2.0) * dr); // high
                                                       // normal branching here:
        for (std::size_t l = 1; l < prev_node_size - 1; ++l) {
          mean = std::log(prev_rate_states.at(l)) + theta;
          rate_states[l - 1] = std::exp(mean - dr); // down
          rate_states[l] = std::exp(mean);          // mid
          rate_states[l + 1] = std::exp(mean + dr); // high
        }
        // Branching downward here:
        mean = std::log(prev_rate_states.at(prev_node_size - 1)) + theta;
        rate_states[prev_node_size - 3] =
            std::exp(mean - static_cast<T>(2.0) * dr);         // down
        rate_states[prev_node_size - 2] = std::exp(mean - dr); // mid
        rate_states[prev_node_size - 1] = std::exp(mean);      // high

        T sum{};
        std::size_t const states_size = arrow_debreu_states.size();
        for (std::size_t i = 0; i < states_size; ++i) {
          sum += arrow_debreu_states.at(i) * dcf(rate_states.at(i), dt);
        }
        return ((sum - market_discount) * (sum - market_discount));
      };
    }
  }

  // Calibration forward function:
  auto calibration_forward_generator() const {
    T const sig = params_.volatility;
    T const a = params_.reversion_speed;

    return [=](T theta, T value, T dt, std::size_t leaf_idx,
               bool is_mean_reverting = false) -> std::tuple<T, T, T> {
      T const sqrtdt = std::sqrt(static_cast<T>(3.0) * dt);
      T const dr = sig * sqrtdt;
      T const mean = std::log(value) + theta;
      return _branching(mean, dr, leaf_idx, is_mean_reverting);
    };
  }
};

} // namespace lattice
#endif ///_LATTICE_BLACK_KARASINSKI_HPP_