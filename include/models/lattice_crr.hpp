/**
 * @file lattice_crr.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief  CRR and modified CRR model
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_CRR_HPP_)
#define _LATTICE_CRR_HPP_

#include "../utilities/lattice_enums.hpp"
#include "lattice_binomial_model.hpp"
#include "lattice_model_params.hpp"
#include <tuple>

namespace lattice {

/**
 * @brief Cox-Rubinstein-Ross model object
 *
 * @tparam T
 */
template <typename T = double>
class cox_rubinstein_ross : public binomial_model<1, T> {
private:
  const T chalf_{0.5};
  T prob_;
  model_params<1, asset_class::Equity, T> params_;

public:
  explicit cox_rubinstein_ross(
      model_params<1, asset_class::Equity, T> const &params)
      : params_{params}, prob_{chalf_} {}

  // Forward generator
  std::tuple<T, T> operator()(T value, T dt, std::size_t leaf_idx,
                              std::size_t time_idx,
                              bool is_mean_reverting) const override {
    T const q = params_.dividend_rate;
    T const sig = params_.volatility;
    T const r = params_.risk_free_rate;
    T const r1 = (r - q - chalf_ * sig * sig) * dt;
    T const r2 = sig * std::sqrt(dt);
    T const up = std::exp(r1 + r2);
    T const down = std::exp(r1 - r2);
    return std::make_tuple(down * value, up * value);
  }

  // Backward generator
  T operator()(T curr_value, T down_value, T up_value, T dt) const override {
    T const r = params_.risk_free_rate;
    T const disc = std::exp(static_cast<T>(-1.0) * r * dt);
    return (disc *
            (prob_ * up_value + (static_cast<T>(1.0) - prob_) * down_value));
  }

  static std::string const name() {
    return std::string{"Cox-Rubinstein-Ross model"};
  }

  static constexpr asset_class class_of_asset() { return asset_class::Equity; }
};

/**
 * @brief Modified Cox-Rubinstein-Ross model object
 *
 * @tparam T
 */
template <typename T = double>
class cox_rubinstein_ross_modified : public binomial_model<1, T> {
private:
  model_params<1, asset_class::Equity, T> params_;
  std::size_t n_;

public:
  explicit cox_rubinstein_ross_modified(
      model_params<1, asset_class::Equity, T> const &params,
      std::size_t number_periods)
      : params_{params}, n_{number_periods} {}

  // Forward generator:
  std::tuple<T, T> operator()(T value, T dt, std::size_t leaf_idx,
                              std::size_t time_idx,
                              bool is_mean_reverting) const override {
    T const s = params_.spot;
    T const sig = params_.volatility;
    T const k = params_.strike;
    T const K_n = (std::log(k / s) / static_cast<T>(n_));
    T const V_n = sig * std::sqrt(dt);
    T const up = std::exp(K_n + V_n);
    T const down = std::exp(K_n - V_n);
    return std::make_tuple(down * value, up * value);
  }

  // Backward generator:
  T operator()(T curr_value, T down_value, T up_value, T dt) const override {
    T const q = params_.dividend_rate;
    T const s = params_.spot;
    T const sig = params_.volatility;
    T const r = params_.risk_free_rate;
    T const k = params_.strike;
    T const K_n = (std::log(k / s) / static_cast<T>(n_));
    T const V_n = sig * std::sqrt(dt);
    T const up = std::exp(K_n + V_n);
    T const down = std::exp(K_n - V_n);
    T const disc = std::exp(static_cast<T>(-1.0) * r * dt);
    T const prob = (std::exp((r - q) * dt) - down) / (up - down);
    return (disc *
            (prob * up_value + (static_cast<T>(1.0) - prob) * down_value));
  }

  static std::string const name() {
    return std::string{"Modified Cox-Rubinstein-Ross model"};
  }

  static constexpr asset_class class_of_asset() { return asset_class::Equity; }
};

/**
 * @brief 2-factor Cox-Rubinstein-Ross model object
 *
 * @tparam T
 */
template <typename T = double>
class cox_rubinstein_ross_2_factor : public binomial_model<2, T> {
private:
  T rho_;
  model_params<2, asset_class::Equity, T> params_;

public:
  explicit cox_rubinstein_ross_2_factor(
      model_params<2, asset_class::Equity, T> const &params)
      : params_(params) {}

  // Returns risk-neutral probability:
  std::tuple<T, T, T, T> node_probability(T dt) const {
    T const r = params_.risk_free_rate;
    T const q1 = params_.dividend_rate_1;
    T const q2 = params_.dividend_rate_2;
    T const sig1 = params_.volatility_1;
    T const sig2 = params_.volatility_2;
    T const nu1 = r - q1 - static_cast<T>(0.5) * sig1 * sig1;
    T const nu2 = r - q2 - static_cast<T>(0.5) * sig2 * sig2;
    T const dx1 = sig1 * std::sqrt(dt);
    T const dx2 = sig2 * std::sqrt(dt);
    T const mix = dx1 * dx2;
    T const corr = params_.correlation * sig1 * sig2;

    T const pdd =
        (mix + (static_cast<T>(-1.0) * dx2 * nu1 - dx1 * nu2 + corr) * dt) /
        (static_cast<T>(4.0) * mix);
    T const pdu =
        (mix + (static_cast<T>(-1.0) * dx2 * nu1 + dx1 * nu2 - corr) * dt) /
        (static_cast<T>(4.0) * mix);
    T const pud = (mix + (dx2 * nu1 - dx1 * nu2 - corr) * dt) /
                  (static_cast<T>(4.0) * mix);
    T const puu = (mix + (dx2 * nu1 + dx1 * nu2 + corr) * dt) /
                  (static_cast<T>(4.0) * mix);
    return std::make_tuple(pdd, pdu, pud, puu);
  }

  // Forward generators:
  std::pair<forward_generator_t<T, T, T>, forward_generator_t<T, T, T>>
  forward_generator() const override {
    model_params<1, asset_class::Equity, T> params1;
    params1.risk_free_rate = params_.risk_free_rate;
    params1.dividend_rate = params_.dividend_rate_1;
    params1.volatility = params_.volatility_1;
    params1.spot = params_.spot_1;
    params1.strike = params_.strike;
    model_params<1, asset_class::Equity, T> params2;
    params2.risk_free_rate = params_.risk_free_rate;
    params2.dividend_rate = params_.dividend_rate_2;
    params2.volatility = params_.volatility_2;
    params2.spot = params_.spot_2;
    params2.strike = params_.strike;
    cox_rubinstein_ross_model<T> factor1{params1};
    cox_rubinstein_ross_model<T> factor2{params2};
    forward_generator_t<T, T, T> first =
        [=](T value, T dt, std::size_t leaf_idx, std::size_t time_idx,
            bool is_mean_reverting) -> std::tuple<T, T> {
      return factor1(value, dt, leaf_idx, time_idx, is_mean_reverting);
    };
    forward_generator_t<T, T, T> second =
        [=](T value, T dt, std::size_t leaf_idx, std::size_t time_idx,
            bool is_mean_reverting) -> std::tuple<T, T> {
      return factor2(value, dt, leaf_idx, time_idx, is_mean_reverting);
    };
    return std::make_pair(first, second);
  }

  // Backward generator:
  T operator()(T curr_value, T down_down_value, T down_up_value,
               T up_down_value, T up_up_value, T dt) {
    // taking risk-free rate for discounting from first factor data:
    T const r = params_.risk_free_rate;
    T const disc = std::exp(static_cast<T>(-1.0) * r * dt);
    std::tuple<T, T, T, T> const prob = node_risk_neutral_prob(dt);
    T const value =
        (std::get<0>(prob) * down_down_value +
         std::get<1>(prob) * down_up_value + std::get<2>(prob) * up_down_value +
         std::get<3>(prob) * up_up_value);
    return (disc * value);
  }

  static std::string const name() {
    return std::string{"Cox-Rubinstein-Ross 2-factor model"};
  }

  static constexpr asset_class class_of_asset() { return asset_class::Equity; }
};

} // namespace lattice
#endif ///_LATTICE_CRR_HPP_