/**
 * @file lattice_model_params.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief model parameters
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_MODEL_PARAMS_HPP_)
#define _LATTICE_MODEL_PARAMS_HPP_

#include "../utilities/lattice_enums.hpp"
#include <typeinfo>

namespace lattice {

template <std::size_t factor_count, asset_class asset, typename T>
struct model_params {};

template <typename T> struct model_params<1, asset_class::InterestRate, T> {
  T reversion_speed;
  T volatility;
};

template <typename T> struct model_params<1, asset_class::Equity, T> {
  T risk_free_rate;
  T volatility;
  T dividend_rate;
  T spot;
  T strike;
  T barrier;
};

template <typename T> struct model_params<2, asset_class::Equity, T> {
  T risk_free_rate;
  T dividend_rate_1;
  T dividend_rate_2;
  T volatility_1;
  T volatility_2;
  T spo_1;
  T spot_2;
  T strike;
  T correlation;
};
} // namespace lattice
#endif ///_LATTICE_MODEL_PARAMS_HPP_