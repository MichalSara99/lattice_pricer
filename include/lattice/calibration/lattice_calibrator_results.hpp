/**
 * @file lattice_calibrator_results.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief calibrator results
 * @version 0.1
 * @date 2022-03-21
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_CALIBRATOR_RESULTS_HPP_)
#define _LATTICE_CALIBRATOR_RESULTS_HPP_

#include "../utilities/lattice_enums.hpp"
#include "../utilities/lattice_macros.hpp"
#include <tuple>
#include <vector>

namespace lattice {

template <asset_class a_class, typename lattice_object, typename... others>
struct calibrator_results {};

/**
 * @brief Partial specialisation for Interest rate asset class
 *
 * @tparam lattice_object
 */
template <typename lattice_object>
struct calibrator_results<asset_class::InterestRate, lattice_object> {
  typedef typename lattice_object::node_t node;

  std::vector<std::tuple<node, node, node, std::size_t>> theta_optimizers;

  lattice_object arrow_debreu_lattice;

  explicit calibrator_results(
      lattice_object const &arrow_debreu,
      std::vector<std::tuple<node, node, node, std::size_t>> const &theta)
      : theta_optimizers{theta}, arrow_debreu_lattice{arrow_debreu} {}
};

/**
 * @brief Partial specialisation for equity asset class
 *
 * @tparam lattice_object
 * @tparam probability_holder
 */
template <typename lattice_object, typename probability_holder>
struct calibrator_results<asset_class::Equity, lattice_object,
                          probability_holder> {

  lattice_object result_lattice;

  std::vector<std::vector<probability_holder>> implied_probabilities;

  explicit calibrator_results(lattice_object const &results)
      : result_lattice{results} {}

  explicit calibrator_results(
      lattice_object const &results,
      std::vector<std::vector<probability_holder>> const &implied_probs)
      : result_lattice{results}, implied_probabilities(implied_probs) {}
};

/**
 * @brief probability pair type
 *
 * @tparam lattice_object
 */
template <typename lattice_object>
using probability_pair_t =
    std::pair<typename lattice_object::node_t, typename lattice_object::node_t>;

/**
 * @brief probability triplet type
 *
 * @tparam lattice_object
 */
template <typename lattice_object>
using probability_triplet_t =
    std::tuple<typename lattice_object::node_t, typename lattice_object::node_t,
               typename lattice_object::node_t>;

template <typename lattice_object>
using calibrator_ir_results_t =
    calibrator_results<asset_class::InterestRate, lattice_object>;

template <typename lattice_object>
using calibrator_ir_results_ptr =
    std::shared_ptr<calibrator_ir_results_t<lattice_object>>;

template <typename lattice_object>
using calibrator_trinomial_equity_results_t =
    calibrator_results<asset_class::Equity, lattice_object,
                       probability_triplet_t<lattice_object>>;

template <typename lattice_object>
using calibrator_binomial_equity_results_t =
    calibrator_results<asset_class::Equity, lattice_object,
                       probability_pair_t<lattice_object>>;

template <typename lattice_object>
using calibrator_trinomial_equity_results_ptr =
    std::shared_ptr<calibrator_trinomial_equity_results_t<lattice_object>>;

template <typename lattice_object>
using calibrator_binomial_equity_results_ptr =
    std::shared_ptr<calibrator_binomial_equity_results_t<lattice_object>>;

} // namespace lattice

#endif ///_LATTICE_CALIBRATOR_RESULTS_HPP_