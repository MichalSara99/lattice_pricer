#if !defined(_LATTICE_BACKWARD_HPP_)
#define _LATTICE_BACKWARD_HPP_

#include "../barriers/lattice_adjusters.hpp"
#include "../calibration/lattice_calibrator_results.hpp"
#include "../utilities/lattice_component_holders.hpp"
#include "../utilities/lattice_enums.hpp"
#include "../utilities/lattice_macros.hpp"

namespace lattice {

template <lattice_type type, typename time_axis, typename delta_time>
struct traverse_backward {};

/**
 * @brief Specialization for binomial lattice
 *
 * @tparam time_axis
 * @tparam delta_time
 */
template <typename time_axis, typename delta_time>
struct traverse_backward<lattice_type::Binomial, time_axis, delta_time> {
private:
  template <typename lattice_object, typename generator, typename payoff>
  static void _back_traverse(lattice_object &lattice, generator &&gen,
                             delta_time const &dt, payoff &&pay);

  template <typename lattice_object, typename generator, typename payoff>
  static void
  _back_traverse_barrier(lattice_object &lattice, generator &&gen,
                         delta_time const &dt, payoff &&pay,
                         barrier_type barrier,
                         typename lattice_object::node_t const &barrier_val,
                         typename lattice_object::node_t const &rebate_val);

  //	This one is for payoffadjusted
  template <typename lattice_object, typename generator, typename payoff,
            typename payoff_adjuster>
  static void _back_traverse(lattice_object &lattice, generator &&gen,
                             delta_time const &dt, payoff &&pay,
                             payoff_adjuster &&payoff_adj);

  template <typename lattice_object, typename generator, typename payoff,
            typename payoff_adjuster>
  static void
  _back_traverse_barrier(lattice_object &lattice, generator &&gen,
                         delta_time const &dt, payoff &&pay,
                         payoff_adjuster &&payoff_adj, barrier_type barrier,
                         typename lattice_object::node_t const &barrier_val,
                         typename lattice_object::node_t const &rebate_val);

public:
  template <typename lattice_object, typename generator, typename payoff>
  static void traverse(lattice_object &lattice, generator &&gen,
                       delta_time const &dt, payoff &&pay) {
    _back_traverse(lattice, std::forward<generator>(gen), dt,
                   std::forward<payoff>(pay));
  }

  template <typename lattice_object, typename generator, typename payoff>
  static void
  traverse_barrier(lattice_object &lattice, generator &&gen,
                   delta_time const &dt, payoff &&pay, barrier_type barrier,
                   typename lattice_object::node_t const &barrier_val,
                   typename lattice_object::node_t const &rebate_val) {
    _back_traverse_barrier(lattice, std::forward<generator>(gen), dt,
                           std::forward<payoff>(pay), barrier, barrier_val,
                           rebate_val);
  }

  template <typename lattice_object, typename generator, typename payoff,
            typename payoff_adjuster>
  static void traverse(lattice_object &lattice, generator &&gen,
                       delta_time const &dt, payoff &&pay,
                       payoff_adjuster &&payoff_adj) {
    _back_traverse(lattice, std::forward<generator>(gen), dt,
                   std::forward<payoff>(pay),
                   std::forward<payoff_adjuster>(payoff_adj));
  }

  template <typename lattice_object, typename generator, typename payoff,
            typename payoff_adjuster>
  static void
  traverse_barrier(lattice_object &lattice, generator &&gen,
                   delta_time const &dt, payoff &&pay,
                   payoff_adjuster &&payoff_adj, barrier_type barrier,
                   typename lattice_object::node_t const &barrier_val,
                   typename lattice_object::node_t const &rebate_val) {
    _back_traverse_barrier(lattice, std::forward<generator>(gen), dt,
                           std::forward<payoff>(pay),
                           std::forward<payoff_adjuster>(payoff_adj), barrier,
                           barrier_val, rebate_val);
  }
};

/**
 * @brief Specialisation for 2 variable binomial lattice
 *
 * @tparam time_axis
 * @tparam delta_time
 */
template <typename time_axis, typename delta_time>
struct traverse_backward<lattice_type::TwoVariableBinomial, time_axis,
                         delta_time> {
private:
  template <typename lattice_object, typename multidim_lattice_object,
            typename generator, typename payoff>
  static void _back_traverse(lattice_object &price_lattice,
                             multidim_lattice_object const &lattice,
                             generator &&gen, delta_time const &dt,
                             payoff &&pay);

  //	This one is for payoffadjusted
  template <typename lattice_object, typename multidim_lattice_object,
            typename generator, typename payoff, typename payoff_adjuster>
  static void _back_traverse(lattice_object &priceLattice,
                             multidim_lattice_object const &lattice,
                             generator &&gen, delta_time const &dt,
                             payoff &&payoff, payoff_adjuster &&payoff_adj);

public:
  template <typename lattice_object, typename multidim_lattice_object,
            typename generator, typename payoff>
  static void traverse(lattice_object &price_lattice,
                       multidim_lattice_object const &lattice, generator &&gen,
                       delta_time const &dt, payoff &&pay) {
    _back_traverse(price_lattice, lattice, std::forward<generator>(gen), dt,
                   std::forward<payoff>(pay));
  }

  template <typename lattice_object, typename multidim_lattice_object,
            typename generator, typename payoff, typename payoff_adjuster>
  static void traverse(lattice_object &price_lattice,
                       multidim_lattice_object const &lattice, generator &&gen,
                       delta_time const &dt, payoff &&pay,
                       payoff_adjuster &&payoff_adj) {
    _back_traverse(price_lattice, lattice, std::forward<generator>(gen), dt,
                   std::forward<payoff>(pay),
                   std::forward<payoff_adjuster>(payoff_adj));
  }
};

/**
 * @brief Specialisation for trinomial lattice
 *
 * @tparam TimeAxis
 * @tparam DeltaTime
 */
template <typename time_axis, typename delta_time>
struct traverse_backward<lattice_type::Trinomial, time_axis, delta_time> {
private:
  template <typename lattice_object, typename generator>
  static void _back_traverse_normal(std::size_t time_idx,
                                    lattice_object &lattice, generator &&gen,
                                    delta_time const &dt);

  template <typename lattice_object, typename generator>
  static void _back_traverse_reverting(std::size_t time_idx,
                                       lattice_object &lattice, generator &&gen,
                                       delta_time const &dt);

  template <typename lattice_object, typename generator, typename payoff>
  static void _back_traverse(lattice_object &lattice, generator &&gen,
                             delta_time const &dt, payoff &&pay);

  template <typename lattice_object, typename generator>
  static void _back_traverse_normal_barrier(
      std::size_t time_idx, lattice_object &lattice, generator &&gen,
      delta_time const &dt, barrier_type barrier,
      typename lattice_object::node_t const &barrier_val,
      typename lattice_object::node_t const &rebate_val);

  template <typename lattice_object, typename generator>
  static void _back_traverse_reverting_barrier(
      std::size_t time_idx, lattice_object &lattice, generator &&gen,
      delta_time const &dt, barrier_type barrier,
      typename lattice_object::node_t const &barrier_val,
      typename lattice_object::node_t const &rebate_val);

  template <typename lattice_object, typename generator, typename payoff>
  static void
  _back_traverse_barrier(lattice_object &lattice, generator &&gen,
                         delta_time const &dt, payoff &&pay,
                         barrier_type barrier,
                         typename lattice_object::node_t const &barrier_val,
                         typename lattice_object::node_t const &rebate_val);

  template <typename lattice_object, typename generator,
            typename payoff_adjuster>
  static void _back_traverse_normal(std::size_t time_idx,
                                    lattice_object &lattice, generator &&gen,
                                    delta_time const &dt,
                                    payoff_adjuster &&payoff_adj);

  template <typename lattice_object, typename generator,
            typename payoff_adjuster>
  static void _back_traverse_reverting(std::size_t time_idx,
                                       lattice_object &lattice, generator &&gen,
                                       delta_time const &dt,
                                       payoff_adjuster &&payoff_adj);

  template <typename lattice_object, typename generator, typename payoff,
            typename payoff_adjuster>
  static void _back_traverse(lattice_object &lattice, generator &&gen,
                             delta_time const &dt, payoff &&pay,
                             payoff_adjuster &&payoff_adj);

  template <typename lattice_object, typename generator,
            typename payoff_adjuster>
  static void _back_traverse_normal_barrier(
      std::size_t time_idx, lattice_object &lattice, generator &&gen,
      delta_time const &dt, payoff_adjuster &&payoff_adj, barrier_type barrier,
      typename lattice_object::node_t const &barrier_val,
      typename lattice_object::node_t const &rebate_val);

  template <typename lattice_object, typename generator,
            typename payoff_adjuster>
  static void _back_traverse_reverting_barrier(
      std::size_t time_idx, lattice_object &lattice, generator &&gen,
      delta_time const &dt, payoff_adjuster &&payoff_adj, barrier_type barrier,
      typename lattice_object::node_t const &barrier_val,
      typename lattice_object::node_t const &rebate_val);

  template <typename lattice_object, typename generator, typename payoff,
            typename payoff_adjuster>
  static void
  _back_traverse_barrier(lattice_object &lattice, generator &&gen,
                         delta_time const &dt, payoff &&pay,
                         payoff_adjuster &&payoff_adj, barrier_type barrier,
                         typename lattice_object::node_t const &barrier_val,
                         typename lattice_object::node_t const &rebate_val);

public:
  template <typename lattice_object, typename generator, typename payoff>
  static void traverse(lattice_object &lattice, generator &&gen,
                       delta_time const &dt, payoff &&pay) {
    _back_traverse(lattice, std::forward<Generator>(gen), dt,
                   std::forward<payoff>(pay));
  }

  template <typename lattice_object, typename generator, typename payoff,
            typename payoff_adjuster>
  static void traverse(lattice_object &lattice, generator &&gen,
                       delta_time const &dt, payoff &&pay,
                       payoff_adjuster &&payoff_adj) {
    _back_traverse(lattice, std::forward<generator>(gen), dt,
                   std::forward<payoff>(pay),
                   std::forward<payoff_adjuster>(payoff_adj));
  }

  template <typename lattice_object, typename generator, typename payoff>
  static void
  traverse_barrier(lattice_object &lattice, generator &&gen,
                   delta_time const &dt, payoff &&pay, barrier_type barrier,
                   typename lattice_object::node_t const &barrier_val,
                   typename lattice_object::node_t const &rebate_val) {
    _back_traverse_barrier(lattice, std::forward<generator>(gen), dt,
                           std::forward<payoff>(pay), barrier, barrier_val,
                           rebate_val);
  }

  template <typename lattice_object, typename generator, typename payoff,
            typename payoff_adjuster>
  static void
  traverse_barrier(lattice_object &lattice, generator &&gen,
                   delta_time const &dt, payoff &&pay,
                   payoff_adjuster &&payoff_adj, barrier_type barrier,
                   typename lattice_object::node_t const &barrier_val,
                   typename lattice_object::node_t const &rebate_val) {
    _back_traverse_barrier(lattice, std::forward<generator>(gen), dt,
                           std::forward<payoff>(pay),
                           std::forward<payoff_adjuster>(payoff_adj), barrier,
                           barrier_val, rebate_val);
  }
};

// ======================================================================================
// ==================== Implied Backward Traversal Algorithms
// ===========================
// ======================================================================================

template <lattice_type type, typename delta_time, typename risk_free_rate>
struct traverse_implied_backward {};

template <typename delta_time, typename risk_free_rate>
struct traverse_implied_backward<lattice_type::Trinomial, delta_time,
                                 risk_free_rate> {
private:
  template <typename lattice_object, typename payoff>
  static void _back_traverse(
      lattice_object &option_lattice, lattice_object const &spot_lattice,
      calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
      delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
      discounting_style style);

  template <typename lattice_object, typename payoff, typename payoff_adjuster>
  static void _back_traverse(
      lattice_object &option_lattice, lattice_object const &spot_lattice,
      calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
      delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
      payoff_adjuster &&payoff_adj, discounting_style style);

  template <typename lattice_object, typename payoff>
  static void _back_traverse_barrier(
      lattice_object &option_lattice, lattice_object const &spot_lattice,
      calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
      delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
      barrier_type barrier, typename lattice_object::node_t const &barrier_val,
      typename lattice_object::node_t const &rebate_val,
      discounting_style style);

  template <typename lattice_object, typename payoff, typename payoff_adjuster>
  static void _back_traverse_barrier(
      lattice_object &option_lattice, lattice_object const &spot_lattice,
      calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
      delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
      payoff_adjuster &&payoff_adj, barrier_type barrier,
      typename lattice_object::node_t const &barrier_val,
      typename lattice_object::node_t const &rebate_val,
      discounting_style style);

  template <typename lattice_object, typename payoff>
  static void _back_traverse_barrier_dke_adjustment(
      lattice_object &option_lattice, lattice_object const &spot_lattice,
      calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
      delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
      barrier_type barrier, typename lattice_object::node_t const &barrier_val,
      typename lattice_object::node_t const &rebate_val,
      discounting_style style);

  template <typename lattice_object, typename payoff, typename payoff_adjuster>
  static void _back_traverse_barrier_dke_adjustment(
      lattice_object &option_lattice, lattice_object const &spot_lattice,
      calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
      delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
      payoff_adjuster &&payoff_adj, barrier_type barrier,
      typename lattice_object::node_t const &barrier_val,
      typename lattice_object::node_t const &rebate_val,
      discounting_style style);

public:
  template <typename lattice_object, typename payoff>
  static void traverse(
      lattice_object &option_lattice, lattice_object const &spot_lattice,
      calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
      delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
      discounting_style style) {
    _back_traverse(option_lattice, spot_lattice, results, dt, rfr,
                   std::forward<payoff>(pay), style);
  }

  template <typename lattice_object, typename payoff, typename payoff_adjuster>
  static void traverse(
      lattice_object &option_lattice, lattice_object const &spot_lattice,
      calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
      delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
      payoff_adjuster &&payoff_adj, discounting_style style) {
    _back_traverse(option_lattice, spot_lattice, results, dt, rfr,
                   std::forward<payoff>(pay),
                   std::forward<payoff_adjuster>(payoff_adj), style);
  }

  template <typename lattice_object, typename payoff>
  static void traverse_barrier(
      lattice_object &option_lattice, lattice_object const &spot_lattice,
      calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
      delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
      barrier_type barrier, typename lattice_object::node_t const &barrier_val,
      typename lattice_object::node_t const &rebate_val,
      bool derman_kani_ergener_adujstment, discounting_style style) {
    if (derman_kani_ergener_adujstment)
      _back_traverse_barrier_dke_adjustment(
          option_lattice, spot_lattice, results, dt, rfr,
          std::forward<payoff>(pay), barrier, barrier_val, rebate_val, style);
    else
      _back_traverse_barrier(option_lattice, spot_lattice, results, dt, rfr,
                             std::forward<payoff>(pay), barrier, barrier_val,
                             rebate_val, style);
  }

  template <typename lattice_object, typename payoff, typename payoff_adjuster>
  static void traverse_barrier(
      lattice_object &option_lattice, lattice_object const &spot_lattice,
      calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
      delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
      payoff_adjuster &&payoff_adj, barrier_type barrier,
      typename lattice_object::node_t const &barrier_val,
      typename lattice_object::node_t const &rebate_val,
      bool derman_kani_ergener_adujstment, discounting_style style) {
    if (derman_kani_ergener_adujstment)
      _back_traverse_barrier_dke_adjustment(
          option_lattice, spot_lattice, results, dt, rfr,
          std::forward<payoff>(pay), std::forward<payoff_adjuster>(payoff_adj),
          barrier, barrier_val, rebate_val, style);
    else
      _back_traverse_barrier(option_lattice, spot_lattice, results, dt, rfr,
                             std::forward<payoff>(pay),
                             std::forward<payoff_adjuster>(payoff_adj), barrier,
                             barrier_val, rebate_val, style);
  }
};

/// =============== IMPLEMENTATION =================

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator, typename payoff>
void traverse_backward<lattice_type::Binomial, time_axis,
                       delta_time>::_back_traverse(lattice_object &lattice,
                                                   generator &&gen,
                                                   delta_time const &dt,
                                                   payoff &&pay) {

  typedef delta_time_holder<delta_time> DT;
  const std::size_t last_idx = lattice.time_dimension() - 1;
  const std::size_t last_nodes_size = lattice.nodes_at_idx(last_idx).size();
  typename lattice_object::node_t delta{};

  for (auto i = 0; i < last_nodes_size; ++i) {
    lattice(last_idx, i) = pay(lattice(last_idx, i));
  }
  std::size_t nodes_size{0};
  for (auto n = last_idx - 1; n > 0; --n) {
    nodes_size = lattice.nodes_at_idx(n).size();
    delta = DT::delta_time(n, dt);
    for (auto i = 0; i < nodes_size; ++i) {
      lattice(n, i) =
          gen(lattice(n, i), lattice(n + 1, i), lattice(n + 1, i + 1), delta);
    }
  }
  delta = DT::delta_time(0, dt);
  lattice(0, 0) = gen(lattice(0, 0), lattice(1, 0), lattice(1, 1), dt);
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator, typename payoff>
void traverse_backward<lattice_type::Binomial, time_axis, delta_time>::
    _back_traverse_barrier(lattice_object &lattice, generator &&gen,
                           delta_time const &dt, payoff &&pay,
                           barrier_type barrier,
                           typename lattice_object::node_t const &barrier_val,
                           typename lattice_object::node_t const &rebate_val) {

  typedef delta_time_holder<delta_time> DT;
  typedef typename lattice_object::node_t node;
  typedef barrier_comparer<node> BC;
  auto cmp = BC::comparer(barrier);
  const std::size_t last_idx = lattice.time_dimension() - 1;
  const std::size_t last_nodes_size = lattice.nodes_at_idx(last_idx).size();
  node delta{};

  for (auto i = 0; i < last_nodes_size; ++i) {
    if (cmp(lattice(last_idx, i), barrier_val)) {
      lattice(last_idx, i) = pay(lattice(last_idx, i));
    } else {
      lattice(last_idx, i) = node{};
    }
  }
  std::size_t nodes_size{0};
  for (auto n = last_idx - 1; n > 0; --n) {
    nodes_size = lattice.nodes_at_idx(n).size();
    delta = DT::delta_time(n, dt);
    for (auto i = 0; i < nodes_size; ++i) {
      if (cmp(lattice(n, i), barrier_val)) {
        lattice(n, i) =
            gen(lattice(n, i), lattice(n + 1, i), lattice(n + 1, i + 1), delta);
      } else {
        lattice(n, i) = rebate;
      }
    }
  }
  delta = DT::delta_time(0, dt);
  if (cmp(lattice(0, 0), barrier_val)) {
    lattice(0, 0) = gen(lattice(0, 0), lattice(1, 0), lattice(1, 1), delta);
  } else {
    lattice(0, 0) = rebate_val;
  }
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator, typename payoff,
          typename payoff_adjuster>
void traverse_backward<lattice_type::Binomial, time_axis, delta_time>::
    _back_traverse(lattice_object &lattice, generator &&gen,
                   delta_time const &dt, payoff &&pay,
                   payoff_adjuster &&payoff_adj) {

  typedef delta_time_holder<delta_time> DT;
  const std::size_t last_idx = lattice.time_dimension() - 1;
  const std::size_t last_nodes_size = lattice.nodes_at_idx(last_idx).size();
  typename lattice_object::node_t delta{};
  typename lattice_object::node_t value{};

  for (auto i = 0; i < last_nodes_size; ++i) {
    lattice(last_idx, i) = pay(lattice(last_idx, i));
  }

  std::size_t nodes_size{0};
  for (auto n = last_idx - 1; n > 0; --n) {
    nodes_size = lattice.nodes_at_idx(n).size();
    delta = DT::delta_time(n, dt);
    for (auto i = 0; i < nodes_size; ++i) {
      value =
          gen(lattice(n, i), lattice(n + 1, i), lattice(n + 1, i + 1), delta);
      payoff_adj(value, lattice(n, i));
      lattice(n, i) = value;
    }
  }
  delta = DT::delta_time(0, dt);
  value = gen(lattice(0, 0), lattice(1, 0), lattice(1, 1), delta);
  payoff_adj(value, lattice(0, 0));
  lattice(0, 0) = value;
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator, typename payoff,
          typename payoff_adjuster>
void traverse_backward<lattice_type::Binomial, time_axis, delta_time>::
    _back_traverse_barrier(lattice_object &lattice, generator &&gen,
                           delta_time const &dt, payoff &&pay,
                           payoff_adjuster &&payoff_adj, barrier_type barrier,
                           typename lattice_object::node_t const &barrier_val,
                           typename lattice_object::node_t const &rebate_val) {

  typedef delta_time_holder<delta_time> DT;
  typedef typename lattice_object::node_t node;
  typedef barrier_comparer<node> BC;
  auto cmp = BC::comparer(barrier);
  const std::size_t last_idx = lattice.time_dimension() - 1;
  const std::size_t last_nodes_size = lattice.nodes_at_idx(last_idx).size();
  node delta{};
  node value{};

  for (auto i = 0; i < last_nodes_size; ++i) {
    if (cmp(lattice(last_idx, i), barrier_val)) {
      lattice(last_idx, i) = pay(lattice(last_idx, i));
    } else {
      lattice(last_idx, i) = node{};
    }
  }

  std::size_t nodes_size{0};
  for (auto n = last_idx - 1; n > 0; --n) {
    nodes_size = lattice.nodesAtIdx(n).size();
    delta = DT::delta_time(n, dt);
    for (auto i = 0; i < nodes_size; ++i) {
      if (cmp(lattice(n, i), barrier_val)) {
        value =
            gen(lattice(n, i), lattice(n + 1, i), lattice(n + 1, i + 1), delta);
        payoff_adj(value, lattice(n, i));
        lattice(n, i) = value;
      } else {
        lattice(n, i) = rebate_val;
      }
    }
  }
  delta = DT::delta_time(0, dt);
  if (cmp(lattice(0, 0), barrier_val)) {
    value = gen(lattice(0, 0), lattice(1, 0), lattice(1, 1), delta);
    payoff_adj(value, lattice(0, 0));
    lattice(0, 0) = value;
  } else {
    lattice(0, 0) = rebate_val;
  }
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename multidim_lattice_object,
          typename generator, typename payoff>
void traverse_backward<
    lattice_type::TwoVariableBinomial, time_axis,
    delta_time>::_back_traverse(lattice_object &price_lattice,
                                multidim_lattice_object const &lattice,
                                generator &&gen, delta_time const &dt,
                                payoff &&pay) {

  typedef delta_time_holder<delta_time> DT;
  typename lattice_object::node_t delta{};

  auto tree1 = lattice.get_factor(0);
  auto tree2 = lattice.get_factor(1);

  std::size_t const tree_size = price_lattice.time_dimension();
  std::size_t const last_idx = tree_size - 1;
  std::size_t const last_nodes_size =
      price_lattice.nodes_at_idx(last_idx).size();
  std::size_t factor_last_nodes_size = tree1.nodes_at_idx(last_idx).size();
  std::size_t col{0}, row{0};

  for (auto i = 0; i < last_nodes_size; ++i) {
    col = i % factor_last_nodes_size;
    if ((i > 0) && ((i % factor_last_nodes_size) == 0))
      row++;
    price_lattice(last_idx, i) =
        pay(tree1(last_idx, row), tree2(last_idx, col));
  }

  std::size_t nodes_size{0};
  std::size_t factor_nodes_size{0};
  for (auto n = last_idx - 1; n > 0; --n) {
    nodes_size = price_lattice.nodes_at_idx(n).size();
    factor_nodes_size = tree1.nodes_at_idx(n).size();
    delta = DT::delta_time(n, dt);
    row = 0;
    for (auto i = 0; i < nodes_size; ++i) {
      if ((i > 0) && ((i % factor_nodes_size) == 0))
        row++;
      price_lattice(n, i) = gen(
          price_lattice(n, i), /*down-down*/ price_lattice(n + 1, i + row),
          /*down-up*/ price_lattice(n + 1, i + row + 1),
          /*up-down*/ price_lattice(n + 1, i + row + factor_last_nodes_size),
          /*up-up*/ price_lattice(n + 1, i + row + factor_last_nodes_size + 1),
          delta);
    }
    factor_last_nodes_size = factor_nodes_size;
  }
  delta = DT::delta_time(0, dt);
  price_lattice(0, 0) =
      gen(price_lattice(0, 0), price_lattice(1, 0), price_lattice(1, 1),
          price_lattice(1, 2), price_lattice(1, 3), delta);
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename multidim_lattice_object,
          typename generator, typename payoff, typename payoff_adjuster>
void traverse_backward<
    lattice_type::TwoVariableBinomial, time_axis,
    delta_time>::_back_traverse(lattice_object &price_lattice,
                                multidim_lattice_object const &lattice,
                                generator &&gen, delta_time const &dt,
                                payoff &&pay, payoff_adjuster &&payoff_adj) {

  typedef delta_time_holder<delta_time> DT;
  typename lattice_object::node_t delta{};
  typename lattice_object::node_t value{};

  auto tree1 = lattice.get_factor(0);
  auto tree2 = lattice.get_factor(1);

  std::size_t const tree_size = price_lattice.time_dimension();
  std::size_t const last_idx = tree_size - 1;
  std::size_t const last_nodes_size =
      price_lattice.nodes_at_idx(last_idx).size();
  std::size_t factor_last_nodes_size = tree1.nodes_at_idx(last_idx).size();
  std::size_t col{0}, row{0};

  for (auto i = 0; i < last_nodes_size; ++i) {
    col = i % factor_last_nodes_size;
    if ((i > 0) && ((i % factor_last_nodes_size) == 0))
      row++;
    price_lattice(last_idx, i) =
        pay(tree1(last_idx, row), tree2(last_idx, col));
  }

  std::size_t nodes_size{0};
  std::size_t factor_nodes_size{0};
  for (auto n = last_idx - 1; n > 0; --n) {
    nodes_size = price_lattice.nodesAtIdx(n).size();
    factor_nodes_size = tree1.nodesAtIdx(n).size();
    delta = DT::delta_time(n, dt);
    col = 0;
    row = 0;
    for (auto i = 0; i < nodes_size; ++i) {
      col = i % factor_nodes_size;
      if ((i > 0) && ((i % factor_nodes_size) == 0))
        row++;
      value = gen(
          price_lattice(n, i), /*down-down*/ price_lattice(n + 1, i + row),
          /*down-up*/ price_lattice(n + 1, i + row + 1),
          /*up-down*/ price_lattice(n + 1, i + row + factor_last_nodes_size),
          /*up-up*/ price_lattice(n + 1, i + row + factor_last_nodes_size + 1),
          delta);
      payoffAdjuster(value, tree1(n, row), tree2(n, col));
      price_lattice(n, i) = value;
    }
    factor_last_nodes_size = factor_nodes_size;
  }
  delta = DT::delta_time(0, dt);
  value = gen(price_lattice(0, 0), price_lattice(1, 0), price_lattice(1, 1),
              price_lattice(1, 2), price_lattice(1, 3), delta);
  payoff_adj(value, tree1(0, 0), tree2(0, 0));
  price_lattice(0, 0) = value;
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator>
void traverse_backward<lattice_type::Trinomial, time_axis, delta_time>::
    _back_traverse_normal(std::size_t time_idx, lattice_object &lattice,
                          generator &&gen, delta_time const &dt) {

  typedef delta_time_holder<delta_time> DT;
  typename lattice_object::node_t delta{};

  std::size_t const revert_branches_size = time_idx;

  std::size_t nodes_size{0};
  for (auto n = time_idx - 1; n > 0; --n) {
    delta = DT::delta_time(n, dt);
    nodes_size = lattice.nodes_at_idx(n).size();
    for (auto i = 0; i < nodes_size; ++i) {
      lattice(n, i) = gen(lattice(n, i), lattice(n + 1, i),
                          lattice(n + 1, i + 1), lattice(n + 1, i + 2), delta,
                          revert_branches_size, nodes_size, i);
    }
  }
  delta = DT::delta_time(0, dt);
  lattice(0, 0) = gen(lattice(0, 0), lattice(1, 0), lattice(1, 1),
                      lattice(1, 2), delta, revert_branches_size, 1, 0);
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator>
void traverse_backward<lattice_type::Trinomial, time_axis, delta_time>::
    _back_traverse_reverting(std::size_t time_idx, lattice_object &lattice,
                             generator &&gen, delta_time const &dt) {

  typedef delta_time_holder<delta_time> DT;
  typename lattice_object::node_t delta{};
  const std::size_t last_idx = lattice.time_dimension() - 1;

  std::size_t const revert_branches_size = time_idx;

  std::size_t nodes_size{0};
  for (auto n = last_idx - 1; n >= time_idx; --n) {
    delta = DT::delta_time(n, dt);
    nodes_size = lattice.nodes_at_idx(n).size();
    lattice(n, 0) =
        gen(lattice(n, 0), lattice(n + 1, 0), lattice(n + 1, 1),
            lattice(n + 1, 2), delta, revert_branches_size, nodes_size, 0);
    for (auto i = 1; i < nodes_size - 1; ++i) {
      lattice(n, i) = gen(lattice(n, i), lattice(n + 1, i - 1),
                          lattice(n + 1, i), lattice(n + 1, i + 1), delta,
                          revert_branches_size, nodes_size, i);
    }
    lattice(n, nodes_size - 1) =
        gen(lattice(n, nodes_size - 1), lattice(n + 1, nodes_size - 3),
            lattice(n + 1, nodes_size - 2), lattice(n + 1, nodes_size - 1),
            delta, revert_branches_size, nodes_size, nodes_size - 1);
  }
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator, typename payoff>
void traverse_backward<lattice_type::Trinomial, time_axis,
                       delta_time>::_back_traverse(lattice_object &lattice,
                                                   generator &&gen,
                                                   delta_time const &dt,
                                                   payoff &&pay) {

  const std::size_t first_revert_idx = lattice.first_reverting_idx();
  const std::size_t tree_size = lattice.time_dimension();

  const std::size_t last_idx = tree_size - 1;
  const std::size_t last_nodes_size = lattice.nodes_at_idx(last_idx).size();

  for (auto i = 0; i < last_nodes_size; ++i) {
    lattice(last_idx, i) = pay(lattice(last_idx, i));
  }

  if (first_revert_idx == 0) {
    // This trinomial tree does not have reverting property:
    _back_traverse_normal(last_idx, lattice, std::forward<generator>(gen), dt);
  } else {
    // This trinomial tree does have reverting property:
    _back_traverse_reverting(first_revert_idx - 1, lattice,
                             std::forward<generator>(gen), dt);
    _back_traverse_normal(first_revert_idx - 1, lattice,
                          std::forward<generator>(gen), dt);
  }
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator>
void traverse_backward<lattice_type::Trinomial, time_axis, delta_time>::
    _back_traverse_normal_barrier(
        std::size_t time_idx, lattice_object &lattice, generator &&gen,
        delta_time const &dt, barrier_type barrier,
        typename lattice_object::node_t const &barrier_val,
        typename lattice_object::node_t const &rebate_val) {

  typedef delta_time_holder<delta_time> DT;
  typedef typename lattice_object::node_t node;
  typedef barrier_comparer<node> BC;
  std::size_t const revert_branches_size = time_idx;
  auto cmp = BC::comparer(barrier);

  node delta{};
  std::size_t nodes_size{0};
  for (auto n = time_idx - 1; n > 0; --n) {
    delta = DT::delta_time(n, dt);
    nodes_size = lattice.nodes_at_idx(n).size();
    for (auto i = 0; i < nodes_size; ++i) {
      if (cmp(lattice(n, i), barrier_val)) {
        lattice(n, i) = gen(lattice(n, i), lattice(n + 1, i),
                            lattice(n + 1, i + 1), lattice(n + 1, i + 2), delta,
                            revert_branches_size, nodes_size, i);
      } else {
        lattice(n, i) = rebate_val;
      }
    }
  }
  delta = DT::delta_time(0, dt);
  if (cmp(lattice(0, 0), barrier_val)) {
    lattice(0, 0) = gen(lattice(0, 0), lattice(1, 0), lattice(1, 1),
                        lattice(1, 2), delta, revert_branches_size, 1, 0);
  } else {
    lattice(0, 0) = rebate_val;
  }
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator>
void traverse_backward<lattice_type::Trinomial, time_axis, delta_time>::
    _back_traverse_reverting_barrier(
        std::size_t time_idx, lattice_object &lattice, generator &&gen,
        delta_time const &dt, barrier_type barrier,
        typename lattice_object::node_t const &barrier_val,
        typename lattice_object::node_t const &rebate_val) {

  typedef delta_time_holder<delta_time> DT;
  typedef typename lattice_object::node_t node;
  typedef barrier_comparer<node> BC;
  auto cmp = BC::comparer(barrier);
  const std::size_t last_idx = lattice.time_dimension() - 1;
  std::size_t const revert_branches_size = time_idx;

  node delta{};
  std::size_t nodes_size{0};
  for (auto n = last_idx - 1; n >= time_idx; --n) {
    delta = DT::delta_time(n, dt);
    nodes_size = lattice.nodes_at_idx(n).size();
    if (cmp(lattice(n, 0), barrier_val)) {
      lattice(n, 0) =
          gen(lattice(n, 0), lattice(n + 1, 0), lattice(n + 1, 1),
              lattice(n + 1, 2), delta, revert_branches_size, nodes_size, 0);
    } else {
      lattice(n, 0) = rebate_val;
    }
    for (auto i = 1; i < nodes_size - 1; ++i) {
      if (cmp(lattice(n, i), barrier_val)) {
        lattice(n, i) = gen(lattice(n, i), lattice(n + 1, i - 1),
                            lattice(n + 1, i), lattice(n + 1, i + 1), delta,
                            revert_branches_size, nodes_size, i);
      } else {
        lattice(n, i) = rebate_val;
      }
    }
    if (cmp(lattice(n, nodes_size - 1), barrier_val)) {
      lattice(n, nodes_size - 1) =
          gen(lattice(n, nodes_size - 1), lattice(n + 1, nodes_size - 3),
              lattice(n + 1, nodes_size - 2), lattice(n + 1, nodes_size - 1),
              delta, revert_branches_size, nodes_size, nodes_size - 1);
    } else {
      lattice(n, nodes_size - 1) = rebate_val;
    }
  }
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator, typename payoff>
void traverse_backward<lattice_type::Trinomial, time_axis, delta_time>::
    _back_traverse_barrier(lattice_object &lattice, generator &&gen,
                           delta_time const &dt, payoff &&pay,
                           barrier_type barrier,
                           typename lattice_object::node_t const &barrier_val,
                           typename lattice_object::node_t const &rebate_val) {

  typedef typename lattice_object::node_t node;
  typedef barrier_comparer<node> BC;
  auto cmp = BC::comparer(barrier);
  const std::size_t first_revert_idx = lattice.first_reverting_idx();
  const std::size_t tree_size = lattice.time_dimension();

  const std::size_t last_idx = tree_size - 1;
  const std::size_t last_nodes_size = lattice.nodes_at_idx(last_idx).size();

  for (auto i = 0; i < last_nodes_size; ++i) {
    if (cmp(lattice(last_idx, i), barrier_val)) {
      lattice(last_idx, i) = pay(lattice(last_idx, i));
    } else {
      lattice(last_idx, i) = node{};
    }
  }

  if (first_revert_idx == 0) {
    // This trinomial tree does not have reverting property:
    _back_traverse_normal_barrier(last_idx, lattice,
                                  std::forward<generator>(gen), dt, barrier,
                                  barrier_val, rebate_val);
  } else {
    // This trinomial tree does have reverting property:
    _back_traverse_reverting_barrier(first_revert_idx - 1, lattice,
                                     std::forward<generator>(gen), dt, barrier,
                                     barrier_val, rebate_val);
    _back_traverse_normal_barrier(first_revert_idx - 1, lattice,
                                  std::forward<generator>(gen), dt, barrier,
                                  barrier_val, rebate_val);
  }
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator, typename payoff_adjuster>
void traverse_backward<lattice_type::Trinomial, time_axis, delta_time>::
    _back_traverse_normal(std::size_t time_idx, lattice_object &lattice,
                          generator &&gen, delta_time const &dt,
                          payoff_adjuster &&payoff_adj) {

  typedef delta_time_holder<delta_time> DT;
  std::size_t const revert_branches_size = time_idx;
  typename lattice_object::node_t delta{};
  typename lattice_object::node_t value{};

  std::size_t nodes_size{0};
  for (auto n = time_idx - 1; n > 0; --n) {
    delta = DT::delta_time(n, dt);
    nodes_size = lattice.nodes_at_idx(n).size();
    for (auto i = 0; i < nodes_size; ++i) {
      value = gen(lattice(n, i), lattice(n + 1, i), lattice(n + 1, i + 1),
                  lattice(n + 1, i + 2), delta, revert_branches_size,
                  nodes_size, i);
      payoff_adj(value, lattice(n, i));
      lattice(n, i) = value;
    }
  }
  delta = DT::delta_time(0, dt);
  value = gen(lattice(0, 0), lattice(1, 0), lattice(1, 1), lattice(1, 2), delta,
              revert_branches_size, 1, 0);
  payoff_adj(value, lattice(0, 0));
  lattice(0, 0) = value;
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator, typename payoff_adjuster>
void traverse_backward<lattice_type::Trinomial, time_axis, delta_time>::
    _back_traverse_reverting(std::size_t time_idx, lattice_object &lattice,
                             generator &&gen, delta_time const &dt,
                             payoff_adjuster &&payoff_adj) {

  typedef delta_time_holder<delta_time> DT;
  typename lattice_object::node_t delta{};
  typename lattice_object::node_t value{};
  const std::size_t last_idx = lattice.time_dimension() - 1;
  std::size_t const revert_branches_size = time_idx;

  std::size_t nodes_size{0};
  for (auto n = last_idx - 1; n >= time_idx; --n) {
    delta = DT::delta_time(n, dt);
    nodes_size = lattice.nodes_at_idx(n).size();

    value = gen(lattice(n, 0), lattice(n + 1, 0), lattice(n + 1, 1),
                lattice(n + 1, 2), delta, revert_branches_size, nodes_size, 0);
    payoff_adj(value, lattice(n, 0));
    lattice(n, 0) = value;

    for (auto i = 1; i < nodes_size - 1; ++i) {
      value = gen(lattice(n, i), lattice(n + 1, i - 1), lattice(n + 1, i),
                  lattice(n + 1, i + 1), delta, revert_branches_size,
                  nodes_size, i);
      payoff_adj(value, lattice(n, i));
      lattice(n, i) = value;
    }

    value = gen(lattice(n, nodes_size - 1), lattice(n + 1, nodes_size - 3),
                lattice(n + 1, nodes_size - 2), lattice(n + 1, nodes_size - 1),
                delta, revert_branches_size, nodes_size, nodes_size - 1);
    payoff_adj(value, lattice(n, nodes_size - 1));
    lattice(n, nodes_size - 1) = value;
  }
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator, typename payoff,
          typename payoff_adjuster>
void traverse_backward<lattice_type::Trinomial, time_axis, delta_time>::
    _back_traverse(lattice_object &lattice, generator &&gen,
                   delta_time const &dt, payoff &&pay,
                   payoff_adjuster &&payoff_adj) {

  const std::size_t first_revert_idx = lattice.first_reverting_idx();
  const std::size_t tree_size = lattice.time_dimension();

  const std::size_t last_idx = tree_size - 1;
  const std::size_t last_nodes_size = lattice.nodes_at_idx(last_idx).size();

  for (auto i = 0; i < last_nodes_size; ++i) {
    lattice(last_idx, i) = pay(lattice(last_idx, i));
  }

  if (first_revert_idx == 0) {
    // This trinomial tree does not have reverting property:
    _back_traverse_normal(last_idx, lattice, std::forward<generator>(gen), dt,
                          std::forward<payoff_adjuster>(payoff_adj));
  } else {
    // This trinomial tree does have reverting property:
    _back_traverse_reverting(first_revert_idx - 1, lattice,
                             std::forward<generator>(gen), dt,
                             std::forward<payoff_adjuster>(payoff_adj));
    _back_traverse_normal(first_revert_idx - 1, lattice,
                          std::forward<generator>(gen), dt,
                          std::forward<payoff_adjuster>(payoff_adj));
  }
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator, typename payoff_adjuster>
void traverse_backward<lattice_type::Trinomial, time_axis, delta_time>::
    _back_traverse_normal_barrier(
        std::size_t time_idx, lattice_object &lattice, generator &&gen,
        delta_time const &dt, payoff_adjuster &&payoff_adj,
        barrier_type barrier,
        typename lattice_object::node_t const &barrier_val,
        typename lattice_object::node_t const &rebate_val) {

  typedef delta_time_holder<delta_time> DT;
  typedef typename lattice_object::node_t node;
  typedef barrier_comparer<node> BC;
  auto cmp = BC::comparer(barrier);
  std::size_t const revert_branches_size = time_idx;

  node delta{};
  node value{};
  std::size_t nodes_size{0};
  for (auto n = time_idx - 1; n > 0; --n) {
    delta = DT::delta_time(n, dt);
    nodes_size = lattice.nodes_at_idx(n).size();
    for (auto i = 0; i < nodes_size; ++i) {
      if (cmp(lattice(n, i), barrier_val)) {
        value = gen(lattice(n, i), lattice(n + 1, i), lattice(n + 1, i + 1),
                    lattice(n + 1, i + 2), delta, revert_branches_size,
                    nodes_size, i);
        payoff_adj(value, lattice(n, i));
        lattice(n, i) = value;
      } else {
        lattice(n, i) = rebate_val;
      }
    }
  }
  delta = DT::delta_time(0, dt);
  if (cmp(lattice(0, 0), barrier_val)) {
    value = gen(lattice(0, 0), lattice(1, 0), lattice(1, 1), lattice(1, 2),
                delta, revert_branches_size, 1, 0);
    payoff_adj(value, lattice(0, 0));
    lattice(0, 0) = value;
  } else {
    lattice(0, 0) = rebate_val;
  }
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator, typename payoff_adjuster>
void traverse_backward<lattice_type::Trinomial, time_axis, delta_time>::
    _back_traverse_reverting_barrier(
        std::size_t time_idx, lattice_object &lattice, generator &&gen,
        delta_time const &dt, payoff_adjuster &&payoff_adj,
        barrier_type barrier,
        typename lattice_object::node_t const &barrier_val,
        typename lattice_object::node_t const &rebate_val) {

  typedef delta_time_holder<delta_time> DT;
  typedef typename lattice_object::node_t node;
  typedef barrier_comparer<node> BC;
  auto cmp = BC::comparer(barrier);
  const std::size_t last_idx = lattice.time_dimension() - 1;
  std::size_t const revert_branches_size = time_idx;

  node delta{};
  node value{};
  std::size_t nodes_size{0};
  for (auto n = last_idx - 1; n >= time_idx; --n) {
    delta = DT::delta_time(n, dt);
    nodes_size = lattice.nodes_at_idx(n).size();

    if (cmp(lattice(n, 0), barrier_val)) {
      value =
          gen(lattice(n, 0), lattice(n + 1, 0), lattice(n + 1, 1),
              lattice(n + 1, 2), delta, revert_branches_size, nodes_size, 0);
      payoff_adj(value, lattice(n, 0));
      lattice(n, 0) = value;
    } else {
      lattice(n, 0) = rebate_val;
    }

    for (auto i = 1; i < nodes_size - 1; ++i) {
      if (cmp(lattice(n, i), barrier_val)) {
        value = gen(lattice(n, i), lattice(n + 1, i - 1), lattice(n + 1, i),
                    lattice(n + 1, i + 1), delta, revert_branches_size,
                    nodes_size, i);
        payoff_adj(value, lattice(n, i));
        lattice(n, i) = value;
      } else {
        lattice(n, i) = rebate_val;
      }
    }
    if (cmp(lattice(n, nodes_size - 1), barrier_val)) {
      value =
          gen(lattice(n, nodes_size - 1), lattice(n + 1, nodes_size - 3),
              lattice(n + 1, nodes_size - 2), lattice(n + 1, nodes_size - 1),
              delta, revert_branches_size, nodes_size, nodes_size - 1);
      payoff_adj(value, lattice(n, nodes_size - 1));
      lattice(n, nodes_size - 1) = value;
    } else {
      lattice(n, nodes_size - 1) = rebate_val;
    }
  }
}

template <typename time_axis, typename delta_time>
template <typename lattice_object, typename generator, typename payoff,
          typename payoff_adjuster>
void traverse_backward<lattice_type::Trinomial, time_axis, delta_time>::
    _back_traverse_barrier(lattice_object &lattice, generator &&gen,
                           delta_time const &dt, payoff &&pay,
                           payoff_adjuster &&payoff_adj, barrier_type barrier,
                           typename lattice_object::node_t const &barrier_val,
                           typename lattice_object::node_t const &rebate_val) {

  typedef typename lattice_object::node_t node;
  typedef barrier_comparer<node> BC;

  auto cmp = BC::comparer(barrier);
  const std::size_t first_revert_idx = lattice.first_reverting_idx();
  const std::size_t tree_size = lattice.time_dimension();

  const std::size_t last_idx = tree_size - 1;
  const std::size_t last_nodes_size = lattice.nodes_at_idx(last_idx).size();

  for (auto i = 0; i < last_nodes_size; ++i) {
    if (cmp(lattice(last_idx, i), barrier_val)) {
      lattice(last_idx, i) = pay(lattice(last_idx, i));
    } else {
      lattice(last_idx, i) = node{};
    }
  }

  if (first_revert_idx == 0) {
    // This trinomial tree does not have reverting property:
    _back_traverse_normal_barrier(last_idx, lattice,
                                  std::forward<generator>(gen), dt,
                                  std::forward<payoff_adjuster>(payoff_adj),
                                  barrier, barrier_val, rebate_val);
  } else {
    // This trinomial tree does have reverting property:
    _back_traverse_reverting_barrier(first_revert_idx - 1, lattice,
                                     std::forward<generator>(gen), dt,
                                     std::forward<payoff_adjuster>(payoff_adj),
                                     barrier, barrier_val, rebate_val);
    _back_traverse_normal_barrier(first_revert_idx - 1, lattice,
                                  std::forward<generator>(gen), dt,
                                  std::forward<payoff_adjuster>(payoff_adj),
                                  barrier, barrier_val, rebate_val);
  }
}

template <typename delta_time, typename risk_free_rate>
template <typename lattice_object, typename payoff>
void traverse_implied_backward<lattice_type::Trinomial, delta_time,
                               risk_free_rate>::
    _back_traverse(
        lattice_object &option_lattice, lattice_object const &spot_lattice,
        calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
        delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
        discounting_style style) {

  // typedefs:
  typedef typename lattice_object::node_t node;
  // typedef RiskFreeRateHolder:
  typedef risk_free_rate_holder<risk_free_rate> RFR;
  // typedef DeltaTimeHolder:
  typedef delta_time_holder<delta_time> DT;
  // typedef discounting factor:
  typedef discounting_factor<node> DCF;

  // get correct discounting factor style:
  auto dcf = DCF::function(style);

  // get the size of optionLattice object:
  std::size_t const tree_size = option_lattice.time_dimension();
  std::size_t const last_idx = tree_size - 1;
  std::size_t const last_nodes_size =
      option_lattice.nodes_at_idx(last_idx).size();

  // unpack the impled probabilities from calibrationResults:
  auto implied_probs = results->implied_probabilities;
  LASSERT(implied_probs.size() >= (tree_size - 1),
          "Must have enough implied volatilities for pricing.");

  // last payoff first:
  for (std::size_t i = 0; i < last_nodes_size; ++i) {
    option_lattice(last_idx, i) = pay(spot_lattice(last_idx, i));
  }

  std::size_t nodes_size{0};
  node df{};
  std::vector<std::tuple<node, node, node>> probs;
  std::tuple<node, node, node> tpl;
  for (auto t = last_idx - 1; t > 0; --t) {
    df = dcf(RFR::rate(t, rfr), DT::delta_time(t, dt));
    nodes_size = option_lattice.nodes_at_idx(t).size();
    probs = implied_probs.at(t);
    for (auto i = 0; i < nodes_size; ++i) {
      tpl = probs.at(i);
      option_lattice(t, i) =
          df * (std::get<0>(tpl) * option_lattice(t + 1, i) +
                std::get<1>(tpl) * option_lattice(t + 1, i + 1) +
                std::get<2>(tpl) * option_lattice(t + 1, i + 2));
    }
  }
  df = dcf(RFR::rate(0, rfr), DT::delta_time(0, dt));
  probs = implied_probs.at(0);
  tpl = probs.at(0);
  option_lattice(0, 0) = df * (std::get<0>(tpl) * option_lattice(1, 0) +
                               std::get<1>(tpl) * option_lattice(1, 1) +
                               std::get<2>(tpl) * option_lattice(1, 2));
}

template <typename delta_time, typename risk_free_rate>
template <typename lattice_object, typename payoff, typename payoff_adjuster>
void traverse_implied_backward<lattice_type::Trinomial, delta_time,
                               risk_free_rate>::
    _back_traverse(
        lattice_object &option_lattice, lattice_object const &spot_lattice,
        calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
        delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
        payoff_adjuster &&payoff_adj, discounting_style style) {

  // typedefs:
  typedef typename lattice_object::node_t node;
  // typedef RiskFreeRateHolder:
  typedef risk_free_rate_holder<risk_free_rate> RFR;
  // typedef DeltaTimeHolder:
  typedef delta_time_holder<delta_time> DT;
  // typedef discounting factor:
  typedef discounting_factor<node> DCF;

  // get correct discounting factor style:
  auto dcf = DCF::function(style);

  // get the size of optionLattice object:
  std::size_t const tree_size = option_lattice.time_dimension();
  std::size_t const last_idx = treeSize - 1;
  std::size_t const last_nodes_size =
      option_lattice.nodes_at_idx(last_idx).size();

  // unpack the impled probabilities from calibrationResults:
  auto implied_probs = results->implied_probabilities;
  LASSERT(implied_probs.size() >= (tree_size - 1),
          "Must have enough implied volatilities for pricing.");

  // last payoff first:
  for (std::size_t i = 0; i < last_nodes_size; ++i) {
    option_lattice(last_idx, i) = pay(spot_lattice(last_idx, i));
  }

  std::size_t nodes_size{0};
  node df{};
  node value{};
  std::vector<std::tuple<node, node, node>> probs;
  std::tuple<node, node, node> tpl;
  for (auto t = last_idx - 1; t > 0; --t) {
    df = dcf(RFR::rate(t, rfr), DT::delta_time(t, dt));
    nodes_size = option_lattice.nodes_at_idx(t).size();
    probs = implied_probs.at(t);
    for (auto i = 0; i < nodes_size; ++i) {
      tpl = probs.at(i);
      value = df * (std::get<0>(tpl) * option_lattice(t + 1, i) +
                    std::get<1>(tpl) * option_lattice(t + 1, i + 1) +
                    std::get<2>(tpl) * option_lattice(t + 1, i + 2));
      payoff_adj(value, spot_lattice(t, i));
      option_lattice(t, i) = value;
    }
  }
  df = dcf(RFR::rate(0, rfr), DT::delta_time(0, dt));
  probs = implied_probs.at(0);
  tpl = probs.at(0);
  value = df * (std::get<0>(tpl) * option_lattice(1, 0) +
                std::get<1>(tpl) * option_lattice(1, 1) +
                std::get<2>(tpl) * option_lattice(1, 2));
  payoff_adj(value, spot_lattice(0, 0));
  option_lattice(0, 0) = value;
}

template <typename delta_time, typename risk_free_rate>
template <typename lattice_object, typename payoff>
void traverse_implied_backward<lattice_type::Trinomial, delta_time,
                               risk_free_rate>::
    _back_traverse_barrier(
        lattice_object &option_lattice, lattice_object const &spot_lattice,
        calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
        delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
        barrier_type barrier,
        typename lattice_object::node_t const &barrier_val,
        typename lattice_object::node_t const &rebate_val,
        discounting_style style) {

  // typedefs:
  typedef typename lattice_object::node_t node;
  // typedef RiskFreeRateHolder:
  typedef risk_free_rate_holder<risk_free_rate> RFR;
  // typedef DeltaTimeHolder:
  typedef delta_time_holder<delta_time> DT;
  // typedef discounting factor:
  typedef discounting_factor<node> DCF;
  // typedef BarrierComparer:
  typedef barrier_comparer<node> BC;

  // get correct comparer function:
  auto cmp = BC::comparer(barrier);
  // get correct discounting factor style:
  auto dcf = DCF::function(style);

  // get the size of optionLattice object:
  std::size_t const tree_size = option_lattice.time_dimension();
  std::size_t const last_idx = treeSize - 1;
  std::size_t const last_nodes_size =
      option_lattice.nodes_at_idx(last_idx).size();

  // unpack the impled probabilities from calibrationResults:
  auto implied_probs = results->implied_probabilities;
  LASSERT(implied_probs.size() >= (tree_size - 1),
          "Must have enough implied volatilities for pricing.");

  // last payoff first:
  for (std::size_t i = 0; i < last_nodes_size; ++i) {
    if (cmp(spot_lattice(last_idx, i), barrier_val)) {
      option_lattice(last_idx, i) = payoff(spot_lattice(last_idx, i));
    } else {
      option_lattice(last_idx, i) = 0.0;
    }
  }

  std::size_t nodes_size{0};
  node df{};
  std::vector<std::tuple<node, node, node>> probs;
  std::tuple<node, node, node> tpl;
  for (auto t = last_idx - 1; t > 0; --t) {
    df = dcf(RFR::rate(t, rfr), DT::delta_time(t, dt));
    nodes_size = option_lattice.nodes_at_idx(t).size();
    probs = implied_probs.at(t);
    for (auto i = 0; i < nodes_size; ++i) {
      if (cmp(spot_lattice(t, i), barrier_val)) {
        tpl = probs.at(i);
        option_lattice(t, i) =
            df * (std::get<0>(tpl) * option_lattice(t + 1, i) +
                  std::get<1>(tpl) * option_lattice(t + 1, i + 1) +
                  std::get<2>(tpl) * option_lattice(t + 1, i + 2));
      } else {
        option_lattice(t, i) = rebate_val;
      }
    }
  }

  df = dcf(RFR::rate(0, rfr), DT::delta_time(0, dt));
  probs = implied_probs.at(0);
  if (cmp(spot_lattice(0, 0), barrier_val)) {
    tpl = probs.at(0);
    option_lattice(0, 0) = df * (std::get<0>(tpl) * option_lattice(1, 0) +
                                 std::get<1>(tpl) * option_lattice(1, 1) +
                                 std::get<2>(tpl) * option_lattice(1, 2));
  } else {
    option_lattice(0, 0) = rebate_val;
  }
}

template <typename delta_time, typename risk_free_rate>
template <typename lattice_object, typename payoff, typename payoff_adjuster>
void traverse_implied_backward<lattice_type::Trinomial, delta_time,
                               risk_free_rate>::
    _back_traverse_barrier(
        lattice_object &option_lattice, lattice_object const &spot_lattice,
        calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
        delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
        payoff_adjuster &&payoff_adj, barrier_type barrier,
        typename lattice_object::node_t const &barrier_val,
        typename lattice_object::node_t const &rebate_val,
        discounting_style style) {

  // typedefs:
  typedef typename lattice_object::node_t node;
  // typedef RiskFreeRateHolder:
  typedef risk_free_rate_holder<risk_free_rate> RFR;
  // typedef DeltaTimeHolder:
  typedef delta_time_holder<delta_time> DT;
  // typedef discounting factor:
  typedef discounting_factor<node> DCF;
  // typedef BarrierComparer:
  typedef barrier_comparer<node> BC;

  // get correct comparer function:
  auto cmp = BC::comparer(barrier);
  // get correct discounting factor style:
  auto dcf = DCF::function(style);

  // get the size of optionLattice object:
  std::size_t const tree_size = option_lattice.time_dimension();
  std::size_t const last_idx = tree_size - 1;
  std::size_t const last_nodes_size =
      option_lattice.nodes_at_idx(last_idx).size();

  // unpack the impled probabilities from calibrationResults:
  auto implied_probs = results->implied_probabilities;
  LASSERT(implied_probs.size() >= (treeSize - 1),
          "Must have enough implied volatilities for pricing.");

  // last payoff first:
  for (std::size_t i = 0; i < last_nodes_size; ++i) {
    if (cmp(spot_lattice(last_idx, i), barrier_val)) {
      option_lattice(last_idx, i) = pay(spot_lattice(last_idx, i));
    } else {
      option_lattice(last_idx, i) = 0.0;
    }
  }

  std::size_t nodes_size{0};
  node df{};
  node value{};
  std::vector<std::tuple<node, node, node>> probs;
  std::tuple<node, node, node> tpl;
  for (auto t = last_idx - 1; t > 0; --t) {
    df = dcf(RFR::rate(t, rfr), DT::delta_time(t, dt));
    nodes_size = option_lattice.nodes_at_idx(t).size();
    probs = implied_probs.at(t);
    for (auto i = 0; i < nodes_size; ++i) {
      if (cmp(spot_lattice(t, i), barrier_val)) {
        tpl = probs.at(i);
        value = df * (std::get<0>(tpl) * option_lattice(t + 1, i) +
                      std::get<1>(tpl) * option_lattice(t + 1, i + 1) +
                      std::get<2>(tpl) * option_lattice(t + 1, i + 2));
        payoff_adj(value, spotLattice(t, i));
        option_lattice(t, i) = value;

      } else {
        option_lattice(t, i) = rebate_val;
      }
    }
  }

  df = dcf(RFR::rate(0, rfr), DT::delta_time(0, dt));
  probs = implied_probs.at(0);
  if (cmp(spot_lattice(0, 0), barrier_val)) {
    tpl = probs.at(0);
    value = df * (std::get<0>(tpl) * option_lattice(1, 0) +
                  std::get<1>(tpl) * option_lattice(1, 1) +
                  std::get<2>(tpl) * option_lattice(1, 2));
    payoff_adj(value, spot_lattice(0, 0));
    option_lattice(0, 0) = value;
  } else {
    option_lattice(0, 0) = rebate_val;
  }
}

template <typename delta_time, typename risk_free_rate>
template <typename lattice_object, typename payoff>
void traverse_implied_backward<lattice_type::Trinomial, delta_time,
                               risk_free_rate>::
    _back_traverse_barrier_dke_adjustment(
        lattice_object &option_lattice, lattice_object const &spot_lattice,
        calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
        delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
        barrier_type barrier,
        typename lattice_object::node_t const &barrier_val,
        typename lattice_object::node_t const &rebate_val,
        discounting_style style) {

  // typedefs:
  typedef typename lattice_object::node_t node;
  // typedef RiskFreeRateHolder:
  typedef risk_free_rate_holder<risk_free_rate> RFR;
  // typedef DeltaTimeHolder:
  typedef delta_time_holder<delta_time> DT;
  // typedef discounting factor:
  typedef discounting_factor<node> DCF;
  // typedef BarrierComparer:
  typedef barrier_comparer<node> BC;
  // typedef DermanKaniErgenerAdjuster:
  typedef derman_kani_ergener_adjuster<node> DKEA;

  // get correct comparer function:
  auto cmp = BC::comparer(barrier);
  // get correct adjuster:
  auto adjust_pair = DKEA::adjuster(barrier);
  // unpack checker and adjuster:
  auto dke_checker = adjust_pair.first;
  auto dke_adjuster = adjust_pair.second;
  // get correct discounting factor style:
  auto dcf = DCF::function(style);

  // get the size of optionLattice object:
  std::size_t const tree_size = option_lattice.time_dimension();
  std::size_t const last_idx = tree_size - 1;
  std::size_t const last_nodes_size =
      option_lattice.nodes_at_idx(last_idx).size();

  // unpack the impled probabilities from calibrationResults:
  auto implied_probs = results->implied_probabilities;
  LASSERT(implied_probs.size() >= (tree_size - 1),
          "Must have enough implied volatilities for pricing.");

  // last payoff first:
  for (std::size_t i = 0; i < last_nodes_size; ++i) {
    if (cmp(spot_lattice(last_idx, i), barrier_val)) {
      option_lattice(last_idx, i) = payoff(spot_lattice(last_idx, i));
    } else {
      option_lattice(last_idx, i) = 0.0;
    }
  }

  std::size_t nodes_size{0};
  node df{};
  std::vector<std::tuple<node, node, node>> probs;
  std::tuple<node, node, node> tpl;
  for (auto t = last_idx - 1; t > 0; --t) {
    df = dcf(RFR::rate(t, rfr), DT::delta_time(t, dt));
    nodes_size = option_lattice.nodes_at_idx(t).size();
    probs = implied_probs.at(t);
    for (auto i = 0; i < nodes_size; ++i) {
      if (cmp(spot_lattice(t, i), barrier_val)) {
        tpl = probs.at(i);
        option_lattice(t, i) =
            df * (std::get<0>(tpl) * option_lattice(t + 1, i) +
                  std::get<1>(tpl) * option_lattice(t + 1, i + 1) +
                  std::get<2>(tpl) * option_lattice(t + 1, i + 2));
      } else {
        option_lattice(t, i) = rebate_val;
      }
      // Derman-Kani-Ergener adjustment:
      if (((i > 0) && (i < nodes_size - 1)) &&
          dke_checker(spot_lattice(t, i), spot_lattice(t, i - 1),
                      spot_lattice(t, i + 1), barrier_val)) {
        option_lattice(t, i) = dke_adjuster(
            spot_lattice(t, i), spot_lattice(t, i - 1), spot_lattice(t, i + 1),
            barrier_val, rebate_val, option_lattice(t, i));
      }
    }
  }

  df = dcf(RFR::rate(0, rfr), DT::delta_time(0, dt));
  probs = implied_probs.at(0);
  if (cmp(spot_lattice(0, 0), barrier_val)) {
    tpl = probs.at(0);
    option_lattice(0, 0) = df * (std::get<0>(tpl) * option_lattice(1, 0) +
                                 std::get<1>(tpl) * option_lattice(1, 1) +
                                 std::get<2>(tpl) * option_lattice(1, 2));
  } else {
    option_lattice(0, 0) = rebate_val;
  }
}

template <typename delta_time, typename risk_free_rate>
template <typename lattice_object, typename payoff, typename payoff_adjuster>
void traverse_implied_backward<lattice_type::Trinomial, delta_time,
                               risk_free_rate>::
    _back_traverse_barrier_dke_adjustment(
        lattice_object &option_lattice, lattice_object const &spot_lattice,
        calibrator_trinomial_equity_results_ptr<lattice_object> const &results,
        delta_time const &dt, risk_free_rate const &rfr, payoff &&pay,
        payoff_adjuster &&payoff_adj, barrier_type barrier,
        typename lattice_object::node_t const &barrier_val,
        typename lattice_object::node_t const &rebate_val,
        discounting_style style) {

  // typedefs:
  typedef typename lattice_object::node_t node;
  // typedef RiskFreeRateHolder:
  typedef risk_free_rate_holder<risk_free_rate> RFR;
  // typedef DeltaTimeHolder:
  typedef delta_time_holder<delta_time> DT;
  // typedef discounting factor:
  typedef discounting_factor<node> DCF;
  // typedef BarrierComparer:
  typedef barrier_comparer<node> BC;
  // typedef DermanKaniErgenerAdjuster:
  typedef derman_kani_ergener_adjuster<node> DKEA;

  // get correct comparer function:
  auto cmp = BC::comparer(barrier);
  // get correct adjuster:
  auto adjust_pair = DKEA::adjuster(barrier);
  // unpack checker and adjuster:
  auto dke_checker = adjust_pair.first;
  auto dke_adjuster = adjust_pair.second;
  // get correct discounting factor style:
  auto dcf = DCF::function(style);

  // get the size of optionLattice object:
  std::size_t const tree_size = option_lattice.time_dimension();
  std::size_t const last_idx = tree_size - 1;
  std::size_t const last_nodes_size =
      option_lattice.nodes_at_idx(last_idx).size();

  // unpack the impled probabilities from calibrationResults:
  auto implied_probs = results->implied_probabilities;
  LASSERT(implied_probs.size() >= (tree_size - 1),
          "Must have enough implied volatilities for pricing.");

  // last payoff first:
  for (std::size_t i = 0; i < last_nodes_size; ++i) {
    if (cmp(spot_lattice(last_idx, i), barrier_val)) {
      option_lattice(last_idx, i) = payoff(spot_lattice(last_idx, i));
    } else {
      option_lattice(last_idx, i) = 0.0;
    }
  }

  std::size_t nodes_size{0};
  node df{};
  node value{};
  std::vector<std::tuple<node, node, node>> probs;
  std::tuple<node, node, node> tpl;
  for (auto t = last_idx - 1; t > 0; --t) {
    df = dcf(RFR::rate(t, rfr), DT::delta_time(t, dt));
    nodes_size = option_lattice.nodes_at_idx(t).size();
    probs = implied_probs.at(t);
    for (auto i = 0; i < nodes_size; ++i) {
      if (cmp(spot_lattice(t, i), barrier_val)) {
        tpl = probs.at(i);
        value = df * (std::get<0>(tpl) * option_lattice(t + 1, i) +
                      std::get<1>(tpl) * option_lattice(t + 1, i + 1) +
                      std::get<2>(tpl) * option_lattice(t + 1, i + 2));
        payoff_adj(value, spot_lattice(t, i));
        option_lattice(t, i) = value;

      } else {
        option_lattice(t, i) = rebate_val;
      }
      // Derman-Kani-Ergener adjustment:
      if (((i > 0) && (i < nodes_size - 1)) &&
          dke_checker(spot_lattice(t, i), spot_lattice(t, i - 1),
                      spot_lattice(t, i + 1), barrier_val)) {
        option_lattice(t, i) = dke_adjuster(
            spot_lattice(t, i), spot_lattice(t, i - 1), spot_lattice(t, i + 1),
            barrier_val, rebate_val, option_lattice(t, i));
      }
    }
  }

  df = dcf(RFR::rate(0, rfr), DT::deltaTime(0, dt));
  probs = implied_probs.at(0);
  if (cmp(spot_lattice(0, 0), barrier_val)) {
    tpl = probs.at(0);
    value = df * (std::get<0>(tpl) * option_lattice(1, 0) +
                  std::get<1>(tpl) * option_lattice(1, 1) +
                  std::get<2>(tpl) * option_lattice(1, 2));
    payoff_adj(value, spot_lattice(0, 0));
    option_lattice(0, 0) = value;
  } else {
    option_lattice(0, 0) = rebate_val;
  }
}

} // namespace lattice

#endif ///_LATTICE_BACKWARD_TRAVERSALS