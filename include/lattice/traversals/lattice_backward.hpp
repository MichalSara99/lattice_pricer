#if !defined(_LATTICE_BACKWARD_HPP_)
#define _LATTICE_BACKWARD_HPP_

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

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator, typename Payoff>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Binomial, TimeAxis,
    DeltaTime>::_backTraverse(LatticeObject &lattice, Generator &&generator,
                              DeltaTime const &deltaTime, Payoff &&payoff) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  const std::size_t lastIdx = lattice.timeDimension() - 1;
  const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();
  typename LatticeObject::Node_type dt{};

  for (auto i = 0; i < lastNodesSize; ++i) {
    lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
  }
  std::size_t nodesSize{0};
  for (auto n = lastIdx - 1; n > 0; --n) {
    nodesSize = lattice.nodesAtIdx(n).size();
    dt = DT::deltaTime(n, deltaTime);
    for (auto i = 0; i < nodesSize; ++i) {
      lattice(n, i) = generator(lattice(n, i), lattice(n + 1, i),
                                lattice(n + 1, i + 1), dt);
    }
  }
  dt = DT::deltaTime(0, deltaTime);
  lattice(0, 0) = generator(lattice(0, 0), lattice(1, 0), lattice(1, 1), dt);
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator, typename Payoff>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime>::
    _backTraverseBarrier(LatticeObject &lattice, Generator &&generator,
                         DeltaTime const &deltaTime, Payoff &&payoff,
                         BarrierType barrierType,
                         typename LatticeObject::Node_type const &barrier,
                         typename LatticeObject::Node_type const &rebate) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  typedef typename LatticeObject::Node_type Node;
  typedef BarrierComparer<Node> BC;
  auto cmp = BC::comparer(barrierType);
  const std::size_t lastIdx = lattice.timeDimension() - 1;
  const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();
  Node dt{};

  for (auto i = 0; i < lastNodesSize; ++i) {
    if (cmp(lattice(lastIdx, i), barrier)) {
      lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
    } else {
      lattice(lastIdx, i) = Node{};
    }
  }
  std::size_t nodesSize{0};
  for (auto n = lastIdx - 1; n > 0; --n) {
    nodesSize = lattice.nodesAtIdx(n).size();
    dt = DT::deltaTime(n, deltaTime);
    for (auto i = 0; i < nodesSize; ++i) {
      if (cmp(lattice(n, i), barrier)) {
        lattice(n, i) = generator(lattice(n, i), lattice(n + 1, i),
                                  lattice(n + 1, i + 1), dt);
      } else {
        lattice(n, i) = rebate;
      }
    }
  }
  dt = DT::deltaTime(0, deltaTime);
  if (cmp(lattice(0, 0), barrier)) {
    lattice(0, 0) = generator(lattice(0, 0), lattice(1, 0), lattice(1, 1), dt);
  } else {
    lattice(0, 0) = rebate;
  }
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator, typename Payoff,
          typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Binomial, TimeAxis,
    DeltaTime>::_backTraverse(LatticeObject &lattice, Generator &&generator,
                              DeltaTime const &deltaTime, Payoff &&payoff,
                              PayoffAdjuster &&payoffAdjuster) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  const std::size_t lastIdx = lattice.timeDimension() - 1;
  const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();
  typename LatticeObject::Node_type dt{};
  typename LatticeObject::Node_type value{};

  for (auto i = 0; i < lastNodesSize; ++i) {
    lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
  }

  std::size_t nodesSize{0};
  for (auto n = lastIdx - 1; n > 0; --n) {
    nodesSize = lattice.nodesAtIdx(n).size();
    dt = DT::deltaTime(n, deltaTime);
    for (auto i = 0; i < nodesSize; ++i) {
      value = generator(lattice(n, i), lattice(n + 1, i), lattice(n + 1, i + 1),
                        dt);
      payoffAdjuster(value, lattice(n, i));
      lattice(n, i) = value;
    }
  }
  dt = DT::deltaTime(0, deltaTime);
  value = generator(lattice(0, 0), lattice(1, 0), lattice(1, 1), dt);
  payoffAdjuster(value, lattice(0, 0));
  lattice(0, 0) = value;
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator, typename Payoff,
          typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime>::
    _backTraverseBarrier(LatticeObject &lattice, Generator &&generator,
                         DeltaTime const &deltaTime, Payoff &&payoff,
                         PayoffAdjuster &&payoffAdjuster,
                         BarrierType barrierType,
                         typename LatticeObject::Node_type const &barrier,
                         typename LatticeObject::Node_type const &rebate) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  typedef typename LatticeObject::Node_type Node;
  typedef BarrierComparer<Node> BC;
  const std::size_t lastIdx = lattice.timeDimension() - 1;
  const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();
  Node dt{};
  Node value{};

  for (auto i = 0; i < lastNodesSize; ++i) {
    if (cmp(lattice(lastIdx, i), barrier)) {
      lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
    } else {
      lattice(lastIdx, i) = Node{};
    }
  }

  std::size_t nodesSize{0};
  for (auto n = lastIdx - 1; n > 0; --n) {
    nodesSize = lattice.nodesAtIdx(n).size();
    dt = DT::deltaTime(n, deltaTime);
    for (auto i = 0; i < nodesSize; ++i) {
      if (cmp(lattice(n, i), barrier)) {
        value = generator(lattice(n, i), lattice(n + 1, i),
                          lattice(n + 1, i + 1), dt);
        payoffAdjuster(value, lattice(n, i));
        lattice(n, i) = value;
      } else {
        lattice(n, i) = rebate;
      }
    }
  }
  dt = DT::deltaTime(0, deltaTime);
  if (cmp(lattice(0, 0), barrier)) {
    value = generator(lattice(0, 0), lattice(1, 0), lattice(1, 1), dt);
    payoffAdjuster(value, lattice(0, 0));
    lattice(0, 0) = value;
  } else {
    lattice(0, 0) = rebate;
  }
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename MultidimLatticeObject,
          typename Generator, typename Payoff>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::TwoVariableBinomial, TimeAxis,
    DeltaTime>::_backTraverse(LatticeObject &priceLattice,
                              MultidimLatticeObject const &lattice,
                              Generator &&generator, DeltaTime const &deltaTime,
                              Payoff &&payoff) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  typename LatticeObject::Node_type dt{};

  auto tree1 = lattice.getFactor(0);
  auto tree2 = lattice.getFactor(1);

  std::size_t const treeSize = priceLattice.timeDimension();
  std::size_t const lastIdx = treeSize - 1;
  std::size_t const lastNodesSize = priceLattice.nodesAtIdx(lastIdx).size();
  std::size_t factorLastNodesSize = tree1.nodesAtIdx(lastIdx).size();
  std::size_t col{0}, row{0};

  for (auto i = 0; i < lastNodesSize; ++i) {
    col = i % factorLastNodesSize;
    if ((i > 0) && ((i % factorLastNodesSize) == 0))
      row++;
    priceLattice(lastIdx, i) = payoff(tree1(lastIdx, row), tree2(lastIdx, col));
  }

  std::size_t nodesSize{0};
  std::size_t factorNodesSize{0};
  for (auto n = lastIdx - 1; n > 0; --n) {
    nodesSize = priceLattice.nodesAtIdx(n).size();
    factorNodesSize = tree1.nodesAtIdx(n).size();
    dt = DT::deltaTime(n, deltaTime);
    row = 0;
    for (auto i = 0; i < nodesSize; ++i) {
      if ((i > 0) && ((i % factorNodesSize) == 0))
        row++;
      priceLattice(n, i) = generator(
          priceLattice(n, i), /*down-down*/ priceLattice(n + 1, i + row),
          /*down-up*/ priceLattice(n + 1, i + row + 1),
          /*up-down*/ priceLattice(n + 1, i + row + factorLastNodesSize),
          /*up-up*/ priceLattice(n + 1, i + row + factorLastNodesSize + 1), dt);
    }
    factorLastNodesSize = factorNodesSize;
  }
  dt = DT::deltaTime(0, deltaTime);
  priceLattice(0, 0) =
      generator(priceLattice(0, 0), priceLattice(1, 0), priceLattice(1, 1),
                priceLattice(1, 2), priceLattice(1, 3), dt);
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename MultidimLatticeObject,
          typename Generator, typename Payoff, typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::TwoVariableBinomial, TimeAxis,
    DeltaTime>::_backTraverse(LatticeObject &priceLattice,
                              MultidimLatticeObject const &lattice,
                              Generator &&generator, DeltaTime const &deltaTime,
                              Payoff &&payoff,
                              PayoffAdjuster &&payoffAdjuster) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  typename LatticeObject::Node_type dt{};
  typename LatticeObject::Node_type value{};

  auto tree1 = lattice.getFactor(0);
  auto tree2 = lattice.getFactor(1);

  std::size_t const treeSize = priceLattice.timeDimension();
  std::size_t const lastIdx = treeSize - 1;
  std::size_t const lastNodesSize = priceLattice.nodesAtIdx(lastIdx).size();
  std::size_t factorLastNodesSize = tree1.nodesAtIdx(lastIdx).size();
  std::size_t col{0}, row{0};

  for (auto i = 0; i < lastNodesSize; ++i) {
    col = i % factorLastNodesSize;
    if ((i > 0) && ((i % factorLastNodesSize) == 0))
      row++;
    priceLattice(lastIdx, i) = payoff(tree1(lastIdx, row), tree2(lastIdx, col));
  }

  std::size_t nodesSize{0};
  std::size_t factorNodesSize{0};
  for (auto n = lastIdx - 1; n > 0; --n) {
    nodesSize = priceLattice.nodesAtIdx(n).size();
    factorNodesSize = tree1.nodesAtIdx(n).size();
    dt = DT::deltaTime(n, deltaTime);
    col = 0;
    row = 0;
    for (auto i = 0; i < nodesSize; ++i) {
      col = i % factorNodesSize;
      if ((i > 0) && ((i % factorNodesSize) == 0))
        row++;
      value = generator(
          priceLattice(n, i), /*down-down*/ priceLattice(n + 1, i + row),
          /*down-up*/ priceLattice(n + 1, i + row + 1),
          /*up-down*/ priceLattice(n + 1, i + row + factorLastNodesSize),
          /*up-up*/ priceLattice(n + 1, i + row + factorLastNodesSize + 1), dt);
      payoffAdjuster(value, tree1(n, row), tree2(n, col));
      priceLattice(n, i) = value;
    }
    factorLastNodesSize = factorNodesSize;
  }
  dt = DT::deltaTime(0, deltaTime);
  value = generator(priceLattice(0, 0), priceLattice(1, 0), priceLattice(1, 1),
                    priceLattice(1, 2), priceLattice(1, 3), dt);
  payoffAdjuster(value, tree1(0, 0), tree2(0, 0));
  priceLattice(0, 0) = value;
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Trinomial, TimeAxis,
    DeltaTime>::_backTraverseNormal(std::size_t timeIdx, LatticeObject &lattice,
                                    Generator &&generator,
                                    DeltaTime const &deltaTime) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  typename LatticeObject::Node_type dt{};

  std::size_t const revertBranchesSize = timeIdx;

  std::size_t nodesSize{0};
  for (auto n = timeIdx - 1; n > 0; --n) {
    dt = DT::deltaTime(n, deltaTime);
    nodesSize = lattice.nodesAtIdx(n).size();
    for (auto i = 0; i < nodesSize; ++i) {
      lattice(n, i) = generator(lattice(n, i), lattice(n + 1, i),
                                lattice(n + 1, i + 1), lattice(n + 1, i + 2),
                                dt, revertBranchesSize, nodesSize, i);
    }
  }
  dt = DT::deltaTime(0, deltaTime);
  lattice(0, 0) = generator(lattice(0, 0), lattice(1, 0), lattice(1, 1),
                            lattice(1, 2), dt, revertBranchesSize, 1, 0);
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Trinomial, TimeAxis,
    DeltaTime>::_backTraverseReverting(std::size_t timeIdx,
                                       LatticeObject &lattice,
                                       Generator &&generator,
                                       DeltaTime const &deltaTime) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  typename LatticeObject::Node_type dt{};
  const std::size_t lastIdx = lattice.timeDimension() - 1;

  std::size_t const revertBranchesSize = timeIdx;

  std::size_t nodesSize{0};
  for (auto n = lastIdx - 1; n >= timeIdx; --n) {
    dt = DT::deltaTime(n, deltaTime);
    nodesSize = lattice.nodesAtIdx(n).size();
    lattice(n, 0) =
        generator(lattice(n, 0), lattice(n + 1, 0), lattice(n + 1, 1),
                  lattice(n + 1, 2), dt, revertBranchesSize, nodesSize, 0);
    for (auto i = 1; i < nodesSize - 1; ++i) {
      lattice(n, i) = generator(lattice(n, i), lattice(n + 1, i - 1),
                                lattice(n + 1, i), lattice(n + 1, i + 1), dt,
                                revertBranchesSize, nodesSize, i);
    }
    lattice(n, nodesSize - 1) =
        generator(lattice(n, nodesSize - 1), lattice(n + 1, nodesSize - 3),
                  lattice(n + 1, nodesSize - 2), lattice(n + 1, nodesSize - 1),
                  dt, revertBranchesSize, nodesSize, nodesSize - 1);
  }
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator, typename Payoff>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Trinomial, TimeAxis,
    DeltaTime>::_backTraverse(LatticeObject &lattice, Generator &&generator,
                              DeltaTime const &deltaTime, Payoff &&payoff) {

  const std::size_t firstRevertIdx = lattice.firstRevertingIdx();
  const std::size_t treeSize = lattice.timeDimension();

  const std::size_t lastIdx = treeSize - 1;
  const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

  for (auto i = 0; i < lastNodesSize; ++i) {
    lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
  }

  if (firstRevertIdx == 0) {
    // This trinomial tree does not have reverting property:
    _backTraverseNormal(lastIdx, lattice, std::forward<Generator>(generator),
                        deltaTime);
  } else {
    // This trinomial tree does have reverting property:
    _backTraverseReverting(firstRevertIdx - 1, lattice,
                           std::forward<Generator>(generator), deltaTime);
    _backTraverseNormal(firstRevertIdx - 1, lattice,
                        std::forward<Generator>(generator), deltaTime);
  }
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
    _backTraverseNormalBarrier(
        std::size_t timeIdx, LatticeObject &lattice, Generator &&generator,
        DeltaTime const &deltaTime, BarrierType barrierType,
        typename LatticeObject::Node_type const &barrier,
        typename LatticeObject::Node_type const &rebate) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  typedef typename LatticeObject::Node_type Node;
  typedef BarrierComparer<Node> BC;
  std::size_t const revertBranchesSize = timeIdx;
  auto cmp = BC::comparer(barrierType);

  Node dt{};
  std::size_t nodesSize{0};
  for (auto n = timeIdx - 1; n > 0; --n) {
    dt = DT::deltaTime(n, deltaTime);
    nodesSize = lattice.nodesAtIdx(n).size();
    for (auto i = 0; i < nodesSize; ++i) {
      if (cmp(lattice(n, i), barrier)) {
        lattice(n, i) = generator(lattice(n, i), lattice(n + 1, i),
                                  lattice(n + 1, i + 1), lattice(n + 1, i + 2),
                                  dt, revertBranchesSize, nodesSize, i);
      } else {
        lattice(n, i) = rebate;
      }
    }
  }
  dt = DT::deltaTime(0, deltaTime);
  if (cmp(lattice(0, 0), barrier)) {
    lattice(0, 0) = generator(lattice(0, 0), lattice(1, 0), lattice(1, 1),
                              lattice(1, 2), dt, revertBranchesSize, 1, 0);
  } else {
    lattice(0, 0) = rebate;
  }
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
    _backTraverseRevertingBarrier(
        std::size_t timeIdx, LatticeObject &lattice, Generator &&generator,
        DeltaTime const &deltaTime, BarrierType barrierType,
        typename LatticeObject::Node_type const &barrier,
        typename LatticeObject::Node_type const &rebate) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  typedef typename LatticeObject::Node_type Node;
  typedef BarrierComparer<Node> BC;
  auto cmp = BC::comparer(barrierType);
  const std::size_t lastIdx = lattice.timeDimension() - 1;
  std::size_t const revertBranchesSize = timeIdx;

  Node dt{};
  std::size_t nodesSize{0};
  for (auto n = lastIdx - 1; n >= timeIdx; --n) {
    dt = DT::deltaTime(n, deltaTime);
    nodesSize = lattice.nodesAtIdx(n).size();
    if (cmp(lattice(n, 0), barrier)) {
      lattice(n, 0) =
          generator(lattice(n, 0), lattice(n + 1, 0), lattice(n + 1, 1),
                    lattice(n + 1, 2), dt, revertBranchesSize, nodesSize, 0);
    } else {
      lattice(n, 0) = rebate;
    }
    for (auto i = 1; i < nodesSize - 1; ++i) {
      if (cmp(lattice(n, i), barrier)) {
        lattice(n, i) = generator(lattice(n, i), lattice(n + 1, i - 1),
                                  lattice(n + 1, i), lattice(n + 1, i + 1), dt,
                                  revertBranchesSize, nodesSize, i);
      } else {
        lattice(n, i) = rebate;
      }
    }
    if (cmp(lattice(n, nodesSize - 1), barrier)) {
      lattice(n, nodesSize - 1) = generator(
          lattice(n, nodesSize - 1), lattice(n + 1, nodesSize - 3),
          lattice(n + 1, nodesSize - 2), lattice(n + 1, nodesSize - 1), dt,
          revertBranchesSize, nodesSize, nodesSize - 1);
    } else {
      lattice(n, nodesSize - 1) = rebate;
    }
  }
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator, typename Payoff>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
    _backTraverseBarrier(LatticeObject &lattice, Generator &&generator,
                         DeltaTime const &deltaTime, Payoff &&payoff,
                         BarrierType barrierType,
                         typename LatticeObject::Node_type const &barrier,
                         typename LatticeObject::Node_type const &rebate) {

  typedef typename LatticeObject::Node_type Node;
  typedef BarrierComparer<Node> BC;
  auto cmp = BC::comparer(barrierType);
  const std::size_t firstRevertIdx = lattice.firstRevertingIdx();
  const std::size_t treeSize = lattice.timeDimension();

  const std::size_t lastIdx = treeSize - 1;
  const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

  for (auto i = 0; i < lastNodesSize; ++i) {
    if (cmp(lattice(lastIdx, i), barrier)) {
      lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
    } else {
      lattice(lastIdx, i) = Node{};
    }
  }

  if (firstRevertIdx == 0) {
    // This trinomial tree does not have reverting property:
    _backTraverseNormalBarrier(lastIdx, lattice,
                               std::forward<Generator>(generator), deltaTime,
                               barrierType, barrier, rebate);
  } else {
    // This trinomial tree does have reverting property:
    _backTraverseRevertingBarrier(firstRevertIdx - 1, lattice,
                                  std::forward<Generator>(generator), deltaTime,
                                  barrierType, barrier, rebate);
    _backTraverseNormalBarrier(firstRevertIdx - 1, lattice,
                               std::forward<Generator>(generator), deltaTime,
                               barrierType, barrier, rebate);
  }
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator, typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Trinomial, TimeAxis,
    DeltaTime>::_backTraverseNormal(std::size_t timeIdx, LatticeObject &lattice,
                                    Generator &&generator,
                                    DeltaTime const &deltaTime,
                                    PayoffAdjuster &&payoffAdjuster) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  std::size_t const revertBranchesSize = timeIdx;
  typename LatticeObject::Node_type dt{};
  typename LatticeObject::Node_type value{};

  std::size_t nodesSize{0};
  for (auto n = timeIdx - 1; n > 0; --n) {
    dt = DT::deltaTime(n, deltaTime);
    nodesSize = lattice.nodesAtIdx(n).size();
    for (auto i = 0; i < nodesSize; ++i) {
      value = generator(lattice(n, i), lattice(n + 1, i), lattice(n + 1, i + 1),
                        lattice(n + 1, i + 2), dt, revertBranchesSize,
                        nodesSize, i);
      payoffAdjuster(value, lattice(n, i));
      lattice(n, i) = value;
    }
  }
  dt = DT::deltaTime(0, deltaTime);
  value = generator(lattice(0, 0), lattice(1, 0), lattice(1, 1), lattice(1, 2),
                    dt, revertBranchesSize, 1, 0);
  payoffAdjuster(value, lattice(0, 0));
  lattice(0, 0) = value;
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator, typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Trinomial, TimeAxis,
    DeltaTime>::_backTraverseReverting(std::size_t timeIdx,
                                       LatticeObject &lattice,
                                       Generator &&generator,
                                       DeltaTime const &deltaTime,
                                       PayoffAdjuster &&payoffAdjuster) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  typename LatticeObject::Node_type dt{};
  typename LatticeObject::Node_type value{};
  const std::size_t lastIdx = lattice.timeDimension() - 1;
  std::size_t const revertBranchesSize = timeIdx;

  std::size_t nodesSize{0};
  for (auto n = lastIdx - 1; n >= timeIdx; --n) {
    dt = DT::deltaTime(n, deltaTime);
    nodesSize = lattice.nodesAtIdx(n).size();

    value = generator(lattice(n, 0), lattice(n + 1, 0), lattice(n + 1, 1),
                      lattice(n + 1, 2), dt, revertBranchesSize, nodesSize, 0);
    payoffAdjuster(value, lattice(n, 0));
    lattice(n, 0) = value;

    for (auto i = 1; i < nodesSize - 1; ++i) {
      value = generator(lattice(n, i), lattice(n + 1, i - 1), lattice(n + 1, i),
                        lattice(n + 1, i + 1), dt, revertBranchesSize,
                        nodesSize, i);
      payoffAdjuster(value, lattice(n, i));
      lattice(n, i) = value;
    }

    value =
        generator(lattice(n, nodesSize - 1), lattice(n + 1, nodesSize - 3),
                  lattice(n + 1, nodesSize - 2), lattice(n + 1, nodesSize - 1),
                  dt, revertBranchesSize, nodesSize, nodesSize - 1);
    payoffAdjuster(value, lattice(n, nodesSize - 1));
    lattice(n, nodesSize - 1) = value;
  }
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator, typename Payoff,
          typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Trinomial, TimeAxis,
    DeltaTime>::_backTraverse(LatticeObject &lattice, Generator &&generator,
                              DeltaTime const &deltaTime, Payoff &&payoff,
                              PayoffAdjuster &&payoffAdjuster) {

  const std::size_t firstRevertIdx = lattice.firstRevertingIdx();
  const std::size_t treeSize = lattice.timeDimension();

  const std::size_t lastIdx = treeSize - 1;
  const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

  for (auto i = 0; i < lastNodesSize; ++i) {
    lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
  }

  if (firstRevertIdx == 0) {
    // This trinomial tree does not have reverting property:
    _backTraverseNormal(lastIdx, lattice, std::forward<Generator>(generator),
                        deltaTime,
                        std::forward<PayoffAdjuster>(payoffAdjuster));
  } else {
    // This trinomial tree does have reverting property:
    _backTraverseReverting(firstRevertIdx - 1, lattice,
                           std::forward<Generator>(generator), deltaTime,
                           std::forward<PayoffAdjuster>(payoffAdjuster));
    _backTraverseNormal(firstRevertIdx - 1, lattice,
                        std::forward<Generator>(generator), deltaTime,
                        std::forward<PayoffAdjuster>(payoffAdjuster));
  }
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator, typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
    _backTraverseNormalBarrier(
        std::size_t timeIdx, LatticeObject &lattice, Generator &&generator,
        DeltaTime const &deltaTime, PayoffAdjuster &&payoffAdjuster,
        BarrierType barrierType,
        typename LatticeObject::Node_type const &barrier,
        typename LatticeObject::Node_type const &rebate) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  typedef typename LatticeObject::Node_type Node;
  typedef BarrierComparer<Node> BC;
  auto cmp = BC::comparer(barrierType);
  std::size_t const revertBranchesSize = timeIdx;

  Node dt{};
  Node value{};
  std::size_t nodesSize{0};
  for (auto n = timeIdx - 1; n > 0; --n) {
    dt = DT::deltaTime(n, deltaTime);
    nodesSize = lattice.nodesAtIdx(n).size();
    for (auto i = 0; i < nodesSize; ++i) {
      if (cmp(lattice(n, i), barrier)) {
        value = generator(lattice(n, i), lattice(n + 1, i),
                          lattice(n + 1, i + 1), lattice(n + 1, i + 2), dt,
                          revertBranchesSize, nodesSize, i);
        payoffAdjuster(value, lattice(n, i));
        lattice(n, i) = value;
      } else {
        lattice(n, i) = rebate;
      }
    }
  }
  dt = DT::deltaTime(0, deltaTime);
  if (cmp(lattice(0, 0), barrier)) {
    value = generator(lattice(0, 0), lattice(1, 0), lattice(1, 1),
                      lattice(1, 2), dt, revertBranchesSize, 1, 0);
    payoffAdjuster(value, lattice(0, 0));
    lattice(0, 0) = value;
  } else {
    lattice(0, 0) = rebate;
  }
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator, typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
    _backTraverseRevertingBarrier(
        std::size_t timeIdx, LatticeObject &lattice, Generator &&generator,
        DeltaTime const &deltaTime, PayoffAdjuster &&payoffAdjuster,
        BarrierType barrierType,
        typename LatticeObject::Node_type const &barrier,
        typename LatticeObject::Node_type const &rebate) {

  typedef DeltaTimeHolder<DeltaTime> DT;
  typedef typename LatticeObject::Node_type Node;
  typedef BarrierComparer<Node> BC;
  auto cmp = BC::comparer(barrierType);
  const std::size_t lastIdx = lattice.timeDimension() - 1;
  std::size_t const revertBranchesSize = timeIdx;

  Node dt{};
  Node value{};
  std::size_t nodesSize{0};
  for (auto n = lastIdx - 1; n >= timeIdx; --n) {
    dt = DT::deltaTime(n, deltaTime);
    nodesSize = lattice.nodesAtIdx(n).size();

    if (cmp(lattice(n, 0), barrier)) {
      value =
          generator(lattice(n, 0), lattice(n + 1, 0), lattice(n + 1, 1),
                    lattice(n + 1, 2), dt, revertBranchesSize, nodesSize, 0);
      payoffAdjuster(value, lattice(n, 0));
      lattice(n, 0) = value;
    } else {
      lattice(n, 0) = rebate;
    }

    for (auto i = 1; i < nodesSize - 1; ++i) {
      if (cmp(lattice(n, i), barrier)) {
        value = generator(lattice(n, i), lattice(n + 1, i - 1),
                          lattice(n + 1, i), lattice(n + 1, i + 1), dt,
                          revertBranchesSize, nodesSize, i);
        payoffAdjuster(value, lattice(n, i));
        lattice(n, i) = value;
      } else {
        lattice(n, i) = rebate;
      }
    }
    if (cmp(lattice(n, nodesSize - 1), barrier)) {
      value = generator(
          lattice(n, nodesSize - 1), lattice(n + 1, nodesSize - 3),
          lattice(n + 1, nodesSize - 2), lattice(n + 1, nodesSize - 1), dt,
          revertBranchesSize, nodesSize, nodesSize - 1);
      payoffAdjuster(value, lattice(n, nodesSize - 1));
      lattice(n, nodesSize - 1) = value;
    } else {
      lattice(n, nodesSize - 1) = rebate;
    }
  }
}

template <typename TimeAxis, typename DeltaTime>
template <typename LatticeObject, typename Generator, typename Payoff,
          typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<
    lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
    _backTraverseBarrier(LatticeObject &lattice, Generator &&generator,
                         DeltaTime const &deltaTime, Payoff &&payoff,
                         PayoffAdjuster &&payoffAdjuster,
                         BarrierType barrierType,
                         typename LatticeObject::Node_type const &barrier,
                         typename LatticeObject::Node_type const &rebate) {

  typedef typename LatticeObject::Node_type Node;
  typedef BarrierComparer<Node> BC;

  auto cmp = BC::comparer(barrierType);
  const std::size_t firstRevertIdx = lattice.firstRevertingIdx();
  const std::size_t treeSize = lattice.timeDimension();

  const std::size_t lastIdx = treeSize - 1;
  const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

  for (auto i = 0; i < lastNodesSize; ++i) {
    if (cmp(lattice(lastIdx, i), barrier)) {
      lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
    } else {
      lattice(lastIdx, i) = Node{};
    }
  }

  if (firstRevertIdx == 0) {
    // This trinomial tree does not have reverting property:
    _backTraverseNormalBarrier(lastIdx, lattice,
                               std::forward<Generator>(generator), deltaTime,
                               std::forward<PayoffAdjuster>(payoffAdjuster),
                               barrierType, barrier, rebate);
  } else {
    // This trinomial tree does have reverting property:
    _backTraverseRevertingBarrier(firstRevertIdx - 1, lattice,
                                  std::forward<Generator>(generator), deltaTime,
                                  std::forward<PayoffAdjuster>(payoffAdjuster),
                                  barrierType, barrier, rebate);
    _backTraverseNormalBarrier(firstRevertIdx - 1, lattice,
                               std::forward<Generator>(generator), deltaTime,
                               std::forward<PayoffAdjuster>(payoffAdjuster),
                               barrierType, barrier, rebate);
  }
}

template <typename DeltaTime, typename RiskFreeRate>
template <typename LatticeObject, typename Payoff>
void lattice_backward_traversals::ImpliedBackwardTraversal<
    lattice_types::LatticeType::Trinomial, DeltaTime,
    RiskFreeRate>::_backTraverse(LatticeObject &optionLattice,
                                 LatticeObject const &spotLattice,
                                 CalibratorTrinomialEquityResultsPtr<
                                     LatticeObject> const &calibrationResults,
                                 DeltaTime const &deltaTime,
                                 RiskFreeRate const &riskFreeRate,
                                 Payoff &&payoff, DiscountingStyle style) {

  // typedefs:
  typedef typename LatticeObject::Node_type Node;
  // typedef RiskFreeRateHolder:
  typedef RiskFreeRateHolder<RiskFreeRate> RFR;
  // typedef DeltaTimeHolder:
  typedef DeltaTimeHolder<DeltaTime> DT;
  // typedef discounting factor:
  typedef DiscountingFactor<Node> DCF;

  // get correct discounting factor style:
  auto dcf = DCF::function(style);

  // get the size of optionLattice object:
  std::size_t const treeSize = optionLattice.timeDimension();
  std::size_t const lastIdx = treeSize - 1;
  std::size_t const lastNodesSize = optionLattice.nodesAtIdx(lastIdx).size();

  // unpack the impled probabilities from calibrationResults:
  auto impliedProbs = calibrationResults->impliedProbabilities;
  LASSERT(impliedProbs.size() >= (treeSize - 1),
          "Must have enough implied volatilities for pricing.");

  // last payoff first:
  for (std::size_t i = 0; i < lastNodesSize; ++i) {
    optionLattice(lastIdx, i) = payoff(spotLattice(lastIdx, i));
  }

  std::size_t nodesSize{0};
  Node df{};
  std::vector<std::tuple<Node, Node, Node>> probs;
  std::tuple<Node, Node, Node> tpl;
  for (auto t = lastIdx - 1; t > 0; --t) {
    df = dcf(RFR::rate(t, riskFreeRate), DT::deltaTime(t, deltaTime));
    nodesSize = optionLattice.nodesAtIdx(t).size();
    probs = impliedProbs.at(t);
    for (auto i = 0; i < nodesSize; ++i) {
      tpl = probs.at(i);
      optionLattice(t, i) =
          df * (std::get<0>(tpl) * optionLattice(t + 1, i) +
                std::get<1>(tpl) * optionLattice(t + 1, i + 1) +
                std::get<2>(tpl) * optionLattice(t + 1, i + 2));
    }
  }
  df = dcf(RFR::rate(0, riskFreeRate), DT::deltaTime(0, deltaTime));
  probs = impliedProbs.at(0);
  tpl = probs.at(0);
  optionLattice(0, 0) = df * (std::get<0>(tpl) * optionLattice(1, 0) +
                              std::get<1>(tpl) * optionLattice(1, 1) +
                              std::get<2>(tpl) * optionLattice(1, 2));
}

template <typename DeltaTime, typename RiskFreeRate>
template <typename LatticeObject, typename Payoff, typename PayoffAdjuster>
void lattice_backward_traversals::ImpliedBackwardTraversal<
    lattice_types::LatticeType::Trinomial, DeltaTime, RiskFreeRate>::
    _backTraverse(LatticeObject &optionLattice,
                  LatticeObject const &spotLattice,
                  CalibratorTrinomialEquityResultsPtr<LatticeObject> const
                      &calibrationResults,
                  DeltaTime const &deltaTime, RiskFreeRate const &riskFreeRate,
                  Payoff &&payoff, PayoffAdjuster &&payoffAdjuster,
                  DiscountingStyle style) {

  // typedefs:
  typedef typename LatticeObject::Node_type Node;
  // typedef RiskFreeRateHolder:
  typedef RiskFreeRateHolder<RiskFreeRate> RFR;
  // typedef DeltaTimeHolder:
  typedef DeltaTimeHolder<DeltaTime> DT;
  // typedef discounting factor:
  typedef DiscountingFactor<Node> DCF;

  // get correct discounting factor style:
  auto dcf = DCF::function(style);

  // get the size of optionLattice object:
  std::size_t const treeSize = optionLattice.timeDimension();
  std::size_t const lastIdx = treeSize - 1;
  std::size_t const lastNodesSize = optionLattice.nodesAtIdx(lastIdx).size();

  // unpack the impled probabilities from calibrationResults:
  auto impliedProbs = calibrationResults->impliedProbabilities;
  LASSERT(impliedProbs.size() >= (treeSize - 1),
          "Must have enough implied volatilities for pricing.");

  // last payoff first:
  for (std::size_t i = 0; i < lastNodesSize; ++i) {
    optionLattice(lastIdx, i) = payoff(spotLattice(lastIdx, i));
  }

  std::size_t nodesSize{0};
  Node df{};
  Node value{};
  std::vector<std::tuple<Node, Node, Node>> probs;
  std::tuple<Node, Node, Node> tpl;
  for (auto t = lastIdx - 1; t > 0; --t) {
    df = dcf(RFR::rate(t, riskFreeRate), DT::deltaTime(t, deltaTime));
    nodesSize = optionLattice.nodesAtIdx(t).size();
    probs = impliedProbs.at(t);
    for (auto i = 0; i < nodesSize; ++i) {
      tpl = probs.at(i);
      value = df * (std::get<0>(tpl) * optionLattice(t + 1, i) +
                    std::get<1>(tpl) * optionLattice(t + 1, i + 1) +
                    std::get<2>(tpl) * optionLattice(t + 1, i + 2));
      payoffAdjuster(value, spotLattice(t, i));
      optionLattice(t, i) = value;
    }
  }
  df = dcf(RFR::rate(0, riskFreeRate), DT::deltaTime(0, deltaTime));
  probs = impliedProbs.at(0);
  tpl = probs.at(0);
  value = df * (std::get<0>(tpl) * optionLattice(1, 0) +
                std::get<1>(tpl) * optionLattice(1, 1) +
                std::get<2>(tpl) * optionLattice(1, 2));
  payoffAdjuster(value, spotLattice(0, 0));
  optionLattice(0, 0) = value;
}

template <typename DeltaTime, typename RiskFreeRate>
template <typename LatticeObject, typename Payoff>
void lattice_backward_traversals::ImpliedBackwardTraversal<
    lattice_types::LatticeType::Trinomial, DeltaTime, RiskFreeRate>::
    _backTraverseBarrier(
        LatticeObject &optionLattice, LatticeObject const &spotLattice,
        CalibratorTrinomialEquityResultsPtr<LatticeObject> const
            &calibrationResults,
        DeltaTime const &deltaTime, RiskFreeRate const &riskFreeRate,
        Payoff &&payoff, BarrierType barrierType,
        typename LatticeObject::Node_type const &barrier,
        typename LatticeObject::Node_type const &rebate,
        DiscountingStyle style) {

  // typedefs:
  typedef typename LatticeObject::Node_type Node;
  // typedef RiskFreeRateHolder:
  typedef RiskFreeRateHolder<RiskFreeRate> RFR;
  // typedef DeltaTimeHolder:
  typedef DeltaTimeHolder<DeltaTime> DT;
  // typedef discounting factor:
  typedef DiscountingFactor<Node> DCF;
  // typedef BarrierComparer:
  typedef BarrierComparer<Node> BC;

  // get correct comparer function:
  auto cmp = BC::comparer(barrierType);
  // get correct discounting factor style:
  auto dcf = DCF::function(style);

  // get the size of optionLattice object:
  std::size_t const treeSize = optionLattice.timeDimension();
  std::size_t const lastIdx = treeSize - 1;
  std::size_t const lastNodesSize = optionLattice.nodesAtIdx(lastIdx).size();

  // unpack the impled probabilities from calibrationResults:
  auto impliedProbs = calibrationResults->impliedProbabilities;
  LASSERT(impliedProbs.size() >= (treeSize - 1),
          "Must have enough implied volatilities for pricing.");

  // last payoff first:
  for (std::size_t i = 0; i < lastNodesSize; ++i) {
    if (cmp(spotLattice(lastIdx, i), barrier)) {
      optionLattice(lastIdx, i) = payoff(spotLattice(lastIdx, i));
    } else {
      optionLattice(lastIdx, i) = 0.0;
    }
  }

  std::size_t nodesSize{0};
  Node df{};
  std::vector<std::tuple<Node, Node, Node>> probs;
  std::tuple<Node, Node, Node> tpl;
  for (auto t = lastIdx - 1; t > 0; --t) {
    df = dcf(RFR::rate(t, riskFreeRate), DT::deltaTime(t, deltaTime));
    nodesSize = optionLattice.nodesAtIdx(t).size();
    probs = impliedProbs.at(t);
    for (auto i = 0; i < nodesSize; ++i) {
      if (cmp(spotLattice(t, i), barrier)) {
        tpl = probs.at(i);
        optionLattice(t, i) =
            df * (std::get<0>(tpl) * optionLattice(t + 1, i) +
                  std::get<1>(tpl) * optionLattice(t + 1, i + 1) +
                  std::get<2>(tpl) * optionLattice(t + 1, i + 2));
      } else {
        optionLattice(t, i) = rebate;
      }
    }
  }

  df = dcf(RFR::rate(0, riskFreeRate), DT::deltaTime(0, deltaTime));
  probs = impliedProbs.at(0);
  if (cmp(spotLattice(0, 0), barrier)) {
    tpl = probs.at(0);
    optionLattice(0, 0) = df * (std::get<0>(tpl) * optionLattice(1, 0) +
                                std::get<1>(tpl) * optionLattice(1, 1) +
                                std::get<2>(tpl) * optionLattice(1, 2));
  } else {
    optionLattice(0, 0) = rebate;
  }
}

template <typename DeltaTime, typename RiskFreeRate>
template <typename LatticeObject, typename Payoff, typename PayoffAdjuster>
void lattice_backward_traversals::ImpliedBackwardTraversal<
    lattice_types::LatticeType::Trinomial, DeltaTime, RiskFreeRate>::
    _backTraverseBarrier(
        LatticeObject &optionLattice, LatticeObject const &spotLattice,
        CalibratorTrinomialEquityResultsPtr<LatticeObject> const
            &calibrationResults,
        DeltaTime const &deltaTime, RiskFreeRate const &riskFreeRate,
        Payoff &&payoff, PayoffAdjuster &&payoffAdjuster,
        BarrierType barrierType,
        typename LatticeObject::Node_type const &barrier,
        typename LatticeObject::Node_type const &rebate,
        DiscountingStyle style) {

  // typedefs:
  typedef typename LatticeObject::Node_type Node;
  // typedef RiskFreeRateHolder:
  typedef RiskFreeRateHolder<RiskFreeRate> RFR;
  // typedef DeltaTimeHolder:
  typedef DeltaTimeHolder<DeltaTime> DT;
  // typedef discounting factor:
  typedef DiscountingFactor<Node> DCF;
  // typedef BarrierComparer:
  typedef BarrierComparer<Node> BC;

  // get correct comparer function:
  auto cmp = BC::comparer(barrierType);
  // get correct discounting factor style:
  auto dcf = DCF::function(style);

  // get the size of optionLattice object:
  std::size_t const treeSize = optionLattice.timeDimension();
  std::size_t const lastIdx = treeSize - 1;
  std::size_t const lastNodesSize = optionLattice.nodesAtIdx(lastIdx).size();

  // unpack the impled probabilities from calibrationResults:
  auto impliedProbs = calibrationResults->impliedProbabilities;
  LASSERT(impliedProbs.size() >= (treeSize - 1),
          "Must have enough implied volatilities for pricing.");

  // last payoff first:
  for (std::size_t i = 0; i < lastNodesSize; ++i) {
    if (cmp(spotLattice(lastIdx, i), barrier)) {
      optionLattice(lastIdx, i) = payoff(spotLattice(lastIdx, i));
    } else {
      optionLattice(lastIdx, i) = 0.0;
    }
  }

  std::size_t nodesSize{0};
  Node df{};
  Node value{};
  std::vector<std::tuple<Node, Node, Node>> probs;
  std::tuple<Node, Node, Node> tpl;
  for (auto t = lastIdx - 1; t > 0; --t) {
    df = dcf(RFR::rate(t, riskFreeRate), DT::deltaTime(t, deltaTime));
    nodesSize = optionLattice.nodesAtIdx(t).size();
    probs = impliedProbs.at(t);
    for (auto i = 0; i < nodesSize; ++i) {
      if (cmp(spotLattice(t, i), barrier)) {
        tpl = probs.at(i);
        value = df * (std::get<0>(tpl) * optionLattice(t + 1, i) +
                      std::get<1>(tpl) * optionLattice(t + 1, i + 1) +
                      std::get<2>(tpl) * optionLattice(t + 1, i + 2));
        payoffAdjuster(value, spotLattice(t, i));
        optionLattice(t, i) = value;

      } else {
        optionLattice(t, i) = rebate;
      }
    }
  }

  df = dcf(RFR::rate(0, riskFreeRate), DT::deltaTime(0, deltaTime));
  probs = impliedProbs.at(0);
  if (cmp(spotLattice(0, 0), barrier)) {
    tpl = probs.at(0);
    value = df * (std::get<0>(tpl) * optionLattice(1, 0) +
                  std::get<1>(tpl) * optionLattice(1, 1) +
                  std::get<2>(tpl) * optionLattice(1, 2));
    payoffAdjuster(value, spotLattice(0, 0));
    optionLattice(0, 0) = value;
  } else {
    optionLattice(0, 0) = rebate;
  }
}

template <typename DeltaTime, typename RiskFreeRate>
template <typename LatticeObject, typename Payoff>
void lattice_backward_traversals::ImpliedBackwardTraversal<
    lattice_types::LatticeType::Trinomial, DeltaTime, RiskFreeRate>::
    _backTraverseBarrierDKEAdjustment(
        LatticeObject &optionLattice, LatticeObject const &spotLattice,
        CalibratorTrinomialEquityResultsPtr<LatticeObject> const
            &calibrationResults,
        DeltaTime const &deltaTime, RiskFreeRate const &riskFreeRate,
        Payoff &&payoff, BarrierType barrierType,
        typename LatticeObject::Node_type const &barrier,
        typename LatticeObject::Node_type const &rebate,
        DiscountingStyle style) {

  // typedefs:
  typedef typename LatticeObject::Node_type Node;
  // typedef RiskFreeRateHolder:
  typedef RiskFreeRateHolder<RiskFreeRate> RFR;
  // typedef DeltaTimeHolder:
  typedef DeltaTimeHolder<DeltaTime> DT;
  // typedef discounting factor:
  typedef DiscountingFactor<Node> DCF;
  // typedef BarrierComparer:
  typedef BarrierComparer<Node> BC;
  // typedef DermanKaniErgenerAdjuster:
  typedef DermanKaniErgenerAdjuster<Node> DKEA;

  // get correct comparer function:
  auto cmp = BC::comparer(barrierType);
  // get correct adjuster:
  auto adjustPair = DKEA::adjuster(barrierType);
  // unpack checker and adjuster:
  auto dkeChecker = adjustPair.first;
  auto dkeAdjuster = adjustPair.second;
  // get correct discounting factor style:
  auto dcf = DCF::function(style);

  // get the size of optionLattice object:
  std::size_t const treeSize = optionLattice.timeDimension();
  std::size_t const lastIdx = treeSize - 1;
  std::size_t const lastNodesSize = optionLattice.nodesAtIdx(lastIdx).size();

  // unpack the impled probabilities from calibrationResults:
  auto impliedProbs = calibrationResults->impliedProbabilities;
  LASSERT(impliedProbs.size() >= (treeSize - 1),
          "Must have enough implied volatilities for pricing.");

  // last payoff first:
  for (std::size_t i = 0; i < lastNodesSize; ++i) {
    if (cmp(spotLattice(lastIdx, i), barrier)) {
      optionLattice(lastIdx, i) = payoff(spotLattice(lastIdx, i));
    } else {
      optionLattice(lastIdx, i) = 0.0;
    }
  }

  std::size_t nodesSize{0};
  Node df{};
  std::vector<std::tuple<Node, Node, Node>> probs;
  std::tuple<Node, Node, Node> tpl;
  for (auto t = lastIdx - 1; t > 0; --t) {
    df = dcf(RFR::rate(t, riskFreeRate), DT::deltaTime(t, deltaTime));
    nodesSize = optionLattice.nodesAtIdx(t).size();
    probs = impliedProbs.at(t);
    for (auto i = 0; i < nodesSize; ++i) {
      if (cmp(spotLattice(t, i), barrier)) {
        tpl = probs.at(i);
        optionLattice(t, i) =
            df * (std::get<0>(tpl) * optionLattice(t + 1, i) +
                  std::get<1>(tpl) * optionLattice(t + 1, i + 1) +
                  std::get<2>(tpl) * optionLattice(t + 1, i + 2));
      } else {
        optionLattice(t, i) = rebate;
      }
      // Derman-Kani-Ergener adjustment:
      if (((i > 0) && (i < nodesSize - 1)) &&
          dkeChecker(spotLattice(t, i), spotLattice(t, i - 1),
                     spotLattice(t, i + 1), barrier)) {
        optionLattice(t, i) = dkeAdjuster(
            spotLattice(t, i), spotLattice(t, i - 1), spotLattice(t, i + 1),
            barrier, rebate, optionLattice(t, i));
      }
    }
  }

  df = dcf(RFR::rate(0, riskFreeRate), DT::deltaTime(0, deltaTime));
  probs = impliedProbs.at(0);
  if (cmp(spotLattice(0, 0), barrier)) {
    tpl = probs.at(0);
    optionLattice(0, 0) = df * (std::get<0>(tpl) * optionLattice(1, 0) +
                                std::get<1>(tpl) * optionLattice(1, 1) +
                                std::get<2>(tpl) * optionLattice(1, 2));
  } else {
    optionLattice(0, 0) = rebate;
  }
}

template <typename DeltaTime, typename RiskFreeRate>
template <typename LatticeObject, typename Payoff, typename PayoffAdjuster>
void lattice_backward_traversals::ImpliedBackwardTraversal<
    lattice_types::LatticeType::Trinomial, DeltaTime, RiskFreeRate>::
    _backTraverseBarrierDKEAdjustment(
        LatticeObject &optionLattice, LatticeObject const &spotLattice,
        CalibratorTrinomialEquityResultsPtr<LatticeObject> const
            &calibrationResults,
        DeltaTime const &deltaTime, RiskFreeRate const &riskFreeRate,
        Payoff &&payoff, PayoffAdjuster &&payoffAdjuster,
        BarrierType barrierType,
        typename LatticeObject::Node_type const &barrier,
        typename LatticeObject::Node_type const &rebate,
        DiscountingStyle style) {

  // typedefs:
  typedef typename LatticeObject::Node_type Node;
  // typedef RiskFreeRateHolder:
  typedef RiskFreeRateHolder<RiskFreeRate> RFR;
  // typedef DeltaTimeHolder:
  typedef DeltaTimeHolder<DeltaTime> DT;
  // typedef discounting factor:
  typedef DiscountingFactor<Node> DCF;
  // typedef BarrierComparer:
  typedef BarrierComparer<Node> BC;
  // typedef DermanKaniErgenerAdjuster:
  typedef DermanKaniErgenerAdjuster<Node> DKEA;

  // get correct comparer function:
  auto cmp = BC::comparer(barrierType);
  // get correct adjuster:
  auto adjustPair = DKEA::adjuster(barrierType);
  // unpack checker and adjuster:
  auto dkeChecker = adjustPair.first;
  auto dkeAdjuster = adjustPair.second;
  // get correct discounting factor style:
  auto dcf = DCF::function(style);

  // get the size of optionLattice object:
  std::size_t const treeSize = optionLattice.timeDimension();
  std::size_t const lastIdx = treeSize - 1;
  std::size_t const lastNodesSize = optionLattice.nodesAtIdx(lastIdx).size();

  // unpack the impled probabilities from calibrationResults:
  auto impliedProbs = calibrationResults->impliedProbabilities;
  LASSERT(impliedProbs.size() >= (treeSize - 1),
          "Must have enough implied volatilities for pricing.");

  // last payoff first:
  for (std::size_t i = 0; i < lastNodesSize; ++i) {
    if (cmp(spotLattice(lastIdx, i), barrier)) {
      optionLattice(lastIdx, i) = payoff(spotLattice(lastIdx, i));
    } else {
      optionLattice(lastIdx, i) = 0.0;
    }
  }

  std::size_t nodesSize{0};
  Node df{};
  Node value{};
  std::vector<std::tuple<Node, Node, Node>> probs;
  std::tuple<Node, Node, Node> tpl;
  for (auto t = lastIdx - 1; t > 0; --t) {
    df = dcf(RFR::rate(t, riskFreeRate), DT::deltaTime(t, deltaTime));
    nodesSize = optionLattice.nodesAtIdx(t).size();
    probs = impliedProbs.at(t);
    for (auto i = 0; i < nodesSize; ++i) {
      if (cmp(spotLattice(t, i), barrier)) {
        tpl = probs.at(i);
        value = df * (std::get<0>(tpl) * optionLattice(t + 1, i) +
                      std::get<1>(tpl) * optionLattice(t + 1, i + 1) +
                      std::get<2>(tpl) * optionLattice(t + 1, i + 2));
        payoffAdjuster(value, spotLattice(t, i));
        optionLattice(t, i) = value;

      } else {
        optionLattice(t, i) = rebate;
      }
      // Derman-Kani-Ergener adjustment:
      if (((i > 0) && (i < nodesSize - 1)) &&
          dkeChecker(spotLattice(t, i), spotLattice(t, i - 1),
                     spotLattice(t, i + 1), barrier)) {
        optionLattice(t, i) = dkeAdjuster(
            spotLattice(t, i), spotLattice(t, i - 1), spotLattice(t, i + 1),
            barrier, rebate, optionLattice(t, i));
      }
    }
  }

  df = dcf(RFR::rate(0, riskFreeRate), DT::deltaTime(0, deltaTime));
  probs = impliedProbs.at(0);
  if (cmp(spotLattice(0, 0), barrier)) {
    tpl = probs.at(0);
    value = df * (std::get<0>(tpl) * optionLattice(1, 0) +
                  std::get<1>(tpl) * optionLattice(1, 1) +
                  std::get<2>(tpl) * optionLattice(1, 2));
    payoffAdjuster(value, spotLattice(0, 0));
    optionLattice(0, 0) = value;
  } else {
    optionLattice(0, 0) = rebate;
  }
}

} // namespace lattice

#endif ///_LATTICE_BACKWARD_TRAVERSALS