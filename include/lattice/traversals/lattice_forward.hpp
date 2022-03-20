/**
 * @file lattice_forward.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief  Forward lattice traversals
 * @version 0.1
 * @date 2022-03-20
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_FORWARD_HPP_)
#define _LATTICE_FORWARD_HPP_

#include "../utilities/lattice_component_holders.hpp"
#include "../utilities/lattice_enums.hpp"
#include "../utilities/lattice_macros.hpp"
#include <tuple>

namespace lattice {

template <lattice_type type, typename time_axis, typename delta_time,
          typename node>
struct traverse_forward {};

//	LatticeType::Binomial specialization
template <typename time_axis, typename delta_time, typename node>
struct traverse_forward<lattice_type::Binomial, time_axis, delta_time, node> {
private:
  //	This one is for compound DeltaTime object
  template <typename lattice_object, typename generator>
  static void _traverse(lattice_object &lattice, generator &&gen,
                        delta_time const &dt, node apex);

  //	This one is for compound DeltaTime object
  template <typename lattice_object, typename generator>
  static void _traverse(lattice_object &lattice, generator &&gen,
                        delta_time const &dt, node apex,
                        std::map<time_axis, node> const &dividend_data);

public:
  template <typename lattice_object, typename generator>
  static void traverse(lattice_object &lattice, generator &&gen,
                       delta_time const &dt, node apex) {
    _traverse(lattice, std::forward<generator>(gen), dt, apex);
  }

  template <typename lattice_object, typename generator>
  static void traverse(lattice_object &lattice, generator &&gen,
                       delta_time const &dt, node apex,
                       std::map<time_axis, node> const &dividend_data) {
    _traverse(lattice, std::forward<generator>(gen), dt, apex, dividend_data);
  }
};

//	LatticeType::TwoVariableBinomial specialization
template <typename time_axis, typename delta_time, typename node>
struct traverse_forward<lattice_type::TwoVariableBinomial, time_axis,
                        delta_time, node> {
private:
  //	This one is for compound DeltaTime object
  template <typename multidim_lattice_object, typename generator>
  static void _traverse(multidim_lattice_object &lattice, generator &&gen,
                        delta_time const &dt,
                        std::pair<node, node> const &apex);

  //	This one is for compound DeltaTime object
  template <typename multidim_lattice_object, typename generator>
  static void
  _traverse(multidim_lattice_object &lattice, generator &&gen,
            delta_time const &dt, std::pair<node, node> const &apex,
            std::pair<std::map<time_axis, node>,
                      std::map<time_axis, node>> const &dividend_data);

public:
  template <typename multidim_lattice_object, typename generator>
  static void traverse(multidim_lattice_object &lattice, generator &&gen,
                       delta_time const &dt,
                       std::pair<node, node> const &apex) {
    _traverse(lattice, std::forward<generator>(gen), dt, apex);
  }

  template <typename multidim_lattice_object, typename generator>
  static void
  traverse(multidim_lattice_object &lattice, generator &&gen,
           delta_time const &dt, std::pair<node, node> const &apex,
           std::pair<std::map<time_axis, node>, std::map<time_axis, node>> const
               &dividend_data) {
    _traverse(lattice, std::forward<generator>(gen), dt, apex, dividend_data);
  }
};

//	LatticeType::Trinomial specialization
template <typename time_axis, typename delta_time, typename node>
struct traverse_forward<lattice_type::Trinomial, time_axis, delta_time, node> {
private:
  //	This one is for compound DeltaTime object
  template <typename lattice_object, typename generator>
  static void _traverse_normal(std::size_t time_idx, lattice_object &lattice,
                               generator &&gen, delta_time const &dt,
                               node apex);

  //	This one is for compound DeltaTime object
  template <typename lattice_object, typename generator>
  static void _traverse_reverting(std::size_t time_idx, lattice_object &lattice,
                                  generator &&gen, delta_time const &dt,
                                  node apex);

  template <typename lattice_object, typename generator>
  static void _traverse(lattice_object &lattice, generator &&gen,
                        delta_time const &dt, node apex);

  //	This one is for compound DeltaTime object
  template <typename lattice_object, typename generator>
  static void _traverse_normal(std::size_t time_idx, lattice_object &lattice,
                               generator &&gen, delta_time const &dt, node apex,
                               std::map<time_axis, node> const &dividend_data);

  //	This one is for compound DeltaTime object
  template <typename lattice_object, typename generator>
  static void
  _traverse_reverting(std::size_t time_idx, lattice_object &lattice,
                      generator &&gen, delta_time const &dt, node apex,
                      std::map<time_axis, node> const &dividend_data);

  template <typename lattice_object, typename generator>
  static void _traverse(lattice_object &lattice, generator &&gen,
                        delta_time const &dt, node apex,
                        std::map<time_axis, node> const &dividend_data);

public:
  template <typename lattice_object, typename generator>
  static void traverse(lattice_object &lattice, generator &&gen,
                       delta_time const &dt, node apex) {
    _traverse(lattice, std::forward<generator>(gen), dt, apex);
  }

  template <typename lattice_object, typename generator>
  static void traverse(lattice_object &lattice, generator &&gen,
                       delta_time const &dt, node apex,
                       std::map<time_axis, node> const &dividend_data) {
    _traverse(lattice, std::forward<generator>(gen), dt, apex, dividend_data);
  }
};

// ============ IMPLEMENTATIONS ============

template <typename time_axis, typename delta_time, typename node>
template <typename multidim_lattice_object, typename generator>
void traverse_forward<lattice_type::TwoVariableBinomial, time_axis, delta_time,
                      node>::_traverse(multidim_lattice_object &lattice,
                                       generator &&gen, delta_time const &dt,
                                       std::pair<node, node> const &apex) {

  typedef delta_time_holder<delta_time> DT;
  node delta{};
  std::size_t nodes_size{};
  std::size_t const tree_size = lattice.get_factor(0).time_timension();
  auto gens = gen.forward_generator();
  auto generator0 = gens.first;
  auto generator1 = gens.second;

  lattice(0, 0, 0) = apex.first;
  lattice(1, 0, 0) = apex.second;
  std::tuple<node, node> tuple;

  for (std::size_t t = 1; t < tree_size; ++t) {
    delta = DT::delta_time(t - 1, dt);
    nodes_size = lattice.get_factor(0).nodes_at_idx(t - 1).size();
    for (std::size_t l = 0; l < nodes_size; ++l) {
      tuple = generator0(lattice(0, t - 1, l), delta, l, t, false);
      lattice(0, t, l) = std::get<0>(tuple);
      lattice(0, t, l + 1) = std::get<1>(tuple);
      tuple = generator1(lattice(1, t - 1, l), delta, l, t, false);
      lattice(1, t, l) = std::get<0>(tuple);
      lattice(1, t, l + 1) = std::get<1>(tuple);
    }
  }
}

template <typename time_axis, typename delta_time, typename node>
template <typename multidim_lattice_object, typename generator>
void traverse_forward<lattice_type::TwoVariableBinomial, time_axis, delta_time,
                      node>::
    _traverse(multidim_lattice_object &lattice, generator &&gen,
              delta_time const &dt, std::pair<node, node> const &apex,
              std::pair<std::map<time_axis, node>,
                        std::map<time_axis, node>> const &dividend_data) {

  typedef delta_time_holder<delta_time> DT;
  node delta{};
  std::size_t nodes_size{};
  std::size_t const tree_size = lattice.get_factor(0).time_dimension();
  auto gens = gen.forward_generator();
  auto generator0 = gens.first;
  auto generator1 = gens.second;

  lattice(0, 0, 0) = apex.first;
  lattice(1, 0, 0) = apex.second;
  std::tuple<node, node> tuple;

  for (std::size_t t = 1; t < tree_size; ++t) {
    delta = DT::delta_time(t - 1, dt);
    nodes_size = lattice.get_factor(0).nodes_at_idx(t - 1).size();
    for (std::size_t l = 0; l < nodes_size; ++l) {
      tuple = generator0(lattice(0, t - 1, l), delta, l, t, false);
      lattice(0, t, l) = std::get<0>(tuple);
      lattice(0, t, l + 1) = std::get<1>(tuple);
      tuple = generator1(lattice(1, t - 1, l), delta, l, t, false);
      lattice(1, t, l) = std::get<0>(tuple);
      lattice(1, t, l + 1) = std::get<1>(tuple);
    }
  }
  // Adjust the tree for discretely paying divdiends:
  auto div_data = dividend_data.first;
  auto curr = div_data.begin();
  auto last = div_data.end();
  if (curr == last)
    return;
  std::size_t first_ex_idx = lattice.get_factor(0).index_of(curr->first);
  if (first_ex_idx >= lattice.get_factor(0).time_dimension())
    return;

  node factor{1.0};
  for (std::size_t t = first_ex_idx; t < tree_size; ++t) {
    auto next_ex_itr = div_data.find(lattice.get_factor(0).time_at(t));
    factor = 1.0;
    if (next_ex_itr != last) {
      factor *= (factor - next_ex_itr->second);
    }
    delta = DT::delta_time(t - 1, dt);
    nodes_size = lattice.get_factor(0).nodes_at_idx(t - 1).size();
    for (std::size_t l = 0; l < nodes_size; ++l) {
      tuple = generator0(lattice(0, t - 1, l), delta, l, t, false);
      lattice(0, t, l) = factor * std::get<0>(tuple);
      lattice(0, t, l + 1) = factor * std::get<1>(tuple);
    }
  }

  // Adjust the tree for discretely paying divdiends:
  div_data = dividend_data.second;
  curr = div_data.begin();
  last = div_data.end();
  if (curr == last)
    return;
  first_ex_idx = lattice.get_factor(1).index_of(curr->first);
  if (first_ex_idx >= lattice.get_factor(1).time_dimension())
    return;

  factor = 1.0;
  for (std::size_t t = first_ex_idx; t < tree_size; ++t) {
    auto next_ex_itr = div_data.find(lattice.get_factor(1).time_at(t));
    factor = 1.0;
    if (next_ex_itr != last) {
      factor *= (factor - next_ex_itr->second);
    }
    delta = DT::delta_time(t - 1, dt);
    nodes_size = lattice.get_factor(1).nodes_at_idx(t - 1).size();
    for (std::size_t l = 0; l < nodes_size; ++l) {
      tuple = generator1(lattice(1, t - 1, l), delta, l, t, false);
      lattice(1, t, l) = factor * std::get<0>(tuple);
      lattice(1, t, l + 1) = factor * std::get<1>(tuple);
    }
  }
}

template <typename time_axis, typename delta_time, typename node>
template <typename lattice_object, typename generator>
void traverse_forward<lattice_type::Binomial, time_axis, delta_time,
                      node>::_traverse(lattice_object &lattice, generator &&gen,
                                       delta_time const &dt, node apex) {

  typedef delta_time_holder<delta_time> DT;
  node delta{};

  lattice(0, 0) = apex;
  std::tuple<node, node> tuple;
  for (std::size_t t = 1; t < lattice.time_dimension(); ++t) {
    delta = DT::delta_time(t - 1, dt);
    for (std::size_t l = 0; l < lattice.nodes_at_idx(t - 1).size(); ++l) {
      tuple = gen(lattice(t - 1, l), delta, l, t, false);
      lattice(t, l) = std::get<0>(tuple);
      lattice(t, l + 1) = std::get<1>(tuple);
    }
  }
}

template <typename time_axis, typename delta_time, typename node>
template <typename lattice_object, typename generator>
void traverse_forward<lattice_type::Binomial, time_axis, delta_time, node>::
    _traverse(lattice_object &lattice, generator &&gen, delta_time const &dt,
              node apex, std::map<time_axis, node> const &dividend_data) {

  typedef delta_time_holder<delta_time> DT;
  node delta{};

  lattice(0, 0) = apex;
  std::tuple<node, node> tuple;
  for (std::size_t t = 1; t < lattice.time_dimension(); ++t) {
    delta = DT::delta_time(t - 1, dt);
    for (std::size_t l = 0; l < lattice.nodes_at_idx(t - 1).size(); ++l) {
      tuple = gen(lattice(t - 1, l), delta, l, t, false);
      lattice(t, l) = std::get<0>(tuple);
      lattice(t, l + 1) = std::get<1>(tuple);
    }
  }
  // Adjust the tree for discretely paying divdiends:
  auto curr = dividend_data.begin();
  auto last = dividend_data.end();
  if (curr == last)
    return;
  std::size_t first_ex_idx = lattice.index_of(curr->first);
  if (first_ex_idx >= lattice.time_dimension())
    return;

  node factor{1.0};
  for (std::size_t t = first_ex_idx; t < lattice.time_dimension(); ++t) {
    auto next_ex_itr = dividend_data.find(lattice.time_at(t));
    factor = 1.0;
    if (next_ex_itr != last) {
      factor *= (factor - next_ex_itr->second);
    }
    delta = DT::delta_time(t - 1, dt);
    for (std::size_t l = 0; l < lattice.nodes_at_idx(t - 1).size(); ++l) {
      tuple = gen(lattice(t - 1, l), delta, l, t, false);
      lattice(t, l) = factor * std::get<0>(tuple);
      lattice(t, l + 1) = factor * std::get<1>(tuple);
    }
  }
}

template <typename time_axis, typename delta_time, typename node>
template <typename lattice_object, typename generator>
void traverse_forward<lattice_type::Trinomial, time_axis, delta_time,
                      node>::_traverse_normal(std::size_t time_idx,
                                              lattice_object &lattice,
                                              generator &&gen,
                                              delta_time const &dt, node apex) {

  typedef delta_time_holder<delta_time> DT;
  node delta{};

  lattice(0, 0) = apex;
  std::size_t nodes_size{0};
  std::tuple<node, node, node> tuple;
  for (std::size_t t = 1; t < time_idx; ++t) {
    nodes_size = lattice.nodes_at_idx(t - 1).size();
    delta = DT::delta_time(t - 1, dt);
    for (std::size_t l = 0; l < nodes_size; ++l) {
      tuple = gen(lattice(t - 1, l), delta, l, t, false);
      lattice(t, l) = std::get<0>(tuple);
      lattice(t, l + 1) = std::get<1>(tuple);
      lattice(t, l + 2) = std::get<2>(tuple);
    }
  }
}

template <typename time_axis, typename delta_time, typename node>
template <typename lattice_object, typename generator>
void traverse_forward<lattice_type::Trinomial, time_axis, delta_time,
                      node>::_traverse_reverting(std::size_t time_idx,
                                                 lattice_object &lattice,
                                                 generator &&gen,
                                                 delta_time const &dt,
                                                 node apex) {

  typedef delta_time_holder<delta_time> DT;
  node delta{};

  std::size_t const tree_size = lattice.time_dimension();
  std::size_t nodes_size{0};
  std::tuple<node, node, node> tuple;
  for (std::size_t t = time_idx; t < tree_size; ++t) {
    delta = DT::delta_time(t - 1, dt);
    nodes_size = lattice.nodes_at_idx(t - 1).size();
    // for lowest node:
    tuple = gen(lattice(t - 1, 0), delta, 0, t, true);
    lattice(t, 0) = std::get<0>(tuple); // low
    lattice(t, 1) = std::get<1>(tuple); // mid
    lattice(t, 2) = std::get<2>(tuple); // high

    for (std::size_t l = 1; l < nodes_size - 1; ++l) {
      tuple = gen(lattice(t - 1, l), delta, l, t, false);
      lattice(t, l - 1) = std::get<0>(tuple);
      lattice(t, l) = std::get<1>(tuple);
      lattice(t, l + 1) = std::get<2>(tuple);
    }
    // for highest node:
    tuple = gen(lattice(t - 1, nodes_size - 1), delta, nodes_size - 1, t, true);
    lattice(t, nodes_size - 1 - 2) = std::get<0>(tuple);
    lattice(t, nodes_size - 1 - 1) = std::get<1>(tuple);
    lattice(t, nodes_size - 1) = std::get<2>(tuple);
  }
}

template <typename time_axis, typename delta_time, typename node>
template <typename lattice_object, typename generator>
void traverse_forward<lattice_type::Trinomial, time_axis, delta_time,
                      node>::_traverse(lattice_object &lattice, generator &&gen,
                                       delta_time const &dt, node apex) {

  const std::size_t first_revert_idx = lattice.first_reverting_idx();
  const std::size_t tree_size = lattice.time_dimension();

  if (first_revert_idx == 0) {
    // This trinomial tree does not have reverting property:
    _traverse_normal(tree_size, lattice, std::forward<generator>(gen), dt,
                     apex);
  } else {
    // This trinomial tree does have reverting property:
    _traverse_normal(first_revert_idx, lattice, std::forward<generator>(gen),
                     dt, apex);
    _traverse_reverting(first_revert_idx, lattice, std::forward<generator>(gen),
                        dt, apex);
  }
}

template <typename time_axis, typename delta_time, typename node>
template <typename lattice_object, typename generator>
void traverse_forward<lattice_type::Trinomial, time_axis, delta_time, node>::
    _traverse_normal(std::size_t time_idx, lattice_object &lattice,
                     generator &&gen, delta_time const &dt, node apex,
                     std::map<time_axis, node> const &dividend_data) {

  typedef delta_time_holder<delta_time> DT;
  node delta{};
  std::size_t nodes_size{0};

  lattice(0, 0) = apex;
  std::tuple<node, node, node> tuple;
  for (std::size_t t = 1; t < time_idx; ++t) {
    delta = DT::delta_time(t - 1, dt);
    nodes_size = lattice.nodes_at_idx(t - 1).size();
    for (std::size_t l = 0; l < nodes_size; ++l) {
      tuple = gen(lattice(t - 1, l), delta, l, t, false);
      lattice(t, l) = std::get<0>(tuple);
      lattice(t, l + 1) = std::get<1>(tuple);
      lattice(t, l + 2) = std::get<2>(tuple);
    }
  }
  // Adjust the tree for discretely paying divdiends:
  auto curr = dividend_data.begin();
  auto last = dividend_data.end();
  if (curr == last)
    return;
  std::size_t first_ex_idx = lattice.index_of(curr->first);
  if (first_ex_idx >= time_idx)
    return;

  node factor{1.0};
  for (std::size_t t = first_ex_idx; t < time_idx; ++t) {
    auto next_ex_itr = dividend_data.find(lattice.time_at(t));
    factor = 1.0;
    if (next_ex_itr != last) {
      factor *= (factor - next_ex_itr->second);
    }
    nodes_size = lattice.nodes_at_idx(t - 1).size();
    delta = DT::delta_time(t - 1, dt);
    for (std::size_t l = 0; l < nodes_size; ++l) {
      tuple = gen(lattice(t - 1, l), delta, l, t, false);
      lattice(t, l) = factor * std::get<0>(tuple);
      lattice(t, l + 1) = factor * std::get<1>(tuple);
      lattice(t, l + 2) = factor * std::get<2>(tuple);
    }
  }
}

template <typename time_axis, typename delta_time, typename node>
template <typename lattice_object, typename generator>
void traverse_forward<lattice_type::Trinomial, time_axis, delta_time, node>::
    _traverse_reverting(std::size_t time_idx, lattice_object &lattice,
                        generator &&gen, delta_time const &dt, node apex,
                        std::map<time_axis, node> const &dividend_data) {

  typedef delta_time_holder<delta_time> DT;
  node delta{};
  std::size_t const tree_size = lattice.time_dimension();
  std::size_t nodes_size{0};

  lattice(0, 0) = apex;
  std::tuple<node, node, node> tuple;
  for (std::size_t t = time_idx; t < tree_size; ++t) {
    delta = DT::delta_time(t - 1, dt);
    nodes_size = lattice.nodes_at_idx(t - 1).size();
    // for lowest node:
    tuple = gen(lattice(t - 1, 0), delta, 0, t, true);
    lattice(t, 0) = std::get<0>(tuple); // low
    lattice(t, 1) = std::get<1>(tuple); // mid
    lattice(t, 2) = std::get<2>(tuple); // high

    for (std::size_t l = 1; l < nodes_size - 1; ++l) {
      tuple = gen(lattice(t - 1, l), delta, l, t, false);
      lattice(t, l - 1) = std::get<0>(tuple);
      lattice(t, l) = std::get<1>(tuple);
      lattice(t, l + 1) = std::get<2>(tuple);
    }
    // for highest node:
    tuple = gen(lattice(t - 1, nodes_size - 1), delta, nodes_size - 1, t, true);
    lattice(t, nodes_size - 1 - 2) = std::get<0>(tuple);
    lattice(t, nodes_size - 1 - 1) = std::get<1>(tuple);
    lattice(t, nodes_size - 1) = std::get<2>(tuple);
  }

  // Adjust the tree for discretely paying divdiends:
  auto curr = dividend_data.begin();
  auto last = dividend_data.end();
  if (curr == last)
    return;
  std::size_t first_ex_idx = lattice.index_of(curr->first);
  if ((first_ex_idx >= tree_size) || (first_ex_idx < time_idx))
    return;

  node factor{1.0};
  for (std::size_t t = first_ex_idx; t < tree_size; ++t) {
    auto next_ex_itr = dividend_data.find(lattice.time_at(t));
    factor = 1.0;
    if (next_ex_itr != last) {
      factor *= (factor - next_ex_itr->second);
    }
    nodes_size = lattice.nodes_at_idx(t - 1).size();
    delta = DT::delta_time(t - 1, dt);
    // for lowest node:
    tuple = gen(lattice(t - 1, 0), delta, 0, t, true);
    lattice(t, 0) = factor * std::get<0>(tuple); // low
    lattice(t, 1) = factor * std::get<1>(tuple); // mid
    lattice(t, 2) = factor * std::get<2>(tuple); // high

    for (std::size_t l = 1; l < nodes_size - 1; ++l) {
      tuple = gen(lattice(t - 1, l), delta, l, t, false);
      lattice(t, l - 1) = factor * std::get<0>(tuple);
      lattice(t, l) = factor * std::get<1>(tuple);
      lattice(t, l + 1) = factor * std::get<2>(tuple);
    }
    // for highest node:
    tuple = gen(lattice(t - 1, nodes_size - 1), delta, nodes_size - 1, t, true);
    lattice(t, nodes_size - 1 - 2) = factor * std::get<0>(tuple);
    lattice(t, nodes_size - 1 - 1) = factor * std::get<1>(tuple);
    lattice(t, nodes_size - 1) = factor * std::get<2>(tuple);
  }
}

template <typename time_axis, typename delta_time, typename node>
template <typename lattice_object, typename generator>
void traverse_forward<lattice_type::Trinomial, time_axis, delta_time, node>::
    _traverse(lattice_object &lattice, generator &&gen, delta_time const &dt,
              node apex, std::map<time_axis, node> const &dividend_data) {

  const std::size_t first_revert_idx = lattice.first_reverting_idx();
  const std::size_t tree_size = lattice.time_dimension();

  if (first_revert_idx == 0) {
    // This trinomial tree does not have reverting property:
    _traverse_normal(tree_size, lattice, std::forward<generator>(gen), dt, apex,
                     dividend_data);
  } else {
    // This trinomial tree does have reverting property:
    _traverse_normal(first_revert_idx, lattice, std::forward<generator>(gen),
                     dt, apex, dividend_data);
    _traverse_reverting(first_revert_idx, lattice, std::forward<generator>(gen),
                        dt, apex, dividend_data);
  }
}

} // namespace lattice

#endif ///_LATTICE_FORWARD_HPP_