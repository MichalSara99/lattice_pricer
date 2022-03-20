/**
 * @file lattice_multidim.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief Multidimensional extension of lattice
 * @version 0.1
 * @date 2022-03-20
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_MULTIDIM_HPP_)
#define _LATTICE_MULTIDIM_HPP_

#include "lattice.hpp"
#include "utilities/lattice_enums.hpp"
#include "utilities/lattice_macros.hpp"
#include <array>
#include <type_traits>

namespace lattice {

/**
 * @brief Multidimensional general lattice object
 *
 * @tparam dimension
 * @tparam type
 * @tparam node
 * @tparam time_axis
 * @tparam node_container_type
 * @tparam lattice_object
 * @tparam std::enable_if<(dimension > 1)
 */
template <std::size_t dimension, lattice_type type, typename node,
          typename time_axis, typename node_container_type,
          template <lattice_type, typename...> typename lattice_object,
          typename = typename std::enable_if<(dimension > 1)>::type>
class multidim_lattice {
public:
  virtual ~multidim_lattice() {}

  std::size_t constexpr factors() const { return dimension; }
};

/**
 * @brief Multidimensional mean-reverting general lattice object
 *
 * @tparam dimension
 * @tparam node
 * @tparam TimeAxis
 * @tparam node_container_type
 * @tparam lattice_object
 * @tparam std::enable_if<(dimension > 1)
 */
template <std::size_t dimension, typename node, typename time_axis,
          typename node_container_type,
          template <typename...> typename lattice_object,
          typename = typename std::enable_if<(dimension > 1)>::type>
class multidim_mean_reverting_lattice {
public:
  virtual ~multidim_mean_reverting_lattice() {}

  std::size_t constexpr factors() const { return dimension; }
};

/**
 * @brief Multidimensional indexed lattice (indexed by integer)
 *
 * @tparam dimension
 * @tparam type
 * @tparam node
 */
template <std::size_t dimension, lattice_type type, typename node>
class multidim_ilattice
    : public multidim_lattice<dimension, type, node, std::size_t,
                              std::vector<node>, ilattice> {
protected:
  std::vector<ilattice<type, node>> multi_tree_;

public:
  explicit multidim_ilattice(std::size_t number_periods) {
    for (std::size_t t = 0; t < dimension; ++t)
      multi_tree_.emplace_back(std::move(ilattice<type, node>{number_periods}));
  }

  ilattice<type, node> &get_factor(std::size_t factor_idx) {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx];
  }

  ilattice<type, node> const &get_factor(std::size_t factor_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx];
  }

  node const &operator()(std::size_t factor_idx, std::size_t time_idx,
                         std::size_t leaf_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx](time_idx, leaf_idx);
  }

  node &operator()(std::size_t factor_idx, std::size_t time_idx,
                   std::size_t leaf_idx) {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx](time_idx, leaf_idx);
  }

  node const &at(std::size_t factor_idx, std::size_t time_idx,
                 std::size_t leaf_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx](time_idx, leaf_idx);
  }

  node apex(std::size_t factor_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx].apex();
  }
};

/**
 * @brief Multidimensional glattice object
 *
 * @tparam dimension
 * @tparam type
 * @tparam node
 * @tparam time_axis
 */
template <std::size_t dimension, lattice_type type, typename node,
          typename time_axis>
class multidim_glattice
    : public multidim_lattice<dimension, type, node, time_axis,
                              std::vector<node>, glattice> {
protected:
  std::vector<glattice<type, node, time_axis>> multi_tree_;

public:
  explicit multidim_glattice(std::set<time_axis> const &fixing_dates_set) {
    for (std::size_t t = 0; t < dimension; ++t)
      multi_tree_.emplace_back(
          std::move(glattice<type, node, time_axis>{fixing_dates_set}));
  }

  explicit multidim_glattice(
      std::initializer_list<time_axis> const &fixing_dates) {
    for (std::size_t t = 0; t < dimension; ++t)
      multi_tree_.emplace_back(
          std::move(glattice<type, node, time_axis>{fixing_dates}));
  }

  glattice<type, node, time_axis> &get_factor(std::size_t factor_idx) {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx];
  }

  glattice<type, node, time_axis> const &
  get_factor(std::size_t factor_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx];
  }

  node const &operator()(std::size_t factor_idx, std::size_t time_idx,
                         std::size_t leaf_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx](time_idx, leaf_idx);
  }

  node &operator()(std::size_t factor_idx, std::size_t time_idx,
                   std::size_t leaf_idx) {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx](time_idx, leaf_idx);
  }

  node const &at(std::size_t factor_idx, std::size_t time_idx,
                 std::size_t leaf_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx](time_idx, leaf_idx);
  }

  node apex(std::size_t factor_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx].apex();
  }
};

/**
 * @brief Multidimensional mean-reverting indexed lattice object (indexed by
 * integer)
 *
 * @tparam dimension
 * @tparam node
 */
template <std::size_t dimension, typename node>
class multidim_mean_reverting_ilattice
    : public multidim_mean_reverting_lattice<dimension, node, std::size_t,
                                             std::vector<node>,
                                             mean_reverting_ilattice> {
protected:
  std::vector<mean_reverting_ilattice<node>> multi_tree_;

public:
  template <typename DeltaTime>
  explicit multidim_mean_reverting_ilattice(std::size_t number_periods,
                                            node reversion_speed,
                                            delta_time const &dt) {
    for (std::size_t t = 0; t < dimension; ++t)
      multi_tree_.emplace_back(std::move(
          mean_reverting_ilattice<node>{number_periods, reversion_speed, dt}));
  }

  mean_reverting_ilattice<Node> &get_factor(std::size_t factor_idx) {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx];
  }

  mean_reverting_ilattice<node> const &
  get_factor(std::size_t factor_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx];
  }

  node const &operator()(std::size_t factor_idx, std::size_t time_idx,
                         std::size_t leaf_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx](time_idx, leaf_idx);
  }

  node &operator()(std::size_t factor_idx, std::size_t time_idx,
                   std::size_t leaf_idx) {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx](time_idx, leaf_idx);
  }

  node const &at(std::size_t factor_idx, std::size_t time_idx,
                 std::size_t leaf_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx](time_idx, leaf_idx);
  }

  node apex(std::size_t factor_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx].apex();
  }
};

/**
 * @brief Multidimensional mean-reverting glattice object
 *
 * @tparam dimension
 * @tparam node
 * @tparam time_axis
 */
template <std::size_t dimension, typename node, typename time_axis>
class multidim_mean_reverting_glattice
    : public multidim_mean_reverting_lattice<dimension, node, time_axis,
                                             std::vector<node>,
                                             mean_reverting_glattice> {
protected:
  std::vector<mean_reverting_glattice<node, time_axis>> multi_tree_;

public:
  template <typename delta_time>
  explicit multidim_mean_reverting_glattice(
      std::set<time_axis> const &fixing_dates_set, node reversion_speed,
      delta_time const &dt) {
    for (std::size_t t = 0; t < dimension; ++t)
      multi_tree_.emplace_back(
          std::move(mean_reverting_glattice<node, time_axis>{
              fixing_dates_set, reversion_speed, dt}));
  }

  mean_reverting_glattice<node, time_axis> &get_factor(std::size_t factor_idx) {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx];
  }

  mean_reverting_glattice<node, time_axis> const &
  get_factor(std::size_t factor_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx];
  }

  node const &operator()(std::size_t factor_idx, std::size_t time_idx,
                         std::size_t leaf_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx](time_idx, leaf_idx);
  }

  node &operator()(std::size_t factor_idx, std::size_t time_idx,
                   std::size_t leaf_idx) {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx](time_idx, leaf_idx);
  }

  node const &at(std::size_t factor_idx, std::size_t time_idx,
                 std::size_t leaf_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx](time_idx, leaf_idx);
  }

  node apex(std::size_t factor_idx) const {
    LASSERT(factor_idx < dimension, "factor index is out of range");
    return multi_tree_[factor_idx].apex();
  }
};

} // namespace lattice

#endif ///_LATTICE_MULTIDIM