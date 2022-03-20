/**
 * @file lattice.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief lattice structure
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_HPP_)
#define _LATTICE_HPP_

#include "utilities/lattice_enums.hpp"
#include "utilities/lattice_macros.hpp"
#include "utilities/lattice_miscellaneous.hpp"
#include <algorithm>
#include <initializer_list>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

namespace lattice {

/**
 * @brief General lattice structure object
 *
 * @tparam type lattice_enums::lattice_type
 * @tparam node
 * @tparam time_axis
 * @tparam node_container_t
 */
template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
class general_lattice {
private:
  std::size_t first_revert_idx_;

  template <typename arg>
  void build_tree_impl(arg const &argument, std::true_type);

  template <typename arg>
  void build_tree_impl(arg const &argument, std::false_type);

  template <typename arg, typename delta_time>
  void build_reverting_tree_impl(arg const &argument, node reversion_speed,
                                 delta_time const &dt, std::true_type);

  template <typename arg, typename delta_time>
  void build_reverting_tree_impl(arg const &argument, node reversion_speed,
                                 delta_time const &dt, std::false_type);

  std::size_t index_of_impl(time_axis time, std::true_type) const;

  std::size_t index_of_impl(time_axis time, std::false_type) const;

  node apex_impl(std::true_type) const;

  node apex_impl(std::false_type) const;

  node const &at_impl(time_axis time_idx, std::size_t leaf_idx,
                      std::true_type) const;

  node const &at_impl(time_axis time_idx, std::size_t leaf_idx,
                      std::false_type) const;

  node const &operator_const_impl(std::size_t time_idx, std::size_t leaf_idx,
                                  std::true_type) const;

  node const &operator_const_impl(std::size_t time_idx, std::size_t leaf_idx,
                                  std::false_type) const;

  node &operator_impl(std::size_t time_idx, std::size_t leaf_idx,
                      std::true_type);

  node &operator_impl(std::size_t time_idx, std::size_t leaf_idx,
                      std::false_type);

  node_container_t nodes_at_impl(time_axis time_idx, std::true_type) const;

  node_container_t nodes_at_impl(time_axis time_idx, std::false_type) const;

  node_container_t nodes_at_idx_impl(std::size_t idx, std::true_type) const;

  node_container_t nodes_at_idx_impl(std::size_t idx, std::false_type) const;

  time_axis const &time_at_impl(std::size_t idx, std::true_type) const;

  time_axis const &time_at_impl(std::size_t idx, std::false_type) const;

  template <typename delta_time>
  typename delta_time::value_type const
  get_min_delta_time_impl(delta_time const &dt, std::true_type) const {
    auto itr = std::min_element(dt.cbegin(), dt.cend());
    return *itr;
  }

  template <typename delta_time>
  delta_time const get_min_delta_time_impl(delta_time const &dt,
                                           std::false_type) const {
    return dt;
  }

public:
  typedef typename std::conditional<
      std::is_integral<time_axis>::value, std::vector<node_container_t>,
      std::map<time_axis, node_container_t>>::type tree_t;

  typedef typename tree_t::iterator iterator_t;
  typedef typename tree_t::const_iterator const_iterator_t;
  typedef node node_t;
  typedef time_axis time_axis_t;
  typedef node_container_t node_container_t_t;

protected:
  tree_t tree_;

  void _first_reverting_idx();

  template <typename delta_time>
  auto const get_min_delta_time(delta_time const &dt) const {
    return get_min_delta_time_impl(dt, std::is_compound<delta_time>());
  }

  template <typename arg> void build_tree(arg const &a);

  template <typename arg, typename delta_time>
  void build_reverting_tree(arg const &a, node reversion_speed,
                            delta_time const &dt);

  std::size_t number_nodes(std::size_t time_idx) const;

public:
  tree_t const &tree() const { return this->tree_; }

  static constexpr lattice_type type_of_lattice() { return type; }

  std::size_t constexpr factors() const { return 1; }

  constexpr std::size_t time_dimension() const {
    return std::distance(tree_.begin(), tree_.end());
  }

  node const &at(time_axis time_idx, std::size_t leaf_idx) const {
    return at_impl(time_idx, leaf_idx, is_map<tree_t>());
  }

  node &at(time_axis time_idx, std::size_t leaf_idx) {
    return (this->tree_[time_idx].at(leaf_idx));
  }

  node const &operator()(std::size_t time_idx, std::size_t leaf_idx) const {
    return operator_const_impl(time_idx, leaf_idx, is_map<tree_t>());
  }

  node &operator()(std::size_t time_idx, std::size_t leaf_idx) {
    return operator_impl(time_idx, leaf_idx, is_map<tree_t>());
  }

  node_container_t nodes_at(time_axis time_idx) const {
    return nodes_at_impl(time_idx, is_map<tree_t>());
  }

  node_container_t &nodes_at(time_axis time_idx) {
    return this->tree_[time_idx];
  }

  node_container_t nodes_at_idx(std::size_t idx) const {
    return nodes_at_idx_impl(idx, is_map<tree_t>());
  }

  std::size_t index_of(time_axis time) const {
    return index_of_impl(time, is_map<tree_t>());
  }

  time_axis const &time_at(std::size_t idx) const {
    return time_at_impl(idx, is_map<tree_t>());
  }

  node apex() const;

  const_iterator_t cbegin() const noexcept { return this->tree_.cbegin(); }

  const_iterator_t cend() const noexcept { return this->tree_.cend(); }

  iterator_t begin() { return this->tree_.begin(); }

  iterator_t end() { return this->tree_.end(); }

  std::size_t const first_reverting_idx() const {
    return this->first_revert_idx_;
  }

  bool is_first_reverting(std::size_t time_idx) {
    return (time_idx <= 0) ? false : time_idx == first_revert_idx_;
  }
};

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
void general_lattice<type, node, time_axis,
                     node_container_t>::_first_reverting_idx() {
  std::vector<int> states;
  const std::size_t tree_size = time_dimension();
  for (std::size_t t = 0; t < tree_size; ++t) {
    states.push_back(nodes_at_idx(t).size());
  }
  std::adjacent_difference(states.begin(), states.end(), states.begin());
  auto zeroItr = std::find_if(states.begin(), states.end(),
                              [](int val) { return (val == 0); });
  if (zeroItr == states.end()) {
    this->first_revert_idx_ = 0;
  } else {
    this->first_revert_idx_ = std::distance(states.begin(), zeroItr);
  }
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
node_container_t
general_lattice<type, node, time_axis, node_container_t>::nodes_at_idx_impl(
    std::size_t idx, std::true_type) const {
  LASSERT(idx >= 0, "Index must be nonegative");
  auto first = this->tree_.cbegin();
  if (idx == 0) {
    return first->second;
  } else {
    auto itr = std::next(first, idx);
    return itr->second;
  }
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
node_container_t
general_lattice<type, node, time_axis, node_container_t>::nodes_at_idx_impl(
    std::size_t idx, std::false_type) const {
  LASSERT(idx >= 0, "Index must be nonegative");
  return this->tree_[idx];
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
time_axis const &
general_lattice<type, node, time_axis, node_container_t>::time_at_impl(
    std::size_t idx, std::true_type) const {
  LASSERT(idx >= 0, "Index must be nonegative");
  auto first = this->tree_.cbegin();
  if (idx == 0) {
    return first->first;
  } else {
    auto itr = std::next(first, idx);
    return itr->first;
  }
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
time_axis const &
general_lattice<type, node, time_axis, node_container_t>::time_at_impl(
    std::size_t idx, std::false_type) const {
  return idx;
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
std::size_t
general_lattice<type, node, time_axis, node_container_t>::index_of_impl(
    time_axis time_idx, std::true_type) const {
  typename tree_t::const_iterator citer(this->tree_.find(time_idx));
  return ((citer != this->tree_.cend())
              ? (std::distance(this->tree_.cbegin(), citer))
              : throw std::out_of_range("Error: time_idx out of range.\n"));
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
std::size_t
general_lattice<type, node, time_axis, node_container_t>::index_of_impl(
    time_axis time_idx, std::false_type) const {
  return time_idx;
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
node const &general_lattice<type, node, time_axis, node_container_t>::at_impl(
    time_axis time_idx, std::size_t leaf_idx, std::true_type) const {
  typename tree_t::const_iterator citer(this->tree_.find(time_idx));
  return ((citer != this->tree_.end())
              ? (citer->second.at(leaf_idx))
              : throw std::out_of_range("Error: time_idx out of range.\n"));
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
node const &general_lattice<type, node, time_axis, node_container_t>::at_impl(
    time_axis time_idx, std::size_t leaf_idx, std::false_type) const {
  return (this->tree_[time_idx].at(leaf_idx));
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
node &general_lattice<type, node, time_axis, node_container_t>::operator_impl(
    std::size_t time_idx, std::size_t leaf_idx, std::true_type) {
  LASSERT(time_idx >= 0, "Index must be nonegative");
  auto first = this->tree_.begin();
  if (time_idx == 0) {
    return first->second.at(leaf_idx);
  } else {
    auto itr = std::next(first, time_idx);
    return itr->second.at(leaf_idx);
  }
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
node &general_lattice<type, node, time_axis, node_container_t>::operator_impl(
    std::size_t time_idx, std::size_t leaf_idx, std::false_type) {
  LASSERT(time_idx >= 0, "Index must be nonegative");
  auto first = this->tree_.begin();
  if (time_idx == 0) {
    return first->at(leaf_idx);
  } else {
    auto itr = std::next(first, time_idx);
    return itr->at(leaf_idx);
  }
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
node const &
general_lattice<type, node, time_axis, node_container_t>::operator_const_impl(
    std::size_t time_idx, std::size_t leaf_idx, std::true_type) const {
  LASSERT(time_idx >= 0, "Index must be nonegative");
  auto first = this->tree_.cbegin();
  if (time_idx == 0) {
    return first->second.at(leaf_idx);
  } else {
    auto itr = std::next(first, time_idx);
    return itr->second.at(leaf_idx);
  }
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
node const &
general_lattice<type, node, time_axis, node_container_t>::operator_const_impl(
    std::size_t time_idx, std::size_t leaf_idx, std::false_type) const {
  LASSERT(time_idx >= 0, "Index must be nonegative");
  auto first = this->tree_.cbegin();
  if (time_idx == 0) {
    return first->at(leaf_idx);
  } else {
    auto itr = std::next(first, time_idx);
    return itr->at(leaf_idx);
  }
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
node_container_t
general_lattice<type, node, time_axis, node_container_t>::nodes_at_impl(
    time_axis time_idx, std::true_type) const {
  typename tree_t::const_iterator citer(this->tree_.find(time_idx));
  return ((citer != this->tree_.end())
              ? (citer->second)
              : throw std::out_of_range("Error: time_idx out of range.\n"));
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
node_container_t
general_lattice<type, node, time_axis, node_container_t>::nodes_at_impl(
    time_axis time_idx, std::false_type) const {
  return (this->tree_[time_idx]);
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
std::size_t
general_lattice<type, node, time_axis, node_container_t>::number_nodes(
    std::size_t time_idx) const {
  if (type == lattice_type::TwoVariableBinomial) {
    return ((time_idx + 1) * (time_idx + 1));
  }
  auto intRep = static_cast<std::underlying_type<lattice_type>::type>(type);
  return ((intRep + 1) * (time_idx + 1) - (intRep));
};

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
template <typename arg>
void general_lattice<type, node, time_axis, node_container_t>::build_tree_impl(
    arg const &a, std::true_type) {
  tree_.clear();
  for (std::size_t t = 0; t < a.size(); ++t) {
    tree_[a[t]] = std::move(node_container_t(number_nodes(t)));
  }
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
template <typename arg>
void general_lattice<type, node, time_axis, node_container_t>::build_tree_impl(
    arg const &a, std::false_type) {
  tree_.clear();
  tree_.reserve(a);
  for (std::size_t t = 0; t <= a; ++t) {
    tree_.emplace_back(node_container_t(number_nodes(t)));
  }
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
template <typename arg>
void general_lattice<type, node, time_axis, node_container_t>::build_tree(
    arg const &a) {
  build_tree_impl(a, std::is_compound<arg>());
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
template <typename arg, typename delta_time>
void general_lattice<type, node, time_axis, node_container_t>::
    build_reverting_tree_impl(arg const &a, node reversion_speed,
                              delta_time const &dt, std::true_type) {
  auto const d = get_min_delta_time(dt);
  std::size_t const maxStatesUp =
      static_cast<std::size_t>(
          std::floor(static_cast<node>(1.0) /
                     (static_cast<node>(2.0) * reversion_speed * d))) +
      1;
  auto const intRep = static_cast<std::underlying_type<lattice_type>::type>(
      lattice_type::Trinomial);

  auto number_nodes = [&](std::size_t time_idx) {
    std::size_t normalnodesSize = ((intRep + 1) * (time_idx + 1) - (intRep));
    return std::min(normalnodesSize, 2 * maxStatesUp + 1);
  };
  for (std::size_t t = 0; t < a.size(); ++t) {
    tree_[a[t]] = std::move(node_container_t(number_nodes(t)));
  }
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
template <typename arg, typename delta_time>
void general_lattice<type, node, time_axis, node_container_t>::
    build_reverting_tree_impl(arg const &a, node reversion_speed,
                              delta_time const &dt, std::false_type) {

  auto const d = get_min_delta_time(dt);
  std::size_t const maxStatesUp =
      static_cast<std::size_t>(
          std::floor(static_cast<node>(1) /
                     (static_cast<node>(2.0) * reversion_speed * d))) +
      1;
  auto const intRep = static_cast<std::underlying_type<lattice_type>::type>(
      lattice_type::Trinomial);

  auto number_nodes = [&](std::size_t time_idx) {
    std::size_t normalnodesSize = ((intRep + 1) * (time_idx + 1) - (intRep));
    return std::min(normalnodesSize, 2 * maxStatesUp + 1);
  };
  tree_.reserve(a);
  for (std::size_t t = 0; t <= a; ++t) {
    tree_.emplace_back(node_container_t(number_nodes(t)));
  }
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
template <typename arg, typename delta_time>
void general_lattice<type, node, time_axis, node_container_t>::
    build_reverting_tree(arg const &a, node reversion_speed,
                         delta_time const &dt) {
  build_reverting_tree_impl(a, reversion_speed, dt, std::is_compound<arg>());
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
node general_lattice<type, node, time_axis, node_container_t>::apex_impl(
    std::false_type) const {
  return *((*this->tree_.begin()).second.begin());
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
node general_lattice<type, node, time_axis, node_container_t>::apex_impl(
    std::true_type) const {
  return *((*this->tree_.begin()).begin());
}

template <lattice_type type, typename node, typename time_axis,
          typename node_container_t>
node general_lattice<type, node, time_axis, node_container_t>::apex() const {
  return apex_impl(std::is_integral<time_axis>());
}

/**
 * @brief Indexed lattice object (time is indexed by integers)
 *
 * @tparam type
 * @tparam node
 */
template <lattice_type type, typename node>
class ilattice
    : public general_lattice<type, node, std::size_t, std::vector<node>> {
private:
  std::size_t min_index_{0};
  std::size_t max_index_;

public:
  explicit ilattice(std::size_t number_periods) : max_index_{number_periods} {
    this->build_tree(number_periods);
  }

  virtual ~ilattice() {}

  ilattice(ilattice<type, node> const &cpy)
      : min_index_{cpy.min_index_}, max_index_{cpy.max_index_} {
    this->tree_ = cpy.tree_;
  }

  ilattice(ilattice<type, node> &&other) noexcept
      : min_index_{std::move(other.min_index_)}, max_index_{std::move(
                                                     other.max_index_)} {
    this->tree_ = std::move(other.tree_);
  }

  ilattice &operator=(ilattice<type, node> const &cpy) {
    if (this != &cpy) {
      min_index_ = cpy.min_index_;
      max_index_ = cpy.max_index_;
      this->tree_ = cpy.tree_;
    }
    return *this;
  }

  ilattice &operator=(ilattice<type, node> &&other) noexcept {
    if (this != &other) {
      min_index_ = std::move(other.min_index_);
      max_index_ = std::move(other.max_index_);
      this->tree_ = std::move(other.tree_);
    }
    return *this;
  }

  static lattice_class const class_of_lattice() {
    return lattice_class::Normal;
  }

  std::size_t min_index() const { return min_index_; }
  std::size_t max_index() const { return max_index_; }
};

/**
 * @brief Lattice object with any time axis
 *
 * @tparam type
 * @tparam node
 * @tparam time_axis
 */
template <lattice_type type, typename node, typename time_axis>
class glattice
    : public general_lattice<type, node, time_axis, std::vector<node>> {
private:
  std::set<time_axis> fixing_dates_set_;
  std::vector<time_axis> fixing_dates_;

public:
  explicit glattice(std::set<time_axis> const &fixing_dates_set)
      : fixing_dates_set_{fixing_dates_set} {
    fixing_dates_.reserve(fixing_dates_set_.size());
    for (auto e : fixing_dates_set_) {
      fixing_dates_.emplace_back(std::move(e));
    }
    this->build_tree(fixing_dates_);
  }

  glattice(std::initializer_list<time_axis> const &fixing_dates)
      : fixing_dates_set_{fixing_dates} {
    fixing_dates_.reserve(fixing_dates_set_.size());
    for (auto e : fixing_dates_set_) {
      fixing_dates_.emplace_back(std::move(e));
    }
    this->build_tree(fixing_dates_);
  }

  virtual ~glattice() {}

  glattice(glattice<type, node, time_axis> const &cpy)
      : fixing_dates_{cpy.fixing_dates_} {
    this->tree_ = cpy.tree_;
  }

  glattice(glattice<type, node, time_axis> &&other) noexcept
      : fixing_dates_{std::move(other.fixing_dates_)} {
    this->tree_ = std::move(other.tree_);
  }

  glattice &operator=(glattice<type, node, time_axis> const &cpy) {
    if (this != &cpy) {
      fixing_dates_ = cpy.fixing_dates_;
      this->tree = cpy.tree_;
    }
    return *this;
  }

  glattice &operator=(glattice<type, node, time_axis> &&other) noexcept {
    if (this != &other) {
      fixing_dates_ = std::move(other.fixing_dates_);
      this->tree = std::move(other.tree_);
    }
    return *this;
  }

  static lattice_class const class_of_lattice() {
    return lattice_class::Normal;
  }

  std::vector<time_axis> fixing_dates() const { return fixing_dates_; }
};

/**
 * @brief Mean-reverting indexed lattice for trinomial trees
 *
 * @tparam node
 */
template <typename node>
class mean_reverting_ilattice
    : public general_lattice<lattice_type::Trinomial, node, std::size_t,
                             std::vector<node>> {
private:
  std::size_t min_index_{0};
  std::size_t max_index_;

public:
  template <typename delta_time>
  explicit mean_reverting_ilattice(std::size_t number_periods,
                                   node reversion_speed, delta_time const &dt)
      : max_index_{number_periods} {
    this->build_reverting_tree(number_periods, reversion_speed, dt);
    this->_first_reverting_idx();
  }

  virtual ~mean_reverting_ilattice() {}

  mean_reverting_ilattice(mean_reverting_ilattice<node> const &cpy)
      : min_index_{cpy.min_index_}, max_index_{cpy.max_index_} {
    this->tree_ = cpy.tree_;
  }

  mean_reverting_ilattice(mean_reverting_ilattice<node> &&other) noexcept
      : min_index_{std::move(other.min_index_)}, max_index_{std::move(
                                                     other.max_index_)} {
    this->tree_ = std::move(other.tree_);
  }

  mean_reverting_ilattice &operator=(mean_reverting_ilattice<node> const &cpy) {
    if (this != &cpy) {
      min_index_ = cpy.min_index_;
      max_index_ = cpy.max_index_;
      this->tree_ = cpy.tree_;
    }
    return *this;
  }

  mean_reverting_ilattice &
  operator=(mean_reverting_ilattice<node> &&other) noexcept {
    if (this != &other) {
      min_index_ = std::move(other.min_index_);
      max_index_ = std::move(other.max_index_);
      this->tree_ = std::move(other.tree_);
    }
    return *this;
  }

  static lattice_class const class_of_lattice() {
    return lattice_class::MeanReverting;
  }

  std::size_t min_index() const { return min_index_; }
  std::size_t max_index() const { return max_index_; }
};

/**
 * @brief Mean-reverting lattice with any time axis
 *
 * @tparam node
 * @tparam time_axis
 */
template <typename node, typename time_axis>
class mean_reverting_glattice
    : public general_lattice<lattice_type::Trinomial, node, time_axis,
                             std::vector<node>> {

private:
  std::set<time_axis> fixing_dates_set_;
  std::vector<time_axis> fixing_dates_;

public:
  template <typename delta_time>
  explicit mean_reverting_glattice(std::set<time_axis> const &fixing_dates_set,
                                   node reversion_speed, delta_time const &dt)
      : fixing_dates_set_{fixing_dates_set} {
    fixing_dates_.reserve(fixing_dates_set_.size());
    for (auto e : fixing_dates_set_) {
      fixing_dates_.emplace_back(std::move(e));
    }
    this->build_reverting_tree(fixing_dates_, reversion_speed, dt);
    this->_first_reverting_idx();
  }

  virtual ~mean_reverting_glattice() {}

  mean_reverting_glattice(mean_reverting_glattice<node, time_axis> const &cpy)
      : fixing_dates_{cpy.fixing_dates_} {
    this->tree_ = cpy.tree_;
  }

  mean_reverting_glattice(
      mean_reverting_glattice<node, time_axis> &&other) noexcept
      : fixing_dates_{std::move(other.fixing_dates_)} {
    this->tree_ = std::move(other.tree_);
  }

  mean_reverting_glattice &
  operator=(mean_reverting_glattice<node, time_axis> const &cpy) {
    if (this != &cpy) {
      fixing_dates_ = cpy.fixing_dates_;
      this->tree = cpy.tree_;
    }
    return *this;
  }

  mean_reverting_glattice &
  operator=(mean_reverting_glattice<node, time_axis> &&other) noexcept {
    if (this != &other) {
      fixing_dates_ = std::move(other.fixing_dates_);
      this->tree = std::move(other.tree_);
    }
    return *this;
  }

  static lattice_class const class_of_lattice() {
    return lattice_class::MeanReverting;
  }

  std::vector<time_axis> fixing_dates() const { return fixing_dates_; }
};

} // namespace lattice

#endif //_LATTICE_HPP_