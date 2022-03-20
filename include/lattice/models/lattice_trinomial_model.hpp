/**
 * @file lattice_trinomial_model.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief  Base trinomial trees
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_TRINOMIAL_MODEL_HPP_)
#define _LATTICE_TRINOMIAL_MODEL_HPP_

#include "../utilities/lattice_enums.hpp"
#include "../utilities/lattice_leaf_generators.hpp"
#include "lattice_model_base.hpp"
#include <typeinfo>

namespace lattice {

template <typename T> class trinomial_model<1, T> {
public:
  // Forward generator:
  virtual std::tuple<T, T, T> operator()(T value, T dt, std::size_t leaf_idx,
                                         std::size_t time_idx,
                                         bool is_mean_reverting) const = 0;

  // Backward generator:
  virtual T operator()(T curr_value, T up_value, T mid_value, T down_value,
                       T dt, std::size_t revert_branches_size,
                       std::size_t nodes_size, std::size_t leaf_idx) const = 0;

  // Factor count:
  enum { factor_count = 1 };

  // LatticeType:
  lattice_type type_of_lattice() const { return lattice_type::Trinomial; }
};

template <typename T> class trinomial_model<2, T> {
public:
  // Forward generators:
  virtual std::pair<forward_generator_t<T, T, T, T>,
                    forward_generator_t<T, T, T, T>>
  forward_generator() const = 0;

  // Forward generator 2:
  virtual std::pair<backward_generator_t<T, T, T, T, T, T>,
                    backward_generator_t<T, T, T, T, T, T>>
  backward_generator() const = 0;

  // Factor count:
  enum { factor_count = 2 };

  // LatticeType:
  lattice_type type_of_lattice() const { return lattice_type::Trinomial; }
};
} // namespace lattice
#endif //_LATTICE_TRINOMIAL_MODEL_HPP_