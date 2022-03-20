/**
 * @file lattice_binomial_model.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief  Base binomial models
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_BINOMIAL_MODEL_HPP_)
#define _LATTICE_BINOMIAL_MODEL_HPP_

#include "../utilities/lattice_enums.hpp"
#include "../utilities/lattice_leaf_generators.hpp"
#include "lattice_model_base.hpp"
#include <typeinfo>

namespace lattice {

template <typename T> class binomial_model<1, T> {
public:
  // Forward generator:
  virtual std::tuple<T, T> operator()(T value, T dt, std::size_t leaf_idx,
                                      std::size_t time_idx,
                                      bool is_mean_reverting) const = 0;

  // Backward generator:
  virtual T operator()(T curr_value, T up_value, T down_value, T dt) const = 0;

  // Factor count:
  enum { factor_count = 1 };

  // LatticeType:
  lattice_type type_of_lattice() const { return lattice_type::Binomial; }
};

template <typename T> class binomial_model<2, T> {
public:
  // Forward generators:
  virtual std::pair<forward_generator_t<T, T, T>, forward_generator_t<T, T, T>>
  forward_generator() const = 0;

  // Factor count:
  enum { factor_count = 2 };

  // LatticeType:
  lattice_type type_of_lattice() const { return lattice_type::Binomial; }
};

} // namespace lattice

#endif ///_LATTICE_BINOMIAL_MODEL_HPP_