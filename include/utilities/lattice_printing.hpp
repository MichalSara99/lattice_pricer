/**
 * @file lattice_printing.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief Printing facilities
 * @version 0.1
 * @date 2022-03-18
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#if !defined(_LATTICE_PRINTING_HPP_)
#define _LATTICE_PRINTING_HPP_

#include "lattice_enums.hpp"
#include "lattice_traits.hpp"
#include <iostream>

namespace lattice {

template <typename lattice_object, typename citer,
          typename traits = print_traits<typename lattice_object::node_t,
                                         lattice_object::type_of_lattice()>>
void _print_impl(lattice_object const &lattice, citer cbegin, citer cend,
                 std::ostream &out, std::true_type) {

  for (auto itr = cbegin; itr != cend; ++itr) {
    out << "[" << traits::print_line(*itr);
    out << "]\n";
  }
}

template <typename lattice_object, typename citer,
          typename traits = print_traits<typename lattice_object::node_t,
                                         lattice_object::type_of_lattice()>>
void _print_impl(lattice_object const &lattice, citer cbegin, citer cend,
                 std::ostream &out, std::false_type) {
  std::string bl{""};
  if (lattice_object::type_of_lattice() == lattice_type::TwoVariableBinomial)
    bl = "\n";
  for (auto itr = cbegin; itr != cend; ++itr) {
    out << "(" << (*itr).first << "):" << bl << "["
        << traits::print_line((*itr).second);
    out << "]\n";
  }
}

/**
 * @brief Global function for printing lattices
 *
 * @tparam lattice_object
 * @tparam citer
 * @tparam print_traits<typename lattice_object::node_type,
 * lattice_object::type()>
 * @param lattice lattice object
 * @param cbegin constant begin iterator
 * @param cend constant end iterator
 * @param out output stream
 */
template <typename lattice_object, typename citer,
          typename traits = print_traits<typename lattice_object::node_t,
                                         lattice_object::type_of_lattice()>>
void print(lattice_object const &lattice, citer cbegin, citer cend,
           std::ostream &out = std::cout) {
  _print_impl(lattice, cbegin, cend, out,
              std::is_integral<typename lattice_object::time_axis_t>());
}

/**
 * @brief Global function for printing lattices
 *
 * @tparam lattice_object
 * @tparam print_traits<typename lattice_object::node_t,
 * lattice_object::type_of_lattice()>
 * @param lattice lattice object
 * @param out output stream
 */
template <typename lattice_object,
          typename traits = print_traits<typename lattice_object::node_t,
                                         lattice_object::type_of_lattice()>>
void print(lattice_object const &lattice, std::ostream &out = std::cout) {
  _print_impl(lattice, lattice.cbegin(), lattice.cend(), out,
              std::is_integral<typename lattice_object::time_axis_t>());
}

} // namespace lattice

#endif ///_LATTICE_PRINTING_HPP_