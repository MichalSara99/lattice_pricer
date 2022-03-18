/**
 * @file lattice_leaf_generators.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief Lattice leaf generators
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_LEAF_GENERATORS_HPP_)
#define _LATTICE_LEAF_GENERATORS_HPP_

#include <functional>
#include <tuple>

namespace lattice {

template <typename Node, typename... Nodes>
using forward_generator_t = std::function<std::tuple<Nodes...>(
    Node, Node, std::size_t, std::size_t, bool)>;

template <typename Node, typename... Nodes>
using backward_generator_t = std::function<Node(Nodes...)>;

} // namespace lattice

#endif ///_LATTICE_LEAF_GENERATORS_HPP_
