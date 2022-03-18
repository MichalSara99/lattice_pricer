/**
 * @file lattice_traits.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief Traits for merging and printing
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_TRAITS_HPP_)
#define _LATTICE_TRAITS_HPP_

#include "lattice_enums.hpp"
#include "lattice_miscellaneous.hpp"
#include <array>
#include <iostream>
#include <sstream>
#include <tuple>

namespace lattice {

/**
 * @brief Traits for merging nodes on the tree
 *
 * @tparam n
 * @tparam node
 */
template <std::size_t n, typename node> struct merge_traits {
  using node_holder = std::array<node, n>;

  template <typename T, typename... Ts>
  static std::array<node, n> holder(T t, Ts... ts) {
    return std::array<node, n>{t, ts...};
  }
};

template <typename node> struct merge_traits<1, node> {
  using node_holder = std::tuple<node, node>;

  template <typename T, typename... Ts>
  static std::tuple<node, node> holder(T t, Ts... ts) {
    return std::make_tuple(t, ts...);
  }
};
template <typename node> struct merge_traits<2, node> {
  using node_holder = std::tuple<node, node, node>;

  template <typename T, typename... Ts>
  static std::tuple<node, node, node> holder(T t, Ts... ts) {
    return std::make_tuple(t, ts...);
  }
};
template <typename node> struct merge_traits<3, node> {
  using node_holder = std::tuple<node, node, node, node>;

  template <typename T, typename... Ts>
  static std::tuple<node, node, node, node> holder(T t, Ts... ts) {
    return std::make_tuple(t, ts...);
  }
};

/**
 * @brief
 *
 * @tparam node
 * @tparam type
 */

template <typename node, lattice_type type> struct print_traits {
  static const std::size_t size = 1;

  template <typename container> static std::string print_line(container &cont) {
    std::ostringstream ss;
    for (auto e : cont) {
      ss << e << ", ";
    }
    return clip_comma(std::move(ss.str()));
  }
};

template <typename node>
struct print_traits<node, lattice_type::TwoVariableBinomial> {
  static const std::size_t size = 1;

  template <typename container> static std::string print_line(container &cont) {
    std::ostringstream ss;
    std::size_t const side = static_cast<std::size_t>(std::sqrt(cont.size()));
    for (std::size_t i = 0; i < cont.size(); ++i) {
      if ((i > 0) && (i % side == 0))
        ss << "\n ";
      ss << cont[i] << ", ";
    }
    return clip_comma(std::move(ss.str()));
  }
};

template <typename node, lattice_type Type>
struct print_traits<std::tuple<node, node>, Type> {
  static const std::size_t size = 2;

  template <typename container> static std::string print_line(container &cont) {
    std::ostringstream ss;
    for (auto e : cont) {
      ss << "(" << std::get<0>(e) << "," << std::get<1>(e) << "), ";
    }
    return clip_comma(std::move(ss.str()));
  }
};

template <typename node>
struct print_traits<std::tuple<node, node>, lattice_type::TwoVariableBinomial> {
  static const std::size_t size = 1;

  template <typename container> static std::string print_line(container &cont) {
    std::ostringstream ss;
    std::size_t const side = static_cast<std::size_t>(std::sqrt(cont.size()));
    std::tuple<node, node> val{};
    for (std::size_t i = 0; i < cont.size(); ++i) {
      if ((i > 0) && (i % side == 0))
        ss << "\n ";
      val = cont[i];
      ss << "(" << std::get<0>(val) << "," << std::get<1>(val) << "), ";
    }
    return clip_comma(std::move(ss.str()));
  }
};

template <typename node, lattice_type Type>
struct print_traits<std::tuple<node, node, node>, Type> {
  static const std::size_t size = 3;

  template <typename container> static std::string print_line(container &cont) {
    std::ostringstream ss;
    for (auto e : cont) {
      ss << "(" << std::get<0>(e) << "," << std::get<1>(e) << ","
         << std::get<2>(e) << "), ";
    }
    return clip_comma(std::move(ss.str()));
  }
};

template <typename node>
struct print_traits<std::tuple<node, node, node>,
                    lattice_type::TwoVariableBinomial> {
  static const std::size_t size = 1;

  template <typename container> static std::string print_line(container &cont) {
    std::ostringstream ss;
    std::size_t const side = static_cast<std::size_t>(std::sqrt(cont.size()));
    std::tuple<node, node, node> val{};
    for (std::size_t i = 0; i < cont.size(); ++i) {
      if ((i > 0) && (i % side == 0))
        ss << "\n ";
      val = cont[i];
      ss << "(" << std::get<0>(val) << "," << std::get<1>(val) << ","
         << std::get<2>(val) << "), ";
    }
    return clip_comma(std::move(ss.str()));
  }
};

template <typename node, lattice_type Type>
struct print_traits<std::tuple<node, node, node, node>, Type> {
  static const std::size_t size = 4;

  template <typename container> static std::string print_line(container &cont) {
    std::ostringstream ss;
    for (auto e : cont) {
      ss << "(" << std::get<0>(e) << "," << std::get<1>(e) << ","
         << std::get<2>(e) << "," << std::get<3>(e) << "), ";
    }
    return clip_comma(std::move(ss.str()));
  }
};

template <typename node>
struct print_traits<std::tuple<node, node, node, node>,
                    lattice_type::TwoVariableBinomial> {
  static const std::size_t size = 1;

  template <typename container> static std::string print_line(container &cont) {
    std::ostringstream ss;
    std::size_t const side = static_cast<std::size_t>(std::sqrt(cont.size()));
    std::tuple<node, node, node, node> val{};
    for (std::size_t i = 0; i < cont.size(); ++i) {
      if ((i > 0) && (i % side == 0))
        ss << "\n ";
      val = cont[i];
      ss << "(" << std::get<0>(val) << "," << std::get<1>(val) << ","
         << std::get<2>(val) << "," << std::get<3>(val) << "), ";
    }
    return clip_comma(std::move(ss.str()));
  }
};

template <typename node, lattice_type Type, std::size_t n>
struct print_traits<std::array<node, n>, Type> {
  static const std::size_t size = n;

  template <typename container> static std::string print_line(container &cont) {
    std::ostringstream ss;
    for (auto e : cont) {
      ss << "(";
      for (auto i = 0; i < n; ++i) {
        ss << e[i] << ",";
      }
      ss << "), ";
    }
    return clip_comma(std::move(ss.str()));
  }
};

template <typename node, std::size_t n>
struct print_traits<std::array<node, n>, lattice_type::TwoVariableBinomial> {
  static const std::size_t size = n;

  template <typename container> static std::string print_line(container &cont) {
    std::ostringstream ss;
    std::size_t const side = static_cast<std::size_t>(std::sqrt(cont.size()));
    std::array<node, n> val{};
    for (std::size_t i = 0; i < cont.size(); ++i) {
      if ((i > 0) && (i % side == 0))
        ss << "\n ";
      val = cont[i];
      ss << "(";
      for (std::size_t j = 0; j < n; ++j) {
        ss << val[j] << ",";
      }
      ss << "), ";
    }
    return clip_comma(std::move(ss.str()));
  }
};

} // namespace lattice

#endif ///_LATTICE_TRAITS_HPP_