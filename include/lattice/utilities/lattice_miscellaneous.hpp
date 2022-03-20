/**
 * @file lattice_miscellaneous.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief Miscellaneous stuff
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_MISCELLANEOUS_HPP_)
#define _LATTICE_MISCELLANEOUS_HPP_

#include "lattice_logging.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

namespace lattice {

/**
 * @brief Function that clips comma in the passed string str
 *
 * @param str
 * @return std::string const
 */
static std::string const clip_comma(std::string const &str) {
  return str.substr(0, str.size() - 2);
}

/**
 * @brief Signum function
 *
 * @tparam T
 * @param x
 * @return T
 */
template <typename T> T sign(T x) {
  return x < T{0} ? T{-1} : (x > T{0} ? T{1} : T{0});
}

/**
 * @brief Function that linearly interpolates between two
 * given values y0 and y1
 *
 * @tparam T
 * @param y0
 * @param y1
 * @param t
 * @return T
 */
template <typename T> T lerp(T y0, T y1, T t) {
  return std::fma(t, (y1 - y0), y0);
}

/**
 * @brief Scoped thread that destroys itself on the exit of scope
 *
 */
class scoped_thread {
private:
  std::thread t_;

public:
  explicit scoped_thread(std::thread thread) : t_{std::move(thread)} {
    if (!t_.joinable()) {
      throw std::logic_error("No thread\n");
    }
  }

  ~scoped_thread() { t_.join(); }

  scoped_thread(scoped_thread const &other) = delete;

  scoped_thread &operator=(scoped_thread const &other) = delete;
};

/**
 * @brief Function that keeps the probability capped
 * in the usual interval <0,1>. Any over/underflows
 * are logged
 *
 * @tparam T
 * @param x
 * @return T
 */
template <typename T> T probability_cap(T x) {
  if (x < 0.0 || x > 1.0) {
    std::stringstream ss;
    ss << "Probability value " << x << " has been modified to fit <0,1>\n";
    logger::get().critical(std::cout, ss.str());
  }
  return std::min(1.0, std::max(0.0, x));
}

namespace map_traits {

template <typename...> struct voider { using type = void; };

template <typename... T> using void_t = typename voider<T...>::type;

template <typename T, typename U = void>
struct is_map_impl : std::false_type {};

template <typename T>
struct is_map_impl<T, void_t<typename T::key_type, typename T::mapped_type,
                             decltype(std::declval<T &>()[std::declval<
                                 const typename T::key_type &>()])>>
    : std::true_type {};
} // namespace map_traits

template <typename T> struct is_map : map_traits::is_map_impl<T>::type {};

} // namespace lattice

#endif //_LATTICE_MISCELLANEOUS_HPP_