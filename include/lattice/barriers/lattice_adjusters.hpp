/**
 * @file lattice_adjusters.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief Adjusters for the lattice
 * @version 0.1
 * @date 2022-03-20
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_ADJUSTERS_HPP_)
#define _LATTICE_ADJUSTERS_HPP_

#include "../utilities/lattice_enums.hpp"
#include <functional>

namespace lattice {

template <typename T>
using check_adjuster_pair_t = std::pair<std::function<bool(T, T, T, T)>,
                                        std::function<T(T, T, T, T, T, T)>>;

/**
 * @brief Barrier comparer struct for adjusters below
 *
 * @tparam T
 */
template <typename T> struct barrier_comparer {
  static std::function<bool(T, T)> const comparer(barrier_type b_type) {
    switch (b_type) {
    case lattice_types::barrier_type::DownAndOut:
      return [](T stock, T barrier) -> bool { return (stock > barrier); };
      break;
    case lattice_types::barrier_type::UpAndIn:
      return [](T stock, T barrier) -> bool { return (stock >= barrier); };
      break;
    case lattice_types::barrier_type::DownAndIn:
      return [](T stock, T barrier) -> bool { return (stock <= barrier); };
      break;
    case lattice_types::barrier_type::UpAndOut:
      return [](T stock, T barrier) -> bool { return (stock < barrier); };
      break;
    }
  }
};

/**
 * @brief Derman-Kani-Ergener adjuster struct
 *
 * @tparam T
 */
template <typename T> struct derman_kani_ergener_adjuster {
  static check_adjuster_pair_t<T> const
  adjuster(lattice_types::barrier_type b_type) {
    auto cmp = barrier_comparer<T>::comparer(b_type);
    switch (b_type) {
    case lattice_types::barrier_type::DownAndOut:
    case lattice_types::barrier_type::UpAndIn: {
      auto checker = [=](T stock, T stockDown, T stockUp, T barrier) -> bool {
        return (cmp(stock, barrier) && (stockDown <= barrier));
      };

      auto adjuster = [](T stock, T stockDown, T stockUp, T barrier, T rebate,
                         T optionPrice) -> T {
        return (std::max(rebate - optionPrice, optionPrice - rebate) /
                (stockDown - stock)) *
               (barrier - stock);
      };
      return std::make_pair(checker, adjuster);
    } break;
    case lattice_types::barrier_type::DownAndIn:
    case lattice_types::barrier_type::UpAndOut: {
      auto checker = [=](T stock, T stockDown, T stockUp, T barrier) -> bool {
        return (cmp(stock, barrier) && (stockUp >= barrier));
      };

      auto adjuster = [](T stock, T stockDown, T stockUp, T barrier, T rebate,
                         T optionPrice) -> T {
        return (std::max(rebate - optionPrice, optionPrice - rebate) /
                (stockUp - stock)) *
               (barrier - stock);
      };
      return std::make_pair(checker, adjuster);
    } break;
    }
  }
};

} // namespace lattice

#endif ///_LATTICE_ADJUSTERS_HPP_