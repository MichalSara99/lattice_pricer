/**
 * @file lattice_enums.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief Enumerations
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_ENUMS_HPP_)
#define _LATTICE_ENUMS_HPP_

namespace lattice {

enum class launch_policy { Sequential, Parallel };

enum class lattice_type { Binomial, Trinomial, TwoVariableBinomial };

enum class asset_class { InterestRate, Equity };

enum class discounting_style { Continuous, Discrete };

enum class minimizer_method { Analytic, Numeric };

enum class lattice_class { Normal, MeanReverting };

enum class branching_style { Normal, Reverting };

enum class barrier_type { DownAndOut, DownAndIn, UpAndOut, UpAndIn };

} // namespace lattice

#endif ///_LATTICE_ENUMS_HPP_