/**
 * @file lattice_model_base.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief  base models
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_MODEL_BASE_HPP_)
#define _LATTICE_MODEL_BASE_HPP_

#include <typeinfo>

namespace lattice {

template <std::size_t factor_count, typename T> class binomial_model {};

template <std::size_t factor_count, typename T> class trinomial_model {};

} // namespace lattice

#endif ///_LATTICE_MODEL_BASE_HPP_