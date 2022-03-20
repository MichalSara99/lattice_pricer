/**
 * @file lattice_version.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief version of lattice library
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_VERSION_HPP_)
#define _LATTICE_VERSION_HPP_

#include <string>

namespace lattice {

std::string const version() { return std::string{"v0.0.0-beta"}; }

} // namespace lattice
#endif //_LATTICE_VERSION_HPP_