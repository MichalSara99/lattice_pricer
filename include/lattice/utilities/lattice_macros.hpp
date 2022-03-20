/**
 * @file lattice_macros.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief  Some macros
 * @version 0.1
 * @date 2022-03-18
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#if !defined(_LATTICE_MACROS_HPP_)
#define _LATTICE_MACROS_HPP_

#define LASSERT(condition, message)                                            \
  do {                                                                         \
    if (!(condition)) {                                                        \
      std::cerr << "Assertion `" #condition "` failed in " << __FILE__         \
                << " line " << __LINE__ << ": " << message << std::endl;       \
      std::terminate();                                                        \
    }                                                                          \
  } while (false)

#endif ///_LATTICE_MACROS_HPP_