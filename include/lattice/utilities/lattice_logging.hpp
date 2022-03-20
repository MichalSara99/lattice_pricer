/**
 * @file lattice_logging.hpp
 * @author Michal Sara (michal.sara99@gmail.com)
 * @brief  Logger
 * @version 0.1
 * @date 2022-03-18
 *
 * @copyright Copyright (c) 2022
 *
 */
#if !defined(_LATTICE_LOGGING_HPP_)
#define _LATTICE_LOGGING_HPP_

#include <iostream>
#include <string>

namespace lattice {

/**
 * @brief logger object
 *
 */
class logger {
protected:
  explicit logger(){};

  void inline what(std::ostream &out, std::string text) const { out << text; }

public:
  virtual ~logger() {}

  logger(logger const &) = delete;
  logger(logger &&) = delete;
  logger &operator=(logger const &) = delete;
  logger &operator=(logger &&) = delete;

  static inline logger &get() {
    static logger log;
    return log;
  }

  void inline warning(std::ostream &out, std::string text) const {
    this->what(out, std::move(std::string{"WARNING: " + text}));
  }

  void inline critical(std::ostream &out, std::string text) const {
    this->what(out, std::move(std::string{"CRITICAL: " + text}));
  }
};

} // namespace lattice

#endif ///_LATTICE_LOGGING_HPP_