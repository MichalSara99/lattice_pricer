#define BOOST_TEST_MODULE

#include "lattice/lattice.hpp"
#include "lattice/utilities/lattice_enums.hpp"
#include "lattice/utilities/lattice_printing.hpp"
#include <fstream>
#include <iostream>

#define BOOST_TEST_DYN_LINK
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/gregorian/gregorian_types.hpp>
#include <boost/test/unit_test.hpp>

using namespace boost::gregorian;
using namespace lattice;

BOOST_AUTO_TEST_SUITE()

BOOST_AUTO_TEST_CASE(create_indexed_bi_lattice) {
  // creating 2-period binomial double tree:
  std::size_t period = 2;
  ilattice<lattice_type::Binomial, double> il(period);
#if defined __SEEOUTPUT_
  std::cout << "type of tree: "
            << typeid(ilattice<lattice_type::Binomial, double>::tree_t).name()
            << "\n";
#endif //__SEEOUTPUT_
  print(il);
  // print to file:
  std::ofstream file("tree_bi_indexed.txt", std::ios::out);
  print(il, file);
  file.close();

  print(il, il.cbegin() + 1, il.cend());

#if defined __SEEOUTPUT_
  std::cout << il(0, 0) << "\n";
#endif //__SEEOUTPUT_
  auto first = il(0, 0);
  auto pi = 3.1415;
  il(0, 0) = pi;
  BOOST_CHECK_EQUAL(il(0, 0), pi);
#if defined __SEEOUTPUT_
  std::cout << il(0, 0) << "\n";
#endif //__SEEOUTPUT_

#if defined __SEEOUTPUT_
  std::cout << "copy assignment: \n";
#endif //__SEEOUTPUT_
  auto il_copy = il;
  print(il_copy, il_copy.begin(), il_copy.end());
  BOOST_CHECK_EQUAL(il_copy(0, 0), pi);
  BOOST_CHECK_THROW(il_copy(2, 3), std::out_of_range);
}

BOOST_AUTO_TEST_SUITE_END()