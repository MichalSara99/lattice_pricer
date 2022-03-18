#define BOOST_TEST_MODULE

#include "lattice/lattice.hpp"
#include "utilities/lattice_enums.hpp"
#include "utilities/lattice_printing.hpp"
#include <fstream>
#include <iostream>

#define BOOST_TEST_DYN_LINK
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/gregorian/gregorian_types.hpp>
#include <boost/test/unit_test.hpp>

using namespace boost::gregorian;
using namespace lattice;

BOOST_AUTO_TEST_SUITE()

BOOST_AUTO_TEST_CASE(print_tree_bi_indexed) {
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

BOOST_AUTO_TEST_CASE(print_tree_tri_indexed) {
  // creating 2-period binomial double tree:
  std::size_t period = 2;
  ilattice<lattice_type::Trinomial, double> il(period);
#if defined __SEEOUTPUT_
  std::cout << "type of tree: "
            << typeid(ilattice<lattice_type::Trinomial, double>::tree_t).name()
            << "\n";
#endif //__SEEOUTPUT_

  print(il, il.begin() + 1, il.end());

#if defined __SEEOUTPUT_
  std::cout << il(0, 0) << "\n";
#endif //__SEEOUTPUT_
  auto first = il.at(0, 0);
  auto pi = 3.1415;
  il.at(0, 0) = pi;
  BOOST_CHECK_EQUAL(il.at(0, 0), pi);
#if defined __SEEOUTPUT_
  std::cout << il(0, 0) << "\n";
#endif //__SEEOUTPUT_

#if defined __SEEOUTPUT_
  std::cout << "copy assignment: \n";
#endif //__SEEOUTPUT_
  auto il_copy = il;
  print(il_copy, il_copy.begin(), il_copy.end());
  BOOST_CHECK_EQUAL(il_copy.at(0, 0), pi);
  BOOST_CHECK_THROW(il_copy.at(2, 5), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(print_tree_tri_dated) {
  // creating 3-period binomial double tree:
  auto today = date(day_clock::local_day());
  auto dd = date_duration{1};
  std::set<date> fixingDates;

  glattice<lattice_type::Trinomial, double, date> la = {
      today, today + date_duration(2), today, today + date_duration(1)};

  print(la);
  // print to file:
  std::ofstream file("tree_tri_dated.txt", std::ios::out);
  print(la, file);
  file.close();
  print(la, std::next(la.cbegin(), 1), la.cend());

  auto pi = 3.1415;
  la(0, 0) = pi;
  BOOST_CHECK_EQUAL(la(0, 0), pi);
  auto const la_copy = la;
  BOOST_CHECK_THROW(la_copy(2, 5), std::out_of_range);

  auto first = la.begin();

  print(la, first, std::next(first, 1));
}

BOOST_AUTO_TEST_SUITE_END()