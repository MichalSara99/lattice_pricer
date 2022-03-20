#define BOOST_TEST_MODULE

#include "lattice/lattice.hpp"
#include "lattice/models/lattice_black_derman_toy.hpp"
#include "lattice/models/lattice_black_karasinski.hpp"
#include "lattice/models/lattice_model_params.hpp"
#include "lattice/traversals/lattice_forward.hpp"
#include "lattice/utilities/lattice_enums.hpp"
#include "lattice/utilities/lattice_printing.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#define BOOST_TEST_DYN_LINK
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/gregorian/gregorian_types.hpp>
#include <boost/test/unit_test.hpp>

using namespace boost::gregorian;
using namespace lattice;

namespace {

std::uniform_real_distribution<> uni(0.23, 0.35);
std::default_random_engine eng;

void uniform_vector(std::vector<double> &rnd_vector) {
  std::generate(rnd_vector.begin(), rnd_vector.end(),
                [&]() { return uni(eng); });
}

} // namespace

BOOST_AUTO_TEST_SUITE()

BOOST_AUTO_TEST_CASE(indexed_forward_bdt) {
  const std::size_t periods = 10;
  std::vector<double> theta(periods);
  uniform_vector(theta);
  ilattice<lattice_type::Binomial, double> bdt_lattice(periods);
  model_params<1, asset_class::InterestRate, double> params;
  params.reversion_speed = 0.15;
  params.volatility = 0.05;
  black_derman_toy<double> bdt(params, theta);
  auto const dt = 0.05;
  auto const apex = 0.3;
  traverse_forward<lattice_type::Binomial, std::size_t, double,
                   double>::traverse(bdt_lattice, bdt, dt, apex);

#if defined __SEEOUTPUT_
  std::cout << "type of tree: "
            << typeid(ilattice<lattice_type::Binomial, double>::tree_t).name()
            << "\n";
#endif //__SEEOUTPUT_
  print(bdt_lattice);
  // print to file:
  std::ofstream file("bdt_indexed.txt", std::ios::out);
  print(bdt_lattice, file);
  file.close();

  print(bdt_lattice, bdt_lattice.cbegin() + 1, bdt_lattice.cend());

#if defined __SEEOUTPUT_
  std::cout << il(0, 0) << "\n";
#endif //__SEEOUTPUT_
  auto first = bdt_lattice(0, 0);
  bdt_lattice(0, 0) = apex;
  BOOST_CHECK_EQUAL(bdt_lattice(0, 0), apex);
#if defined __SEEOUTPUT_
  std::cout << il(0, 0) << "\n";
#endif //__SEEOUTPUT_

#if defined __SEEOUTPUT_
  std::cout << "copy assignment: \n";
#endif //__SEEOUTPUT_
  auto il_copy = bdt_lattice;
  print(il_copy, il_copy.begin(), il_copy.end());
  BOOST_CHECK_EQUAL(il_copy(0, 0), apex);
  BOOST_CHECK_THROW(il_copy(10, 11), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(indexed_forward_bk) {
  const std::size_t periods = 9;
  std::vector<double> theta(periods);
  uniform_vector(theta);
  ilattice<lattice_type::Trinomial, double> bk_lattice(periods);
  model_params<1, asset_class::InterestRate, double> params;
  params.reversion_speed = 0.15;
  params.volatility = 0.05;
  black_karasinski<double> bk(params, theta);
  auto const dt = 0.05;
  auto const apex = 0.3;
  traverse_forward<lattice_type::Trinomial, std::size_t, double,
                   double>::traverse(bk_lattice, bk, dt, apex);

#if defined __SEEOUTPUT_
  std::cout << "type of tree: "
            << typeid(ilattice<lattice_type::Binomial, double>::tree_t).name()
            << "\n";
#endif //__SEEOUTPUT_
  print(bk_lattice);
  // print to file:
  std::ofstream file("bk_indexed.txt", std::ios::out);
  print(bk_lattice, file);
  file.close();

  print(bk_lattice, bk_lattice.cbegin() + 1, bk_lattice.cend());

#if defined __SEEOUTPUT_
  std::cout << bk_lattice(0, 0) << "\n";
#endif //__SEEOUTPUT_
  auto first = bk_lattice(0, 0);
  bk_lattice(0, 0) = apex;
  BOOST_CHECK_EQUAL(bk_lattice(0, 0), apex);
#if defined __SEEOUTPUT_
  std::cout << bk_lattice(0, 0) << "\n";
#endif //__SEEOUTPUT_

#if defined __SEEOUTPUT_
  std::cout << "copy assignment: \n";
#endif //__SEEOUTPUT_
  auto il_copy = bk_lattice;
  print(il_copy, il_copy.begin(), il_copy.end());
  BOOST_CHECK_EQUAL(il_copy(0, 0), apex);
  BOOST_CHECK_THROW(il_copy(1, 11), std::out_of_range);
}

BOOST_AUTO_TEST_SUITE_END()