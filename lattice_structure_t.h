#pragma once
#if !defined(_LATTICE_STRUCTURE_T)
#define _LATTICE_STRUCTURE_T

#include"lattice_structure.h"
#include<iostream>
#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>

using namespace boost::gregorian;


void createBinomialIndexedLattice1() {

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(2);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << il(0, 0) << "\n";
	auto first = il(0, 0);
	il(0, 0) = 3.1415;
	std::cout << il(0, 0) << "\n";

	std::cout << "copy assignment: \n";
	auto il_copy = il;
	for (auto const &v : il_copy.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

}



void createTrinomialIndexedLattice1() {

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(3);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << il(0, 0) << "\n";
	auto first = il(0, 0);
	il(0, 0) = 3.1415;
	std::cout << il(0, 0) << "\n";

	std::cout << "copy assignment: \n";
	auto il_copy = il;
	for (auto const &v : il_copy.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

}

void createBinomialLattice1() {

	auto today = date(day_clock::local_day());
	auto dd = date_duration{ 1 };
	std::set<date> fixingDates;

	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date>
		la = { today,today + date_duration(2),today,today + date_duration(1) };

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";


	std::cout << la(date(day_clock::local_day()), 0) << "\n";
	auto first_node = la(today, 0);
	la(today, 0) = 3.1415;
	std::cout << la(today, 0) << "\n";
}

void createTrinomialLattice1() {

	auto today = date(day_clock::local_day());
	auto dd = date_duration{ 1 };
	std::set<date> fixingDates;

	lattice_structure::Lattice<lattice_types::LatticeType::Trinomial, double, date>
		la = { today,today + date_duration(2),today,today + date_duration(1) };

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";


	std::cout << la(date(day_clock::local_day()), 0) << "\n";
	auto first_node = la(today, 0);
	la(today, 0) = 3.1415;
	std::cout << la(today, 0) << "\n";
}




#endif //_LATTICE_STRUCTURE_T