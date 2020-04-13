#pragma once
#if !defined(_LATTICE_STRUCTURE_T)
#define _LATTICE_STRUCTURE_T

#include"lattice_structure.h"
#include"lattice_utility.h"
#include<iostream>
#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>

using namespace boost::gregorian;


void createBinomialIndexedLattice() {

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(5);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";

	lattice_utility::print(il, il.begin(), il.end());

	std::cout << il(0, 0) << "\n";
	auto first = il(0, 0);
	il(0, 0) = 3.1415;
	std::cout << il(0, 0) << "\n";

	std::cout << "copy assignment: \n";
	auto il_copy = il;
	lattice_utility::print(il_copy, il_copy.begin(), il_copy.end());

}



void createTrinomialIndexedLattice() {

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(5);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	lattice_utility::print(il, il.begin(), il.end());

	std::cout << il(0, 0) << "\n";
	auto first = il(0, 0);
	il(0, 0) = 3.1415;
	std::cout << il(0, 0) << "\n";

	std::cout << "copy assignment: \n";
	auto il_copy = il;
	lattice_utility::print(il_copy, il_copy.begin(), il_copy.end());

}


void testIndexedLatticeCreation() {
	std::cout << "=======================================================\n";
	std::cout << "========== Indexed Lattices Creation - TEST ===========\n";
	std::cout << "=======================================================\n";

	createBinomialIndexedLattice();
	createTrinomialIndexedLattice();

	std::cout << "=======================================================\n";

}


void createBinomialLattice() {

	auto today = date(day_clock::local_day());
	auto dd = date_duration{ 1 };
	std::set<date> fixingDates;

	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date>
		la = { today,today + date_duration(2),today,today + date_duration(1),today + date_duration(3) };

	lattice_utility::print(la, la.begin(), la.end());


	std::cout << la.at(date(day_clock::local_day()), 0) << "\n";
	auto first_node = la.at(today, 0);
	la.at(today, 0) = 3.1415;
	std::cout << la.at(today, 0) << "\n";

	std::cout << "copy assignment: \n";
	auto la_copy = la;
	lattice_utility::print(la_copy, la_copy.begin(), la_copy.end());
}

void createTrinomialLattice() {

	auto today = date(day_clock::local_day());
	auto dd = date_duration{ 1 };
	std::set<date> fixingDates;

	lattice_structure::Lattice<lattice_types::LatticeType::Trinomial, double, date>
		la = { today,today + date_duration(2),today,today + date_duration(1),today + date_duration(3) };

	lattice_utility::print(la, la.begin(), la.end());


	std::cout << la.at(date(day_clock::local_day()), 0) << "\n";
	auto first_node = la.at(today, 0);
	la.at(today, 0) = 3.1415;
	std::cout << la.at(today, 0) << "\n";
	std::cout << "copy assignment: \n";
	auto la_copy = la;
	lattice_utility::print(la_copy, la_copy.begin(), la_copy.end());
}


void testLatticeCreation() {
	std::cout << "=======================================================\n";
	std::cout << "=============== Lattices Creation - TEST ==============\n";
	std::cout << "=======================================================\n";

	createBinomialLattice();
	createTrinomialLattice();

	std::cout << "=======================================================\n";

}


#endif //_LATTICE_STRUCTURE_T