#pragma once
#if !defined(_LATTICE_UTILITY_T)
#define _LATTICE_UTILITY_T

#include<iostream>
#include"lattice_utility.h"

void utilityTest1() {

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(2);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";


	lattice_utility::print(il, il.begin()+1, il.end());

	std::cout << il(0, 0) << "\n";
	auto first = il(0, 0);
	il(0, 0) = 3.1415;
	std::cout << il(0, 0) << "\n";

	

	std::cout << "copy assignment: \n";
	auto il_copy = il;
	lattice_utility::print(il_copy, il_copy.begin(), il_copy.end());
}

void utilityTest2() {

	auto today = date(day_clock::local_day());
	auto dd = date_duration{ 1 };
	std::set<date> fixingDates;

	lattice_structure::Lattice<lattice_types::LatticeType::Trinomial, double, date>
		la = { today,today + date_duration(2),today,today + date_duration(1) };

	lattice_utility::print(la, la.begin(), la.end());

	auto first = la.begin();

	lattice_utility::print(la, first, std::next(first, 1));


}



#endif //_LATTICE_UTILITY_T