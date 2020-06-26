#pragma once
#if !defined(_LATTICE_MULTIDIMENSIONAL_T)
#define _LATTICE_MULTIDIMENSIONAL_T


#include"lattice_types.h"
#include"lattice_structure.h"
#include"lattice_multidimensional.h"
#include"lattice_utility.h"

void testCreate2DIndexedLattice(){


	using lattice_multidimensional::MultidimIndexedLattice;
	using lattice_structure::IndexedLattice;
	using lattice_types::LatticeType;
	using lattice_utility::print;

	std::size_t periods = 5;

	MultidimIndexedLattice<2,LatticeType::Binomial, double> index2Dtree(periods);

	std::cout << "Number of factors: " << index2Dtree.factors() << "\n";
	index2Dtree(0, 0, 0) = 3.1415;
	index2Dtree(1, 1, 0) = std::exp(1.0);

	auto firstFactor = index2Dtree.getFactor(0);
	auto secondFactor = index2Dtree.getFactor(1);

	std::cout << "Print first factor: \n";
	print(firstFactor, firstFactor.begin(), firstFactor.end());
	std::cout << "\nPrint second factor: \n";
	print(secondFactor, secondFactor.begin(), secondFactor.end());


	std::cout << "My apex: \n";
	std::cout << "1: " << index2Dtree.apex(0) << "\n";
	std::cout << "2: " << index2Dtree.apex(1) << "\n";

}











#endif ///_LATTICE_MULTIDIMENSIONAL_T