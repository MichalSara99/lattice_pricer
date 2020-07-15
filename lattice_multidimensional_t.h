#pragma once
#if !defined(_LATTICE_MULTIDIMENSIONAL_T)
#define _LATTICE_MULTIDIMENSIONAL_T


#include"lattice_types.h"
#include"lattice_structure.h"
#include"lattice_multidimensional.h"
#include"lattice_utility.h"


#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>

using namespace boost::gregorian;

void testCreate2DIndexedLattice(){


	std::cout << "This is 2D indexed-based lattice\n";

	using lattice_multidimensional::MultidimIndexedLattice;
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

void testCreate3DLattice() {

	std::cout << "This is 3D date-based lattice\n";

	using lattice_multidimensional::MultidimLattice;
	using lattice_types::LatticeType;
	using lattice_utility::print;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 10 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	MultidimLattice<3, LatticeType::Binomial, double, date> tree3D(fixingDates);

	std::cout << "Number of factors: " << tree3D.factors() << "\n";
	tree3D(0, 0, 0) = 3.1415;
	tree3D(1, 1, 0) = std::exp(1.0);

	auto firstFactor = tree3D.getFactor(0);
	auto secondFactor = tree3D.getFactor(1);

	std::cout << "Print first factor: \n";
	print(firstFactor, firstFactor.begin(), firstFactor.end());
	std::cout << "\nPrint second factor: \n";
	print(firstFactor, secondFactor.begin(), secondFactor.end());


	std::cout << "My apex: \n";
	std::cout << "1: " << tree3D.apex(0) << "\n";
	std::cout << "2: " << tree3D.apex(1) << "\n";

}


void testCreateMultidimensionalLattice() {
	std::cout << "=======================================================\n";
	std::cout << "====== Create Multidimensional Lattice - TEST =========\n";
	std::cout << "=======================================================\n";

	testCreate2DIndexedLattice();
	testCreate3DLattice();

	std::cout << "=======================================================\n";
}

void testCreate2DMeanrevertingIndexedLattice() {


	std::cout << "This is 2D indexed-based mean-reverting lattice\n";

	using lattice_multidimensional::MultidimMeanRevertingIndexedLattice;
	using lattice_types::LatticeType;
	using lattice_utility::print;


	double reversionSpeed = 0.25;
	std::size_t periods = 10;
	double dt{ 0.5 };

	MultidimMeanRevertingIndexedLattice<2, double> index2DMRtree(periods, reversionSpeed, dt);

	std::cout << "Number of factors: " << index2DMRtree.factors() << "\n";
	index2DMRtree(0, 0, 0) = 3.1415;
	index2DMRtree(1, 1, 0) = std::exp(1.0);

	auto firstFactor = index2DMRtree.getFactor(0);
	auto secondFactor = index2DMRtree.getFactor(1);

	std::cout << "Print first factor: \n";
	print(firstFactor, firstFactor.begin(), firstFactor.end());
	std::cout << "\nPrint second factor: \n";
	print(secondFactor, secondFactor.begin(), secondFactor.end());


	std::cout << "My apex: \n";
	std::cout << "1: " << index2DMRtree.apex(0) << "\n";
	std::cout << "2: " << index2DMRtree.apex(1) << "\n";

}

void testCreate3DMeanRevertingLattice() {

	std::cout << "This is 3D date-based mean-reverting lattice\n";

	using lattice_multidimensional::MultidimMeanRevertingLattice;
	using lattice_types::LatticeType;
	using lattice_utility::print;

	double reversionSpeed = 0.25;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 10 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	double dt{ 0.5 };
	MultidimMeanRevertingLattice<3, double, date> tree3DMR(fixingDates, reversionSpeed, dt);

	std::cout << "Number of factors: " << tree3DMR.factors() << "\n";
	tree3DMR(0, 0, 0) = 3.1415;
	tree3DMR(1, 1, 0) = std::exp(1.0);

	auto firstFactor = tree3DMR.getFactor(0);
	auto secondFactor = tree3DMR.getFactor(1);

	std::cout << "Print first factor: \n";
	print(firstFactor, firstFactor.begin(), firstFactor.end());
	std::cout << "\nPrint second factor: \n";
	print(firstFactor, secondFactor.begin(), secondFactor.end());


	std::cout << "My apex: \n";
	std::cout << "1: " << tree3DMR.apex(0) << "\n";
	std::cout << "2: " << tree3DMR.apex(1) << "\n";

}


void testCreateMultidimensionalMRLattice() {
	std::cout << "=========================================================\n";
	std::cout << "= Create Multidimensional Mean Reverting Lattice - TEST =\n";
	std::cout << "=========================================================\n";

	testCreate2DMeanrevertingIndexedLattice();
	testCreate3DMeanRevertingLattice();

	std::cout << "==========================================================\n";
}





#endif ///_LATTICE_MULTIDIMENSIONAL_T