#pragma once
#if !defined(_LATTICE_STRUCTURE_T)
#define _LATTICE_STRUCTURE_T

#include"lattice_structure.h"
#include"lattice_miscellaneous.h"
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



void createMeanRevertingTrinomialIndexedLattice() {

	using lattice_miscellaneous::MeanRevertingParams;

	MeanRevertingParams<double> revertingParams;
	revertingParams.ReversionSpeed = 0.25;
	double dt{ 0.5 };

	lattice_structure::MeanRevertingIndexedLattice<double> il(15, revertingParams, 0.5);
	std::cout << "type of tree: " << typeid(lattice_structure::MeanRevertingIndexedLattice<double>::TreeType).name() << "\n";

	lattice_utility::print(il, il.begin(), il.end());

	std::cout << il(0, 0) << "\n";
	auto first = il(0, 0);
	il(0, 0) = 3.1415;
	std::cout << il(0, 0) << "\n";

	std::cout << "copy assignment: \n";
	auto il_copy = il;
	lattice_utility::print(il_copy, il_copy.begin(), il_copy.end());

}



void createMeanRevertingTrinomialLattice() {
	using lattice_miscellaneous::MeanRevertingParams;

	MeanRevertingParams<float> revertingParams;
	revertingParams.ReversionSpeed = 0.25;

	std::size_t daysInhalfYear{ 125 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	for (std::size_t t = 1; t <= 15; ++t) {
		fixingDates.emplace_back(today + date_duration(daysInhalfYear*t));
		fixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	float daysInYear{ 365.0 };
	std::vector<float> timeDeltas(fixingDates.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fixingDates[i + 1] - fixingDates[i]).days() / daysInYear);
	}

	lattice_structure::MeanRevertingLattice<float, date>
		la(fixingDatesSet, revertingParams, timeDeltas);

	lattice_utility::print(la, la.begin(), la.end());


	std::cout << la.at(date(day_clock::local_day()), 0) << "\n";
	auto first_node = la.at(today, 0);
	la.at(today, 0) = 3.1415;
	std::cout << la.at(today, 0) << "\n";
	std::cout << "copy assignment: \n";
	auto la_copy = la;
	lattice_utility::print(la_copy, la_copy.begin(), la_copy.end());
}


void createMeanRevertingTrinomialLattice1() {
	using lattice_miscellaneous::MeanRevertingParams;

	MeanRevertingParams<float> revertingParams;
	revertingParams.ReversionSpeed = 0.25;

	std::size_t daysInYear{ 365 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	fixingDatesSet.emplace(today + date_duration(daysInYear * 1));
	fixingDates.emplace_back(today + date_duration(daysInYear*1));
	fixingDatesSet.emplace(today + date_duration(daysInYear * 2));
	fixingDates.emplace_back(today + date_duration(daysInYear * 2));
	fixingDatesSet.emplace(today + date_duration(daysInYear * 3));
	fixingDates.emplace_back(today + date_duration(daysInYear * 3));
	fixingDatesSet.emplace(today + date_duration(daysInYear * 4));
	fixingDates.emplace_back(today + date_duration(daysInYear * 4));
	fixingDatesSet.emplace(today + date_duration(daysInYear * 5));
	fixingDates.emplace_back(today + date_duration(daysInYear * 5));
	fixingDatesSet.emplace(today + date_duration(daysInYear * 7));
	fixingDates.emplace_back(today + date_duration(daysInYear * 7));
	fixingDatesSet.emplace(today + date_duration(daysInYear * 10));
	fixingDates.emplace_back(today + date_duration(daysInYear * 10));
	fixingDatesSet.emplace(today + date_duration(daysInYear * 15));
	fixingDates.emplace_back(today + date_duration(daysInYear * 15));
	fixingDatesSet.emplace(today + date_duration(daysInYear * 20));
	fixingDates.emplace_back(today + date_duration(daysInYear * 20));
	fixingDatesSet.emplace(today + date_duration(daysInYear * 25));
	fixingDates.emplace_back(today + date_duration(daysInYear * 25));
	fixingDatesSet.emplace(today + date_duration(daysInYear * 30));
	fixingDates.emplace_back(today + date_duration(daysInYear * 30));



	float daysInYearFloat{ 365.0 };
	std::vector<float> timeDeltas(fixingDates.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fixingDates[i + 1] - fixingDates[i]).days() / daysInYearFloat);
	}

	lattice_structure::MeanRevertingLattice<float, date>
		la(fixingDatesSet, revertingParams, timeDeltas);

	lattice_utility::print(la, la.begin(), la.end());

	auto idx = la.firstRevertingIdx();
	std::cout << "First mean-reverting timeidx:" << idx << "\n";

	std::cout << "\n\n";
	std::cout << la.at(date(day_clock::local_day()), 0) << "\n";
	auto first_node = la.at(today, 0);
	la.at(today, 0) = 3.1415;
	std::cout << la.at(today, 0) << "\n";
	std::cout << "copy assignment: \n";
	auto la_copy = la;
	lattice_utility::print(la_copy, la_copy.begin(), la_copy.end());
}

void testMeanRevertingLatticeCreation() {
	std::cout << "=======================================================\n";
	std::cout << "==== Mean-Reverting Lattices Creation - TEST ==========\n";
	std::cout << "=======================================================\n";

	createMeanRevertingTrinomialIndexedLattice();
	createMeanRevertingTrinomialLattice();
	createMeanRevertingTrinomialLattice1();

	std::cout << "=======================================================\n";

}

#endif //_LATTICE_STRUCTURE_T