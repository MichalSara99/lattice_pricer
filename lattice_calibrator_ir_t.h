#pragma once
#if! defined(_LATTICE_CALIBRATOR_IR_T)
#define _LATTICE_CALIBRATOR_IR_T

#include"lattice_structure.h"
#include"lattice_model.h"
#include"lattice_miscellaneous.h"
#include"lattice_utility.h"
#include"lattice_miscellaneous.h"
#include"lattice_calibrator.h"

#include<iostream>
#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>

using namespace boost::gregorian;

void testIndexedBDTCalibration() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

	lattice_miscellaneous::OptionData<double> option;
	option.Volatility = 0.005;

	double dt{ 0.5 };
	std::vector<double>  discount_curve = {
		1.00000,0.97584,0.95223,0.92914,0.90712,
		0.88629,0.86643,0.84724,0.82856,0.81032,
		0.79250, 0.77506,0.75799,0.74127,0.72489,
		0.70884,0.69312,0.67770,0.66258,0.64776,
		0.63322,0.61896,0.60497,0.59125,0.57780,
		0.56460,0.55165,0.53895,0.52649,0.51427,
		0.50229,0.49055,0.47903,0.46774,0.45668,
		0.44584,0.43523,0.42483,0.41464,0.40467,
		0.39492,0.38537,0.37604,0.36690,0.35798,
		0.34925,0.34073, 0.33241,0.32428,0.31635,
		0.30862,0.30107,0.29372,0.28655,0.27957,
		0.27277, 0.26615,0.25971,0.25345,0.24736,
		0.24145
	};

	std::size_t periods{ discount_curve.size() - 1 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> rateTree(periods);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(option);

	std::cout << "\nModel name: " << decltype(bdt)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Binomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> bdt_calibrator;
	// Instantiate Calibrator and launch it:
	bdt_calibrator calibrator(discount_curve);
	auto result = calibrator(rateTree, bdt, dt);

	// Print the part of generated lattice:
	auto first = rateTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(rateTree, first, last);

	// Print theta:
	std::cout << "\n";
	std::vector<double> theta;
	for (auto tpl : result->thetaOptimizers) {
		std::cout << "minimiser: " << std::get<0>(tpl) << "|";
		std::cout << "fun value: " << std::get<1>(tpl) << "|";
		std::cout << "iterations: " << std::get<3>(tpl) << "\n";
	}

}

void testIndexedBDTSanityCheck() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

	lattice_miscellaneous::OptionData<double> option;
	option.Volatility = 0.005;

	double dt{ 0.5 };
	std::vector<double>  discount_curve = {
		1.00000,0.97584,0.95223,0.92914,0.90712,
		0.88629,0.86643,0.84724,0.82856,0.81032,
		0.79250, 0.77506,0.75799,0.74127,0.72489,
		0.70884,0.69312,0.67770,0.66258,0.64776,
		0.63322,0.61896,0.60497,0.59125,0.57780,
		0.56460,0.55165,0.53895,0.52649,0.51427,
		0.50229,0.49055,0.47903,0.46774,0.45668,
		0.44584,0.43523,0.42483,0.41464,0.40467,
		0.39492,0.38537,0.37604,0.36690,0.35798,
		0.34925,0.34073, 0.33241,0.32428,0.31635,
		0.30862,0.30107,0.29372,0.28655,0.27957,
		0.27277, 0.26615,0.25971,0.25345,0.24736,
		0.24145
	};

	std::size_t periods{ discount_curve.size() - 1 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(periods);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(option);

	std::cout << "\nModel name: " << decltype(bdt)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Binomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> bdt_calibrator;
	// Instantiate Calibrator and launch it:
	bdt_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bdt, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// collect the theta into vector:
	std::vector<double> theta;
	for (auto const &t : result->thetaOptimizers) {
		theta.push_back(std::get<0>(t));
	}

	std::cout << "Forward induction after calibration:\n";
	// construct tree from calibrated theta:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> rateTree(periods);
	// setting theta for the model:
	bdt.setTheta(theta);

	//forward traversal algorithm:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, std::size_t, double, double>
		forward_induction;

	forward_induction fwd_induction;

	// this rate I picked up from previous calibartion procedure
	auto r{ calibratedTree.apex() };
	fwd_induction(rateTree, bdt, dt, r);

	// Print the part of generated lattice:
	first = rateTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(rateTree, first, last);

}

void testIndexedBDT() {
	std::cout << "=======================================================\n";
	std::cout << "====== Indexed Black-Derman-Toy Model - TEST ==========\n";
	std::cout << "=======================================================\n";

	testIndexedBDTCalibration();
	testIndexedBDTSanityCheck();

	std::cout << "=======================================================\n";
}





void testBDTCalibration() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

	lattice_miscellaneous::OptionData<double> option;
	option.Volatility = 0.005;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	std::vector<double>  discount_curve = {
		1.00000,0.97584,0.95223,0.92914,0.90712,
		0.88629,0.86643,0.84724,0.82856,0.81032,
		0.79250, 0.77506,0.75799,0.74127,0.72489,
		0.70884,0.69312,0.67770,0.66258,0.64776,
		0.63322,0.61896,0.60497,0.59125,0.57780,
		0.56460,0.55165,0.53895,0.52649,0.51427,
		0.50229,0.49055,0.47903,0.46774,0.45668,
		0.44584,0.43523,0.42483,0.41464,0.40467,
		0.39492,0.38537,0.37604,0.36690,0.35798,
		0.34925,0.34073, 0.33241,0.32428,0.31635,
		0.30862,0.30107,0.29372,0.28655,0.27957,
		0.27277, 0.26615,0.25971,0.25345,0.24736,
		0.24145
	};

	// 
	std::size_t year{ 180 };
	std::size_t periods{ discount_curve.size() - 1 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}


	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> rateTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(option);

	std::cout << "\nModel name: " << decltype(bdt)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Binomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> bdt_calibrator;
	// Instantiate Calibrator and launch it:
	bdt_calibrator calibrator(discount_curve);

	// Get time deltas via boost date
	double daysInYear{ 365.0 };
	auto fd = rateTree.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	auto result = calibrator(rateTree, bdt, timeDeltas);

	// Print the part of generated lattice:
	auto first = rateTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(rateTree, first, last);

	// Print theta:
	std::cout << "\n";
	std::vector<double> theta;
	for (auto tpl : result->thetaOptimizers) {
		std::cout << "minimiser: " << std::get<0>(tpl) << "|";
		std::cout << "fun value: " << std::get<1>(tpl) << "|";
		std::cout << "iterations: " << std::get<3>(tpl) << "\n";
	}

}


void testBDTSanityCheck() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

	lattice_miscellaneous::OptionData<double> option;
	option.Volatility = 0.005;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	std::vector<double>  discount_curve = {
		1.00000,0.97584,0.95223,0.92914,0.90712,
		0.88629,0.86643,0.84724,0.82856,0.81032,
		0.79250, 0.77506,0.75799,0.74127,0.72489,
		0.70884,0.69312,0.67770,0.66258,0.64776,
		0.63322,0.61896,0.60497,0.59125,0.57780,
		0.56460,0.55165,0.53895,0.52649,0.51427,
		0.50229,0.49055,0.47903,0.46774,0.45668,
		0.44584,0.43523,0.42483,0.41464,0.40467,
		0.39492,0.38537,0.37604,0.36690,0.35798,
		0.34925,0.34073, 0.33241,0.32428,0.31635,
		0.30862,0.30107,0.29372,0.28655,0.27957,
		0.27277, 0.26615,0.25971,0.25345,0.24736,
		0.24145
	};

	// 
	std::size_t year{ 180 };
	std::size_t periods{ discount_curve.size() - 1 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}


	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(option);

	std::cout << "\nModel name: " << decltype(bdt)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Binomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> bdt_calibrator;
	// Instantiate Calibrator and launch it:
	bdt_calibrator calibrator(discount_curve);

	// Get time deltas via boost date
	double daysInYear{ 365.0 };
	auto fd = calibratedTree.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	auto result = calibrator(calibratedTree, bdt, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// collect the theta into vector:
	std::vector<double> theta;
	for (auto const &t : result->thetaOptimizers) {
		theta.push_back(std::get<0>(t));
	}

	std::cout << "Forward induction after calibration:\n";
	// construct tree from calibrated theta:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> rateTree(fixingDates);
	// setting theta for the model:
	bdt.setTheta(theta);

	//forward traversal algorithm:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, date, std::vector<double>, double>
		forward_induction;

	forward_induction fwd_induction;

	// this rate I picked up from previous calibartion procedure
	auto r{ calibratedTree.apex() };
	fwd_induction(rateTree, bdt, timeDeltas, r);

	// Print the part of generated lattice:
	first = rateTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(rateTree, first, last);

}


void testBDT() {
	std::cout << "=======================================================\n";
	std::cout << "============= Black-Derman-Toy Model - TEST ==========\n";
	std::cout << "=======================================================\n";

	testBDTCalibration();
	testBDTSanityCheck();

	std::cout << "=======================================================\n";
}


void testIndexedHLCalibration() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

	lattice_miscellaneous::OptionData<double> option;
	option.Volatility = 0.005;

	double dt{ 0.5 };
	std::vector<double>  discount_curve = {
		1.00000,0.97584,0.95223,0.92914,0.90712,
		0.88629,0.86643,0.84724,0.82856,0.81032,
		0.79250, 0.77506,0.75799,0.74127,0.72489,
		0.70884,0.69312,0.67770,0.66258,0.64776,
		0.63322,0.61896,0.60497,0.59125,0.57780,
		0.56460,0.55165,0.53895,0.52649,0.51427,
		0.50229,0.49055,0.47903,0.46774,0.45668,
		0.44584,0.43523,0.42483,0.41464,0.40467,
		0.39492,0.38537,0.37604,0.36690,0.35798,
		0.34925,0.34073, 0.33241,0.32428,0.31635,
		0.30862,0.30107,0.29372,0.28655,0.27957,
		0.27277, 0.26615,0.25971,0.25345,0.24736,
		0.24145
	};

	std::size_t periods{ discount_curve.size() - 1 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> rateTree(periods);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(option);

	std::cout << "\nModel name: " << decltype(hlm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Binomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> bdt_calibrator;
	// Instantiate Calibrator and launch it:
	bdt_calibrator calibrator(discount_curve);
	auto result = calibrator(rateTree, hlm, dt);

	// Print the part of generated lattice:
	auto first = rateTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(rateTree, first, last);

	// Print theta:
	std::cout << "\n";
	std::vector<double> theta;
	for (auto tpl : result->thetaOptimizers) {
		std::cout << "minimiser: " << std::get<0>(tpl) << "|";
		std::cout << "fun value: " << std::get<1>(tpl) << "|";
		std::cout << "iterations: " << std::get<3>(tpl) << "\n";
	}
}

void testIndexedHLSanityCheck() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

	lattice_miscellaneous::OptionData<double> option;
	option.Volatility = 0.005;

	double dt{ 0.5 };
	std::vector<double>  discount_curve = {
		1.00000,0.97584,0.95223,0.92914,0.90712,
		0.88629,0.86643,0.84724,0.82856,0.81032,
		0.79250, 0.77506,0.75799,0.74127,0.72489,
		0.70884,0.69312,0.67770,0.66258,0.64776,
		0.63322,0.61896,0.60497,0.59125,0.57780,
		0.56460,0.55165,0.53895,0.52649,0.51427,
		0.50229,0.49055,0.47903,0.46774,0.45668,
		0.44584,0.43523,0.42483,0.41464,0.40467,
		0.39492,0.38537,0.37604,0.36690,0.35798,
		0.34925,0.34073, 0.33241,0.32428,0.31635,
		0.30862,0.30107,0.29372,0.28655,0.27957,
		0.27277, 0.26615,0.25971,0.25345,0.24736,
		0.24145
	};

	std::size_t periods{ discount_curve.size() - 1 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(periods);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(option);

	std::cout << "\nModel name: " << decltype(hlm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Binomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> bdt_calibrator;
	// Instantiate Calibrator and launch it:
	bdt_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hlm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// collect the theta into vector:
	std::vector<double> theta;
	for (auto const &t : result->thetaOptimizers) {
		theta.push_back(std::get<0>(t));
	}

	std::cout << "Forward induction after calibration:\n";
	// construct tree from calibrated theta:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> rateTree(periods);
	// setting theta for the model:
	hlm.setTheta(theta);

	//forward traversal algorithm:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, std::size_t, double, double>
		forward_induction;

	forward_induction fwd_induction;

	// this rate I picked up from previous calibartion procedure
	auto r{ calibratedTree.apex() };
	fwd_induction(rateTree, hlm, dt, r);

	// Print the part of generated lattice:
	first = rateTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(rateTree, first, last);
}

void testIndexedHL() {
	std::cout << "=======================================================\n";
	std::cout << "============ Indexed Ho-Lee Model - TEST ==============\n";
	std::cout << "=======================================================\n";

	testIndexedHLCalibration();
	testIndexedHLSanityCheck();

	std::cout << "=======================================================\n";
}



void testHLCalibration() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

	lattice_miscellaneous::OptionData<double> option;
	option.Volatility = 0.005;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	std::vector<double>  discount_curve = {
		1.00000,0.97584,0.95223,0.92914,0.90712,
		0.88629,0.86643,0.84724,0.82856,0.81032,
		0.79250, 0.77506,0.75799,0.74127,0.72489,
		0.70884,0.69312,0.67770,0.66258,0.64776,
		0.63322,0.61896,0.60497,0.59125,0.57780,
		0.56460,0.55165,0.53895,0.52649,0.51427,
		0.50229,0.49055,0.47903,0.46774,0.45668,
		0.44584,0.43523,0.42483,0.41464,0.40467,
		0.39492,0.38537,0.37604,0.36690,0.35798,
		0.34925,0.34073, 0.33241,0.32428,0.31635,
		0.30862,0.30107,0.29372,0.28655,0.27957,
		0.27277, 0.26615,0.25971,0.25345,0.24736,
		0.24145
	};

	// 
	std::size_t year{ 180 };
	std::size_t periods{ discount_curve.size() - 1 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double,date> rateTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(option);

	std::cout << "\nModel name: " << decltype(hlm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Binomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> bdt_calibrator;
	// Instantiate Calibrator and launch it:
	bdt_calibrator calibrator(discount_curve);

	// Get time deltas via boost date
	double daysInYear{ 365.0 };
	auto fd = rateTree.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	auto result = calibrator(rateTree, hlm, timeDeltas);

	// Print the part of generated lattice:
	auto first = rateTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(rateTree, first, last);

	// Print theta:
	std::cout << "\n";
	std::vector<double> theta;
	for (auto tpl : result->thetaOptimizers) {
		std::cout << "minimiser: " << std::get<0>(tpl) << "|";
		std::cout << "fun value: " << std::get<1>(tpl) << "|";
		std::cout << "iterations: " << std::get<3>(tpl) << "\n";
	}
}


void testHLSanityCheck() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

	lattice_miscellaneous::OptionData<double> option;
	option.Volatility = 0.005;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	std::vector<double>  discount_curve = {
		1.00000,0.97584,0.95223,0.92914,0.90712,
		0.88629,0.86643,0.84724,0.82856,0.81032,
		0.79250, 0.77506,0.75799,0.74127,0.72489,
		0.70884,0.69312,0.67770,0.66258,0.64776,
		0.63322,0.61896,0.60497,0.59125,0.57780,
		0.56460,0.55165,0.53895,0.52649,0.51427,
		0.50229,0.49055,0.47903,0.46774,0.45668,
		0.44584,0.43523,0.42483,0.41464,0.40467,
		0.39492,0.38537,0.37604,0.36690,0.35798,
		0.34925,0.34073, 0.33241,0.32428,0.31635,
		0.30862,0.30107,0.29372,0.28655,0.27957,
		0.27277, 0.26615,0.25971,0.25345,0.24736,
		0.24145
	};

	// 
	std::size_t year{ 180 };
	std::size_t periods{ discount_curve.size() - 1 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(option);

	std::cout << "\nModel name: " << decltype(hlm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Binomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> bdt_calibrator;
	// Instantiate Calibrator and launch it:
	bdt_calibrator calibrator(discount_curve);

	// Get time deltas via boost date
	double daysInYear{ 365.0 };
	auto fd = calibratedTree.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	auto result = calibrator(calibratedTree, hlm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// collect the theta into vector:
	std::vector<double> theta;
	for (auto const &t : result->thetaOptimizers) {
		theta.push_back(std::get<0>(t));
	}

	std::cout << "Forward induction after calibration:\n";
	// construct tree from calibrated theta:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> rateTree(fixingDates);
	// setting theta for the model:
	hlm.setTheta(theta);

	//forward traversal algorithm:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, date, std::vector<double>, double>
		forward_induction;

	forward_induction fwd_induction;

	// this rate I picked up from previous calibartion procedure
	auto r{ calibratedTree.apex() };
	fwd_induction(rateTree, hlm, timeDeltas, r);

	// Print the part of generated lattice:
	first = rateTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(rateTree, first, last);
}


void testHL() {
	std::cout << "=======================================================\n";
	std::cout << "================ Ho-Lee Model - TEST ==================\n";
	std::cout << "=======================================================\n";

	testHLCalibration();
	testHLSanityCheck();

	std::cout << "=======================================================\n";
}



void testIndexedHWCalibration() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

	lattice_miscellaneous::MeanRevertingParams<double> params;
	lattice_miscellaneous::OptionData<double> option;


	option.ReversionSpeed = 0.25;
	option.Volatility = 0.005;
	params.ReversionSpeed = 0.25;

	double dt{ 0.5 };
	std::vector<double>  discount_curve = {
		1.00000,0.97584,0.95223,0.92914,0.90712,
		0.88629,0.86643,0.84724,0.82856,0.81032,
		0.79250, 0.77506,0.75799,0.74127,0.72489,
		0.70884,0.69312,0.67770,0.66258,0.64776,
		0.63322,0.61896,0.60497,0.59125,0.57780,
		0.56460,0.55165,0.53895,0.52649,0.51427,
		0.50229,0.49055,0.47903,0.46774,0.45668,
		0.44584,0.43523,0.42483,0.41464,0.40467,
		0.39492,0.38537,0.37604,0.36690,0.35798,
		0.34925,0.34073, 0.33241,0.32428,0.31635,
		0.30862,0.30107,0.29372,0.28655,0.27957,
		0.27277, 0.26615,0.25971,0.25345,0.24736,
		0.24145
	};

	std::size_t periods{ discount_curve.size() - 1 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> rateTree(periods,params, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(option);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(rateTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = rateTree.begin();
	//auto last = std::next(first, 17);
	lattice_utility::print(rateTree, first, rateTree.end());

	// Print theta:
	std::cout << "\n";
	std::vector<double> theta;
	for (auto tpl : result->thetaOptimizers) {
		std::cout << "minimiser: " << std::get<0>(tpl) << "|";
		std::cout << "fun value: " << std::get<1>(tpl) << "|";
		std::cout << "iterations: " << std::get<3>(tpl) << "\n";
	}

}


void testIndexedHWSanityCheck() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

	lattice_miscellaneous::MeanRevertingParams<double> params;
	lattice_miscellaneous::OptionData<double> option;


	option.ReversionSpeed = 0.25;
	option.Volatility = 0.005;
	params.ReversionSpeed = 0.25;

	double dt{ 0.5 };
	std::vector<double>  discount_curve = {
		1.00000,0.97584,0.95223,0.92914,0.90712,
		0.88629,0.86643,0.84724,0.82856,0.81032,
		0.79250, 0.77506,0.75799,0.74127,0.72489,
		0.70884,0.69312,0.67770,0.66258,0.64776,
		0.63322,0.61896,0.60497,0.59125,0.57780,
		0.56460,0.55165,0.53895,0.52649,0.51427,
		0.50229,0.49055,0.47903,0.46774,0.45668,
		0.44584,0.43523,0.42483,0.41464,0.40467,
		0.39492,0.38537,0.37604,0.36690,0.35798,
		0.34925,0.34073, 0.33241,0.32428,0.31635,
		0.30862,0.30107,0.29372,0.28655,0.27957,
		0.27277, 0.26615,0.25971,0.25345,0.24736,
		0.24145
	};

	std::size_t periods{ discount_curve.size() - 1 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(option);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = calibratedTree.end();//std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// collect the theta into vector:
	std::vector<double> theta;
	for (auto const &t : result->thetaOptimizers) {
		theta.push_back(std::get<0>(t));
	}

	std::cout << "Forward induction after calibration:\n";
	// construct tree from calibrated theta:
	lattice_structure::MeanRevertingIndexedLattice<double> rateTree(periods, params, dt);
	// setting theta for the model:
	hwm.setTheta(theta);

	//forward traversal algorithm:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial, std::size_t, double, double>
		forward_induction;

	forward_induction fwd_induction;

	// this rate I picked up from previous calibartion procedure
	auto r{ calibratedTree.apex() };
	fwd_induction(rateTree, hwm, dt, r);

	// Print the part of generated lattice:
	first = rateTree.begin();
	last = rateTree.end();
	lattice_utility::print(rateTree, first, last);
}

void testIndexedHW() {
	std::cout << "=======================================================\n";
	std::cout << "======== Indexed Hull-White Model - TEST ==============\n";
	std::cout << "=======================================================\n";

	testIndexedHWCalibration();
	testIndexedHWSanityCheck();

	std::cout << "=======================================================\n";
}


void testHWCalibration() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

	lattice_miscellaneous::MeanRevertingParams<double> params;
	lattice_miscellaneous::OptionData<double> option;


	option.ReversionSpeed = 0.25;
	option.Volatility = 0.005;
	params.ReversionSpeed = 0.25;


	std::vector<double>  discount_curve = {
		1.00000,0.97584,0.95223,0.92914,0.90712,
		0.88629,0.86643,0.84724,0.82856,0.81032,
		0.79250, 0.77506,0.75799,0.74127,0.72489,
		0.70884,0.69312,0.67770,0.66258,0.64776,
		0.63322,0.61896,0.60497,0.59125,0.57780,
		0.56460,0.55165,0.53895,0.52649,0.51427,
		0.50229,0.49055,0.47903,0.46774,0.45668,
		0.44584,0.43523,0.42483,0.41464,0.40467,
		0.39492,0.38537,0.37604,0.36690,0.35798,
		0.34925,0.34073, 0.33241,0.32428,0.31635,
		0.30862,0.30107,0.29372,0.28655,0.27957,
		0.27277, 0.26615,0.25971,0.25345,0.24736,
		0.24145
	};


	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 1 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace_back(today + date_duration(daysInhalfYear*t));
		fixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	float daysInYear{ 365.0 };
	std::vector<double> timeDeltas(fixingDates.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fixingDates[i + 1] - fixingDates[i]).days() / daysInYear);
	}

	// Creating indexed lattice:
	lattice_structure::MeanRevertingLattice<double,date> rateTree(fixingDatesSet, params, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(option);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(rateTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = rateTree.begin();
	//auto last = std::next(first, 17);
	lattice_utility::print(rateTree, first, rateTree.end());

	// Print theta:
	std::cout << "\n";
	std::vector<double> theta;
	for (auto tpl : result->thetaOptimizers) {
		std::cout << "minimiser: " << std::get<0>(tpl) << "|";
		std::cout << "fun value: " << std::get<1>(tpl) << "|";
		std::cout << "iterations: " << std::get<3>(tpl) << "\n";
	}

}


void testHWSanityCheck() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

	lattice_miscellaneous::MeanRevertingParams<double> params;
	lattice_miscellaneous::OptionData<double> option;


	option.ReversionSpeed = 0.25;
	option.Volatility = 0.005;
	params.ReversionSpeed = 0.25;


	std::vector<double>  discount_curve = {
		1.00000,0.97584,0.95223,0.92914,0.90712,
		0.88629,0.86643,0.84724,0.82856,0.81032,
		0.79250, 0.77506,0.75799,0.74127,0.72489,
		0.70884,0.69312,0.67770,0.66258,0.64776,
		0.63322,0.61896,0.60497,0.59125,0.57780,
		0.56460,0.55165,0.53895,0.52649,0.51427,
		0.50229,0.49055,0.47903,0.46774,0.45668,
		0.44584,0.43523,0.42483,0.41464,0.40467,
		0.39492,0.38537,0.37604,0.36690,0.35798,
		0.34925,0.34073, 0.33241,0.32428,0.31635,
		0.30862,0.30107,0.29372,0.28655,0.27957,
		0.27277, 0.26615,0.25971,0.25345,0.24736,
		0.24145
	};


	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 1 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace_back(today + date_duration(daysInhalfYear*t));
		fixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	double daysInYear{ 365.0 };
	std::vector<double> timeDeltas(fixingDates.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fixingDates[i + 1] - fixingDates[i]).days() / daysInYear);
	}

	// Creating indexed lattice:
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(option);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = calibratedTree.end(); //std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// collect the theta into vector:
	std::vector<double> theta;
	for (auto const &t : result->thetaOptimizers) {
		theta.push_back(std::get<0>(t));
	}

	std::cout << "Forward induction after calibration:\n";
	// construct tree from calibrated theta:
	lattice_structure::MeanRevertingLattice<double, date> rateTree(fixingDatesSet, params, timeDeltas);
	// setting theta for the model:
	hwm.setTheta(theta);

	//forward traversal algorithm:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial, date, std::vector<double>, double>
		forward_induction;

	forward_induction fwd_induction;

	// this rate I picked up from previous calibartion procedure
	auto r{ calibratedTree.apex() };
	fwd_induction(rateTree, hwm, timeDeltas, r);

	// Print the part of generated lattice:
	first = rateTree.begin();
	last = rateTree.end();
	lattice_utility::print(rateTree, first, last);


}

void testHW() {
	std::cout << "=======================================================\n";
	std::cout << "============ Hull-White Model - TEST ==================\n";
	std::cout << "=======================================================\n";

	testHWCalibration();
	testHWSanityCheck();

	std::cout << "=======================================================\n";
}

#endif ///_LATTICE_CALIBRATOR_IR_T