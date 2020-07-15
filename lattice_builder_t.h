#pragma once
#if !defined(_LATTICE_BUILDER_T)
#define _LATTICE_BUILDER_T


#include"lattice_builder.h"
#include"lattice_types.h"
#include"lattice_utility.h"
#include"lattice_miscellaneous.h"
#include"lattice_algorithms.h"
#include"lattice_model.h"
#include"lattice_calibrator.h"
#include"lattice_product_builder.h"


#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>

#include<iostream>

using namespace boost::gregorian;


void testCreateIndexBasedOneFactor() {

	using lattice_builder::LatticeBuilder;
	using lattice_types::LatticeType;
	using lattice_utility::print;

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	double maturity = 0.29;
	double dt = (maturity / double(option.periods()));

	// create one-factor index-based lattice:
	auto binLattice = LatticeBuilder<1>::createIndexedBasedLattice<LatticeType::Binomial, double>(option.periods());

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ params };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(binLattice, crr, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = binLattice.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(binLattice, first, last);

}

void testCreateIndexBasedBinomialOneFactor() {

	using lattice_builder::LatticeBuilder;
	using lattice_types::LatticeType;
	using lattice_utility::print;

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	double maturity = 0.29;
	double dt = (maturity / double(option.periods()));

	// create one-factor index-based lattice:
	auto binLattice = LatticeBuilder<1>::createIndexedBasedBinomialLattice<double>(option.periods());

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ params };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(binLattice, crr, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = binLattice.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(binLattice, first, last);

}

void testCreateIndexBasedTrinomialOneFactor() {

	using lattice_builder::LatticeBuilder;
	using lattice_types::LatticeType;
	using lattice_utility::print;

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	double maturity = 0.29;
	double dt = (maturity / double(option.periods()));

	// create one-factor index-based lattice:
	auto triLattice = LatticeBuilder<1>::createIndexedBasedTrinomialLattice<double>(option.periods());

	// Create Boyle model:
	lattice_model::BoyleModel<> bm{ params };

	// Print the model name:
	std::cout << decltype(bm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(triLattice, bm, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = triLattice.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(triLattice, first, last);

}



void testCreateIndexBasedMROneFactor() {

	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_builder::LatticeBuilder;
	using lattice_types::LatticeType;
	using lattice_utility::print;

	ModelParams<1,AssetClass::InterestRate,double> params;
	params.ReversionSpeed = 0.25;
	params.Volatility = 0.005;

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

	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };

	// create one-factor index-based lattice:
	auto triMRLattice = LatticeBuilder<1>::createIndexedBasedMRLattice<double, double>(periods, params.ReversionSpeed, dt);

	// Create Hull-White model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(triMRLattice, hwm, dt);

	// Print the part of generated lattice:
	auto first = triMRLattice.begin();
	//auto last = std::next(first, 17);
	lattice_utility::print(triMRLattice, first, triMRLattice.end());

}


void testCreateIndexBasedLattice() {
	std::cout << "=======================================================\n";
	std::cout << "======== Creation Index-based Lattices - TEST =========\n";
	std::cout << "=======================================================\n";

	testCreateIndexBasedOneFactor();
	testCreateIndexBasedBinomialOneFactor();
	testCreateIndexBasedTrinomialOneFactor();
	testCreateIndexBasedMROneFactor();

	std::cout << "=======================================================\n";
}


void testCreateLatticeOneFactor() {

	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_builder::LatticeBuilder;
	using lattice_types::LatticeType;
	using lattice_utility::print;

	ModelParams<1, AssetClass::Equity, double> params;

	params.Strike = 65.0;
	params.RiskFreeRate = 0.25;
	params.DividendRate = 0.05;
	params.Volatility = 0.3;
	params.Spot = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// create one-factor index-based lattice:
	auto binLattice = LatticeBuilder<1>::createLattice<LatticeType::Binomial, double, date>(fixingDates);

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = binLattice.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ params };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(binLattice, crr, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = binLattice.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(binLattice, first, last);

}

void testCreateBinomialOneFactor() {

	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_builder::LatticeBuilder;
	using lattice_types::LatticeType;
	using lattice_utility::print;

	ModelParams<1, AssetClass::Equity, double> params;

	params.Strike = 65.0;
	params.RiskFreeRate = 0.25;
	params.DividendRate = 0.05;
	params.Volatility = 0.3;
	params.Spot = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// create one-factor index-based lattice:
	auto binLattice = LatticeBuilder<1>::createBinomialLattice<double, date>(fixingDates);

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = binLattice.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::CoxRubinsteinRossModel<> crr{ params };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(binLattice, crr, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = binLattice.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(binLattice, first, last);

}

void testCreateTrinomialOneFactor() {

	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_builder::LatticeBuilder;
	using lattice_types::LatticeType;
	using lattice_utility::print;

	ModelParams<1, AssetClass::Equity, double> params;

	params.Strike = 65.0;
	params.RiskFreeRate = 0.25;
	params.DividendRate = 0.05;
	params.Volatility = 0.3;
	params.Spot = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// create one-factor index-based lattice:
	auto triLattice = LatticeBuilder<1>::createTrinomialLattice<double, date>(fixingDates);

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = triLattice.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}


	// Create Boyle model:
	lattice_model::BoyleModel<> bm{ params };

	// Print the model name:
	std::cout << decltype(bm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>, double> forward_trinomial_induction;

	forward_trinomial_induction fwd_induction;
	fwd_induction(triLattice, bm, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = triLattice.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(triLattice, first, last);

}



void testCreateMROneFactor() {

	using lattice_builder::LatticeBuilder;
	using lattice_types::LatticeType;
	using lattice_utility::print;
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1,AssetClass::InterestRate,double> params;
	params.ReversionSpeed = 0.25;
	params.Volatility = 0.005;

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
	std::vector<date> fixingDates;
	fixingDates.emplace_back(today);

	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace_back(today + date_duration(daysInhalfYear*t));
	}

	float daysInYear{ 365.0 };
	std::vector<double> timeDeltas(fixingDates.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fixingDates[i + 1] - fixingDates[i]).days() / daysInYear);
	}

	// create one-factor index-based lattice:
	auto triMRLattice = LatticeBuilder<1>::createMRLattice<double, date, std::vector<double>>(fixingDates, params.ReversionSpeed, timeDeltas);

	// Create Hull-White model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(triMRLattice, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = triMRLattice.begin();
	//auto last = std::next(first, 17);
	lattice_utility::print(triMRLattice, first, triMRLattice.end());

}


void testCreateLattice() {
	std::cout << "=======================================================\n";
	std::cout << "============= Creation Lattices - TEST ================\n";
	std::cout << "=======================================================\n";

	testCreateLatticeOneFactor();
	testCreateBinomialOneFactor();
	testCreateTrinomialOneFactor();
	testCreateMROneFactor();

	std::cout << "=======================================================\n";
}






#endif ///_LATTICE_BUILDER_T