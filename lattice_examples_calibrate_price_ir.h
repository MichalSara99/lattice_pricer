#pragma once
#if !defined(_LATTICE_EXAMPLES_CALIBRATE_PRICE_IR)
#define _LATTICE_EXAMPLES_CALIBRATE_PRICE_IR


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

// ===========================================================================================
// ============================== coupon bond binomial tree ==================================
// = Convinience function to populate binomial coupon bond lattice
// = To get pure discount bond lattice just pass couponRate with 0.0

using lattice_utility::DiscountingFactor;
using lattice_utility::DeltaTimeHolder;
using lattice_utility::DiscountingStyle;

template<typename T,typename DeltaTime,typename BinomialLattice>
void createBinomialCouponBondLattice(BinomialLattice &bondLattice, 
	BinomialLattice const &calibartedRateLattice,T nominal,T couponRate,
	DeltaTime const &deltaTime,DiscountingStyle style = DiscountingStyle::Continuous) {

	typedef DiscountingFactor<T> DCF;
	typedef DeltaTimeHolder<DeltaTime> DT;
	auto dcf = DCF::function(style);

	LASSERT(bondLattice.timeDimension() == calibartedRateLattice.timeDimension(),
		"Lattices must have the same dimension.");

	const std::size_t lastIdx = bondLattice.timeDimension() - 1;
	const std::size_t lastNodesSize = bondLattice.nodesAtIdx(lastIdx).size();
	T dt{};
	T coupon{};

	for (auto i = 0; i < lastNodesSize; ++i) {
		bondLattice(lastIdx, i) = nominal;
	}

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = bondLattice.nodesAtIdx(n).size();
		dt = DT::deltaTime(n, deltaTime);
		coupon = nominal * couponRate*dt;
		for (auto i = 0; i < nodesSize; ++i) {
			bondLattice(n, i) =
				0.5*(bondLattice(n + 1, i) + bondLattice(n + 1, i + 1) + 2.0 * coupon)*
				dcf(calibartedRateLattice(n, i), dt);
		}
	}

	dt = DT::deltaTime(0, deltaTime);
	bondLattice(0, 0) = 0.5*(bondLattice(1, 0) + bondLattice(1, 1) + 2.0*coupon)*
		dcf(calibartedRateLattice(0, 0), dt);


}


// ===========================================================================================
// ===========================================================================================





void indexedBDTCouponBond() {

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


	double nominal{ 100.0 };
	double couponRate{ 0.05 };


	// Creating indexed bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(periods);

	// populate coupon bond tree:
	createBinomialCouponBondLattice(bondTree, calibratedTree, nominal, couponRate, dt);

	std::cout << "indexed Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Todays price of coupon bond (nominal=100, couponRate=0.05, coupon paid semiannualy):"
		<< bondTree.apex() << "\n";

}



void indexedHLCouponBond() {

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


	double nominal{ 100.0 };
	double couponRate{ 0.05 };

 
	// Creating indexedbond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(periods);

	// populate coupon bond tree:
	createBinomialCouponBondLattice(bondTree, calibratedTree, nominal, couponRate, dt);

	std::cout << "indexed Ho-Lee coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Todays price of coupon bond (nominal=100, couponRate=0.05, coupon paid semiannualy):"
		<< bondTree.apex() << "\n";

}

void testIndexedBinomialCouponBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "========= Indexed Coupon Bond Lattice - TEST ==========\n";
	std::cout << "=======================================================\n";

	indexedBDTCouponBond();
	indexedHLCouponBond();

	std::cout << "=======================================================\n";
}


void BDTCouponBond() {

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

	std::size_t year{ 180 };
	std::size_t periods{ discount_curve.size() - 1 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial,double,date> calibratedTree(fixingDates);

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


	double nominal{ 100.0 };
	double couponRate{ 0.05 };


	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial,double,date> bondTree(fixingDates);

	// populate coupon bond tree:
	createBinomialCouponBondLattice(bondTree, calibratedTree, nominal, couponRate, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Todays price of coupon bond (nominal=100, couponRate=0.05, coupon paid semiannualy):"
		<< bondTree.apex() << "\n";

}



void HLCouponBond() {
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


	double nominal{ 100.0 };
	double couponRate{ 0.05 };


	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// populate coupon bond tree:
	createBinomialCouponBondLattice(bondTree, calibratedTree, nominal, couponRate, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Todays price of coupon bond (nominal=100, couponRate=0.05, coupon paid semiannualy):"
		<< bondTree.apex() << "\n";

}


void testBinomialCouponBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "============= Coupon Bond Lattice - TEST ==============\n";
	std::cout << "=======================================================\n";

	BDTCouponBond();
	HLCouponBond();

	std::cout << "=======================================================\n";
}



#endif ///_LATTICE_EXAMPLES_CALIBRATE_PRICE_IR