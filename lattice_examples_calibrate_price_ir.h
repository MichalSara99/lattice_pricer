#pragma once
#if !defined(_LATTICE_EXAMPLES_CALIBRATE_PRICE_IR)
#define _LATTICE_EXAMPLES_CALIBRATE_PRICE_IR


#include"lattice_structure.h"
#include"lattice_model.h"
#include"lattice_miscellaneous.h"
#include"lattice_utility.h"
#include"lattice_miscellaneous.h"
#include"lattice_calibrator.h"
#include"lattice_bond_builders.h"
#include"lattice_model_params.h"
#include"lattice_product_builder.h"


#include<iostream>
#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>

using namespace boost::gregorian;

// ============================================================================================================
// =============================== Pure Discount Bond Tree ====================================================
// ============================================================================================================


void indexedBDTPureDiscountBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1,AssetClass::InterestRate,double> params;
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
	std::size_t maxPeriods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<size_t> couponDates;

	// Creating indexed bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(maxPeriods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(maxPeriods) << "\n";

}


void indexedBDTPureDiscountBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	// create params fro model:
	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::size_t maxPeriods{ discount_curve.size() - 2 };

	// build pure discount bond:
	auto discoundBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0)	
		.withPeriods(maxPeriods)
		.build();

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(discoundBond.periods());

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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

	// Creating indexed bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(discoundBond.periods());

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate bond tree:
	bb(bondTree, calibratedTree, bdt, discoundBond.nominal(), 0.0, 0.0, std::set<std::size_t>(), dt);

	std::cout << "indexed Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(maxPeriods) << "\n";

}


void indexedHLPureDiscountBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(periods);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<size_t> couponDates;

 
	// Creating indexedbond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(periods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed Ho-Lee coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";

}

void indexedHLPureDiscountBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	// build pure discount bond:
	auto discoundBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0)
		.withPeriods(periods)
		.build();

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(discoundBond.periods());

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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

	// Creating indexedbond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(discoundBond.periods());

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, discoundBond.nominal(), 0.0, 0.0, std::set<size_t>(), dt);

	std::cout << "indexed Ho-Lee coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";

}


void testIndexedBinomialPureDiscountBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "=== Indexed Pure Discount Bond Lattice - TEST =========\n";
	std::cout << "=======================================================\n";

	indexedBDTPureDiscountBond();
	indexedBDTPureDiscountBondWithProduct();
	indexedHLPureDiscountBond();
	indexedHLPureDiscountBondWithProduct();

	std::cout << "=======================================================\n";
}


void indexedBDTPureDiscountBond(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	const std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<size_t> couponDates;

	// Creating indexed bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(bondPeriods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";

}


void indexedBDTPureDiscountBondWithProduct(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	const std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// build pure discount bond product:
	auto discountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0).withPeriods(bondPeriods).build();

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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

	// Creating indexed bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(discountBond.periods());

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, discountBond.nominal(), 0.0,0.0, std::set<std::size_t>(), dt);

	std::cout << "indexed Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";

}


void indexedHLPureDiscountBond(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<size_t> couponDates;


	// Creating indexedbond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(bondPeriods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed Ho-Lee coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";

}


void indexedHLPureDiscountBondWithProduct(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// build pure discount bond product:
	auto discountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0).withPeriods(bondPeriods).build();

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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

	// Creating indexedbond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(discountBond.periods());

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, discountBond.nominal(), 0.0, 0.0, std::set<std::size_t>(), dt);

	std::cout << "indexed Ho-Lee coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";

}


void testIndexedBinomialPureDiscountBondLatticeVariable(std::size_t bondPeriods) {
	std::cout << "=======================================================\n";
	std::cout << "=== Indexed Pure Discount Bond Lattice - TEST =========\n";
	std::cout << "=======================================================\n";

	indexedBDTPureDiscountBond(bondPeriods);
	indexedBDTPureDiscountBondWithProduct(bondPeriods);
	indexedHLPureDiscountBond(bondPeriods);
	indexedHLPureDiscountBondWithProduct(bondPeriods);

	std::cout << "=======================================================\n";
}

void BDTPureDiscountBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial,double,date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<date> couponDates;


	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial,double,date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";

}


void BDTPureDiscountBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// build the product:
	auto discountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0)
		.withPeriods(periods)
		.build();

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, discountBond.nominal(), 0.0, 0.0, std::set<date>(), timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";

}



void HLPureDiscountBond() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<date> couponDates;


	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";

}


void HLPureDiscountBondWithProduct() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// build the product:
	auto biscountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0)
		.withPeriods(periods)
		.build();

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, biscountBond.nominal(), 0.0, 0.0, std::set<date>(), timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";

}


void testBinomialPureDicsountBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "======= Pure Disocunt Bond Lattice - TEST =============\n";
	std::cout << "=======================================================\n";

	BDTPureDiscountBond();
	BDTPureDiscountBondWithProduct();
	HLPureDiscountBond();
	HLPureDiscountBondWithProduct();

	std::cout << "=======================================================\n";
}

void BDTPureDiscountBond(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<date> couponDates;


	// Create fixing dates for bond tree
	std::set<date> bondFixingDates;
	bondFixingDates.emplace(today);
	for (std::size_t t = 1; t <= bondPeriods; ++t) {
		bondFixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(bondFixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";

}


void BDTPureDiscountBondWithProduct(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// create product
	auto discountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0).withPeriods(bondPeriods).build();

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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

	// Create fixing dates for bond tree
	std::set<date> bondFixingDates;
	bondFixingDates.emplace(today);
	for (std::size_t t = 1; t <= discountBond.periods(); ++t) {
		bondFixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(bondFixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, discountBond.nominal(), 0.0, 0.0, std::set<date>(), timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";

}


void HLPureDiscountBond(std::size_t bondPeriods) {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<date> couponDates;

	// Create fixing dates for bond tree
	std::set<date> bondFixingDates;
	bondFixingDates.emplace(today);
	for (std::size_t t = 1; t <= bondPeriods; ++t) {
		bondFixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(bondFixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";

}

void HLPureDiscountBondWithProduct(std::size_t bondPeriods) {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// build pure discount bond:
	auto discountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0).withPeriods(bondPeriods).build();

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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

	// Create fixing dates for bond tree
	std::set<date> bondFixingDates;
	bondFixingDates.emplace(today);
	for (std::size_t t = 1; t <= discountBond.periods(); ++t) {
		bondFixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(bondFixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, discountBond.nominal(), 0.0, 0.0, std::set<date>(), timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";

}

void testBinomialPureDicsountBondLatticeVariable(std::size_t bondPeriods) {
	std::cout << "=======================================================\n";
	std::cout << "======= Pure Disocunt Bond Lattice - TEST =============\n";
	std::cout << "=======================================================\n";

	BDTPureDiscountBond(bondPeriods);
	BDTPureDiscountBondWithProduct(bondPeriods);
	HLPureDiscountBond(bondPeriods);
	HLPureDiscountBondWithProduct(bondPeriods);

	std::cout << "=======================================================\n";
}

void indexedHWPureDiscountBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed,dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<std::size_t> couponDates;


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(periods, params.ReversionSpeed,dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, nominal, couponRate,lastCoupon,couponDates, dt);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";

}


void indexedHWPureDiscountBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// create pure discount bond:
	auto discountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0).withPeriods(periods).build();

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(discountBond.periods(), params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(periods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, discountBond.nominal(), 0.0, 0.0, std::set<std::size_t>(), dt);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";

}


void indexedBKPureDiscountBond() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<std::size_t> couponDates;


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(periods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, nominal, couponRate,lastCoupon,couponDates, dt);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";
}


void indexedBKPureDiscountBondWithProduct() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// create pure dicount bond product:
	auto discountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0).withPeriods(periods).build();

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(discountBond.periods(), params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(discountBond.periods(), params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, discountBond.nominal(), 0.0, 0.0, std::set<std::size_t>() , dt);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";
}


void testIndexedTrinomialPureDiscountBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "==== Indexed Pure Discount Bond Lattice - TEST ========\n";
	std::cout << "=======================================================\n";

	indexedHWPureDiscountBond();
	indexedHWPureDiscountBondWithProduct();
	indexedBKPureDiscountBond();
	indexedBKPureDiscountBondWithProduct();

	std::cout << "=======================================================\n";
}

void indexedHWPureDiscountBond(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(maxPeriods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<std::size_t> couponDates;


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(bondPeriods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";

}


void indexedHWPureDiscountBondWithProduct(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// create pure discount bond product:
	auto discountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0).withPeriods(bondPeriods).build();

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(maxPeriods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(discountBond.periods(), params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, discountBond.nominal(), 0.0, 0.0, std::set<std::size_t>(), dt);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";

}


void indexedBKPureDiscountBond(std::size_t bondPeriods) {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(maxPeriods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<std::size_t> couponDates;


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(bondPeriods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";
}

void indexedBKPureDiscountBondWithProduct(std::size_t bondPeriods) {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// create pure discount bond product:
	auto discountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0).withPeriods(bondPeriods).build();

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(maxPeriods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(discountBond.periods(), params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, discountBond.nominal(), 0.0,0.0, std::set<std::size_t>(), dt);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";
}


void testIndexedTrinomialPureDiscountBondLatticeVariable(std::size_t bondPeriods) {
	std::cout << "=======================================================\n";
	std::cout << "==== Indexed Pure Discount Bond Lattice - TEST ========\n";
	std::cout << "=======================================================\n";

	indexedHWPureDiscountBond(bondPeriods);
	indexedHWPureDiscountBondWithProduct(bondPeriods);
	indexedBKPureDiscountBond(bondPeriods);
	indexedBKPureDiscountBondWithProduct(bondPeriods);

	std::cout << "=======================================================\n";
}


void HWPureDiscountBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<date> couponDates;


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double,date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, nominal, couponRate,lastCoupon,couponDates, timeDeltas);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";

}


void HWPureDiscountBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// create discount bond product:
	auto discountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0).withPeriods(periods).build();

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, discountBond.nominal(), 0.0, 0.0, std::set<date>(), timeDeltas);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";

}


void BKPureDiscountBond() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<date> couponDates;


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, nominal, couponRate,lastCoupon,couponDates, timeDeltas);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";
}


void BKPureDiscountBondWithProduct() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	auto discountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0).withPeriods(periods).build();

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, discountBond.nominal(), 0.0, 0.0, std::set<date>(), timeDeltas);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(periods) << "\n";
}

void testTrinomialPureDiscountBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "========= Pure Discount Bond Lattice - TEST ===========\n";
	std::cout << "=======================================================\n";

	HWPureDiscountBond();
	HWPureDiscountBondWithProduct();
	BKPureDiscountBond();
	BKPureDiscountBondWithProduct();

	std::cout << "=======================================================\n";
}


void HWPureDiscountBond(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace_back(today + date_duration(daysInhalfYear*t));
		fixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	float daysInYear{ 365.0 };
	std::vector<double> timeDeltas(fixingDates.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fixingDates[i + 1] - fixingDates[i]).days() / daysInYear);
	}

	// Creating indexed lattice:
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<date> couponDates;


	// Prepare fixing dates for bond tree
	std::set<date> bondFixingDatesSet;
	bondFixingDatesSet.emplace(today);
	for (std::size_t t = 1; t <= bondPeriods; ++t) {
		bondFixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(bondFixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";

}


void HWPureDiscountBondWithProduct(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace_back(today + date_duration(daysInhalfYear*t));
		fixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	float daysInYear{ 365.0 };
	std::vector<double> timeDeltas(fixingDates.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fixingDates[i + 1] - fixingDates[i]).days() / daysInYear);
	}

	// Creating indexed lattice:
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	// create pure discount bond product:
	auto discountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0).withPeriods(bondPeriods).build();

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// Prepare fixing dates for bond tree
	std::set<date> bondFixingDatesSet;
	bondFixingDatesSet.emplace(today);
	for (std::size_t t = 1; t <= discountBond.periods(); ++t) {
		bondFixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(bondFixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, discountBond.nominal(), 0.0, 0.0, std::set<date>(), timeDeltas);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";

}


void BKPureDiscountBond(std::size_t bondPeriods) {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace_back(today + date_duration(daysInhalfYear*t));
		fixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	float daysInYear{ 365.0 };
	std::vector<double> timeDeltas(fixingDates.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fixingDates[i + 1] - fixingDates[i]).days() / daysInYear);
	}

	// Creating indexed lattice:
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.0 };
	double lastCoupon{ 0.0 };
	std::set<date> couponDates;


	// Prepare fixing dates for bond tree
	std::set<date> bondFixingDatesSet;
	bondFixingDatesSet.emplace(today);
	for (std::size_t t = 1; t <= bondPeriods; ++t) {
		bondFixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(bondFixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";
}

void BKPureDiscountBondWithProduct(std::size_t bondPeriods) {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::PureDiscountBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace_back(today + date_duration(daysInhalfYear*t));
		fixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	float daysInYear{ 365.0 };
	std::vector<double> timeDeltas(fixingDates.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fixingDates[i + 1] - fixingDates[i]).days() / daysInYear);
	}

	// Creating indexed lattice:
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	// create discount bond product:
	auto discountBond = PureDiscountBondBuilder<double>()
		.withNominal(1.0).withPeriods(bondPeriods).build();

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	// Prepare fixing dates for bond tree
	std::set<date> bondFixingDatesSet;
	bondFixingDatesSet.emplace(today);
	for (std::size_t t = 1; t <= discountBond.periods(); ++t) {
		bondFixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(bondFixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, discountBond.nominal(), 0.0, 0.0, std::set<date>(), timeDeltas);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Replicated price of pure discount bond:"
		<< bondTree.apex() << "\n";
	std::cout << "Market price of the pure discount bond:"
		<< discount_curve.at(bondPeriods) << "\n";
}

void testTrinomialPureDiscountBondLatticeVariable(std::size_t bondPeriods) {
	std::cout << "=======================================================\n";
	std::cout << "========= Pure Discount Bond Lattice - TEST ===========\n";
	std::cout << "=======================================================\n";

	HWPureDiscountBond(bondPeriods);
	HWPureDiscountBondWithProduct(bondPeriods);
	BKPureDiscountBond(bondPeriods);
	BKPureDiscountBondWithProduct(bondPeriods);

	std::cout << "=======================================================\n";
}


// ============================================================================================================
// ====================================== Coupon Bond Tree ====================================================
// ============================================================================================================


void indexedBDTCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1,AssetClass::InterestRate,double> params;
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
	std::size_t maxPeriods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < maxPeriods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexed bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(maxPeriods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";
}


void indexedBDTCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::size_t maxPeriods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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

	// build coupon data:
	std::map<std::size_t,double> couponDates;
	for (std::size_t t = 0; t < maxPeriods; ++t) {
		couponDates.emplace(t + 1, 0.05);
	}

	// build coupon bond product:
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withPeriods(maxPeriods)
		.withCouponData(couponDates)
		.withLastCoupon(1.0*0.05*dt)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(maxPeriods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";
}


void indexedHLCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(periods);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < periods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexedbond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(periods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed Ho-Lee coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";

}

void indexedHLCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(periods);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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

	// create coupon data:
	std::map<std::size_t,double> couponData;
	for (std::size_t t = 0; t < periods; ++t) {
		couponData.emplace(t + 1,0.05);
	}

	// build coupon bond product:
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withLastCoupon(1.0*dt*0.05)
		.withCouponData(couponData)
		.withPeriods(periods)
		.build();


	// Creating indexedbond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(periods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed Ho-Lee coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";

}


void testIndexedBinomialCouponBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "========== Indexed Coupon Bond Lattice - TEST =========\n";
	std::cout << "=======================================================\n";

	indexedBDTCouponBond();
	indexedBDTCouponBondWithProduct();
	indexedHLCouponBond();
	indexedHLCouponBondWithProduct();

	std::cout << "=======================================================\n";
}



void indexedBDTCouponBond(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	const std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < bondPeriods; ++t) {
		couponDates.emplace(t + 1);
	}

	// Creating indexed bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(bondPeriods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";

}


void indexedBDTCouponBondWithProduct(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	const std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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

	// construct coupon data:
	std::map<std::size_t,double> couponData;
	for (std::size_t t = 0; t < bondPeriods; ++t) {
		couponData.emplace(t + 1, 0.05);
	}

	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withLastCoupon(1.0*dt*0.05)
		.withPeriods(bondPeriods)
		.withCouponData(couponData)
		.build();


	// Creating indexed bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(bondPeriods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";

}


void indexedHLCouponBond(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < bondPeriods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexedbond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(bondPeriods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed Ho-Lee coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}


void indexedHLCouponBondWithProduct(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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

	std::map<std::size_t,double> couponData;
	for (std::size_t t = 0; t < bondPeriods; ++t) {
		couponData.emplace(t + 1, 0.05);
	}

	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withLastCoupon(1.0*dt*0.05)
		.withPeriods(bondPeriods)
		.withCouponData(couponData)
		.build();


	// Creating indexedbond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(bondPeriods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed Ho-Lee coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}

void testIndexedBinomialCouponBondLatticeVariable(std::size_t bondPeriods) {
	std::cout << "=======================================================\n";
	std::cout << "========== Indexed Coupon Bond Lattice - TEST =========\n";
	std::cout << "=======================================================\n";

	indexedBDTCouponBond(bondPeriods);
	indexedBDTCouponBondWithProduct(bondPeriods);
	indexedHLCouponBond(bondPeriods);
	indexedHLCouponBondWithProduct(bondPeriods);

	std::cout << "=======================================================\n";
}



void BDTCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / daysInYear)*0.05 };
	std::set<date> couponDates;
	couponDates.emplace(today);
	for (std::size_t t = 1; t < periods; ++t) {
		couponDates.emplace(today + date_duration(t*year));
	}


	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";
}


void BDTCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	std::map<date,double> couponData;
	couponData.emplace(today, 0.05);
	for (std::size_t t = 1; t < periods; ++t) {
		couponData.emplace(today + date_duration(t*year), 0.05);
	}

	// build coupon bond product:
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withLastCoupon(1.0*(180. / daysInYear)*0.05)
		.withPeriods(periods)
		.withCouponData(couponData).build();


	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";
}


void HLCouponBond() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / daysInYear)*0.05 };
	std::set<date> couponDates;
	couponDates.emplace(today);
	for (std::size_t t = 1; t < periods; ++t) {
		couponDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}


void HLCouponBondWithProduct() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	std::map<date,double> couponData;
	couponData.emplace(today, 0.05);
	for (std::size_t t = 1; t < periods; ++t) {
		couponData.emplace(today + date_duration(t*year),0.05);
	}

	// creat coupon bond product:
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*(180. / 365.)*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}


void testBinomialCouponBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "============== Coupon Bond Lattice - TEST =============\n";
	std::cout << "=======================================================\n";

	BDTCouponBond();
	BDTCouponBondWithProduct();
	HLCouponBond();
	HLCouponBondWithProduct();

	std::cout << "=======================================================\n";
}

void BDTCouponBond(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	couponDates.emplace(today);
	for (std::size_t t = 1; t < bondPeriods; ++t) {
		couponDates.emplace(today + date_duration(t*year));
	}


	// Create fixing dates for bond tree
	std::set<date> bondFixingDates;
	bondFixingDates.emplace(today);
	for (std::size_t t = 1; t <= bondPeriods; ++t) {
		bondFixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(bondFixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}


void BDTCouponBondWithProduct(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	std::map<date,double> couponData;
	couponData.emplace(today, 0.05);
	for (std::size_t t = 1; t < bondPeriods; ++t) {
		couponData.emplace(today + date_duration(t*year),0.05);
	}

	// build coupon bond product:
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(bondPeriods)
		.withLastCoupon(1.0*(year / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();

	// Create fixing dates for bond tree
	std::set<date> bondFixingDates;
	bondFixingDates.emplace(today);
	for (std::size_t t = 1; t <= bondPeriods; ++t) {
		bondFixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(bondFixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}



void HLCouponBond(std::size_t bondPeriods) {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	couponDates.emplace(today);
	for (std::size_t t = 1; t < bondPeriods; ++t) {
		couponDates.emplace(today + date_duration(t*year));
	}

	// Create fixing dates for bond tree
	std::set<date> bondFixingDates;
	bondFixingDates.emplace(today);
	for (std::size_t t = 1; t <= bondPeriods; ++t) {
		bondFixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(bondFixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}

void HLCouponBondWithProduct(std::size_t bondPeriods) {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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

	std::map<date,double> couponData;
	couponData.emplace(today, 0.05);
	for (std::size_t t = 1; t < bondPeriods; ++t) {
		couponData.emplace(today + date_duration(t*year),0.05);
	}

	// create coupon bond product:
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(bondPeriods)
		.withLastCoupon(1.0*(year / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();

	// Create fixing dates for bond tree
	std::set<date> bondFixingDates;
	bondFixingDates.emplace(today);
	for (std::size_t t = 1; t <= bondPeriods; ++t) {
		bondFixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(bondFixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}


void testBinomialCouponBondLatticeVariable(std::size_t bondPeriods) {
	std::cout << "=======================================================\n";
	std::cout << "=========== Coupon Bond Lattice - TEST ================\n";
	std::cout << "=======================================================\n";

	BDTCouponBond(bondPeriods);
	BDTCouponBondWithProduct(bondPeriods);
	HLCouponBond(bondPeriods);
	HLCouponBondWithProduct(bondPeriods);

	std::cout << "=======================================================\n";
}


void indexedHWCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < periods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(periods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}


void indexedHWCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	std::map<std::size_t,double> couponData;
	for (std::size_t t = 0; t < periods; ++t) {
		couponData.emplace(t + 1, 0.05);
	}

	// build coupon bond product:
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withLastCoupon(1.0*dt*0.05)
		.withPeriods(periods)
		.withCouponData(couponData)
		.build();


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(periods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}

void indexedBKCouponBond() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < periods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(periods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";
}


void indexedBKCouponBondWithProduct() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	std::map<std::size_t,double> couponData;
	for (std::size_t t = 0; t < periods; ++t) {
		couponData.emplace(t + 1,0.05);
	}

	// build coupon bond product:
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withLastCoupon(1.0*dt*0.05)
		.withPeriods(periods)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(periods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";
}

void testIndexedTrinomialCouponBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "========= Indexed Coupon Bond Lattice - TEST ==========\n";
	std::cout << "=======================================================\n";

	indexedHWCouponBond();
	indexedHWCouponBondWithProduct();
	indexedBKCouponBond();
	indexedBKCouponBondWithProduct();

	std::cout << "=======================================================\n";
}


void indexedHWCouponBond(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(maxPeriods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < bondPeriods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(bondPeriods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}

void indexedHWCouponBondWithProduct(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(maxPeriods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	std::map<std::size_t,double> couponData;
	for (std::size_t t = 0; t < bondPeriods; ++t) {
		couponData.emplace(t + 1, 0.05);
	}

	// create coupon bond product:
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withLastCoupon(1.0*dt*0.05)
		.withPeriods(bondPeriods)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(couponBond.periods(), params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}



void indexedBKCouponBond(std::size_t bondPeriods) {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;


	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(maxPeriods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < bondPeriods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(bondPeriods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";
}

void indexedBKCouponBondWithProduct(std::size_t bondPeriods) {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;


	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(maxPeriods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	std::map<std::size_t,double> couponData;
	for (std::size_t t = 0; t < bondPeriods; ++t) {
		couponData.emplace(t + 1,0.05);
	}

	// create coupon bond product:
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withPeriods(bondPeriods)
		.withLastCoupon(1.0*dt*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(couponBond.periods(), params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";
}

void testIndexedTrinomialCouponBondLatticeVariable(std::size_t bondPeriods) {
	std::cout << "=======================================================\n";
	std::cout << "========= Indexed Coupon Bond Lattice - TEST ==========\n";
	std::cout << "=======================================================\n";

	indexedHWCouponBond(bondPeriods);
	indexedHWCouponBondWithProduct(bondPeriods);
	indexedBKCouponBond(bondPeriods);
	indexedBKCouponBondWithProduct(bondPeriods);

	std::cout << "=======================================================\n";
}



void HWCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;


	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	for (std::size_t t = 1; t < periods; ++t) {
		couponDates.emplace(today + date_duration(t*daysInhalfYear));
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}


void HWCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;


	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	std::map<date,double> couponData;
	for (std::size_t t = 1; t < periods; ++t) {
		couponData.emplace(today + date_duration(t*daysInhalfYear), 0.05);
	}

	// create coupn bond product:
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*(daysInhalfYear / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}


void BKCouponBond() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;


	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	for (std::size_t t = 1; t < periods; ++t) {
		couponDates.emplace(today + date_duration(t*daysInhalfYear));
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";
}

void BKCouponBondWithProduct() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;


	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	std::map<date,double> couponData;
	for (std::size_t t = 1; t < periods; ++t) {
		couponData.emplace(today + date_duration(t*daysInhalfYear),0.05);
	}

	// create coupon bond product:
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*(daysInhalfYear / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";
}

void testTrinomialCouponBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "============ Coupon Bond Lattice - TEST ===============\n";
	std::cout << "=======================================================\n";

	HWCouponBond();
	HWCouponBondWithProduct();
	BKCouponBond();
	BKCouponBondWithProduct();

	std::cout << "=======================================================\n";
}


void HWCouponBond(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;


	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace_back(today + date_duration(daysInhalfYear*t));
		fixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	float daysInYear{ 365.0 };
	std::vector<double> timeDeltas(fixingDates.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fixingDates[i + 1] - fixingDates[i]).days() / daysInYear);
	}

	// Creating indexed lattice:
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	for (std::size_t t = 1; t < bondPeriods; ++t) {
		couponDates.emplace(today + date_duration(t*daysInhalfYear));
	}


	// Prepare fixing dates for bond tree
	std::set<date> bondFixingDatesSet;
	bondFixingDatesSet.emplace(today);
	for (std::size_t t = 1; t <= bondPeriods; ++t) {
		bondFixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(bondFixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}


void HWCouponBondWithProduct(std::size_t bondPeriods) {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;


	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace_back(today + date_duration(daysInhalfYear*t));
		fixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	float daysInYear{ 365.0 };
	std::vector<double> timeDeltas(fixingDates.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fixingDates[i + 1] - fixingDates[i]).days() / daysInYear);
	}

	// Creating indexed lattice:
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	std::map<date,double> couponData;
	for (std::size_t t = 1; t < bondPeriods; ++t) {
		couponData.emplace(today + date_duration(t*daysInhalfYear),0.05);
	}

	// create coupon bond:
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(bondPeriods)
		.withLastCoupon(1.0*(daysInhalfYear / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();

	// Prepare fixing dates for bond tree
	std::set<date> bondFixingDatesSet;
	bondFixingDatesSet.emplace(today);
	for (std::size_t t = 1; t <= bondPeriods; ++t) {
		bondFixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(bondFixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

}


void BKCouponBond(std::size_t bondPeriods) {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;


	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace_back(today + date_duration(daysInhalfYear*t));
		fixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	float daysInYear{ 365.0 };
	std::vector<double> timeDeltas(fixingDates.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fixingDates[i + 1] - fixingDates[i]).days() / daysInYear);
	}

	// Creating indexed lattice:
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	for (std::size_t t = 1; t < bondPeriods; ++t) {
		couponDates.emplace(today + date_duration(t*daysInhalfYear));
	}

	// Prepare fixing dates for bond tree
	std::set<date> bondFixingDatesSet;
	bondFixingDatesSet.emplace(today);
	for (std::size_t t = 1; t <= bondPeriods; ++t) {
		bondFixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(bondFixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";
}


void BKCouponBondWithProduct(std::size_t bondPeriods) {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;


	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t daysInhalfYear{ 180 };
	auto today = date(day_clock::local_day());
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t maxPeriods{ discount_curve.size() - 2 };
	LASSERT(bondPeriods <= maxPeriods, "passed periods must be less or equal to maxPeriods.");
	for (std::size_t t = 1; t <= maxPeriods; ++t) {
		fixingDates.emplace_back(today + date_duration(daysInhalfYear*t));
		fixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	float daysInYear{ 365.0 };
	std::vector<double> timeDeltas(fixingDates.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fixingDates[i + 1] - fixingDates[i]).days() / daysInYear);
	}

	// Creating indexed lattice:
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	std::map<date,double> couponData;
	for (std::size_t t = 1; t < bondPeriods; ++t) {
		couponData.emplace(today + date_duration(t*daysInhalfYear),0.05);
	}

	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(bondPeriods)
		.withLastCoupon(1.0*(daysInhalfYear / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();

	// Prepare fixing dates for bond tree
	std::set<date> bondFixingDatesSet;
	bondFixingDatesSet.emplace(today);
	for (std::size_t t = 1; t <= bondPeriods; ++t) {
		bondFixingDatesSet.emplace(today + date_duration(daysInhalfYear*t));
	}

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(bondFixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, bondPeriods + 1);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";
}


void testTrinomialCouponBondLatticeVariable(std::size_t bondPeriods) {
	std::cout << "=======================================================\n";
	std::cout << "============ Coupon Bond Lattice - TEST ===============\n";
	std::cout << "=======================================================\n";

	HWCouponBond(bondPeriods);
	HWCouponBondWithProduct(bondPeriods);
	BKCouponBond(bondPeriods);
	BKCouponBondWithProduct(bondPeriods);

	std::cout << "=======================================================\n";
}


// ============================================================================================================
// ====================================== Option On Bond Tree =================================================
// ============================================================================================================



void indexedBDTEuropeanOptionOnCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1,AssetClass::InterestRate,double> params;
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
	std::size_t maxPeriods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < maxPeriods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexed bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(maxPeriods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> optionTree(maxPeriods);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bdt, dt, call_payoff);

	std::cout << "indexed Black-Derman-Toy option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";
}



void indexedBDTEuropeanOptionOnCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::size_t maxPeriods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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

	std::map<std::size_t, double> couponData;
	for (std::size_t t = 0; t < maxPeriods; ++t) {
		couponData.emplace(t + 1, 0.05);
	}

	// create coupon bond product
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withPeriods(maxPeriods)
		.withNominal(1.0)
		.withLastCoupon(1.0*dt*0.05)
		.withCouponData(couponData)
		.build();


	// Creating indexed bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(couponBond.periods());

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";

	// Create option on coupon bond
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(maxPeriods).withStrike(1.02).build();

	// Creating indexed option on bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> optionTree(option.periods());


	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bdt, dt, call_payoff);

	std::cout << "indexed Black-Derman-Toy option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";
}

void indexedHLEuropeanOptionOnCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(periods);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < periods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexedbond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(periods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed Ho-Lee coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> optionTree(periods);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hlm, dt, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";

}

void indexedHLEuropeanOptionOnCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(periods);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	std::map<std::size_t,double> couponData;
	for (std::size_t t = 0; t < periods; ++t) {
		couponData.emplace(t + 1,0.05);
	}

	// create coupon bond product: 
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*dt*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexedbond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(couponBond.periods());

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed Ho-Lee coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";

	// create option on coupon bond product: 
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(periods).withStrike(1.02).build();


	// Creating indexed option on bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> optionTree(option.periods());

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hlm, dt, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";

}


void testIndexedBinomialEuropeanOptionOnCouponBondLattice() {
	std::cout << "=================================================================\n";
	std::cout << "========== Indexed Option on Coupon Bond Lattice - TEST =========\n";
	std::cout << "=================================================================\n";

	indexedBDTEuropeanOptionOnCouponBond();
	indexedBDTEuropeanOptionOnCouponBondWithProduct();
	indexedHLEuropeanOptionOnCouponBond();
	indexedHLEuropeanOptionOnCouponBondWithProduct();

	std::cout << "=======================================================\n";
}




void indexedBDTAmericanOptionOnCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::size_t maxPeriods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < maxPeriods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexed bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(maxPeriods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> optionTree(maxPeriods);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};
	// create american adjuster:
	auto american_adjuster = [&call_payoff](double &optionPrice,double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bdt, dt, call_payoff, american_adjuster);

	std::cout << "indexed Black-Derman-Toy option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";
}



void indexedBDTAmericanOptionOnCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::size_t maxPeriods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(maxPeriods);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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

	std::map<std::size_t, double> couponData;
	for (std::size_t t = 0; t < maxPeriods; ++t) {
		couponData.emplace(t + 1, 0.05);
	}

	// create coupon bond product
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withPeriods(maxPeriods)
		.withNominal(1.0)
		.withLastCoupon(1.0*dt*0.05)
		.withCouponData(couponData)
		.build();


	// Creating indexed bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(couponBond.periods());

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";

	// Create option on coupon bond
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(maxPeriods).withStrike(1.02).build();

	// Creating indexed option on bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> optionTree(option.periods());


	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// create american adjuster:
	auto american_adjuster = [&call_payoff](double &optionPrice, double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bdt, dt, call_payoff, american_adjuster);

	std::cout << "indexed Black-Derman-Toy option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";
}

void indexedHLAmericanOptionOnCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(periods);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < periods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexedbond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(periods);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed Ho-Lee coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> optionTree(periods);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// create american adjuster:
	auto american_adjuster = [&call_payoff](double &optionPrice, double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hlm, dt, call_payoff, american_adjuster);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";

}

void indexedHLAmericanOptionOnCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> calibratedTree(periods);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	std::map<std::size_t, double> couponData;
	for (std::size_t t = 0; t < periods; ++t) {
		couponData.emplace(t + 1, 0.05);
	}

	// create coupon bond product: 
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*dt*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexedbond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> bondTree(couponBond.periods());

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed Ho-Lee coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modelled price of coupon bond:"
		<< bondTree.apex() << "\n";

	// create option on coupon bond product: 
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(periods).withStrike(1.02).build();


	// Creating indexed option on bond lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> optionTree(option.periods());

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// create american adjuster:
	auto american_adjuster = [&call_payoff](double &optionPrice, double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hlm, dt, call_payoff, american_adjuster);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";

}


void testIndexedBinomialAmericanOptionOnCouponBondLattice() {
	std::cout << "=================================================================\n";
	std::cout << "========== Indexed Option on Coupon Bond Lattice - TEST =========\n";
	std::cout << "=================================================================\n";

	indexedBDTAmericanOptionOnCouponBond();
	indexedBDTAmericanOptionOnCouponBondWithProduct();
	indexedHLAmericanOptionOnCouponBond();
	indexedHLAmericanOptionOnCouponBondWithProduct();

	std::cout << "=======================================================\n";
}


void BDTEuropeanOptionOnCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	for (std::size_t t = 1; t < periods; ++t) {
		couponDates.emplace(today + date_duration(t*year));
	}


	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> optionTree(fixingDates);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bdt, timeDeltas, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";
}


void BDTEuropeanOptionOnCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	std::map<date,double> couponData;
	for (std::size_t t = 1; t < periods; ++t) {
		couponData.emplace(today + date_duration(t*year),0.05);
	}

	// creta coupon bond product
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*(year / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// create option on coupon bond product
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(periods)
		.withStrike(1.02)
		.build();

	// Creating indexed option on bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> optionTree(fixingDates);

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bdt, timeDeltas, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";
}


void HLEuropeanOptionOnCouponBond() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;
	//option.Strike = 1.02;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	for (std::size_t t = 1; t < periods; ++t) {
		couponDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> optionTree(fixingDates);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hlm, timeDeltas, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";
}

void HLEuropeanOptionOnCouponBondWithProduct() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;
	//option.Strike = 1.02;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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

	std::map<date,double> couponData;
	for (std::size_t t = 1; t < periods; ++t) {
		couponData.emplace(today + date_duration(t*year),0.05);
	}

	// create coupon bond product:
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*(year / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();


	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> optionTree(fixingDates);

	// create option on coupon bond product:
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(periods)
		.withStrike(1.02)
		.build();

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hlm, timeDeltas, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";
}


void testBinomialEuropeanOptionOnCouponBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "== European Option on Coupon Bond Lattice - TEST ======\n";
	std::cout << "=======================================================\n";

	BDTEuropeanOptionOnCouponBond();
	BDTEuropeanOptionOnCouponBondWithProduct();
	HLEuropeanOptionOnCouponBond();
	HLEuropeanOptionOnCouponBondWithProduct();

	std::cout << "=======================================================\n";
}

void BDTAmericanOptionOnCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	for (std::size_t t = 1; t < periods; ++t) {
		couponDates.emplace(today + date_duration(t*year));
	}


	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> optionTree(fixingDates);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// create american adjuster:
	auto american_adjuster = [&call_payoff](double &optionPrice, double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bdt, timeDeltas, call_payoff, american_adjuster);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";
}


void BDTAmericanOptionOnCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackDermanToyModel<> bdt(params);

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


	std::map<date, double> couponData;
	for (std::size_t t = 1; t < periods; ++t) {
		couponData.emplace(today + date_duration(t*year), 0.05);
	}

	// creta coupon bond product
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*(year / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bdt, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// create option on coupon bond product
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(periods)
		.withStrike(1.02)
		.build();

	// Creating indexed option on bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> optionTree(fixingDates);

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// create american adjuster:
	auto american_adjuster = [&call_payoff](double &optionPrice, double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bdt, timeDeltas, call_payoff, american_adjuster);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";
}


void HLAmericanOptionOnCouponBond() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;
	//option.Strike = 1.02;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	for (std::size_t t = 1; t < periods; ++t) {
		couponDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> optionTree(fixingDates);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// create american adjuster:
	auto american_adjuster = [&call_payoff](double &optionPrice, double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hlm, timeDeltas, call_payoff, american_adjuster);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";
}

void HLAmericanOptionOnCouponBondWithProduct() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
	params.Volatility = 0.005;
	//option.Strike = 1.02;

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
	// the discount factor 1.0 is not used in calibration therefore
	std::size_t periods{ discount_curve.size() - 2 };
	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t*year));
	}

	// Creating indexed lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> calibratedTree(fixingDates);

	// Create Blac-Derman-Toy model:
	lattice_model::HoLeeModel<> hlm(params);

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

	std::map<date, double> couponData;
	for (std::size_t t = 1; t < periods; ++t) {
		couponData.emplace(today + date_duration(t*year), 0.05);
	}

	// create coupon bond product:
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*(year / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();


	// Creating indexed bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> bondTree(fixingDates);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hlm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "Black-Derman-Toy coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> optionTree(fixingDates);

	// create option on coupon bond product:
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(periods)
		.withStrike(1.02)
		.build();

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// create american adjuster:
	auto american_adjuster = [&call_payoff](double &optionPrice, double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Binomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hlm, timeDeltas, call_payoff, american_adjuster);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";
}




void testBinomialAmericanOptionOnCouponBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "== American Option on Coupon Bond Lattice - TEST ======\n";
	std::cout << "=======================================================\n";

	BDTAmericanOptionOnCouponBond();
	BDTAmericanOptionOnCouponBondWithProduct();
	HLAmericanOptionOnCouponBond();
	HLAmericanOptionOnCouponBondWithProduct();

	std::cout << "=======================================================\n";
}


void indexedHWEuropeanOptionOnCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;


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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < periods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(periods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> optionTree(periods, params.ReversionSpeed, dt);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hwm, dt, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";

}


void indexedHWEuropeanOptionOnCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;


	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	std::map<std::size_t,double> couponData;
	for (std::size_t t = 0; t < periods; ++t) {
		couponData.emplace(t + 1, 0.05);
	}

	// Create coupon bond product:
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*dt*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(couponBond.periods(), params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Create option on coupon bond product:
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(periods)
		.withStrike(1.02)
		.build();

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> optionTree(option.periods(), params.ReversionSpeed, dt);

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hwm, dt, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";

}


void indexedBKEuropeanOptionOnCouponBond() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < periods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(periods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> optionTree(periods, params.ReversionSpeed, dt);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bkm, dt, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";
}

void indexedBKEuropeanOptionOnCouponBondWithProduct() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	std::map<std::size_t,double> couponData;
	for (std::size_t t = 0; t < periods; ++t) {
		couponData.emplace(t + 1, 0.05);
	}

	// create coupon bond product:
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*dt*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(periods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// create option on coupon bond product:
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(periods)
		.withStrike(1.02)
		.build();

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> optionTree(option.periods(), params.ReversionSpeed, dt);

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bkm, dt, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";
}


void testIndexedTrinomialEuropeanOptionOnCouponBondLattice() {
	std::cout << "=========================================================\n";
	std::cout << "= Indexed European Option on Coupon Bond Lattice - TEST =\n";
	std::cout << "=========================================================\n";

	indexedHWEuropeanOptionOnCouponBond();
	indexedHWEuropeanOptionOnCouponBondWithProduct();
	indexedBKEuropeanOptionOnCouponBond();
	indexedBKEuropeanOptionOnCouponBondWithProduct();

	std::cout << "=========================================================\n";
}

void indexedHWAmericanOptionOnCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;


	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < periods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(periods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> optionTree(periods, params.ReversionSpeed, dt);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	auto american_adjuster = [&call_payoff](double &optionPrice,double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hwm, dt, call_payoff, american_adjuster);

	std::cout << "indexed Hull-White option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";

}


void indexedHWAmericanOptionOnCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;


	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	std::map<std::size_t, double> couponData;
	for (std::size_t t = 0; t < periods; ++t) {
		couponData.emplace(t + 1, 0.05);
	}

	// Create coupon bond product:
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*dt*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(couponBond.periods(), params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Create option on coupon bond product:
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(periods)
		.withStrike(1.02)
		.build();

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> optionTree(option.periods(), params.ReversionSpeed, dt);

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};
	auto american_adjuster = [&call_payoff](double &optionPrice, double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hwm, dt, call_payoff, american_adjuster);

	std::cout << "indexed Hull-White option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";

}


void indexedBKAmericanOptionOnCouponBond() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*dt*0.05 };
	std::set<std::size_t> couponDates;
	for (std::size_t t = 0; t < periods; ++t) {
		couponDates.emplace(t + 1);
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(periods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, nominal, couponRate, lastCoupon, couponDates, dt);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> optionTree(periods, params.ReversionSpeed, dt);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	auto american_adjuster = [&call_payoff](double &optionPrice, double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bkm, dt, call_payoff, american_adjuster);

	std::cout << "indexed Black-Karasinski option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";
}

void indexedBKAmericanOptionOnCouponBondWithProduct() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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

	std::size_t periods{ discount_curve.size() - 2 };

	// Creating indexed lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> calibratedTree(periods, params.ReversionSpeed, dt);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		std::size_t, double, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, dt);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	std::map<std::size_t, double> couponData;
	for (std::size_t t = 0; t < periods; ++t) {
		couponData.emplace(t + 1, 0.05);
	}

	// create coupon bond product:
	auto couponBond = CouponBondBuilder<double, std::size_t>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*dt*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> bondTree(periods, params.ReversionSpeed, dt);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, double, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), dt);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// create option on coupon bond product:
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(periods)
		.withStrike(1.02)
		.build();

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingIndexedLattice<double> optionTree(option.periods(), params.ReversionSpeed, dt);

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	auto american_adjuster = [&call_payoff](double &optionPrice, double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, double> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bkm, dt, call_payoff, american_adjuster);

	std::cout << "indexed Black-Karasinski option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";
}


void testIndexedTrinomialAmericanOptionOnCouponBondLattice() {
	std::cout << "=========================================================\n";
	std::cout << "= Indexed American Option on Coupon Bond Lattice - TEST =\n";
	std::cout << "=========================================================\n";

	indexedHWAmericanOptionOnCouponBond();
	indexedHWAmericanOptionOnCouponBondWithProduct();
	indexedBKAmericanOptionOnCouponBond();
	indexedBKAmericanOptionOnCouponBondWithProduct();

	std::cout << "=========================================================\n";
}

void HWEuropeanOptionOnCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	for (std::size_t t = 1; t < periods; ++t) {
		couponDates.emplace(today + date_duration(t*daysInhalfYear));
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingLattice<double,date> optionTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hwm, timeDeltas, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";

}

void HWEuropeanOptionOnCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	std::map<date,double> couponData;
	for (std::size_t t = 1; t < periods; ++t) {
		couponData.emplace(today + date_duration(t*daysInhalfYear),0.05);
	}
	
	// create coupon bond product:
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*(daysInhalfYear / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> optionTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// create option on coupon bond product:
	auto option = OptionOnCouponBondBuilder<double>()
		.withStrike(1.02).withPeriods(periods).build();

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hwm, timeDeltas, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";

}


void BKEuropeanOptionOnCouponBond() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	for (std::size_t t = 1; t < periods; ++t) {
		couponDates.emplace(today + date_duration(t*daysInhalfYear));
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> optionTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bkm, timeDeltas, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";
}

void BKEuropeanOptionOnCouponBondWithProduct() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	std::map<date,double> couponData;
	for (std::size_t t = 1; t < periods; ++t) {
		couponData.emplace(today + date_duration(t*daysInhalfYear),0.05);
	}

	// create coupon bond product:
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*(daysInhalfYear / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> optionTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// create option on coupon bond product:
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(periods)
		.withStrike(1.02)
		.build();

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bkm, timeDeltas, call_payoff);

	std::cout << "indexed Ho-Lee option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of option on coupon bond:"
		<< optionTree.apex() << "\n";
}

void testTrinomialEuropeanOptionOnCouponBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "=== European Option On Coupon Bond Lattice - TEST =====\n";
	std::cout << "=======================================================\n";

	HWEuropeanOptionOnCouponBond();
	HWEuropeanOptionOnCouponBondWithProduct();
	BKEuropeanOptionOnCouponBond();
	BKEuropeanOptionOnCouponBondWithProduct();

	std::cout << "=======================================================\n";
}



void HWAmericanOptionOnCouponBond() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	for (std::size_t t = 1; t < periods; ++t) {
		couponDates.emplace(today + date_duration(t*daysInhalfYear));
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> optionTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	auto american_adjuster = [&call_payoff](double &optionPrice, double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hwm, timeDeltas, call_payoff, american_adjuster);

	std::cout << "indexed Hull-White option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";

}

void HWAmericanOptionOnCouponBondWithProduct() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::HullWhiteModel<> hwm(params);

	std::cout << "\nModel name: " << decltype(hwm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, hwm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	std::map<date, double> couponData;
	for (std::size_t t = 1; t < periods; ++t) {
		couponData.emplace(today + date_duration(t*daysInhalfYear), 0.05);
	}

	// create coupon bond product:
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*(daysInhalfYear / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, hwm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "indexed " << decltype(hwm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> optionTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// create option on coupon bond product:
	auto option = OptionOnCouponBondBuilder<double>()
		.withStrike(1.02).withPeriods(periods).build();

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};
	auto american_adjuster = [&call_payoff](double &optionPrice, double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, hwm, timeDeltas, call_payoff, american_adjuster);

	std::cout << "indexed Hull-White option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";

}


void BKAmericanOptionOnCouponBond() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);


	double nominal{ 1.0 };
	double couponRate{ 0.05 };
	double lastCoupon{ nominal*(180. / 365.)*0.05 };
	std::set<date> couponDates;
	for (std::size_t t = 1; t < periods; ++t) {
		couponDates.emplace(today + date_duration(t*daysInhalfYear));
	}


	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, nominal, couponRate, lastCoupon, couponDates, timeDeltas);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> optionTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create payoff fro the option	
	double K = 1.02;
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	auto american_adjuster = [&call_payoff](double &optionPrice, double bondPrice) {
		optionPrice = std::max(optionPrice, call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bkm, timeDeltas, call_payoff, american_adjuster);

	std::cout << "indexed Black-Karasinski option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";
}

void BKAmericanOptionOnCouponBondWithProduct() {
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_product_builder::CouponBondBuilder;
	using lattice_product_builder::OptionOnCouponBondBuilder;

	ModelParams<1, AssetClass::InterestRate, double> params;
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
	std::set<date> fixingDatesSet;
	std::vector<date> fixingDates;

	fixingDatesSet.emplace(today);
	fixingDates.emplace_back(today);
	std::size_t periods{ discount_curve.size() - 2 };
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
	lattice_structure::MeanRevertingLattice<double, date> calibratedTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// Create Blac-Derman-Toy model:
	lattice_model::BlackKarasinskiModel<> bkm(params);

	std::cout << "\nModel name: " << decltype(bkm)::name() << "\n";

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		date, std::vector<double>, std::vector<double>> hw_calibrator;
	// Instantiate Calibrator and launch it:
	hw_calibrator calibrator(discount_curve);
	auto result = calibrator(calibratedTree, bkm, timeDeltas);

	// Print the part of generated lattice:
	auto first = calibratedTree.begin();
	auto last = std::next(first, 17);
	lattice_utility::print(calibratedTree, first, last);

	std::map<date, double> couponData;
	for (std::size_t t = 1; t < periods; ++t) {
		couponData.emplace(today + date_duration(t*daysInhalfYear), 0.05);
	}

	// create coupon bond product:
	auto couponBond = CouponBondBuilder<double, date>()
		.withNominal(1.0)
		.withPeriods(periods)
		.withLastCoupon(1.0*(daysInhalfYear / daysInYear)*0.05)
		.withCouponData(couponData)
		.build();

	// Creating indexed bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> bondTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// typedef BondBuilder:
	typedef lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>, double> bond_builder;
	bond_builder bb;

	// populate coupon bond tree:
	bb(bondTree, calibratedTree, bkm, couponBond.nominal(), couponBond.lastCoupon(), couponBond.coupons(), timeDeltas);

	std::cout << "indexed " << decltype(bkm)::name() << " coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = bondTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(bondTree, first, last);

	std::cout << "Modeled price of pure discount bond:"
		<< bondTree.apex() << "\n";

	// Creating indexed option on bond lattice:
	lattice_structure::MeanRevertingLattice<double, date> optionTree(fixingDatesSet, params.ReversionSpeed, timeDeltas);

	// create option on coupon bond product:
	auto option = OptionOnCouponBondBuilder<double>()
		.withPeriods(periods)
		.withStrike(1.02)
		.build();

	// Create payoff fro the option	
	double K = option.strike();
	auto call_payoff = [&K](double bondPrice) {
		return std::max(bondPrice - K, 0.0);
	};

	auto american_adjuster = [&call_payoff](double &optionPrice,double bondPrice) {
		optionPrice =  std::max(optionPrice,call_payoff(bondPrice));
	};

	// typedef OptionOnBondBuilder:
	typedef lattice_bond_builders::OptionOnBondBuilder<lattice_types::LatticeType::Trinomial, std::vector<double>> option_builder;
	option_builder ob;
	// populate option tree:
	ob(optionTree, bondTree, calibratedTree, bkm, timeDeltas, call_payoff, american_adjuster);

	std::cout << "indexed Black-Karasinski option on coupon bond lattice:\n";
	// Print the part of generated lattice:
	first = optionTree.begin();
	last = std::next(first, 17);
	lattice_utility::print(optionTree, first, last);

	std::cout << "Modeled price of american call option on coupon bond:"
		<< optionTree.apex() << "\n";
}

void testTrinomialAmericanOptionOnCouponBondLattice() {
	std::cout << "=======================================================\n";
	std::cout << "=== American Option On Coupon Bond Lattice - TEST =====\n";
	std::cout << "=======================================================\n";

	HWAmericanOptionOnCouponBond();
	HWAmericanOptionOnCouponBondWithProduct();
	BKAmericanOptionOnCouponBond();
	BKAmericanOptionOnCouponBondWithProduct();

	std::cout << "=======================================================\n";
}


#endif ///_LATTICE_EXAMPLES_CALIBRATE_PRICE_IR