#pragma once
#if !defined(_LATTICE_MODEL_T)
#define _LATTICE_MODEL_T

#include"lattice_structure.h"
#include"lattice_model.h"
#include"lattice_miscellaneous.h"
#include"lattice_utility.h"
#include"lattice_algorithms.h"
#include"lattice_model_components.h"
#include"lattice_multidimensional.h"
#include"lattice_product_builder.h"

#include<iostream>
#include<numeric>

void crr2FactorIndexedLatticeCall() {

	using lattice_product_builder::SpreadOptionBuilder;

	// build spread option product:
	auto spreadOption = SpreadOptionBuilder<double>()
		.withStrike(1.0).withRate(0.06)
		.withSpot1(100.0).withSpot2(100.0)
		.withDividend1(0.03).withDividend2(0.04)
		.withVolatility1(0.2).withVolatility2(0.3)
		.withCorrelation(0.5).withPeriods(3)
		.build();

	// return model params 
	auto params = spreadOption.modelParams();

	// set maturity with dt:
	double maturity = 1.0;
	double dt = (maturity / double(spreadOption.periods()));


	// Creating indexed lattice:
	lattice_multidimensional::MultidimIndexedLattice<2,lattice_types::LatticeType::Binomial, double> il(spreadOption.periods());

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel2Factor<> crr{ params };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, crr, dt, std::make_pair(params.Spot1, params.Spot2));

	// Print the part of generated lattice:
	std::cout << "First factor: \n";
	auto factor1 = il.getFactor(0);
	auto first = factor1.begin();
	auto last = factor1.end();
	lattice_utility::print(factor1, first, last);
	std::cout << "Second factor: \n";
	auto factor2 = il.getFactor(1);
	first = factor2.begin();
	last = factor2.end();
	lattice_utility::print(factor2, first, last);

	// Backward induction:

	// Prepare payoff for spread call option:
	double K = params.Strike;
	auto call_payoff = [&K](double stock1, double stock2) {return std::max(stock1 - stock2 - K, 0.0); };

	// Creating indexed two-variable binomial lattice:
	lattice_multidimensional::IndexedLattice<lattice_types::LatticeType::TwoVariableBinomial, double> optionTree(spreadOption.periods());


	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
	std::size_t, double> backward_2binomial_induction;
	backward_2binomial_induction brd_induction;

	brd_induction(optionTree, il, crr, dt, call_payoff);

	// Print the part of generated lattice:
	std::cout << "Option price lattice:\n";
	first = optionTree.begin();
	last = optionTree.end();
	lattice_utility::print(optionTree, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << optionTree.apex() << "\n\n";
}

void crr2FactorIndexedLatticePut() {

	using lattice_product_builder::SpreadOptionBuilder;

	// build spread option product:
	auto spreadOption = SpreadOptionBuilder<double>()
		.withStrike(1.0).withRate(0.06)
		.withSpot1(100.0).withSpot2(100.0)
		.withDividend1(0.03).withDividend2(0.04)
		.withVolatility1(0.2).withVolatility2(0.3)
		.withCorrelation(0.5).withPeriods(3)
		.build();

	// extract model params from spread option:
	auto params = spreadOption.modelParams();

	// set maturity with dt:
	double maturity = 1.0;
	double dt = (maturity / double(spreadOption.periods()));


	// Creating indexed lattice:
	lattice_multidimensional::MultidimIndexedLattice<2, lattice_types::LatticeType::Binomial, double> il(spreadOption.periods());

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel2Factor<> crr{ params };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, crr, dt, std::make_pair(params.Spot1, params.Spot2));

	// Print the part of generated lattice:
	std::cout << "First factor: \n";
	auto factor1 = il.getFactor(0);
	auto first = factor1.begin();
	auto last = factor1.end();
	lattice_utility::print(factor1, first, last);
	std::cout << "Second factor: \n";
	auto factor2 = il.getFactor(1);
	first = factor2.begin();
	last = factor2.end();
	lattice_utility::print(factor2, first, last);

	// Backward induction:

	// Prepare payoff for spread put option:
	double K = params.Strike;
	auto put_payoff = [&K](double stock1, double stock2) {return std::max(K - stock1 - stock2, 0.0); };

	// Creating indexed two-variable binomial lattice:
	lattice_multidimensional::IndexedLattice<lattice_types::LatticeType::TwoVariableBinomial, double> optionTree(spreadOption.periods());


	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		std::size_t, double> backward_2binomial_induction;
	backward_2binomial_induction brd_induction;

	brd_induction(optionTree, il, crr, dt, put_payoff);

	// Print the part of generated lattice:
	std::cout << "Option price lattice:\n";
	first = optionTree.begin();
	last = optionTree.end();
	lattice_utility::print(optionTree, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << optionTree.apex() << "\n\n";
}

void crrIndexedLattice() {

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

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(option.periods());

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ params };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;
	
	forward_binomial_induction fwd_induction;
	fwd_induction(il, crr, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	
	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, crr, dt, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}


void mcrrIndexedLattice() {

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

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(option.periods());

	// Create CRR model:
	lattice_model::ModifiedCoxRubinsteinRossModel<> mcrr{ params,option.periods() };

	// Print the model name:
	std::cout << decltype(mcrr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, mcrr, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, mcrr, dt, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}
//
//

void jrIndexedLattice() {
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


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(option.periods());

	// Create CRR model:
	lattice_model::JarrowRuddModel<> jr{ params };

	// Print the model name:
	std::cout << decltype(jr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, jr, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, jr, dt, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}

//
void trimIndexedLattice() {

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


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(option.periods());

	// Create CRR model:
	lattice_model::TrigeorgisModel<> trim{ params };

	// Print the model name:
	std::cout << decltype(trim)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, trim, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, trim, dt, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}
//
//
//
void tmIndexedLattice() {

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


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(option.periods());

	// Create CRR model:
	lattice_model::TianModel<> tm{ params };

	// Print the model name:
	std::cout << decltype(tm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, tm, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, tm, dt, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}


void lrIndexedLattice() {

	using lattice_product_builder::OptionBuilder;
	using lattice_model_components::leisen_reimer_inversion::PeizerPrattSecondInversion;

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


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(option.periods());

	// Prepare inversion formula functor:
	std::size_t numberTimePoints = option.periods() + 1;
	PeizerPrattSecondInversion<> ppi{ numberTimePoints };

	// Create CRR model:
	lattice_model::LeisenReimerModel<> lr{ params,option.periods(),ppi };

	// Print the model name:
	std::cout << decltype(lr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, lr, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, lr, dt, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}

void bmIndexedLattice() {

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


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(option.periods());

	// Create CRR model:
	lattice_model::BoyleModel<> bm{ params };

	// Print the model name:
	std::cout << decltype(bm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, bm, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	// Prepare payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, bm, dt, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}

void testIndexedEuropeanLattices() {
	std::cout << "=======================================================\n";
	std::cout << "========== Indexed European Lattices - TEST ===========\n";
	std::cout << "=======================================================\n";

	crr2FactorIndexedLatticeCall();
	crr2FactorIndexedLatticePut();
	crrIndexedLattice();
	mcrrIndexedLattice();
	jrIndexedLattice();
	trimIndexedLattice();
	tmIndexedLattice();
	lrIndexedLattice();
	bmIndexedLattice();

	std::cout << "=======================================================\n";
}


void crr2FactorIndexedLatticeAmericanCall() {

	using lattice_product_builder::SpreadOptionBuilder;

	// build spread option product:
	auto spreadOption = SpreadOptionBuilder<double>()
		.withStrike(1.0).withRate(0.06)
		.withSpot1(100.0).withSpot2(100.0)
		.withDividend1(0.03).withDividend2(0.04)
		.withVolatility1(0.2).withVolatility2(0.3)
		.withCorrelation(0.5).withPeriods(3)
		.build();

	// extract model params from spread option:
	auto params = spreadOption.modelParams();

	// set maturity with dt:
	double maturity = 1.0;
	double dt = (maturity / double(spreadOption.periods()));

	// Creating indexed lattice:
	lattice_multidimensional::MultidimIndexedLattice<2, lattice_types::LatticeType::Binomial, double> il(spreadOption.periods());

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel2Factor<> crr{ params };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, crr, dt, std::make_pair(params.Spot1, params.Spot2));

	// Print the part of generated lattice:
	std::cout << "First factor: \n";
	auto factor1 = il.getFactor(0);
	auto first = factor1.begin();
	auto last = factor1.end();
	lattice_utility::print(factor1, first, last);
	std::cout << "Second factor: \n";
	auto factor2 = il.getFactor(1);
	first = factor2.begin();
	last = factor2.end();
	lattice_utility::print(factor2, first, last);

	// Backward induction:

	// Prepare spread call option payoff:
	double K = params.Strike;
	auto call_payoff = [&K](double stock1, double stock2) {return std::max(stock1 - stock2 - K, 0.0); };
	// Prepare american spread call option adujster:
	auto american_adjuster = [&call_payoff](double& value, double stock1, double stock2) {
		value = std::max(value, call_payoff(stock1, stock2));
	};

	// Creating indexed two-variable binomial lattice:
	lattice_multidimensional::IndexedLattice<lattice_types::LatticeType::TwoVariableBinomial, double> optionTree(spreadOption.periods());


	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		std::size_t, double> backward_2binomial_induction;
	backward_2binomial_induction brd_induction;

	brd_induction(optionTree, il, crr, dt, call_payoff, american_adjuster);

	// Print the part of generated lattice:
	std::cout << "Option price lattice:\n";
	first = optionTree.begin();
	last = optionTree.end();
	lattice_utility::print(optionTree, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << optionTree.apex() << "\n\n";
}

void crr2FactorIndexedLatticeAmericanPut() {

	using lattice_product_builder::SpreadOptionBuilder;

	// build spread option product:
	auto spreadOption = SpreadOptionBuilder<double>()
		.withStrike(1.0).withRate(0.06)
		.withSpot1(100.0).withSpot2(100.0)
		.withDividend1(0.03).withDividend2(0.04)
		.withVolatility1(0.2).withVolatility2(0.3)
		.withCorrelation(0.5).withPeriods(3)
		.build();

	// extract model params from spread option:
	auto params = spreadOption.modelParams();

	// set maturity with dt:
	double maturity = 1.0;
	double dt = (maturity / double(spreadOption.periods()));


	// Creating indexed lattice:
	lattice_multidimensional::MultidimIndexedLattice<2, lattice_types::LatticeType::Binomial, double> il(spreadOption.periods());

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel2Factor<> crr{ params };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, crr, dt, std::make_pair(params.Spot1, params.Spot2));

	// Print the part of generated lattice:
	std::cout << "First factor: \n";
	auto factor1 = il.getFactor(0);
	auto first = factor1.begin();
	auto last = factor1.end();
	lattice_utility::print(factor1, first, last);
	std::cout << "Second factor: \n";
	auto factor2 = il.getFactor(1);
	first = factor2.begin();
	last = factor2.end();
	lattice_utility::print(factor2, first, last);

	// Backward induction:

	// Prepare spread put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock1, double stock2) {return std::max(K - stock1 - stock2, 0.0); };
	// Prepare american spread put option adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock1, double stock2) {
		value = std::max(value, put_payoff(stock1, stock2));
	};

	// Creating indexed two-variable binomial lattice:
	lattice_multidimensional::IndexedLattice<lattice_types::LatticeType::TwoVariableBinomial, double> optionTree(spreadOption.periods());


	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		std::size_t, double> backward_2binomial_induction;
	backward_2binomial_induction brd_induction;

	brd_induction(optionTree, il, crr, dt, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	std::cout << "Option price lattice:\n";
	first = optionTree.begin();
	last = optionTree.end();
	lattice_utility::print(optionTree, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << optionTree.apex() << "\n\n";
}

void crrIndexedLatticeAmerican() {

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


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(option.periods());

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ params };

	std::cout << decltype(crr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, crr, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american put option adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, crr, dt, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}



void mcrrIndexedLatticeAmerican() {

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


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(option.periods());

	// Create CRR model:
	lattice_model::ModifiedCoxRubinsteinRossModel<> mcrr{ params,option.periods() };

	std::cout << decltype(mcrr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, mcrr, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american put option adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, mcrr, dt, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}


void jrIndexedLatticeAmerican() {

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


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(option.periods());

	// Create CRR model:
	lattice_model::JarrowRuddModel<> jr{ params };

	std::cout << decltype(jr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, jr, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american put option adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, jr, dt, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}


void trimIndexedLatticeAmerican() {

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


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(option.periods());

	// Create CRR model:
	lattice_model::TrigeorgisModel<> trim{ params };

	std::cout << decltype(trim)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, trim, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american put option adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, trim, dt, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}



void tmIndexedLatticeAmerican() {

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


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(option.periods());

	// Create CRR model:
	lattice_model::TianModel<> tm{ params };

	std::cout << decltype(tm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, tm, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american put option adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, tm, dt, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}


void lrIndexedLatticeAmerican() {

	using lattice_model_components::leisen_reimer_inversion::PeizerPrattSecondInversion;
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


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(option.periods());

	// Prepare inversion functor:
	std::size_t numerTimePoints{ option.periods() + 1 };
	PeizerPrattSecondInversion<> ppi{ numerTimePoints };

	// Create CRR model:
	lattice_model::LeisenReimerModel<> lr{ params,option.periods(),ppi };

	std::cout << decltype(lr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, lr, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, lr, dt, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}

void bmIndexedLatticeAmerican() {

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


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(option.periods());

	// Create CRR model:
	lattice_model::BoyleModel<> bm{ params };

	// Print the model name:
	std::cout << decltype(bm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, bm, dt, params.Spot);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american put option adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, bm, dt, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}


void testIndexedAmericanLattices() {
	std::cout << "=======================================================\n";
	std::cout << "========== Indexed American Lattices - TEST ===========\n";
	std::cout << "=======================================================\n";

	crr2FactorIndexedLatticeAmericanCall();
	crr2FactorIndexedLatticeAmericanPut();
	crrIndexedLatticeAmerican();
	mcrrIndexedLatticeAmerican();
	jrIndexedLatticeAmerican();
	trimIndexedLatticeAmerican();
	tmIndexedLatticeAmerican();
	lrIndexedLatticeAmerican();
	bmIndexedLatticeAmerican();

	std::cout << "=======================================================\n";
}

////
//

void crr2FactorLatticeCall() {

	using lattice_product_builder::SpreadOptionBuilder;

	// build spread option product:
	auto spreadOption = SpreadOptionBuilder<double>()
		.withStrike(1.0).withRate(0.06)
		.withSpot1(100.0).withSpot2(100.0)
		.withDividend1(0.03).withDividend2(0.04)
		.withVolatility1(0.2).withVolatility2(0.3)
		.withCorrelation(0.5).withPeriods(3)
		.build();

	// return model params 
	auto params = spreadOption.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= spreadOption.periods(); ++t) {
		fixingDates.emplace(today + date_duration(90*t));
	}

	// Creating indexed lattice:
	lattice_multidimensional::MultidimLattice<2, lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };


	double daysInYear{ 365.0 };
	auto fd = la.getFactor(1).fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel2Factor<> crr{ params };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, crr, timeDeltas, std::make_pair(params.Spot1, params.Spot2));

	// Print the part of generated lattice:
	std::cout << "First factor: \n";
	auto factor1 = la.getFactor(0);
	auto first = factor1.begin();
	auto last = factor1.end();
	lattice_utility::print(factor1, first, last);
	std::cout << "Second factor: \n";
	auto factor2 = la.getFactor(1);
	first = factor2.begin();
	last = factor2.end();
	lattice_utility::print(factor2, first, last);

	// Backward induction:

	// Prepare payoff for spread call option:
	double K = params.Strike;
	auto call_payoff = [&K](double stock1, double stock2) {return std::max(stock1 - stock2 - K, 0.0); };

	// Creating indexed two-variable binomial lattice:
	lattice_multidimensional::Lattice<lattice_types::LatticeType::TwoVariableBinomial, double, date> optionTree{ fixingDates };


	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		date, std::vector<double>> backward_2binomial_induction;
	backward_2binomial_induction brd_induction;

	brd_induction(optionTree, la, crr, timeDeltas, call_payoff);

	// Print the part of generated lattice:
	std::cout << "Option price lattice:\n";
	first = optionTree.begin();
	last = optionTree.end();
	lattice_utility::print(optionTree, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << optionTree.apex() << "\n\n";
}

void crr2FactorLatticePut() {

	using lattice_product_builder::SpreadOptionBuilder;

	// build spread option product:
	auto spreadOption = SpreadOptionBuilder<double>()
		.withStrike(1.0).withRate(0.06)
		.withSpot1(100.0).withSpot2(100.0)
		.withDividend1(0.03).withDividend2(0.04)
		.withVolatility1(0.2).withVolatility2(0.3)
		.withCorrelation(0.5).withPeriods(3)
		.build();

	// extract model params from spread option:
	auto params = spreadOption.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= spreadOption.periods(); ++t) {
		fixingDates.emplace(today + date_duration(90*t));
	}


	// Creating indexed lattice:
	lattice_multidimensional::MultidimLattice<2, lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };

	double daysInYear{ 365.0 };
	auto fd = la.getFactor(1).fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel2Factor<> crr{ params };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, crr, timeDeltas, std::make_pair(params.Spot1, params.Spot2));

	// Print the part of generated lattice:
	std::cout << "First factor: \n";
	auto factor1 = la.getFactor(0);
	auto first = factor1.begin();
	auto last = factor1.end();
	lattice_utility::print(factor1, first, last);
	std::cout << "Second factor: \n";
	auto factor2 = la.getFactor(1);
	first = factor2.begin();
	last = factor2.end();
	lattice_utility::print(factor2, first, last);

	// Backward induction:

	// Prepare payoff for spread put option:
	double K = params.Strike;
	auto put_payoff = [&K](double stock1, double stock2) {return std::max(K - stock1 - stock2, 0.0); };

	// Creating indexed two-variable binomial lattice:
	lattice_multidimensional::Lattice<lattice_types::LatticeType::TwoVariableBinomial, double, date> optionTree{ fixingDates };


	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		date, std::vector<double>> backward_2binomial_induction;
	backward_2binomial_induction brd_induction;

	brd_induction(optionTree, la, crr, timeDeltas, put_payoff);

	// Print the part of generated lattice:
	std::cout << "Option price lattice:\n";
	first = optionTree.begin();
	last = optionTree.end();
	lattice_utility::print(optionTree, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << optionTree.apex() << "\n\n";
}


void crrLattice() {

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::CoxRubinsteinRossModel<> crr{ params };

	// Name of the model:
	std::cout << decltype(crr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date,std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, crr, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la,crr, timeDeltas, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}

void mcrrLattice() {

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la{fixingDates};

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::ModifiedCoxRubinsteinRossModel<> mcrr{ params,option.periods() };

	// Name of the model:
	std::cout << decltype(mcrr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, mcrr, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, mcrr, timeDeltas, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}
//


void jrLattice() {

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::JarrowRuddModel<> jr{ params };

	// Name of the model:
	std::cout << decltype(jr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, jr, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, jr, timeDeltas, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}

void trimLattice() {

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::TrigeorgisModel<> trim{ params };

	// Name of the model:
	std::cout << decltype(trim)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, trim, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, trim, timeDeltas, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}
//

void tmLattice() {

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::TianModel<> tm{ params };

	// Name of the model:
	std::cout << decltype(tm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, tm, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, tm, timeDeltas, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}
//


void lrLattice() {

	using ::lattice_model_components::leisen_reimer_inversion::PeizerPrattSecondInversion;
	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	// Prepare inversion formula:
	std::size_t numberTimePoints{ option.periods() + 1 };
	PeizerPrattSecondInversion<> ppi{ numberTimePoints };


	lattice_model::LeisenReimerModel<> lr{ params,option.periods(),ppi };

	// Name of the model:
	std::cout << decltype(lr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, lr, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, lr, timeDeltas, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}
//


void bmLattice() {

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Trinomial, double, date> la{ fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::BoyleModel<> bm{ params };

	// Name of the model:
	std::cout << decltype(bm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, bm, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, bm, timeDeltas, put_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}





void testEuropeanLattices() {
	std::cout << "=======================================================\n";
	std::cout << "=============== European Lattices - TEST ==============\n";
	std::cout << "=======================================================\n";

	crr2FactorLatticeCall();
	crr2FactorLatticePut();
	crrLattice();
	mcrrLattice();
	jrLattice();
	trimLattice();
	tmLattice();
	lrLattice();
	bmLattice();

	std::cout << "=======================================================\n";
}

void crr2FactorLatticeAmericanCall() {

	using lattice_product_builder::SpreadOptionBuilder;

	// build spread option product:
	auto spreadOption = SpreadOptionBuilder<double>()
		.withStrike(1.0).withRate(0.06)
		.withSpot1(100.0).withSpot2(100.0)
		.withDividend1(0.03).withDividend2(0.04)
		.withVolatility1(0.2).withVolatility2(0.3)
		.withCorrelation(0.5).withPeriods(3)
		.build();

	// return model params 
	auto params = spreadOption.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= spreadOption.periods(); ++t) {
		fixingDates.emplace(today + date_duration(90 * t));
	}

	// Creating indexed lattice:
	lattice_multidimensional::MultidimLattice<2, lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };


	double daysInYear{ 365.0 };
	auto fd = la.getFactor(1).fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel2Factor<> crr{ params };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, crr, timeDeltas, std::make_pair(params.Spot1, params.Spot2));

	// Print the part of generated lattice:
	std::cout << "First factor: \n";
	auto factor1 = la.getFactor(0);
	auto first = factor1.begin();
	auto last = factor1.end();
	lattice_utility::print(factor1, first, last);
	std::cout << "Second factor: \n";
	auto factor2 = la.getFactor(1);
	first = factor2.begin();
	last = factor2.end();
	lattice_utility::print(factor2, first, last);

	// Backward induction:

	// Prepare payoff for spread call option:
	double K = params.Strike;
	auto call_payoff = [&K](double stock1, double stock2) {return std::max(stock1 - stock2 - K, 0.0); };

	// Prepare american call spread option adujster:
	auto american_adjuster = [&call_payoff](double& value, double stock1, double stock2) {
		value = std::max(value, call_payoff(stock1, stock2));
	};

	// Creating indexed two-variable binomial lattice:
	lattice_multidimensional::Lattice<lattice_types::LatticeType::TwoVariableBinomial, double, date> optionTree{ fixingDates };


	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		date, std::vector<double>> backward_2binomial_induction;
	backward_2binomial_induction brd_induction;

	brd_induction(optionTree, la, crr, timeDeltas, call_payoff, american_adjuster);

	// Print the part of generated lattice:
	std::cout << "Option price lattice:\n";
	first = optionTree.begin();
	last = optionTree.end();
	lattice_utility::print(optionTree, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << optionTree.apex() << "\n\n";
}

void crr2FactorLatticeAmericanPut() {

	using lattice_product_builder::SpreadOptionBuilder;

	// build spread option product:
	auto spreadOption = SpreadOptionBuilder<double>()
		.withStrike(1.0).withRate(0.06)
		.withSpot1(100.0).withSpot2(100.0)
		.withDividend1(0.03).withDividend2(0.04)
		.withVolatility1(0.2).withVolatility2(0.3)
		.withCorrelation(0.5).withPeriods(3)
		.build();

	// extract model params from spread option:
	auto params = spreadOption.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= spreadOption.periods(); ++t) {
		fixingDates.emplace(today + date_duration(90 * t));
	}


	// Creating indexed lattice:
	lattice_multidimensional::MultidimLattice<2, lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };

	double daysInYear{ 365.0 };
	auto fd = la.getFactor(1).fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel2Factor<> crr{ params };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, crr, timeDeltas, std::make_pair(params.Spot1, params.Spot2));

	// Print the part of generated lattice:
	std::cout << "First factor: \n";
	auto factor1 = la.getFactor(0);
	auto first = factor1.begin();
	auto last = factor1.end();
	lattice_utility::print(factor1, first, last);
	std::cout << "Second factor: \n";
	auto factor2 = la.getFactor(1);
	first = factor2.begin();
	last = factor2.end();
	lattice_utility::print(factor2, first, last);

	// Backward induction:

	// Prepare payoff for spread put option:
	double K = params.Strike;
	auto put_payoff = [&K](double stock1, double stock2) {return std::max(K - stock1 - stock2, 0.0); };

	// Prepare american put spread option adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock1, double stock2) {
		value = std::max(value, put_payoff(stock1, stock2));
	};

	// Creating indexed two-variable binomial lattice:
	lattice_multidimensional::Lattice<lattice_types::LatticeType::TwoVariableBinomial, double, date> optionTree{ fixingDates };


	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::TwoVariableBinomial,
		date, std::vector<double>> backward_2binomial_induction;
	backward_2binomial_induction brd_induction;

	brd_induction(optionTree, la, crr, timeDeltas, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	std::cout << "Option price lattice:\n";
	first = optionTree.begin();
	last = optionTree.end();
	lattice_utility::print(optionTree, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << optionTree.apex() << "\n\n";
}

void crrLatticeAmerican() {

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}


	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::CoxRubinsteinRossModel<> crr{ params };

	// Name of the model:
	std::cout << decltype(crr)::name() << "\n";


	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, crr, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american put option adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, crr, timeDeltas, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}


void mcrrLatticeAmerican() {

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}


	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::ModifiedCoxRubinsteinRossModel<> mcrr{ params,option.periods() };

	// Name of the model:
	std::cout << decltype(mcrr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, mcrr, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american put option adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, mcrr, timeDeltas, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}


void jrLatticeAmerican() {

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}


	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::JarrowRuddModel<> jr{ params };

	// Name of the model:
	std::cout << decltype(jr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, jr, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american put option adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, jr, timeDeltas, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}


void trimLatticeAmerican() {

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}


	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::TrigeorgisModel<> trim{ params };

	// Name of the model:
	std::cout << decltype(trim)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, trim, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, trim, timeDeltas, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}


void tmLatticeAmerican() {

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}


	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::TianModel<> tm{ params };

	// Name of the model:
	std::cout << decltype(tm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, tm, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american put option adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, tm, timeDeltas, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}

void lrLatticeAmerican() {
	using ::lattice_model_components::leisen_reimer_inversion::PeizerPrattSecondInversion;
	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}


	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la{ fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	// Prepare inversion formula:
	std::size_t numberTimePoints{ option.periods() + 1 };
	PeizerPrattSecondInversion<> ppi{ numberTimePoints };

	lattice_model::LeisenReimerModel<> lr{ params,option.periods(),ppi };

	// Name of the model:
	std::cout << decltype(lr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, lr, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare put option payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american put option adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, lr, timeDeltas, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}


void bmLatticeAmerican() {

	using lattice_product_builder::OptionBuilder;

	// build option product:
	auto option = OptionBuilder<double>()
		.withDividend(0.0).withRate(0.25)
		.withSpot(60.0).withStrike(65.0)
		.withPeriods(100).withVolatility(0.3)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t <= option.periods(); ++t) {
		fixingDates.emplace(today + date_duration(t));
	}


	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Trinomial, double, date> la{ fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::BoyleModel<> bm{ params };

	// Name of the model:
	std::cout << decltype(bm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, bm, timeDeltas, params.Spot);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = params.Strike;
	auto put_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	auto american_adjuster = [&put_payoff](double& value, double stock) {
		value = std::max(value, put_payoff(stock));
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, bm, timeDeltas, put_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << la.apex() << "\n\n";

}

void testAmericanLattices() {
	std::cout << "=======================================================\n";
	std::cout << "=============== American Lattices - TEST ==============\n";
	std::cout << "=======================================================\n";

	crr2FactorLatticeAmericanCall();
	crr2FactorLatticeAmericanPut();
	crrLatticeAmerican();
	mcrrLatticeAmerican();
	jrLatticeAmerican();
	trimLatticeAmerican();
	tmLatticeAmerican();
	lrLatticeAmerican();
	bmLatticeAmerican();

	std::cout << "=======================================================\n";
}



#endif ///_LATTICE_MODEL_T