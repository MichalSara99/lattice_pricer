#pragma once
#ifndef _LATTICE_BARRIER_OPTION_H
#define _LATTICE_BARRIER_OPTION_H


#include"lattice_structure.h"
#include"lattice_model.h"
#include"lattice_miscellaneous.h"
#include"lattice_utility.h"
#include"lattice_algorithms.h"
#include"lattice_model_components.h"
#include"lattice_builder.h"
#include"lattice_product_builder.h"


void crrIndexedBarrierLattice() {

	using lattice_product_builder::BarrierOptionBuilder;
	using lattice_types::BarrierType;

	// build option product:
	auto option = BarrierOptionBuilder<double>()
		.withDividend(0.0).withRate(0.06)
		.withSpot(100.0).withStrike(100.0)
		.withPeriods(4).withVolatility(0.3)
		.withBarrier(110.0).withRebate(1.0)
		.withBarrierType(BarrierType::UpAndOut)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	double maturity = 1.0;
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

	brd_induction(il, crr, dt, put_payoff, option.barrierType(), option.barrier(), option.rebate());

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}


void mcrrIndexedBarrierLattice() {

	using lattice_product_builder::BarrierOptionBuilder;
	using lattice_types::BarrierType;

	// build option product:
	auto option = BarrierOptionBuilder<double>()
		.withDividend(0.0).withRate(0.06)
		.withSpot(100.0).withStrike(100.0)
		.withPeriods(4).withVolatility(0.3)
		.withBarrier(110.0).withRebate(1.0)
		.withBarrierType(BarrierType::UpAndOut)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	double maturity = 1.0;
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

	brd_induction(il, mcrr, dt, put_payoff, option.barrierType(), option.barrier(), option.rebate());

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}
//
//

void jrIndexedBarrierLattice() {
	using lattice_product_builder::BarrierOptionBuilder;
	using lattice_types::BarrierType;

	// build option product:
	auto option = BarrierOptionBuilder<double>()
		.withDividend(0.0).withRate(0.06)
		.withSpot(100.0).withStrike(100.0)
		.withPeriods(4).withVolatility(0.3)
		.withBarrier(110.0).withRebate(1.0)
		.withBarrierType(BarrierType::UpAndOut)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	double maturity = 1.0;
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

	brd_induction(il, jr, dt, put_payoff, option.barrierType(), option.barrier(), option.rebate());

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}

//
void trimIndexedBarrierLattice() {

	using lattice_product_builder::BarrierOptionBuilder;
	using lattice_types::BarrierType;

	// build option product:
	auto option = BarrierOptionBuilder<double>()
		.withDividend(0.0).withRate(0.06)
		.withSpot(100.0).withStrike(100.0)
		.withPeriods(4).withVolatility(0.3)
		.withBarrier(110.0).withRebate(1.0)
		.withBarrierType(BarrierType::UpAndOut)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	double maturity = 1.0;
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

	brd_induction(il, trim, dt, put_payoff, option.barrierType(), option.barrier(), option.rebate());

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}
//
//
//
void tmIndexedBarrierLattice() {

	using lattice_product_builder::BarrierOptionBuilder;
	using lattice_types::BarrierType;

	// build option product:
	auto option = BarrierOptionBuilder<double>()
		.withDividend(0.0).withRate(0.06)
		.withSpot(100.0).withStrike(100.0)
		.withPeriods(4).withVolatility(0.3)
		.withBarrier(110.0).withRebate(1.0)
		.withBarrierType(BarrierType::UpAndOut)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	double maturity = 1.0;
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

	brd_induction(il, tm, dt, put_payoff, option.barrierType(), option.barrier(), option.rebate());

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}


void lrIndexedBarrierLattice() {

	using lattice_model_components::leisen_reimer_inversion::PeizerPrattSecondInversion;
	using lattice_product_builder::BarrierOptionBuilder;
	using lattice_types::BarrierType;

	// build option product:
	auto option = BarrierOptionBuilder<double>()
		.withDividend(0.0).withRate(0.06)
		.withSpot(100.0).withStrike(100.0)
		.withPeriods(4).withVolatility(0.3)
		.withBarrier(110.0).withRebate(1.0)
		.withBarrierType(BarrierType::UpAndOut)
		.build();

	// extract model params from option:
	auto params = option.modelParams();

	double maturity = 1.0;
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

	brd_induction(il, lr, dt, put_payoff, option.barrierType(), option.barrier(), option.rebate());

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}

void bmIndexedBarrierLattice() {

	using lattice_product_builder::BarrierOptionBuilder;
	using lattice_types::BarrierType;

	// build option product:
	auto option = BarrierOptionBuilder<double>()
		.withDividend(0.0).withRate(0.06)
		.withSpot(100.0).withStrike(100.0)
		.withPeriods(4).withVolatility(0.3)
		.withBarrier(110.0).withRebate(1.0)
		.withBarrierType(BarrierType::UpAndOut)
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

	brd_induction(il, bm, dt, put_payoff, option.barrierType(), option.barrier(), option.rebate());

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of put: " << il.apex() << "\n\n";
}

void testIndexedEuropeanBarrierLattices() {
	std::cout << "======================================================================\n";
	std::cout << "========== Indexed European Barrier Option Lattices - TEST ===========\n";
	std::cout << "======================================================================\n";

	crrIndexedBarrierLattice();
	mcrrIndexedBarrierLattice();
	jrIndexedBarrierLattice();
	trimIndexedBarrierLattice();
	tmIndexedBarrierLattice();
	lrIndexedBarrierLattice();
	bmIndexedBarrierLattice();

	std::cout << "======================================================================\n";
}















#endif ///_LATTICE_BARRIER_OPTION_H