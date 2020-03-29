#pragma once
#if !defined(_LATTICE_GREEKS_T)
#define _LATTICE_GREEKS_T

#include"lattice_algorithms.h"
#include"lattice_utility.h"
#include<iostream>
#include"lattice_greeks.h"
#include<iostream>

void indexedLatticeBinomialDelta() {

	const int N = 3;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value, double dt)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	lattice_algorithms::forward_induction<>(il, fwdgen, 1.0, dt);
	std::cout << "After forward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	auto pt = il;

	lattice_types::LeafBackwardGenerator<double, double, double, double> backgen = [](double upValue, double downValue, double dt) {
		double p = 0.5;
		return (p*upValue + (1.0 - p)*downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	lattice_algorithms::backward_induction<>(pt, backgen, payoff, dt);
	std::cout << "After backward induction step:\n";
	for (auto const &v : pt.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << pt.apex() << "\n\n";
	std::cout << "Delta: " << lattice_greeks::delta(il, pt) << "\n";
}

void crrIndexedLatticeGreeks() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 95.0;
	option.RiskFreeRate = 0.025;
	option.DividentRate = 0.0;
	option.Volatility = 0.1;
	option.Underlying = 95.0;

	std::size_t periods{ 100 };
	double maturity = 0.5;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ option };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = crr;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// make a copy for later greeks: 
	auto pt(il);
	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = crr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(pt, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	std::cout << "Delta of call: " << lattice_greeks::delta(il,pt) << "\n";
	std::cout << "Gamma of call: " << lattice_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << lattice_greeks::theta(il, pt, dt) << "\n\n";
}


void mcrrIndexedLatticeGreeks() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 95.0;
	option.RiskFreeRate = 0.025;
	option.DividentRate = 0.0;
	option.Volatility = 0.1;
	option.Underlying = 95.0;

	std::size_t periods{ 100 };
	double maturity = 0.5;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::ModifiedCoxRubinsteinRossModel<> mcrr{ option,periods };

	// Print the model name:
	std::cout << decltype(mcrr)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = mcrr;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// make a copy for later greeks: 
	auto pt(il);
	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = mcrr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(pt, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	std::cout << "Delta of call: " << lattice_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << lattice_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << lattice_greeks::theta(il, pt, dt) << "\n\n";
}


void jrIndexedLatticeGreeks() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 95.0;
	option.RiskFreeRate = 0.025;
	option.DividentRate = 0.0;
	option.Volatility = 0.1;
	option.Underlying = 95.0;

	std::size_t periods{ 100 };
	double maturity = 0.5;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::JarrowRuddModel<> jr{ option };

	// Print the model name:
	std::cout << decltype(jr)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = jr;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// make a copy for later greeks: 
	auto pt(il);
	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);


	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = jr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(pt, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	std::cout << "Delta of call: " << lattice_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << lattice_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << lattice_greeks::theta(il, pt, dt) << "\n\n";
}

void trimIndexedLatticeGreeks() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 95.0;
	option.RiskFreeRate = 0.025;
	option.DividentRate = 0.0;
	option.Volatility = 0.1;
	option.Underlying = 95.0;

	std::size_t periods{ 100 };
	double maturity = 0.5;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::TrigeorgisModel<> trim{ option };

	// Print the model name:
	std::cout << decltype(trim)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = trim;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// make a copy for later greeks: 
	auto pt(il);
	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = trim;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(pt, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	std::cout << "Delta of call: " << lattice_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << lattice_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << lattice_greeks::theta(il, pt, dt) << "\n\n";
}

void tmIndexedLatticeGreeks() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 95.0;
	option.RiskFreeRate = 0.025;
	option.DividentRate = 0.0;
	option.Volatility = 0.1;
	option.Underlying = 95.0;

	std::size_t periods{ 100 };
	double maturity = 0.5;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::TianModel<> tm{ option };

	// Print the model name:
	std::cout << decltype(tm)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = tm;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// make a copy for later greeks: 
	auto pt(il);
	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = tm;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(pt, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	std::cout << "Delta of call: " << lattice_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << lattice_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << lattice_greeks::theta(il, pt, dt) << "\n\n";
}


void lrIndexedLatticeGreeks() {

	lattice_miscellaneous::OptionData<double> option;
	using lattice_model_components::leisen_reimer_inversion::PeizerPrattSecondInversion;


	option.Strike = 95.0;
	option.RiskFreeRate = 0.025;
	option.DividentRate = 0.0;
	option.Volatility = 0.1;
	option.Underlying = 95.0;

	std::size_t periods{ 100 };
	double maturity = 0.5;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Prepare inversion formula functor:
	std::size_t numberTimePoints = periods + 1;
	PeizerPrattSecondInversion<> ppi{ numberTimePoints };

	// Create CRR model:
	lattice_model::LeisenReimerModel<> lr{ option,periods,ppi };

	// Print the model name:
	std::cout << decltype(lr)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = lr;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// make a copy for later greeks: 
	auto pt(il);
	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = lr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(pt, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	std::cout << "Delta of call: " << lattice_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << lattice_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << lattice_greeks::theta(il, pt, dt) << "\n\n";
}

void bmIndexedLatticeGreeks() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 95.0;
	option.RiskFreeRate = 0.025;
	option.DividentRate = 0.0;
	option.Volatility = 0.1;
	option.Underlying = 95.0;

	std::size_t periods{ 100 };
	double maturity = 0.5;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(periods);

	// Create CRR model:
	lattice_model::BoyleModel<> bm{ option };

	// Print the model name:
	std::cout << decltype(bm)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double, double> fwdGen = bm;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// make a copy for later greeks: 
	auto pt(il);
	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double, double> backGen = bm;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(pt, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	std::cout << "Delta of call: " << lattice_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << lattice_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << lattice_greeks::theta(il, pt, dt) << "\n\n";
}


#endif ///_LATTICE_GREEKS_T