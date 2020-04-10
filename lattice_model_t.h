#pragma once
#if !defined(_LATTICE_MODEL_T)
#define _LATTICE_MODEL_T

#include"lattice_structure.h"
#include"lattice_model.h"
#include"lattice_miscellaneous.h"
#include"lattice_utility.h"
#include"lattice_algorithms.h"
#include"lattice_model_components.h"

#include<iostream>
#include<numeric>

void crrIndexedLattice() {
	
	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ option};

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double,double,double> fwdGen = crr;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double,double, double, double> backGen = crr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K-stock, 0.0); };
	lattice_algorithms::backward_induction(il, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

void crrIndexedLatticeNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ option };

	// Print the model name:
	std::cout << decltype(crr)::name() << "\n";


	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;
	
	forward_binomial_induction fwd_induction;
	fwd_induction(il, crr, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = crr;
	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, crr, dt, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}
//
void mcrrIndexedLattice() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
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

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = mcrr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(il, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}



void mcrrIndexedLatticeNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::ModifiedCoxRubinsteinRossModel<> mcrr{ option,periods };

	// Print the model name:
	std::cout << decltype(mcrr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, mcrr, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, mcrr, dt, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

//
void jrIndexedLattice() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
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

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = jr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(il, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}


void jrIndexedLatticeNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::JarrowRuddModel<> jr{ option };

	// Print the model name:
	std::cout << decltype(jr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, jr, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, jr, dt, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}


//
void trimIndexedLattice() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
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

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = trim;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(il, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

void trimIndexedLatticeNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::TrigeorgisModel<> trim{ option };

	// Print the model name:
	std::cout << decltype(trim)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, trim, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, trim, dt, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

//
void tmIndexedLattice() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
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

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = tm;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(il, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}


void tmIndexedLatticeNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::TianModel<> tm{ option };

	// Print the model name:
	std::cout << decltype(tm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, tm, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = tm;
	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, tm, dt, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

//
//
void lrIndexedLattice() {

	lattice_miscellaneous::OptionData<double> option;
	using lattice_model_components::leisen_reimer_inversion::PeizerPrattSecondInversion;


	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
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

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = lr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(il, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}


void lrIndexedLatticeNew() {

	lattice_miscellaneous::OptionData<double> option;
	using lattice_model_components::leisen_reimer_inversion::PeizerPrattSecondInversion;


	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
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
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, lr, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, lr, dt, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

//
void bmIndexedLattice() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(periods);

	// Create CRR model:
	lattice_model::BoyleModel<> bm{ option };

	// Print the model name:
	std::cout << decltype(bm)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double,double, double, double> fwdGen = bm;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double,double, double, double> backGen = bm;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(il, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

void bmIndexedLatticeNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(periods);

	// Create CRR model:
	lattice_model::BoyleModel<> bm{ option };

	// Print the model name:
	std::cout << decltype(bm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, bm, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double, double> backGen = bm;
	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, bm, dt, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}




void crrIndexedLatticeAmerican() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ option };

	std::cout << decltype(crr)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = crr;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = crr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	lattice_types::PayoffAdjuster<double&, double> american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};
	lattice_algorithms::backward_induction(il, backGen, call_payoff, american_adjuster, dt);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}


void crrIndexedLatticeAmericanNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ option };

	std::cout << decltype(crr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, crr, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	auto american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, crr, dt, call_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

//
void mcrrIndexedLatticeAmerican() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::ModifiedCoxRubinsteinRossModel<> mcrr{ option,periods };

	std::cout << decltype(mcrr)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = mcrr;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = mcrr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	lattice_types::PayoffAdjuster<double&, double> american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};
	lattice_algorithms::backward_induction(il, backGen, call_payoff, american_adjuster, dt);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}


void mcrrIndexedLatticeAmericanNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::ModifiedCoxRubinsteinRossModel<> mcrr{ option,periods };

	std::cout << decltype(mcrr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, mcrr, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	auto american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, mcrr, dt, call_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

//
void jrIndexedLatticeAmerican() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::JarrowRuddModel<> jr{ option };

	std::cout << decltype(jr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, jr, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = jr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	lattice_types::PayoffAdjuster<double&, double> american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};
	lattice_algorithms::backward_induction(il, backGen, call_payoff, american_adjuster, dt);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

void jrIndexedLatticeAmericanNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::JarrowRuddModel<> jr{ option };

	std::cout << decltype(jr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, jr, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = jr;
	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	auto american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, jr, dt, call_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

//
void trimIndexedLatticeAmerican() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::TrigeorgisModel<> trim{ option };

	std::cout << decltype(trim)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = trim;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = trim;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	lattice_types::PayoffAdjuster<double&, double> american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};
	lattice_algorithms::backward_induction(il, backGen, call_payoff, american_adjuster, dt);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

void trimIndexedLatticeAmericanNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::TrigeorgisModel<> trim{ option };

	std::cout << decltype(trim)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, trim, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = trim;
	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	auto american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, trim, dt, call_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

//
void tmIndexedLatticeAmerican() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::TianModel<> tm{ option };

	std::cout << decltype(tm)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = tm;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = tm;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	lattice_types::PayoffAdjuster<double&, double> american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};
	lattice_algorithms::backward_induction(il, backGen, call_payoff, american_adjuster, dt);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

void tmIndexedLatticeAmericanNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::TianModel<> tm{ option };

	std::cout << decltype(tm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, tm, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = tm;
	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	auto american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, tm, dt, call_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}
//
void lrIndexedLatticeAmerican() {

	using lattice_model_components::leisen_reimer_inversion::PeizerPrattSecondInversion;
	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Prepare inversion functor:
	std::size_t numerTimePoints{ periods + 1 };
	PeizerPrattSecondInversion<> ppi{ numerTimePoints };

	// Create CRR model:
	lattice_model::LeisenReimerModel<> lr{ option,periods,ppi };

	std::cout << decltype(lr)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = lr;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = lr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	lattice_types::PayoffAdjuster<double&, double> american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};
	lattice_algorithms::backward_induction(il, backGen, call_payoff, american_adjuster, dt);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

void lrIndexedLatticeAmericanNew() {

	using lattice_model_components::leisen_reimer_inversion::PeizerPrattSecondInversion;
	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 100 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Prepare inversion functor:
	std::size_t numerTimePoints{ periods + 1 };
	PeizerPrattSecondInversion<> ppi{ numerTimePoints };

	// Create CRR model:
	lattice_model::LeisenReimerModel<> lr{ option,periods,ppi };

	std::cout << decltype(lr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, lr, dt, option.Underlying);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = lr;
	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	auto american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(il, lr, dt, call_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(il, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}
//
//
void crrLattice() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::CoxRubinsteinRossModel<> crr{ option };

	// Name of the model:
	std::cout << decltype(crr)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = crr;
	lattice_algorithms::forward_induction(la, fwdGen, option.Underlying, timeDeltas);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double,double> backGen = crr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(la, backGen, call_payoff, timeDeltas);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}


void crrLatticeNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::CoxRubinsteinRossModel<> crr{ option };

	// Name of the model:
	std::cout << decltype(crr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date,std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, crr, timeDeltas, option.Underlying);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la,crr, timeDeltas, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}
//
void mcrrLattice() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::ModifiedCoxRubinsteinRossModel<> mcrr{ option,periods };

	// Name of the model:
	std::cout << decltype(mcrr)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = mcrr;
	lattice_algorithms::forward_induction(la, fwdGen, option.Underlying, timeDeltas);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = mcrr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(la, backGen, call_payoff, timeDeltas);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}

void mcrrLatticeNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::ModifiedCoxRubinsteinRossModel<> mcrr{ option,periods };

	// Name of the model:
	std::cout << decltype(mcrr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, mcrr, timeDeltas, option.Underlying);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, mcrr, timeDeltas, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}
//
//
void jrLattice() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::JarrowRuddModel<> jr{ option };

	// Name of the model:
	std::cout << decltype(jr)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = jr;
	lattice_algorithms::forward_induction(la, fwdGen, option.Underlying, timeDeltas);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = jr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(la, backGen, call_payoff, timeDeltas);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}

void jrLatticeNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::JarrowRuddModel<> jr{ option };

	// Name of the model:
	std::cout << decltype(jr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, jr, timeDeltas, option.Underlying);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, jr, timeDeltas, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}
//
void trimLattice() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::TrigeorgisModel<> trim{ option };

	// Name of the model:
	std::cout << decltype(trim)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = trim;
	lattice_algorithms::forward_induction(la, fwdGen, option.Underlying, timeDeltas);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = trim;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(la, backGen, call_payoff, timeDeltas);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}

void trimLatticeNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::TrigeorgisModel<> trim{ option };

	// Name of the model:
	std::cout << decltype(trim)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, trim, timeDeltas, option.Underlying);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, trim, timeDeltas, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}
//
//
void tmLattice() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::TianModel<> tm{ option };

	// Name of the model:
	std::cout << decltype(tm)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = tm;
	lattice_algorithms::forward_induction(la, fwdGen, option.Underlying, timeDeltas);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = tm;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(la, backGen, call_payoff, timeDeltas);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}

void tmLatticeNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::TianModel<> tm{ option };

	// Name of the model:
	std::cout << decltype(tm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, tm, timeDeltas, option.Underlying);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, tm, timeDeltas, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}
//
//
void lrLattice() {

	using ::lattice_model_components::leisen_reimer_inversion::PeizerPrattSecondInversion;
	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	// Prepare inversion formula:
	std::size_t numberTimePoints{ periods + 1 };
	PeizerPrattSecondInversion<> ppi{ numberTimePoints };


	lattice_model::LeisenReimerModel<> lr{ option,periods,ppi };

	// Name of the model:
	std::cout << decltype(lr)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = lr;
	lattice_algorithms::forward_induction(la, fwdGen, option.Underlying, timeDeltas);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = lr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(la, backGen, call_payoff, timeDeltas);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}

void lrLatticeNew() {

	using ::lattice_model_components::leisen_reimer_inversion::PeizerPrattSecondInversion;
	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	// Prepare inversion formula:
	std::size_t numberTimePoints{ periods + 1 };
	PeizerPrattSecondInversion<> ppi{ numberTimePoints };


	lattice_model::LeisenReimerModel<> lr{ option,periods,ppi };

	// Name of the model:
	std::cout << decltype(lr)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, lr, timeDeltas, option.Underlying);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, lr, timeDeltas, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}
//
//
void bmLattice() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Trinomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::BoyleModel<> bm{ option };

	// Name of the model:
	std::cout << decltype(bm)::name() << "\n";

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double, double> fwdGen = bm;
	lattice_algorithms::forward_induction(la, fwdGen, option.Underlying, timeDeltas);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double, double> backGen = bm;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(la, backGen, call_payoff, timeDeltas);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}

void bmLatticeNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.05;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods{ 100 };

	for (std::size_t t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Trinomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::BoyleModel<> bm{ option };

	// Name of the model:
	std::cout << decltype(bm)::name() << "\n";

	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, bm, timeDeltas, option.Underlying);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, bm, timeDeltas, call_payoff);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}
//
//
//
void crrLatticeAmerican() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t < 101; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}


	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::CoxRubinsteinRossModel<> crr{ option };


	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = crr;
	lattice_algorithms::forward_induction(la, fwdGen, option.Underlying, timeDeltas);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = crr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	lattice_types::PayoffAdjuster<double&, double> american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};
	lattice_algorithms::backward_induction(la, backGen, call_payoff,american_adjuster,timeDeltas);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}


void crrLatticeAmericanNew() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);

	for (std::size_t t = 1; t < 101; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}


	// Creating lattice:
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date> la = { fixingDates };

	// Create CRR model:
	double daysInYear{ 365.0 };
	auto fd = la.fixingDates();
	std::vector<double> timeDeltas(fd.size() - 1);
	for (auto i = 0; i < timeDeltas.size(); ++i) {
		timeDeltas[i] = ((fd[i + 1] - fd[i]).days() / daysInYear);
	}

	lattice_model::CoxRubinsteinRossModel<> crr{ option };


	// Forward induction:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, crr, timeDeltas, option.Underlying);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	// Prepare american adujster:
	auto american_adjuster = [&call_payoff](double& value, double stock) {
		value = std::max(value, call_payoff(stock));
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(la, crr, timeDeltas, call_payoff, american_adjuster);

	// Print the part of generated lattice:
	lattice_utility::print(la, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}


#endif ///_LATTICE_MODEL_T