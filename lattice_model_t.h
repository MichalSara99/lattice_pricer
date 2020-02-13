#pragma once
#if !defined(_LATTICE_MODEL_T)
#define _LATTICE_MODEL_T

#include"lattice_structure.h"
#include"lattice_model.h"
#include"lattice_miscellaneous.h"
#include"lattice_utility.h"
#include"lattice_algorithms.h"

#include<iostream>
#include<numeric>

void crrIndexedLattice() {
	
	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 101 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ option};

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

void crrIndexedLatticeMod() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 101 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ option };

	// Forward induction:
	lattice_types::LeafForwardGenerator<double, double, double> fwdGen = crr;
	lattice_algorithms::forward_induction(il, fwdGen, option.Underlying, dt);

	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Make a copy:
	auto il_readonly = il;
	// Backward induction:
	lattice_types::LeafBackwardGenerator<double, double, double, double> backGen = crr;
	// Prepare payoff:
	double K = option.Strike;
	lattice_types::Payoff<double, double> call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	lattice_algorithms::backward_induction(il_readonly, il, backGen, call_payoff, dt);

	// Print the part of generated lattice:
	std::cout << "Option price tree: \n";
	lattice_utility::print(il, first, last);

	// Print the original asset price tree:
	std::cout << "Original asset price tree: \n";
	auto first_c = il_readonly.begin();
	auto last_c = std::next( first_c, 5);
	lattice_utility::print(il_readonly, first_c, last_c);

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

	std::size_t periods{ 101 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ option };

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

void crrIndexedLatticeAmericanMod() {

	lattice_miscellaneous::OptionData<double> option;

	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 60.0;

	std::size_t periods{ 101 };
	double maturity = 0.29;
	double dt = (maturity / double(periods));


	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(periods);

	// Create CRR model:
	lattice_model::CoxRubinsteinRossModel<> crr{ option };

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

	// Make a copy:
	auto il_copy = il;
	lattice_algorithms::backward_induction(il_copy, il, backGen, call_payoff, american_adjuster, dt);

	// Print the part of generated lattice:
	std::cout << "Option price lattice:\n";
	lattice_utility::print(il, first, last);

	// Print the original asset price lattice:
	std::cout << "Option price lattice:\n";
	auto first_c = il_copy.begin();
	auto last_c = std::next(first_c, 5);
	lattice_utility::print(il_copy, first_c, last_c);
	// Print apex: value of option:
	std::cout << "Price of call: " << il.apex() << "\n\n";
}

void crrLattice() {

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


void crrLatticeMod() {

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

	// Make a copy:
	auto la_copy = la;
	lattice_algorithms::backward_induction(la_copy,la, backGen, call_payoff, timeDeltas);


	// Print the part of generated lattice:
	std::cout << "Option price lattice: \n";
	lattice_utility::print(la, first, last);

	// Print the oricinal assetv price lattice:
	std::cout << "Asset price lattice: \n";
	auto first_c = la_copy.begin();
	auto last_c = std::next(first_c, 5);
	lattice_utility::print(la_copy, first_c, last_c);
	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}


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

void crrLatticeAmericanMod() {

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
	// Make a copy:
	auto la_copy = la;

	lattice_algorithms::backward_induction(la_copy,la, backGen, call_payoff, american_adjuster, timeDeltas);

	// Print the part of generated lattice:
	std::cout << "Option price lattice: \n";
	lattice_utility::print(la, first, last);

	std::cout << "Asset price lattice: \n";
	auto first_c = la_copy.begin();
	auto last_c = std::next(first_c, 5);
	lattice_utility::print(la_copy, first_c, last_c);

	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";

}




#endif ///_LATTICE_MODEL_T