#pragma once
#if !defined(_LATTICE_GREEKS_T)
#define _LATTICE_GREEKS_T

#include"lattice_algorithms.h"
#include"lattice_utility.h"
#include<iostream>
#include"lattice_model_components.h"
#include"lattice_greeks.h"
#include"lattice_model.h"
#include<iostream>

#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>
using namespace boost::gregorian;

// ===============================================================================
// ========================== Define some arbitrary model ========================
// ===============================================================================

template<typename T = double>
class ArbitraryBinomialModel :public lattice_model::BinomialModel<1, T> {
public:
	// Forward generator
	std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false) override {
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	}

	// Backward generator
	T operator()(T currValue, T upValue, T downValue, T dt) override {
		double p = 0.5;
		return (p*upValue + (1.0 - p)*downValue);
	}

	static std::string const name() {
		return std::string{ "Arbitrary binomial model" };
	}
};

template<typename T = double>
class ArbitraryTrinomialModel :public lattice_model::TrinomialModel<1, T> {
public:
	// Forward generator
	std::tuple<T, T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false) override {
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	}

	// Backward generator
	T operator()(T currValue, T upValue, T midValue, T downValue, T dt) override {
		double p = (1.0 / 3.0);
		return (p*upValue + p * midValue + p * downValue);
	}

	static std::string const name() {
		return std::string{ "Arbitrary trinomial model" };
	}
};
// ===============================================================================
// ===============================================================================




void indexedLatticeBinomialDelta() {

	const int N = 5;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";

	lattice_utility::print(il, il.cbegin(), il.cend());

	// Instantiate a binomial model:
	ArbitraryBinomialModel<> model;
	std::cout << "Model: " << decltype(model)::name() << "\n";


	// Forward induction algo:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, std::size_t, double, double> forward_binomial_induction;
	forward_binomial_induction fwd_induction;
	fwd_induction(il, model, dt, 1.0);


	std::cout << "After forward induction step:\n";
	lattice_utility::print(il, il.cbegin(), il.cend());
	// Make a copy:
	auto pt = il;

	auto payoff = [](double S) {
		return S;
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;

	backward_binomial_induction bwrd_induction;
	bwrd_induction(pt, model, dt, payoff);

	std::cout << "After backward induction step:\n";
	lattice_utility::print(pt, pt.cbegin(), pt.cend());

	std::cout << "Price: " << pt.apex() << "\n\n";

	typedef lattice_greeks::Greeks<std::size_t> indexed_greeks;

	std::cout << "Delta of call: " << indexed_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << indexed_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << indexed_greeks::theta(il, pt, dt) << "\n";

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
	// Forward induction algo:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, std::size_t, double, double> forward_binomial_induction;
	forward_binomial_induction fwd_induction;

	fwd_induction(il, crr, dt, 1.0);

	// make a copy for later greeks: 
	auto pt(il);
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

	backward_binomial_induction bwrd_induction;
	bwrd_induction(pt, crr, dt, call_payoff);


	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	typedef lattice_greeks::Greeks<std::size_t> indexed_greeks;

	std::cout << "Delta of call: " << indexed_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << indexed_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << indexed_greeks::theta(il, pt, dt) << "\n";

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
	// Forward induction algo:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, std::size_t, double, double> forward_binomial_induction;
	forward_binomial_induction fwd_induction;
	fwd_induction(il, mcrr, dt, 1.0);

	// make a copy for later greeks: 
	auto pt(il);
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

	backward_binomial_induction bwrd_induction;
	bwrd_induction(pt, mcrr, dt, call_payoff);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	typedef lattice_greeks::Greeks<std::size_t> indexed_greeks;

	std::cout << "Delta of call: " << indexed_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << indexed_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << indexed_greeks::theta(il, pt, dt) << "\n";
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
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, std::size_t, double, double> forward_binomial_induction;
	forward_binomial_induction fwd_induction;
	fwd_induction(il, jr, dt, 1.0);

	// make a copy for later greeks: 
	auto pt(il);
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

	backward_binomial_induction bwrd_induction;
	bwrd_induction(pt, jr, dt, call_payoff);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	typedef lattice_greeks::Greeks<std::size_t> indexed_greeks;

	std::cout << "Delta of call: " << indexed_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << indexed_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << indexed_greeks::theta(il, pt, dt) << "\n";
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
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, std::size_t, double, double> forward_binomial_induction;
	forward_binomial_induction fwd_induction;
	fwd_induction(il, trim, dt, 1.0);

	// make a copy for later greeks: 
	auto pt(il);
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

	backward_binomial_induction bwrd_induction;
	bwrd_induction(pt, trim, dt, call_payoff);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	typedef lattice_greeks::Greeks<std::size_t> indexed_greeks;

	std::cout << "Delta of call: " << indexed_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << indexed_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << indexed_greeks::theta(il, pt, dt) << "\n";
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
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, std::size_t, double, double> forward_binomial_induction;
	forward_binomial_induction fwd_induction;
	fwd_induction(il, tm, dt, 1.0);

	// make a copy for later greeks: 
	auto pt(il);
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

	backward_binomial_induction bwrd_induction;
	bwrd_induction(pt, tm, dt, call_payoff);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	typedef lattice_greeks::Greeks<std::size_t> indexed_greeks;

	std::cout << "Delta of call: " << indexed_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << indexed_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << indexed_greeks::theta(il, pt, dt) << "\n";
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
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, std::size_t, double, double> forward_binomial_induction;
	forward_binomial_induction fwd_induction;
	fwd_induction(il, lr, dt, 1.0);

	// make a copy for later greeks: 
	auto pt(il);
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

	backward_binomial_induction bwrd_induction;
	bwrd_induction(pt, lr, dt, call_payoff);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	typedef lattice_greeks::Greeks<std::size_t> indexed_greeks;

	std::cout << "Delta of call: " << indexed_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << indexed_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << indexed_greeks::theta(il, pt, dt) << "\n";
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
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial, std::size_t, double, double> forward_binomial_induction;
	forward_binomial_induction fwd_induction;
	fwd_induction(il, bm, dt, 1.0);

	// make a copy for later greeks: 
	auto pt(il);
	// Print the part of generated lattice:
	auto first = il.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(il, first, last);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };
	
	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double> backward_binomial_induction;

	backward_binomial_induction bwrd_induction;
	bwrd_induction(pt, bm, dt, call_payoff);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	typedef lattice_greeks::Greeks<std::size_t> indexed_greeks;

	std::cout << "Delta of call: " << indexed_greeks::delta(il, pt) << "\n";
	std::cout << "Gamma of call: " << indexed_greeks::gamma(il, pt) << "\n";
	std::cout << "Theta of call: " << indexed_greeks::theta(il, pt, dt) << "\n";
}


void testIndexedLatticesGreeks() {
	std::cout << "=======================================================\n";
	std::cout << "=========== Indexed Lattice Greeks - TEST =============\n";
	std::cout << "=======================================================\n";

	indexedLatticeBinomialDelta();
	crrIndexedLatticeGreeks();
	mcrrIndexedLatticeGreeks();
	jrIndexedLatticeGreeks();
	trimIndexedLatticeGreeks();
	tmIndexedLatticeGreeks();
	lrIndexedLatticeGreeks();
	bmIndexedLatticeGreeks();

	std::cout << "=======================================================\n";
}

void crrLatticeGreeks() {

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
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, crr, timeDeltas, option.Underlying);

	// Print the part of generated lattice:
	auto first = la.begin();
	auto last = std::next(first, 5);
	lattice_utility::print(la, first, last);


	// make a copy for later greeks: 
	auto pt(la);

	// Backward induction:

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(pt, crr, timeDeltas, call_payoff);

	// Print apex: value of option:
	std::cout << "Price of call: " << la.apex() << "\n\n";
	
	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";

	typedef lattice_greeks::Greeks<date> date_greeks;

	std::cout << "Delta of call: " << date_greeks::delta(la, pt) << "\n";
	std::cout << "Gamma of call: " << date_greeks::gamma(la, pt) << "\n";
	std::cout << "Theta of call: " << date_greeks::theta(la, pt, timeDeltas) << "\n";

}

void mcrrLatticeGreeks() {

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
	// make a copy for later greeks: 
	auto pt(la);

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(pt, mcrr, timeDeltas, call_payoff);


	// Print apex: value of option:
	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	typedef lattice_greeks::Greeks<date> date_greeks;

	std::cout << "Delta of call: " << date_greeks::delta(la, pt) << "\n";
	std::cout << "Gamma of call: " << date_greeks::gamma(la, pt) << "\n";
	std::cout << "Theta of call: " << date_greeks::theta(la, pt, timeDeltas) << "\n";
}
//


void jrLatticeGreeks() {

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
	// make a copy for later greeks: 
	auto pt(la);

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(pt, jr, timeDeltas, call_payoff);

	// Print apex: value of option:
	// Print apex: value of option:
	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	typedef lattice_greeks::Greeks<date> date_greeks;

	std::cout << "Delta of call: " << date_greeks::delta(la, pt) << "\n";
	std::cout << "Gamma of call: " << date_greeks::gamma(la, pt) << "\n";
	std::cout << "Theta of call: " << date_greeks::theta(la, pt, timeDeltas) << "\n";
}

void trimLatticeGreeks() {

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
	// make a copy for later greeks: 
	auto pt(la);

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(pt, trim, timeDeltas, call_payoff);

	// Print apex: value of option:
	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	typedef lattice_greeks::Greeks<date> date_greeks;

	std::cout << "Delta of call: " << date_greeks::delta(la, pt) << "\n";
	std::cout << "Gamma of call: " << date_greeks::gamma(la, pt) << "\n";
	std::cout << "Theta of call: " << date_greeks::theta(la, pt, timeDeltas) << "\n";
}
//

void tmLatticeGreeks() {

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

	// make a copy for later greeks: 
	auto pt(la);

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(pt, tm, timeDeltas, call_payoff);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	typedef lattice_greeks::Greeks<date> date_greeks;

	std::cout << "Delta of call: " << date_greeks::delta(la, pt) << "\n";
	std::cout << "Gamma of call: " << date_greeks::gamma(la, pt) << "\n";
	std::cout << "Theta of call: " << date_greeks::theta(la, pt, timeDeltas) << "\n";
}
//


void lrLatticeGreeks() {

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

	// make a copy for later greeks: 
	auto pt(la);

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(pt, lr, timeDeltas, call_payoff);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	typedef lattice_greeks::Greeks<date> date_greeks;

	std::cout << "Delta of call: " << date_greeks::delta(la, pt) << "\n";
	std::cout << "Gamma of call: " << date_greeks::gamma(la, pt) << "\n";
	std::cout << "Theta of call: " << date_greeks::theta(la, pt, timeDeltas) << "\n";
}
//


void bmLatticeGreeks() {

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

	// make a copy for later greeks: 
	auto pt(la);

	// Prepare payoff:
	double K = option.Strike;
	auto call_payoff = [&K](double stock) {return std::max(K - stock, 0.0); };

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>> backward_binomial_induction;
	backward_binomial_induction brd_induction;

	brd_induction(pt, bm, timeDeltas, call_payoff);

	// Print the part of generated lattice:
	first = pt.begin();
	last = std::next(first, 5);
	lattice_utility::print(pt, first, last);
	// Print apex: value of option:
	std::cout << "Price of call: " << pt.apex() << "\n";
	typedef lattice_greeks::Greeks<date> date_greeks;

	std::cout << "Delta of call: " << date_greeks::delta(la, pt) << "\n";
	std::cout << "Gamma of call: " << date_greeks::gamma(la, pt) << "\n";
	std::cout << "Theta of call: " << date_greeks::theta(la, pt, timeDeltas) << "\n";

}





void testLatticesGreeks() {
	std::cout << "=======================================================\n";
	std::cout << "=================== Lattice Greeks - TEST =============\n";
	std::cout << "=======================================================\n";

	crrLatticeGreeks();
	mcrrLatticeGreeks();
	jrLatticeGreeks();
	trimLatticeGreeks();
	tmLatticeGreeks();
	lrLatticeGreeks();
	bmLatticeGreeks();

	std::cout << "=======================================================\n";
}





#endif ///_LATTICE_GREEKS_T