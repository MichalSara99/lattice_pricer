#pragma once
#if !defined(_LATTICE_ALGORITHMS_T)
#define _LATTICE_ALGORITHMS_T

#include"lattice_model.h"
#include"lattice_algorithms.h"
#include"lattice_utility.h"
#include<iostream>
#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>
#include<boost/math/special_functions/binomial.hpp>


using namespace boost::gregorian;

// ===============================================================================
// ========================== Define some arbitrary model ========================
// ===============================================================================

template<typename T = double>
class AbitraryBinomialModel :public lattice_model::BinomialModel<1, T> {
public:
	// Forward generator
	std::tuple<T, T> operator()(T value, T dt,std::size_t leafIdx,std::size_t timeIdx, bool isMeanReverting = false) override {
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
class AbitraryTrinomialModel :public lattice_model::TrinomialModel<1, T> {
public:
	// Forward generator
	std::tuple<T, T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false) override {
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	}

	// Backward generator
	T operator()(T currValue,T downValue,T midValue, T upValue, T dt,
		std::size_t revertBranchesSize, std::size_t nodesSize, std::size_t leafIdx) override {
		double p = (1.0 / 3.0);
		return (p*upValue + p * midValue + p * downValue);
	}

	static std::string const name() {
		return std::string{ "Arbitrary trinomial model" };
	}
};
// ===============================================================================
// ===============================================================================


// =====================================================================================
// ============================= Testing forward-induction algos =======================
// =====================================================================================

void indexedLatticeBinomialForwardInduction() {

	const int N = 12;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";

	lattice_utility::print(il, il.cbegin(), il.cend());

	// Instantiate a binomial model:
	AbitraryBinomialModel<> model;
	std::cout << "Model: " << decltype(model)::name() << "\n";

	// Forward induction algo:
	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, std::size_t, double, double> forward_binomial_induction;
	forward_binomial_induction fwd_induction;
	fwd_induction(il, model, dt,1.0);

	std::cout << "After forward induction step:\n";
	lattice_utility::print(il, il.cbegin(), il.cend());

	std::cout << "Price: " << il.apex() << "\n\n";
}


void indexedLatticeTrinomialForwardInduction() {
	const int N = 12;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	lattice_utility::print(il, il.cbegin(), il.cend());

	// Instantiate a trinomial model:
	AbitraryTrinomialModel<> model;
	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial, std::size_t, double, double> forward_trinomial_induction;

	forward_trinomial_induction fwd_induction;
	fwd_induction(il, model, dt,1.0 );
	std::cout << "After forward induction step:\n";
	lattice_utility::print(il, il.cbegin(), il.cend());

	std::cout << "Price: " << il.apex() << "\n\n";
}


void latticeBinomialForwardInduction() {

	auto today = date(day_clock::local_day());
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date>
		la = { today,today + date_duration(1),
		today + date_duration(2),today + date_duration(3) ,
		today + date_duration(4) ,today + date_duration(5) ,
		today + date_duration(6) ,today + date_duration(7),
		today + date_duration(8) ,today + date_duration(9) ,
		today + date_duration(10) ,today + date_duration(11),
		today + date_duration(12) };

	// get fixing dates and compute deltas:
	auto fixingDates = la.fixingDates();
	double const daysInYear{ 365.0 };
	std::vector<double> deltas(fixingDates.size() - 1);
	for (std::size_t t = 0; t < fixingDates.size() - 1; ++t) {
		deltas[t] = ((fixingDates[t + 1] - fixingDates[t]).days() / daysInYear);
	}

	lattice_utility::print(la, la.cbegin(), la.cend());

	// Instantiate a binomial model:
	AbitraryBinomialModel<> model;
	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, model, deltas, 1.0);

	std::cout << "After forward induction step:\n";
	lattice_utility::print(la, la.cbegin(), la.cend());

	std::cout << "Price: " << la.apex() << "\n\n";
}



void latticeTrinomialForwardInduction() {

	auto today = date(day_clock::local_day());
	lattice_structure::Lattice<lattice_types::LatticeType::Trinomial, double, date>
		la = { today,today + date_duration(1),
		today + date_duration(2),today + date_duration(3) ,
		today + date_duration(4) ,today + date_duration(5) ,
		today + date_duration(6) ,today + date_duration(7),
		today + date_duration(8) ,today + date_duration(9) ,
		today + date_duration(10) ,today + date_duration(11),
		today + date_duration(12) };

	// get fixing dates and compute deltas:
	auto fixingDates = la.fixingDates();
	double const daysInYear{ 365.0 };
	std::vector<double> deltas(fixingDates.size() - 1);
	for (std::size_t t = 0; t < fixingDates.size() - 1; ++t) {
		deltas[t] = ((fixingDates[t + 1] - fixingDates[t]).days() / daysInYear);
	}

	lattice_utility::print(la, la.cbegin(), la.cend());


	// Instantiate a trinomial model:
	AbitraryTrinomialModel<> model;

	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>, double> forward_trinomial_induction;

	forward_trinomial_induction fwd_induction;
	fwd_induction(la, model, deltas, 1.0);

	std::cout << "After forward induction step:\n";
	lattice_utility::print(la, la.cbegin(), la.cend());

	std::cout << "Price: " << la.apex() << "\n\n";
}



void testIndexedForwardInduction() {
	std::cout << "=====================================================\n";
	std::cout << "========= Indexed Forward Induction - TEST ==========\n";
	std::cout << "=====================================================\n";

	indexedLatticeBinomialForwardInduction();
	indexedLatticeTrinomialForwardInduction();


	std::cout << "=====================================================\n";
}

void testForwardInduction() {
	std::cout << "=====================================================\n";
	std::cout << "============ Forward Induction - TEST ===============\n";
	std::cout << "=====================================================\n";

	latticeBinomialForwardInduction();
	latticeTrinomialForwardInduction();


	std::cout << "=====================================================\n";
}


// ==============================================================================================
// ================== Testing forward-induction with discretely paid dividends ==================
// ==============================================================================================



void indexedLatticeBinomialForwardInductionDividend() {

	const int N = 12;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);
	// dicrete dividends:
	std::map<std::size_t, double> dividends = { { 3,0.01 },{ 9,0.015 } };

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";

	lattice_utility::print(il, il.cbegin(), il.cend());

	// Instantiate a trinomial model:
	AbitraryBinomialModel<> model;

	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, model, dt, 1.0, dividends);
	std::cout << "After forward induction step:\n";
	lattice_utility::print(il, il.cbegin(), il.cend());

	std::cout << "Price: " << il.apex() << "\n\n";
}

//


void indexedLatticeTrinomialForwardInductionDividends() {
	const int N = 12;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);
	std::map<std::size_t, double> dividends = { { 3,0.01 },{ 9,0.015 } };

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	lattice_utility::print(il, il.cbegin(), il.cend());

	// Instantiate a trinomial model:
	AbitraryTrinomialModel<> model;

	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double, double> forward_trinomial_induction;

	forward_trinomial_induction fwd_induction;
	fwd_induction(il, model, dt, 1.0, dividends);
	std::cout << "After forward induction step:\n";
	lattice_utility::print(il, il.cbegin(), il.cend());

	std::cout << "Price: " << il.apex() << "\n\n";
}
//
//
void testIndexedForwardInductionDividends() {
	std::cout << "=====================================================\n";
	std::cout << "===== Indexed Forward Induction Dividends - TEST ====\n";
	std::cout << "=====================================================\n";

	indexedLatticeBinomialForwardInductionDividend();
	indexedLatticeTrinomialForwardInductionDividends();


	std::cout << "=====================================================\n";
}
//
//
void latticeBinomialForwardInductionDividends() {

	auto today = date(day_clock::local_day());
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date>
		la = { today,today + date_duration(1),
		today + date_duration(2),today + date_duration(3) ,
		today + date_duration(4) ,today + date_duration(5) ,
		today + date_duration(6) ,today + date_duration(7),
		today + date_duration(8) ,today + date_duration(9) ,
		today + date_duration(10) ,today + date_duration(11),
		today + date_duration(12) };

	std::map<date, double> dividends = { { today + date_duration(3),0.01 },
	{ today + date_duration(9),0.015 } };

	// get fixing dates and compute deltas:
	auto fixingDates = la.fixingDates();
	double const daysInYear{ 365.0 };
	std::vector<double> deltas(fixingDates.size() - 1);
	for (std::size_t t = 0; t < fixingDates.size() - 1; ++t) {
		deltas[t] = ((fixingDates[t + 1] - fixingDates[t]).days() / daysInYear);
	}

	lattice_utility::print(la, la.cbegin(), la.cend());


	// Instantiate a trinomial model:
	AbitraryBinomialModel<> model;

	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date,std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, model, deltas, 1.0, dividends);

	std::cout << "After forward induction step:\n";
	lattice_utility::print(la, la.cbegin(), la.cend());

	std::cout << "Price: " << la.apex() << "\n\n";
}




void latticeTrinomialForwardInductionDividend() {

	auto today = date(day_clock::local_day());
	lattice_structure::Lattice<lattice_types::LatticeType::Trinomial, double, date>
		la = { today,today + date_duration(1),
		today + date_duration(2),today + date_duration(3) ,
		today + date_duration(4) ,today + date_duration(5) ,
		today + date_duration(6) ,today + date_duration(7),
		today + date_duration(8) ,today + date_duration(9) ,
		today + date_duration(10) ,today + date_duration(11),
		today + date_duration(12) };

	std::map<date, double> dividends = { { today + date_duration(3),0.01 },
	{ today + date_duration(9),0.015 } };

	// get fixing dates and compute deltas:
	auto fixingDates = la.fixingDates();
	double const daysInYear{ 365.0 };
	std::vector<double> deltas(fixingDates.size() - 1);
	for (std::size_t t = 0; t < fixingDates.size() - 1; ++t) {
		deltas[t] = ((fixingDates[t + 1] - fixingDates[t]).days() / daysInYear);
	}

	lattice_utility::print(la, la.cbegin(), la.cend());


	// Instantiate a trinomial model:
	AbitraryTrinomialModel<> model;

	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>, double> forward_trinomial_induction;

	forward_trinomial_induction fwd_induction;
	fwd_induction(la, model, deltas, 1.0, dividends);

	std::cout << "After forward induction step:\n";
	lattice_utility::print(la, la.cbegin(), la.cend());

	std::cout << "Price: " << la.apex() << "\n\n";
}

void testForwardInductionDividends() {
	std::cout << "=====================================================\n";
	std::cout << "======= Forward Induction Dividends - TEST ==========\n";
	std::cout << "=====================================================\n";

	latticeBinomialForwardInductionDividends();
	latticeTrinomialForwardInductionDividend();


	std::cout << "=====================================================\n";
}



//// ======================================================================================
//// ============================= Testing backward-induction algos =======================
//// ======================================================================================
//
//
void indexedLatticeBinomialBackwardInduction() {

	const int N = 3;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";

	lattice_utility::print(il, il.cbegin(), il.cend());

	// Instantiate a trinomial model:
	AbitraryBinomialModel<> model;

	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t,double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, model,dt, 1.0);

	std::cout << "After forward induction step:\n";
	lattice_utility::print(il, il.cbegin(), il.cend());

	auto payoff = [](double S) {
		return S;
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;

	backward_binomial_induction bwrd_induction;
	bwrd_induction(il, model, dt, payoff);

	std::cout << "After backward induction step:\n";
	lattice_utility::print(il, il.cbegin(), il.cend());

	std::cout << "Price: " << il.apex() << "\n\n";
}

//
//

void indexedLatticeTrinomialBackwardInduction() {

	const int N = 3;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	lattice_utility::print(il, il.cbegin(), il.cend());

	// Instantiate a trinomial model:
	AbitraryTrinomialModel<> model;

	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double, double> forward_trinomial_induction;

	forward_trinomial_induction fwd_induction;
	fwd_induction(il, model, dt, 1.0);
	std::cout << "After forward induction step:\n";
	lattice_utility::print(il, il.cbegin(), il.cend());

	auto payoff = [](double S) {
		return S;
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double> backward_trinomial_induction;

	backward_trinomial_induction bwrd_induction;
	bwrd_induction(il, model, dt, payoff);

	std::cout << "After backward induction step:\n";
	lattice_utility::print(il, il.cbegin(), il.cend());

	std::cout << "Price: " << il.apex() << "\n\n";

}


void testIndexedBackwardInductionDividends() {
	std::cout << "=====================================================\n";
	std::cout << "======== Indexed Backward Induction  - TEST =========\n";
	std::cout << "=====================================================\n";

	indexedLatticeBinomialBackwardInduction();
	indexedLatticeTrinomialBackwardInduction();


	std::cout << "=====================================================\n";
}


//
void latticeBinomialBackwardInduction() {


	auto today = date(day_clock::local_day());
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date>
		la = { today,today + date_duration(2),today,today + date_duration(1) };

	// get fixing dates and compute deltas:
	auto fixingDates = la.fixingDates();
	double const daysInYear{ 365.0 };
	std::vector<double> deltas(fixingDates.size() - 1);
	for (std::size_t t = 0; t < fixingDates.size() - 1; ++t) {
		deltas[t] = ((fixingDates[t + 1] - fixingDates[t]).days() / daysInYear);
	}

	lattice_utility::print(la, la.cbegin(), la.cend());

	// Instantiate a trinomial model:
	AbitraryBinomialModel<> model;

	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>,double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, model, deltas, 1.0);

	std::cout << "After forward induction step:\n";
	lattice_utility::print(la, la.cbegin(), la.cend());

	auto payoff = [](double S) {
		return S;
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;

	backward_binomial_induction bwrd_induction;
	bwrd_induction(la, model,deltas, payoff);

	std::cout << "After backward induction step:\n";
	lattice_utility::print(la, la.cbegin(), la.cend());


	std::cout << "Price: " << la.apex() << "\n\n";
}



void latticeTrinomialBackwardInduction() {

	auto today = date(day_clock::local_day());
	lattice_structure::Lattice<lattice_types::LatticeType::Trinomial, double, date>
		la = { today,today + date_duration(2),today,today + date_duration(1) };

	// get fixing dates and compute deltas:
	auto fixingDates = la.fixingDates();
	double const daysInYear{ 365.0 };
	std::vector<double> deltas(fixingDates.size() - 1);
	for (std::size_t t = 0; t < fixingDates.size() - 1; ++t) {
		deltas[t] = ((fixingDates[t + 1] - fixingDates[t]).days() / daysInYear);
	}

	lattice_utility::print(la, la.cbegin(), la.cend());

	// Instantiate a trinomial model:
	AbitraryTrinomialModel<> model;

	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>, double> forward_trinomial_induction;

	forward_trinomial_induction fwd_induction;
	fwd_induction(la, model, deltas, 1.0);

	std::cout << "After forward induction step:\n";
	lattice_utility::print(la, la.cbegin(), la.cend());

	auto payoff = [](double S) {
		return S;
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>> backward_trinomial_induction;

	backward_trinomial_induction bwrd_induction;
	bwrd_induction(la, model, deltas, payoff);

	std::cout << "After backward induction step:\n";
	lattice_utility::print(la, la.cbegin(), la.cend());

	std::cout << "Price: " << la.apex() << "\n\n";
}

void testBackwardInductionDividends() {
	std::cout << "=====================================================\n";
	std::cout << "============ Backward Induction  - TEST =============\n";
	std::cout << "=====================================================\n";

	latticeBinomialBackwardInduction();
	latticeTrinomialBackwardInduction();


	std::cout << "=====================================================\n";
}

// ============================================================================================
// =============================== Testing merging algo =======================================
// ============================================================================================	

void mergeIndexedBinomial() {
	const int N = 3;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> stockTree(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";

	lattice_utility::print(stockTree,stockTree.begin(), stockTree.end());

	// Instantiate a trinomial model:
	AbitraryBinomialModel<> model;

	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, std::size_t, double, double> forward_induction;
	forward_induction fwd_induction;
	fwd_induction(stockTree, model,dt, 1.0);

	std::cout << "After forward induction step:\n";
	lattice_utility::print(stockTree,stockTree.begin(), stockTree.end());
	
	// make a copy for merging
	auto optionTree = stockTree;

	auto payoff = [](double S) {
		return S;
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial, std::size_t, double> backward_induction;
	backward_induction bwd_induction;

	bwd_induction(optionTree, model, dt, payoff);
	std::cout << "After backward induction step:\n";
	lattice_utility::print(optionTree,optionTree.begin(), optionTree.end());

	std::cout << "Price: " << optionTree.apex() << "\n\n";

	std::cout << "Merging stockTree with optionTree:\n";
	typedef lattice_algorithms::MergeLattices<std::size_t> merge;
	merge merge2;

	auto merged = merge2(lattice_types::LaunchPolicy::Parallel, stockTree, optionTree);
	lattice_utility::print(merged, merged.cbegin(), merged.cend());
}
//

void mergeIndexedTrinomial() {

	const int N = 3;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> stockTree(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());

	// Instantiate a trinomial model:
	AbitraryTrinomialModel<> model;

	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial, std::size_t, double, double> forward_induction;
	forward_induction fwd_induction;

	fwd_induction(stockTree, model,dt, 1.0);
	std::cout << "After forward induction step:\n";
	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());

	auto payoff = [](double S) {
		return S;
	};
	// Make a copy for merging
	auto optionTree = stockTree;

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial, std::size_t, double> backward_induction;
	backward_induction bwd_induction;

	bwd_induction(optionTree, model,dt, payoff);
	std::cout << "After backward induction step:\n";
	lattice_utility::print(optionTree,optionTree.begin(), optionTree.end());

	std::cout << "Price: " << optionTree.apex() << "\n\n";

	std::cout << "Merged stockTree and optionTree:\n";
	typedef lattice_algorithms::MergeLattices<std::size_t> merge;
	merge merge2;

	auto merged = merge2(lattice_types::LaunchPolicy::Parallel, stockTree, optionTree);
	lattice_utility::print(merged, merged.cbegin(), merged.cend());

}

void testIndexedMerge() {
	std::cout << "=====================================================\n";
	std::cout << "============ Merging Indexed Trees - TEST ===========\n";
	std::cout << "=====================================================\n";

	mergeIndexedBinomial();
	mergeIndexedTrinomial();


	std::cout << "=====================================================\n";
}


////
////
void mergeBinomial() {

	auto today = date(day_clock::local_day());
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date>
		stockTree = { today,today + date_duration(2),today,today + date_duration(1) };

	// get fixing dates and compute deltas:
	auto fixingDates = stockTree.fixingDates();
	double const daysInYear{ 365.0 };
	std::vector<double> deltas(fixingDates.size() - 1);
	for (std::size_t t = 0; t < fixingDates.size() - 1; ++t) {
		deltas[t] = ((fixingDates[t + 1] - fixingDates[t]).days() / daysInYear);
	}

	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());


	// Instantiate a trinomial model:
	AbitraryBinomialModel<> model;

	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, date, std::vector<double>, double> forward_induction;
	forward_induction fwd_induction;

	fwd_induction(stockTree, model, deltas, 1.0);
	std::cout << "After forward induction step:\n";
	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());


	auto payoff = [](double S) {
		return S;
	};

	auto optionTree = stockTree;

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial, date,std::vector<double>> backward_induction;
	backward_induction bwd_induction;

	bwd_induction(optionTree, model, deltas,payoff);
	std::cout << "After backward induction step:\n";
	lattice_utility::print(optionTree, optionTree.begin(), optionTree.end());


	std::cout << "Price: " << optionTree.apex() << "\n\n";

	typedef lattice_algorithms::MergeLattices<date> merge;
	merge merge3;

	auto merged = merge3(lattice_types::LaunchPolicy::Parallel,stockTree, optionTree, stockTree);
	std::cout << "Merged stockTree and optionTree:\n";

	lattice_utility::print(merged, merged.begin(), merged.end());

}
//
//
void mergeTrinomial() {
	auto today = date(day_clock::local_day());
	lattice_structure::Lattice<lattice_types::LatticeType::Trinomial, double, date>
		stockTree = { today,today + date_duration(2),today,today + date_duration(1) };

	// get fixing dates and compute deltas:
	auto fixingDates = stockTree.fixingDates();
	double const daysInYear{ 365.0 };
	std::vector<double> deltas(fixingDates.size() - 1);
	for (std::size_t t = 0; t < fixingDates.size() - 1; ++t) {
		deltas[t] = ((fixingDates[t + 1] - fixingDates[t]).days() / daysInYear);
	}

	lattice_utility::print(stockTree,stockTree.begin(), stockTree.end());


	// Instantiate a trinomial model:
	AbitraryTrinomialModel<> model;

	std::cout << "Model: " << decltype(model)::name() << "\n";

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial, date, std::vector<double>, double> forward_induction;
	forward_induction fwd_induction;

	fwd_induction(stockTree, model, deltas, 1.0);
	std::cout << "After forward induction step:\n";
	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());


	auto payoff = [](double S) {
		return S;
	};

	auto optionTree = stockTree;

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial, date, std::vector<double>> backward_induction;
	backward_induction bwd_induction;

	bwd_induction(optionTree, model,deltas, payoff);
	std::cout << "After backward induction step:\n";
	lattice_utility::print(optionTree, optionTree.begin(), optionTree.end());
	std::cout << "Price: " << optionTree.apex() << "\n\n";

	typedef lattice_algorithms::MergeLattices<date> merge;
	merge merge3;

	auto merged = merge3(lattice_types::LaunchPolicy::Sequential,stockTree, optionTree, stockTree);
	std::cout << "Merged stockTree and optionTree:\n";

	lattice_utility::print(merged, merged.begin(), merged.end());

}

void testMerge() {
	std::cout << "=====================================================\n";
	std::cout << "================= Merging Trees - TEST ==============\n";
	std::cout << "=====================================================\n";

	mergeBinomial();
	mergeTrinomial();


	std::cout << "=====================================================\n";
}

#endif ///_LATTICE_ALGORITHMS_T