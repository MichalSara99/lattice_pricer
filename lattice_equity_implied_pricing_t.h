#pragma once
#if !defined(_LATTICE_EQUITY_IMPLIED_PRICING_T)
#define _LATTICE_EQUITY_IMPLIED_PRICING_T


#include"lattice_structure.h"
#include"lattice_utility.h"
#include"lattice_calibrator.h"
#include"lattice_algorithms.h"



void testIndexedPriceFromImpliedCallOption1() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_algorithms::ImpliedBackwardInduction;

	std::cout << "Replication of call price:\n\n";

	// Firts collect Market data:
	// option strike price container:
	std::vector<double> strikes = { 46.90, 60.37,77.70,100.0,128.71,165.65,213.20 };
	// option maturity container
	std::vector<double> maturities = { 0.0,0.25,0.5,0.75,1.0 };
	// option price surface container of call,put prices
	// The option price surface is of following form
	// ----------------------------------------------------------
	// |    |            |     Maturity							|
	// |---------------------------------------------------------
	// |S	| {call,put} |										|
	// |t	|----------------------------------------------------
	// |r	| {call,put} |										|
	// |i	|----------------------------------------------------
	// |k	|													|
	// |e	| ....												|
	// |	|													|
	// ----------------------------------------------------------
	// Populating optionPrice container of pairs<call,put> is row-major:
	std::vector<std::pair<double, double>>
		optionPrice = { { 53.0966,0.0000 },{ 53.7949,0.0000 },{ 54.4828,0.0000 },{ 55.1604,0.0000 },{ 55.8281,0.0001 },
	{ 39.6325,0.0000 },{ 40.5313,0.0000 },{ 41.4170,0.0003 },{ 42.2922,0.0033 },{ 43.1598,0.0118 },
	{ 22.3035,0.0000 },{ 23.4719,0.0117 },{ 24.7110,0.1112 },{ 25.9912,0.2689 },{ 27.2682,0.4400 },
	{ 0.0000,0.0000 },{ 4.7469,3.2581 },{ 7.1559,4.2004 },{ 9.1754,4.7752 },{ 10.9895,5.1660 },
	{ 0.0000,28.7059 },{ 0.0403,26.8300 },{ 0.4195,25.3216 },{ 1.1159,24.1584 },{ 1.9965,23.2072 },
	{ 0.0000,65.6521 },{ 0.0000,63.1859 },{ 0.0050,60.7613 },{ 0.0412,58.4042 },{ 0.1399,56.1451 },
	{ 0.0000,113.2040 },{ 0.0000,110.0298 },{ 0.0001,106.9030 },{ 0.0011,103.8236 },{ 0.0055,100.7935 } };

	double impliedVol{ 0.29144 };
	double rate{ 0.06 };
	double dt{ 0.25 };
	double initialStock{ 100.0 };
	// Collect the market date into tuple in following order always
	// (strikes,maturities,optionPrices,implied volatilities) 
	auto optionData = std::make_tuple(strikes, maturities, optionPrice, impliedVol);
	std::size_t periods{ maturities.size() - 1 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> stockTree(periods);

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::Equity,
		std::size_t, double, double, decltype(optionData)> spl_calibrator;

	// Instantiate Calibrator for statePrice lattice and launch it:
	spl_calibrator calibrator(optionData, rate);
	auto result = calibrator(stockTree, dt, initialStock, true);
	// get the state prices lattice:
	auto statePriceTree = result->resultLattice;
	// Instantiate Calibrator for implied probabilities and launch it:
	auto impliedProbsResult = calibrator(statePriceTree, stockTree, dt, rate, true);
	
	// Prepare indexed lattice for option price:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> callTree(periods);

	// Declare implied backward alog
	typedef ImpliedBackwardInduction<LatticeType::Trinomial, double, double> implied_induction;
	implied_induction pricer(rate);

	// prepare call payoff:
	auto K = strikes[strikes.size() - 1];
	auto call_payoff = [&K](double stockPrice) {
		return std::max(stockPrice - K, 0.0);
	};

	pricer(callTree, stockTree, impliedProbsResult, dt, call_payoff);

	// print call option lattice:
	std::cout << "Call price lattice: \n";
	lattice_utility::print(callTree, callTree.cbegin(), callTree.cend());
	std::cout << "\n\nGiven call price: " << optionPrice[optionPrice.size() - 1].first << "\n";
	std::cout << "Replicated call price: " << callTree.apex() << "\n";

	std::cout << "==================================================\n";

}


void testIndexedPriceFromImpliedCallOption2() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_algorithms::ImpliedBackwardInduction;

	std::cout << "Replication of call price:\n\n";

	// Firts collect Market data:
	// option strike price container:
	std::vector<double> strikes = { 46.90, 60.37,77.70,100.0,128.71,165.65,213.20 };
	// option maturity container
	std::vector<double> maturities = { 0.0,0.25,0.5,0.75,1.0 };
	// option price surface container of call,put prices
	// The option price surface is of following form
	// ----------------------------------------------------------
	// |    |            |     Maturity							|
	// |---------------------------------------------------------
	// |S	| {call,put} |										|
	// |t	|----------------------------------------------------
	// |r	| {call,put} |										|
	// |i	|----------------------------------------------------
	// |k	|													|
	// |e	| ....												|
	// |	|													|
	// ----------------------------------------------------------
	// Populating optionPrice container of pairs<call,put> is row-major:
	std::vector<std::pair<double, double>>
		optionPrice = { { 53.0966,0.0000 },{ 53.7949,0.0000 },{ 54.4828,0.0000 },{ 55.1604,0.0000 },{ 55.8281,0.0001 },
	{ 39.6325,0.0000 },{ 40.5313,0.0000 },{ 41.4170,0.0003 },{ 42.2922,0.0033 },{ 43.1598,0.0118 },
	{ 22.3035,0.0000 },{ 23.4719,0.0117 },{ 24.7110,0.1112 },{ 25.9912,0.2689 },{ 27.2682,0.4400 },
	{ 0.0000,0.0000 },{ 4.7469,3.2581 },{ 7.1559,4.2004 },{ 9.1754,4.7752 },{ 10.9895,5.1660 },
	{ 0.0000,28.7059 },{ 0.0403,26.8300 },{ 0.4195,25.3216 },{ 1.1159,24.1584 },{ 1.9965,23.2072 },
	{ 0.0000,65.6521 },{ 0.0000,63.1859 },{ 0.0050,60.7613 },{ 0.0412,58.4042 },{ 0.1399,56.1451 },
	{ 0.0000,113.2040 },{ 0.0000,110.0298 },{ 0.0001,106.9030 },{ 0.0011,103.8236 },{ 0.0055,100.7935 } };

	double impliedVol{ 0.29144 };
	double rate{ 0.06 };
	double dt{ 0.25 };
	double initialStock{ 100.0 };
	// Collect the market date into tuple in following order always
	// (strikes,maturities,optionPrices,implied volatilities) 
	auto optionData = std::make_tuple(strikes, maturities, optionPrice, impliedVol);
	std::size_t periods{ maturities.size() - 1 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> stockTree(periods);

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::Equity,
		std::size_t, double, double, decltype(optionData)> spl_calibrator;

	// Instantiate Calibrator for statePrice lattice and launch it:
	spl_calibrator calibrator(optionData, rate);
	auto result = calibrator(stockTree, dt, initialStock, true);
	// get the state prices lattice:
	auto statePriceTree = result->resultLattice;
	// Instantiate Calibrator for implied probabilities and launch it:
	auto impliedProbsResult = calibrator(statePriceTree, stockTree, dt, rate, true);

	// Prepare indexed lattice for option price:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> callTree(2);

	// Declare implied backward alog
	typedef ImpliedBackwardInduction<LatticeType::Trinomial, double, double> implied_induction;
	implied_induction pricer(rate);

	// prepare call payoff:
	auto K = strikes[3];
	auto call_payoff = [&K](double stockPrice) {
		return std::max(stockPrice - K, 0.0);
	};

	pricer(callTree, stockTree, impliedProbsResult, dt, call_payoff);

	// print call option lattice:
	std::cout << "Call price lattice: \n";
	lattice_utility::print(callTree, callTree.cbegin(), callTree.cend());
	std::cout << "\n\nGiven call price: " << optionPrice[17].first << "\n";
	std::cout << "Replicated call price: " << callTree.apex() << "\n";

	std::cout << "==================================================\n";

}

void testIndexedPriceFromImpliedPutOption1() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_algorithms::ImpliedBackwardInduction;

	std::cout << "Replication of put price:\n\n";

	// Firts collect Market data:
	// option strike price container:
	std::vector<double> strikes = { 46.90, 60.37,77.70,100.0,128.71,165.65,213.20 };
	// option maturity container
	std::vector<double> maturities = { 0.0,0.25,0.5,0.75,1.0 };
	// option price surface container of call,put prices
	// The option price surface is of following form
	// ----------------------------------------------------------
	// |    |            |     Maturity							|
	// |---------------------------------------------------------
	// |S	| {call,put} |										|
	// |t	|----------------------------------------------------
	// |r	| {call,put} |										|
	// |i	|----------------------------------------------------
	// |k	|													|
	// |e	| ....												|
	// |	|													|
	// ----------------------------------------------------------
	// Populating optionPrice container of pairs<call,put> is row-major:
	std::vector<std::pair<double, double>>
		optionPrice = { { 53.0966,0.0000 },{ 53.7949,0.0000 },{ 54.4828,0.0000 },{ 55.1604,0.0000 },{ 55.8281,0.0001 },
	{ 39.6325,0.0000 },{ 40.5313,0.0000 },{ 41.4170,0.0003 },{ 42.2922,0.0033 },{ 43.1598,0.0118 },
	{ 22.3035,0.0000 },{ 23.4719,0.0117 },{ 24.7110,0.1112 },{ 25.9912,0.2689 },{ 27.2682,0.4400 },
	{ 0.0000,0.0000 },{ 4.7469,3.2581 },{ 7.1559,4.2004 },{ 9.1754,4.7752 },{ 10.9895,5.1660 },
	{ 0.0000,28.7059 },{ 0.0403,26.8300 },{ 0.4195,25.3216 },{ 1.1159,24.1584 },{ 1.9965,23.2072 },
	{ 0.0000,65.6521 },{ 0.0000,63.1859 },{ 0.0050,60.7613 },{ 0.0412,58.4042 },{ 0.1399,56.1451 },
	{ 0.0000,113.2040 },{ 0.0000,110.0298 },{ 0.0001,106.9030 },{ 0.0011,103.8236 },{ 0.0055,100.7935 } };

	double impliedVol{ 0.29144 };
	double rate{ 0.06 };
	double dt{ 0.25 };
	double initialStock{ 100.0 };
	// Collect the market date into tuple in following order always
	// (strikes,maturities,optionPrices,implied volatilities) 
	auto optionData = std::make_tuple(strikes, maturities, optionPrice, impliedVol);
	std::size_t periods{ maturities.size() - 1 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> stockTree(periods);

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::Equity,
		std::size_t, double, double, decltype(optionData)> spl_calibrator;

	// Instantiate Calibrator for statePrice lattice and launch it:
	spl_calibrator calibrator(optionData, rate);
	auto result = calibrator(stockTree, dt, initialStock, true);
	// get the state prices lattice:
	auto statePriceTree = result->resultLattice;
	// Instantiate Calibrator for implied probabilities and launch it:
	auto impliedProbsResult = calibrator(statePriceTree, stockTree, dt, rate, true);

	// Prepare indexed lattice for option price:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> putTree(periods);

	// Declare implied backward alog
	typedef ImpliedBackwardInduction<LatticeType::Trinomial, double, double> implied_induction;
	implied_induction pricer(rate);

	// prepare call payoff:
	auto K = strikes[strikes.size() - 1];
	auto put_payoff = [&K](double stockPrice) {
		return std::max(K - stockPrice, 0.0);
	};

	pricer(putTree, stockTree, impliedProbsResult, dt, put_payoff);

	// print call option lattice:
	std::cout << "Put price lattice: \n";
	lattice_utility::print(putTree, putTree.cbegin(), putTree.cend());
	std::cout << "\n\nGiven put price: " << optionPrice[optionPrice.size() - 1].second << "\n";
	std::cout << "Replicated put price: " << putTree.apex() << "\n";

	std::cout << "==================================================\n";

}

void testIndexedPriceFromImpliedPutOption2() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_algorithms::ImpliedBackwardInduction;

	std::cout << "Replication of put price:\n\n";

	// Firts collect Market data:
	// option strike price container:
	std::vector<double> strikes = { 46.90, 60.37,77.70,100.0,128.71,165.65,213.20 };
	// option maturity container
	std::vector<double> maturities = { 0.0,0.25,0.5,0.75,1.0 };
	// option price surface container of call,put prices
	// The option price surface is of following form
	// ----------------------------------------------------------
	// |    |            |     Maturity							|
	// |---------------------------------------------------------
	// |S	| {call,put} |										|
	// |t	|----------------------------------------------------
	// |r	| {call,put} |										|
	// |i	|----------------------------------------------------
	// |k	|													|
	// |e	| ....												|
	// |	|													|
	// ----------------------------------------------------------
	// Populating optionPrice container of pairs<call,put> is row-major:
	std::vector<std::pair<double, double>>
		optionPrice = { { 53.0966,0.0000 },{ 53.7949,0.0000 },{ 54.4828,0.0000 },{ 55.1604,0.0000 },{ 55.8281,0.0001 },
	{ 39.6325,0.0000 },{ 40.5313,0.0000 },{ 41.4170,0.0003 },{ 42.2922,0.0033 },{ 43.1598,0.0118 },
	{ 22.3035,0.0000 },{ 23.4719,0.0117 },{ 24.7110,0.1112 },{ 25.9912,0.2689 },{ 27.2682,0.4400 },
	{ 0.0000,0.0000 },{ 4.7469,3.2581 },{ 7.1559,4.2004 },{ 9.1754,4.7752 },{ 10.9895,5.1660 },
	{ 0.0000,28.7059 },{ 0.0403,26.8300 },{ 0.4195,25.3216 },{ 1.1159,24.1584 },{ 1.9965,23.2072 },
	{ 0.0000,65.6521 },{ 0.0000,63.1859 },{ 0.0050,60.7613 },{ 0.0412,58.4042 },{ 0.1399,56.1451 },
	{ 0.0000,113.2040 },{ 0.0000,110.0298 },{ 0.0001,106.9030 },{ 0.0011,103.8236 },{ 0.0055,100.7935 } };

	double impliedVol{ 0.29144 };
	double rate{ 0.06 };
	double dt{ 0.25 };
	double initialStock{ 100.0 };
	// Collect the market date into tuple in following order always
	// (strikes,maturities,optionPrices,implied volatilities) 
	auto optionData = std::make_tuple(strikes, maturities, optionPrice, impliedVol);
	std::size_t periods{ maturities.size() - 1 };

	// Creating indexed lattice:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> stockTree(periods);

	// Declare calibrator:
	typedef lattice_calibrator::Calibrator<LatticeType::Trinomial, AssetClass::Equity,
		std::size_t, double, double, decltype(optionData)> spl_calibrator;

	// Instantiate Calibrator for statePrice lattice and launch it:
	spl_calibrator calibrator(optionData, rate);
	auto result = calibrator(stockTree, dt, initialStock, true);
	// get the state prices lattice:
	auto statePriceTree = result->resultLattice;
	// Instantiate Calibrator for implied probabilities and launch it:
	auto impliedProbsResult = calibrator(statePriceTree, stockTree, dt, rate, true);

	// Prepare indexed lattice for option price:
	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> putTree(2);

	// Declare implied backward alog
	typedef ImpliedBackwardInduction<LatticeType::Trinomial, double, double> implied_induction;
	implied_induction pricer(rate);

	// prepare call payoff:
	auto K = strikes[3];
	auto put_payoff = [&K](double stockPrice) {
		return std::max(K - stockPrice, 0.0);
	};

	pricer(putTree, stockTree, impliedProbsResult, dt, put_payoff);

	// print call option lattice:
	std::cout << "Put price lattice: \n";
	lattice_utility::print(putTree, putTree.cbegin(), putTree.cend());
	std::cout << "\n\nGiven put price: " << optionPrice[17].second << "\n";
	std::cout << "Replicated put price: " << putTree.apex() << "\n";

	std::cout << "==================================================\n";

}


void testIndexedReplicateEuropeanOptionPrice() {
	std::cout << "==================================================================\n";
	std::cout << "====== Indexed Replication of European Option Prices =============\n";
	std::cout << "==================================================================\n";

	testIndexedPriceFromImpliedCallOption1();
	testIndexedPriceFromImpliedCallOption2();
	testIndexedPriceFromImpliedPutOption1();
	testIndexedPriceFromImpliedPutOption2();

	std::cout << "===================================================================\n";
}








#endif ///_LATTICE_EQUITY_IMPLIED_PRICING_T