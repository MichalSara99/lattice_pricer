#pragma once
#if !defined(_LATTICE_CALIBRATOR_EQUITY_T)
#define _LATTICE_CALIBRATOR_EQUITY_T


#include"lattice_structure.h"
#include"lattice_utility.h"
#include"lattice_calibrator.h"


#include<iostream>
#include<vector>
#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>

using namespace boost::gregorian;

// SPL = StatePriceLattice
void testIndexedSPLImpliedCallsLiquid() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

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
		optionPrice = { {53.0966,0.0000},{53.7949,0.0000},{54.4828,0.0000},{55.1604,0.0000},{55.8281,0.0001},
						{39.6325,0.0000},{40.5313,0.0000},{41.4170,0.0003},{42.2922,0.0033},{43.1598,0.0118},
						{22.3035,0.0000},{23.4719,0.0117},{24.7110,0.1112},{25.9912,0.2689},{27.2682,0.4400},
						{0.0000,0.0000},{4.7469,3.2581},{7.1559,4.2004},{9.1754,4.7752},{10.9895,5.1660},
						{0.0000,28.7059},{0.0403,26.8300},{0.4195,25.3216},{1.1159,24.1584},{1.9965,23.2072},
						{0.0000,65.6521},{0.0000,63.1859},{0.0050,60.7613},{0.0412,58.4042},{0.1399,56.1451},
						{0.0000,113.2040},{0.0000,110.0298},{0.0001,106.9030},{0.0011,103.8236},{0.0055,100.7935} };

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
		std::size_t, double,double, decltype(optionData)> spl_calibrator;

	// Instantiate Calibrator and launch it:
	spl_calibrator calibrator(optionData, rate);
	auto result = calibrator(stockTree, dt, initialStock,true);

	// Print the part of generated stocklattice:
	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());

	// Print State price lattice:
	std::cout << "\n";
	std::cout << "Printing state price lattice: \n";
	auto statePriceTree = result->resultLattice;
	lattice_utility::print(statePriceTree, statePriceTree.begin(), statePriceTree.end());
}

void testIndexedSPLImpliedPutsLiquid() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

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

	// Instantiate Calibrator and launch it:
	spl_calibrator calibrator(optionData, rate);
	auto result = calibrator(stockTree, dt, initialStock, false);

	// Print the part of generated stocklattice:
	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());

	// Print State price lattice:
	std::cout << "\n";
	std::cout << "Printing state price lattice: \n";
	auto statePriceTree = result->resultLattice;
	lattice_utility::print(statePriceTree, statePriceTree.begin(), statePriceTree.end());
}


void testIndexedImpliedProbabilitiesCallsLiquid() {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;

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

	// Print the part of generated stocklattice:
	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());

	// Print State price lattice:
	std::cout << "\n";
	std::cout << "Printing state price lattice: \n";
	auto statePriceTree = result->resultLattice;
	lattice_utility::print(statePriceTree, statePriceTree.begin(), statePriceTree.end());

	// Instantiate Calibrator for implied probabilities and launch it:
	auto impliedProbsResult = calibrator(statePriceTree, stockTree, dt, rate);
	auto impliedProbs = impliedProbsResult->impliedProbabilities;

	for (std::size_t t = 0; t < impliedProbs.size(); ++t) {
		std::cout << "time: " << t << "\n";
		auto branch = impliedProbs.at(t);
		for (std::size_t l = 0; l < branch.size(); ++l) {
			std::cout << "(" << std::get<0>(branch.at(l)) << ", "
				<< std::get<1>(branch.at(l)) << ", "
				<< std::get<2>(branch.at(l)) << "), ";

		}
		std::cout << "\n";
	}

}


#endif ///_LATTICE_CALIBRATOR_EQUITY_T