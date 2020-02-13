#pragma once
#if !defined(_LATTICE_EXAMPLES)
#define _LATTICE_EXAMPLES

#include<iostream>
#include<string>
#include<set>
#include<numeric>
#include<thread>
#include<vector>
#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>

#include"lattice.h"

using namespace boost::gregorian;

void binomialLatticeParallelPricing() {

	std::cout << "\nPricing with std::thread: \n";
	using lattice_structure::Lattice;
	using lattice_types::LatticeType;
	using lattice_types::LeafForwardGenerator;
	using lattice_types::LeafBackwardGenerator;
	using lattice_types::Payoff;
	using lattice_types::PayoffAdjuster;
	using lattice_miscellaneous::OptionData;
	using lattice_model::CoxRubinsteinRossModel;
	using lattice_algorithms::forward_induction;
	using lattice_algorithms::backward_induction;
	using lattice_utility::print;


	// Create option data holder:
	OptionData<double> option;
	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 65.0;

	// Generate fixing dates:
	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t N = 365;
	for (auto t = 1; t < N; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Create a binomial lattice:
	Lattice<LatticeType::Binomial, double, date> biLattice{ fixingDates };

	// Create a model:
	CoxRubinsteinRossModel<> crr{ option };

	// Prepare delta times:
	double year{ 365.0 };
	auto fd = biLattice.fixingDates();
	std::vector<double> deltaTimes(fd.size() - 1);
	for (auto i = 0; i < deltaTimes.size(); ++i) {
		deltaTimes[i] = ((fd[i + 1] - fd[i]).days() / year);
	}

	// Build the lattice, i.e. populate it with asset prices from model:
	LeafForwardGenerator<double, double, double> fwdGen = crr;
	forward_induction(biLattice, fwdGen, option.Underlying, deltaTimes);


	// Prepare payoffs:
	double K = option.Strike;
	Payoff<double,double> callPayoff = [&K](double stock){
		return std::max(stock - K, 0.0);
	};
	Payoff<double, double> putPayoff = [&K](double stock) {
		return std::max(K - stock, 0.0);
	};

	// Prepare adjusters form american options:
	PayoffAdjuster<double&, double> callAdjuster = [&callPayoff](double &optionValue,double stock) {
		optionValue = std::max(optionValue, callPayoff(stock));
	};

	PayoffAdjuster<double&, double> putAdjuster = [&putPayoff](double &optionValue,double stock) {
		optionValue = std::max(optionValue, putPayoff(stock));
	};


	// create all the lattices:
	auto euroPutLattice = biLattice;
	auto euroCallLattice = biLattice;
	auto americanCallLattice = biLattice;
	auto americanPutLattice = biLattice;


	// spawn threads to do the forward induction:
	{
		std::thread t1([&]() {
			forward_induction(euroPutLattice, fwdGen, option.Underlying, deltaTimes);
		});

		std::thread t2([&]() {
			forward_induction(euroCallLattice, fwdGen, option.Underlying, deltaTimes);
		});
		std::thread t3([&]() {
			forward_induction(americanCallLattice, fwdGen, option.Underlying, deltaTimes);
		});
		std::thread t4([&]() {
			forward_induction(americanPutLattice, fwdGen, option.Underlying, deltaTimes);
		});

		// wait for the threads to finish:
		t1.join();
		t2.join();
		t3.join();
		t4.join();
	}
	// Build the option price lattices, i.e. populate it with option prices from model:
	LeafBackwardGenerator<double, double, double,double> backGen = crr;

	{
		std::thread t1([&]() {
			backward_induction(biLattice, euroPutLattice, backGen, putPayoff, deltaTimes);
		});

		std::thread t2([&]() {
			backward_induction(biLattice, euroCallLattice, backGen, callPayoff, deltaTimes);
		});
		std::thread t3([&]() {
			backward_induction(biLattice, americanCallLattice, backGen, callPayoff,callAdjuster,deltaTimes);
		});
		std::thread t4([&]() {
			backward_induction(biLattice, americanPutLattice, backGen, putPayoff, putAdjuster, deltaTimes);
		});

		// wait for the threads to finish:
		t1.join();
		t2.join();
		t3.join();
		t4.join();
	}

	// print the asset prices:
	std::cout << "Asset price tree: \n";
	auto first = biLattice.begin();
	auto last = std::next(first, 5);
	print(biLattice, first, last);
	std::cout << "\n";

	// print the option prices:
	std::cout << "Euro call option price tree: \n";
	first = euroCallLattice.begin();
	last = std::next(first, 5);
	print(euroCallLattice, first, last);

	std::cout << "Euro put option price tree: \n";
	first = euroPutLattice.begin();
	last = std::next(first, 5);
	print(euroPutLattice, first, last);

	std::cout << "American call option price tree: \n";
	first = americanCallLattice.begin();
	last = std::next(first, 5);
	print(americanCallLattice, first, last);

	std::cout << "American put option price tree: \n";
	first = americanPutLattice.begin();
	last = std::next(first, 5);
	print(americanPutLattice, first, last);

	std::cout << "European call price: " << euroCallLattice.apex() << "\n";
	std::cout << "European put price: " << euroPutLattice.apex() << "\n";
	std::cout << "American call price: " << americanCallLattice.apex() << "\n";
	std::cout << "American put price: " << americanPutLattice.apex() << "\n";

}


void binomialLatticeParallelPricingScoped() {

	std::cout << "\nPricing with scoped threads:\n";
	using lattice_structure::Lattice;
	using lattice_types::LatticeType;
	using lattice_types::LeafForwardGenerator;
	using lattice_types::LeafBackwardGenerator;
	using lattice_types::Payoff;
	using lattice_types::PayoffAdjuster;
	using lattice_miscellaneous::OptionData;
	using lattice_miscellaneous::scoped_thread;
	using lattice_model::CoxRubinsteinRossModel;
	using lattice_algorithms::forward_induction;
	using lattice_algorithms::backward_induction;
	using lattice_utility::print;


	// Create option data holder:
	OptionData<double> option;
	option.Strike = 65.0;
	option.RiskFreeRate = 0.25;
	option.DividentRate = 0.0;
	option.Volatility = 0.3;
	option.Underlying = 65.0;

	// Generate fixing dates:
	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t N = 365;
	for (auto t = 1; t < N; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Create a binomial lattice:
	Lattice<LatticeType::Binomial, double, date> biLattice{ fixingDates };

	// Create a model:
	CoxRubinsteinRossModel<> crr{ option };

	// Prepare delta times:
	double year{ 365.0 };
	auto fd = biLattice.fixingDates();
	std::vector<double> deltaTimes(fd.size() - 1);
	for (auto i = 0; i < deltaTimes.size(); ++i) {
		deltaTimes[i] = ((fd[i + 1] - fd[i]).days() / year);
	}

	// Build the lattice, i.e. populate it with asset prices from model:
	LeafForwardGenerator<double, double, double> fwdGen = crr;
	forward_induction(biLattice, fwdGen, option.Underlying, deltaTimes);


	// Prepare payoffs:
	double K = option.Strike;
	Payoff<double, double> callPayoff = [&K](double stock) {
		return std::max(stock - K, 0.0);
	};
	Payoff<double, double> putPayoff = [&K](double stock) {
		return std::max(K - stock, 0.0);
	};

	// Prepare adjusters form american options:
	PayoffAdjuster<double&, double> callAdjuster = [&callPayoff](double &optionValue, double stock) {
		optionValue = std::max(optionValue, callPayoff(stock));
	};

	PayoffAdjuster<double&, double> putAdjuster = [&putPayoff](double &optionValue, double stock) {
		optionValue = std::max(optionValue, putPayoff(stock));
	};


	// create all the lattices:
	auto euroPutLattice = biLattice;
	auto euroCallLattice = biLattice;
	auto americanCallLattice = biLattice;
	auto americanPutLattice = biLattice;


	// spawn threads to do the forward induction:
	{
		scoped_thread t1(std::thread([&]() {
			forward_induction(euroPutLattice, fwdGen, option.Underlying, deltaTimes);
		}));
		scoped_thread t2(std::thread([&]() {
			forward_induction(euroCallLattice, fwdGen, option.Underlying, deltaTimes);
		}));
		scoped_thread t3(std::thread([&]() {
			forward_induction(americanCallLattice, fwdGen, option.Underlying, deltaTimes);
		}));
		scoped_thread t4(std::thread([&]() {
			forward_induction(americanPutLattice, fwdGen, option.Underlying, deltaTimes);
		}));
	}
	// Build the option price lattices, i.e. populate it with option prices from model:
	LeafBackwardGenerator<double, double, double, double> backGen = crr;

	{
		scoped_thread t1(std::thread([&]() {
			backward_induction(biLattice, euroPutLattice, backGen, putPayoff, deltaTimes);
		}));
		scoped_thread t2(std::thread([&]() {
			backward_induction(biLattice, euroCallLattice, backGen, callPayoff, deltaTimes);
		}));
		scoped_thread t3(std::thread([&]() {
			backward_induction(biLattice, americanCallLattice, backGen, callPayoff, callAdjuster, deltaTimes);
		}));
		scoped_thread t4(std::thread([&]() {
			backward_induction(biLattice, americanPutLattice, backGen, putPayoff, putAdjuster, deltaTimes);
		}));
	}

	// print the asset prices:
	std::cout << "Asset price tree: \n";
	auto first = biLattice.begin();
	auto last = std::next(first, 5);
	print(biLattice, first, last);
	std::cout << "\n";

	// print the option prices:
	std::cout << "Euro call option price tree: \n";
	first = euroCallLattice.begin();
	last = std::next(first, 5);
	print(euroCallLattice, first, last);

	std::cout << "Euro put option price tree: \n";
	first = euroPutLattice.begin();
	last = std::next(first, 5);
	print(euroPutLattice, first, last);

	std::cout << "American call option price tree: \n";
	first = americanCallLattice.begin();
	last = std::next(first, 5);
	print(americanCallLattice, first, last);

	std::cout << "American put option price tree: \n";
	first = americanPutLattice.begin();
	last = std::next(first, 5);
	print(americanPutLattice, first, last);

	std::cout << "European call price: " << euroCallLattice.apex() << "\n";
	std::cout << "European put price: " << euroPutLattice.apex() << "\n";
	std::cout << "American call price: " << americanCallLattice.apex() << "\n";
	std::cout << "American put price: " << americanPutLattice.apex() << "\n";

}





#endif //_LATTICE_EXAMPLES