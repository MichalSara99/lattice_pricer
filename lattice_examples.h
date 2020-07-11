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



void crrBinomialLatticeParallelPricing() {

	std::cout << "\nPricing with std::thread: \n";
	using lattice_structure::Lattice;
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_model::CoxRubinsteinRossModel;
	using lattice_algorithms::ForwardInduction;
	using lattice_algorithms::BackwardInduction;
	using lattice_utility::print;


	// Create option data holder:
	ModelParams<1,AssetClass::Equity,double> params;
	params.Strike = 65.0;
	params.RiskFreeRate = 0.25;
	params.DividendRate = 0.0;
	params.Volatility = 0.3;
	params.Spot = 65.0;

	// Generate fixing dates:
	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods = 364;
	for (auto t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Create a binomial lattice:
	Lattice<LatticeType::Binomial, double, date> biLattice{ fixingDates };

	// Create a model:
	CoxRubinsteinRossModel<> crr{ params };

	// Prepare delta times:
	double year{ 365.0 };
	auto fd = biLattice.fixingDates();
	std::vector<double> deltaTimes(fd.size() - 1);
	for (auto i = 0; i < deltaTimes.size(); ++i) {
		deltaTimes[i] = ((fd[i + 1] - fd[i]).days() / year);
	}

	// Build the lattice, i.e. populate it with asset prices from model:
	typedef ForwardInduction<LatticeType::Binomial, date, std::vector<double>, double> forward_induction;
	forward_induction fwd_induction;


	// Prepare payoffs:
	double K = params.Strike;
	auto callPayoff = [&K](double stock) {
		return std::max(stock - K, 0.0);
	};
	auto putPayoff = [&K](double stock) {
		return std::max(K - stock, 0.0);
	};

	// Prepare adjusters form american options:
	auto callAdjuster = [&callPayoff](double &optionValue, double stock) {
		optionValue = std::max(optionValue, callPayoff(stock));
	};

	auto putAdjuster = [&putPayoff](double &optionValue, double stock) {
		optionValue = std::max(optionValue, putPayoff(stock));
	};

	// spawn threads to do the forward induction:
	{
		std::thread t0([&]() {
			fwd_induction(biLattice, crr, deltaTimes, params.Spot);
		});
		// wait for the threads to finish:
		t0.join();
	}

	// create all the lattices:
	auto euroPutLattice = biLattice;
	auto euroCallLattice = biLattice;
	auto americanCallLattice = biLattice;
	auto americanPutLattice = biLattice;

	// Build the option price lattices, i.e. populate it with option prices from model:
	typedef BackwardInduction<LatticeType::Binomial, date, std::vector<double>> backward_induction;
	backward_induction bwd_induction;


	{
		std::thread t1([&]() {
			bwd_induction(euroPutLattice, crr,deltaTimes, putPayoff);
		});

		std::thread t2([&]() {
			bwd_induction(euroCallLattice, crr, deltaTimes, callPayoff);
		});
		std::thread t3([&]() {
			bwd_induction(americanCallLattice, crr, deltaTimes, callPayoff, callAdjuster);
		});
		std::thread t4([&]() {
			bwd_induction(americanPutLattice, crr, deltaTimes, putPayoff, putAdjuster);
		});

		// wait for the threads to finish:
		t1.join();
		t2.join();
		t3.join();
		t4.join();
	}

	//// print the asset prices:
	//std::cout << "Asset price tree: \n";
	//auto first = biLattice.begin();
	//auto last = std::next(first, 5);
	//print(biLattice, first, last);
	//std::cout << "\n";

	//// print the option prices:
	//std::cout << "Euro call option price tree: \n";
	//first = euroCallLattice.begin();
	//last = std::next(first, 5);
	//print(euroCallLattice, first, last);

	//std::cout << "Euro put option price tree: \n";
	//first = euroPutLattice.begin();
	//last = std::next(first, 5);
	//print(euroPutLattice, first, last);

	//std::cout << "American call option price tree: \n";
	//first = americanCallLattice.begin();
	//last = std::next(first, 5);
	//print(americanCallLattice, first, last);

	//std::cout << "American put option price tree: \n";
	//first = americanPutLattice.begin();
	//last = std::next(first, 5);
	//print(americanPutLattice, first, last);

	std::cout << "European call price: " << euroCallLattice.apex() << "\n";
	std::cout << "European put price: " << euroPutLattice.apex() << "\n";
	std::cout << "American call price: " << americanCallLattice.apex() << "\n";
	std::cout << "American put price: " << americanPutLattice.apex() << "\n";

}


void crrBinomialLatticeParallelPricingScoped() {

	std::cout << "\nPricing with scoped threads:\n";
	using lattice_structure::Lattice;
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_miscellaneous::scoped_thread;
	using lattice_model::CoxRubinsteinRossModel;
	using lattice_algorithms::ForwardInduction;
	using lattice_algorithms::BackwardInduction;
	using lattice_utility::print;


	// Create option data holder:
	ModelParams<1, AssetClass::Equity, double> params;
	params.Strike = 65.0;
	params.RiskFreeRate = 0.25;
	params.DividendRate = 0.0;
	params.Volatility = 0.3;
	params.Spot = 65.0;

	// Generate fixing dates:
	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods = 364;
	for (auto t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Create a binomial lattice:
	Lattice<LatticeType::Binomial, double, date> biLattice{ fixingDates };

	// Create a model:
	CoxRubinsteinRossModel<> crr{ params };

	// name of the Model:
	std::cout << decltype(crr)::name() << "\n";

	// Prepare delta times:
	double year{ 365.0 };
	auto fd = biLattice.fixingDates();
	std::vector<double> deltaTimes(fd.size() - 1);
	for (auto i = 0; i < deltaTimes.size(); ++i) {
		deltaTimes[i] = ((fd[i + 1] - fd[i]).days() / year);
	}

	// Build the lattice, i.e. populate it with asset prices from model:
	typedef ForwardInduction<LatticeType::Binomial, date, std::vector<double>, double> forward_induction;
	forward_induction fwd_induction;


	// Prepare payoffs:
	double K = params.Strike;
	auto callPayoff = [&K](double stock) {
		return std::max(stock - K, 0.0);
	};
	auto putPayoff = [&K](double stock) {
		return std::max(K - stock, 0.0);
	};

	// Prepare adjusters form american options:
	auto callAdjuster = [&callPayoff](double &optionValue, double stock) {
		optionValue = std::max(optionValue, callPayoff(stock));
	};

	auto putAdjuster = [&putPayoff](double &optionValue, double stock) {
		optionValue = std::max(optionValue, putPayoff(stock));
	};

	// spawn threads to do the forward induction:
	{
		scoped_thread t0(std::thread([&]() {
			fwd_induction(biLattice, crr, deltaTimes , params.Spot);
		}));
	}

	// create all the lattices:
	auto euroPutLattice = biLattice;
	auto euroCallLattice = biLattice;
	auto americanCallLattice = biLattice;
	auto americanPutLattice = biLattice;

	// Build the option price lattices, i.e. populate it with option prices from model:
	typedef BackwardInduction<LatticeType::Binomial, date, std::vector<double>> backward_induction;
	backward_induction bwd_induction;

	{
		scoped_thread t1(std::thread([&]() {
			bwd_induction( euroPutLattice, crr, deltaTimes,putPayoff);
		}));
		scoped_thread t2(std::thread([&]() {
			bwd_induction( euroCallLattice, crr, deltaTimes, callPayoff);
		}));
		scoped_thread t3(std::thread([&]() {
			bwd_induction( americanCallLattice, crr, deltaTimes, callPayoff, callAdjuster);
		}));
		scoped_thread t4(std::thread([&]() {
			bwd_induction( americanPutLattice, crr, deltaTimes, putPayoff, putAdjuster);
		}));
	}

	//// print the asset prices:
	//std::cout << "Asset price tree: \n";
	//auto first = biLattice.begin();
	//auto last = std::next(first, 5);
	//print(biLattice, first, last);
	//std::cout << "\n";

	//// print the option prices:
	//std::cout << "Euro call option price tree: \n";
	//first = euroCallLattice.begin();
	//last = std::next(first, 5);
	//print(euroCallLattice, first, last);

	//std::cout << "Euro put option price tree: \n";
	//first = euroPutLattice.begin();
	//last = std::next(first, 5);
	//print(euroPutLattice, first, last);

	//std::cout << "American call option price tree: \n";
	//first = americanCallLattice.begin();
	//last = std::next(first, 5);
	//print(americanCallLattice, first, last);

	//std::cout << "American put option price tree: \n";
	//first = americanPutLattice.begin();
	//last = std::next(first, 5);
	//print(americanPutLattice, first, last);

	std::cout << "European call price: " << euroCallLattice.apex() << "\n";
	std::cout << "European put price: " << euroPutLattice.apex() << "\n";
	std::cout << "American call price: " << americanCallLattice.apex() << "\n";
	std::cout << "American put price: " << americanPutLattice.apex() << "\n";
	std::cout << "time dimension: " << biLattice.timeDimension() << "\n";
}

void mcrrBinomialLatticeParallelPricingScoped() {

	std::cout << "\nPricing with scoped threads:\n";
	using lattice_structure::Lattice;
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_miscellaneous::scoped_thread;
	using lattice_model::ModifiedCoxRubinsteinRossModel;
	using lattice_algorithms::ForwardInduction;
	using lattice_algorithms::BackwardInduction;
	using lattice_utility::print;


	// Create option data holder:
	ModelParams<1, AssetClass::Equity, double> params;
	params.Strike = 65.0;
	params.RiskFreeRate = 0.25;
	params.DividendRate = 0.0;
	params.Volatility = 0.3;
	params.Spot = 65.0;

	// Generate fixing dates:
	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods = 364;
	for (auto t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Create a binomial lattice:
	Lattice<LatticeType::Binomial, double, date> biLattice{ fixingDates };

	// Create a model:
	ModifiedCoxRubinsteinRossModel<> mcrr{ params,periods };

	// Name of the model:
	std::cout << decltype(mcrr)::name() << "\n";

	// Prepare delta times:
	double year{ 365.0 };
	auto fd = biLattice.fixingDates();
	std::vector<double> deltaTimes(fd.size() - 1);
	for (auto i = 0; i < deltaTimes.size(); ++i) {
		deltaTimes[i] = ((fd[i + 1] - fd[i]).days() / year);
	}

	// Build the lattice, i.e. populate it with asset prices from model:
	typedef ForwardInduction<LatticeType::Binomial, date, std::vector<double>,double> forward_induction;
	forward_induction fwd_induction;


	// Prepare payoffs:
	double K = params.Strike;
	auto callPayoff = [&K](double stock) {
		return std::max(stock - K, 0.0);
	};
	auto putPayoff = [&K](double stock) {
		return std::max(K - stock, 0.0);
	};

	// Prepare adjusters form american options:
	auto callAdjuster = [&callPayoff](double& optionValue, double stock) {
		optionValue = std::max(optionValue, callPayoff(stock));
	};

	auto putAdjuster = [&putPayoff](double& optionValue, double stock) {
		optionValue = std::max(optionValue, putPayoff(stock));
	};

	// spawn threads to do the forward induction:
	{
		scoped_thread t1(std::thread([&]() {
			fwd_induction(biLattice, mcrr, deltaTimes, params.Spot);
		}));
	}

	// create all the lattices:
	auto euroPutLattice = biLattice;
	auto euroCallLattice = biLattice;
	auto americanCallLattice = biLattice;
	auto americanPutLattice = biLattice;
	// Build the option price lattices, i.e. populate it with option prices from model:
	typedef BackwardInduction<LatticeType::Binomial, date, std::vector<double>> backward_induction;
	backward_induction bwd_induction;

	{
		scoped_thread t1(std::thread([&]() {
			bwd_induction(euroPutLattice, mcrr, deltaTimes, putPayoff);
		}));
		scoped_thread t2(std::thread([&]() {
			bwd_induction(euroCallLattice, mcrr, deltaTimes, callPayoff);
		}));
		scoped_thread t3(std::thread([&]() {
			bwd_induction(americanCallLattice, mcrr, deltaTimes, callPayoff, callAdjuster);
		}));
		scoped_thread t4(std::thread([&]() {
			bwd_induction(americanPutLattice, mcrr,deltaTimes, putPayoff, putAdjuster);
		}));
	}

	//// print the asset prices:
	//std::cout << "Asset price tree: \n";
	//auto first = biLattice.begin();
	//auto last = std::next(first, 5);
	//print(biLattice, first, last);
	//std::cout << "\n";

	//// print the option prices:
	//std::cout << "Euro call option price tree: \n";
	//first = euroCallLattice.begin();
	//last = std::next(first, 5);
	//print(euroCallLattice, first, last);

	//std::cout << "Euro put option price tree: \n";
	//first = euroPutLattice.begin();
	//last = std::next(first, 5);
	//print(euroPutLattice, first, last);

	//std::cout << "American call option price tree: \n";
	//first = americanCallLattice.begin();
	//last = std::next(first, 5);
	//print(americanCallLattice, first, last);

	//std::cout << "American put option price tree: \n";
	//first = americanPutLattice.begin();
	//last = std::next(first, 5);
	//print(americanPutLattice, first, last);

	std::cout << "European call price: " << euroCallLattice.apex() << "\n";
	std::cout << "European put price: " << euroPutLattice.apex() << "\n";
	std::cout << "American call price: " << americanCallLattice.apex() << "\n";
	std::cout << "American put price: " << americanPutLattice.apex() << "\n";
	std::cout << "time dimension: " << biLattice.timeDimension() << "\n";
}


void jrBinomialLatticeParallelPricingScoped() {

	std::cout << "\nPricing with scoped threads:\n";
	using lattice_structure::Lattice;
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_miscellaneous::scoped_thread;
	using lattice_model::JarrowRuddModel;
	using lattice_algorithms::ForwardInduction;
	using lattice_algorithms::BackwardInduction;
	using lattice_utility::print;


	// Create option data holder:
	ModelParams<1, AssetClass::Equity, double> params;
	params.Strike = 65.0;
	params.RiskFreeRate = 0.25;
	params.DividendRate = 0.0;
	params.Volatility = 0.3;
	params.Spot = 65.0;

	// Generate fixing dates:
	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods = 364;
	for (auto t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Create a binomial lattice:
	Lattice<LatticeType::Binomial, double, date> biLattice{ fixingDates };

	// Create a model:
	JarrowRuddModel<> jr{ params };

	// Name of the model:
	std::cout << decltype(jr)::name() << "\n";

	// Prepare delta times:
	double year{ 365.0 };
	auto fd = biLattice.fixingDates();
	std::vector<double> deltaTimes(fd.size() - 1);
	for (auto i = 0; i < deltaTimes.size(); ++i) {
		deltaTimes[i] = ((fd[i + 1] - fd[i]).days() / year);
	}

	// Build the lattice, i.e. populate it with asset prices from model:
	typedef ForwardInduction<LatticeType::Binomial, date, std::vector<double>, double> forward_induction;
	forward_induction fwd_induction;


	// Prepare payoffs:
	double K = params.Strike;
	auto callPayoff = [&K](double stock) {
		return std::max(stock - K, 0.0);
	};
	auto putPayoff = [&K](double stock) {
		return std::max(K - stock, 0.0);
	};

	// Prepare adjusters form american options:
	auto callAdjuster = [&callPayoff](double& optionValue, double stock) {
		optionValue = std::max(optionValue, callPayoff(stock));
	};

	auto putAdjuster = [&putPayoff](double& optionValue, double stock) {
		optionValue = std::max(optionValue, putPayoff(stock));
	};


	// spawn threads to do the forward induction:
	{
		scoped_thread t1(std::thread([&]() {
			fwd_induction(biLattice, jr, deltaTimes, params.Spot);
		}));
	}

	// create all the lattices:
	auto euroPutLattice = biLattice;
	auto euroCallLattice = biLattice;
	auto americanCallLattice = biLattice;
	auto americanPutLattice = biLattice;

	// Build the option price lattices, i.e. populate it with option prices from model:
	typedef BackwardInduction<LatticeType::Binomial, date, std::vector<double>> backward_induction;
	backward_induction bwd_induction;

	{
		scoped_thread t1(std::thread([&]() {
			bwd_induction(euroPutLattice, jr, deltaTimes,putPayoff);
		}));
		scoped_thread t2(std::thread([&]() {
			bwd_induction( euroCallLattice, jr, deltaTimes,callPayoff);
		}));
		scoped_thread t3(std::thread([&]() {
			bwd_induction( americanCallLattice, jr, deltaTimes, callPayoff, callAdjuster);
		}));
		scoped_thread t4(std::thread([&]() {
			bwd_induction( americanPutLattice, jr, deltaTimes, putPayoff, putAdjuster);
		}));
	}

	//// print the asset prices:
	//std::cout << "Asset price tree: \n";
	//auto first = biLattice.begin();
	//auto last = std::next(first, 5);
	//print(biLattice, first, last);
	//std::cout << "\n";

	//// print the option prices:
	//std::cout << "Euro call option price tree: \n";
	//first = euroCallLattice.begin();
	//last = std::next(first, 5);
	//print(euroCallLattice, first, last);

	//std::cout << "Euro put option price tree: \n";
	//first = euroPutLattice.begin();
	//last = std::next(first, 5);
	//print(euroPutLattice, first, last);

	//std::cout << "American call option price tree: \n";
	//first = americanCallLattice.begin();
	//last = std::next(first, 5);
	//print(americanCallLattice, first, last);

	//std::cout << "American put option price tree: \n";
	//first = americanPutLattice.begin();
	//last = std::next(first, 5);
	//print(americanPutLattice, first, last);

	std::cout << "European call price: " << euroCallLattice.apex() << "\n";
	std::cout << "European put price: " << euroPutLattice.apex() << "\n";
	std::cout << "American call price: " << americanCallLattice.apex() << "\n";
	std::cout << "American put price: " << americanPutLattice.apex() << "\n";
	std::cout << "time dimension: " << biLattice.timeDimension() << "\n";
}


void trimBinomialLatticeParallelPricingScoped() {

	std::cout << "\nPricing with scoped threads:\n";
	using lattice_structure::Lattice;
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_miscellaneous::scoped_thread;
	using lattice_model::TrigeorgisModel;
	using lattice_algorithms::ForwardInduction;
	using lattice_algorithms::BackwardInduction;
	using lattice_utility::print;


	// Create option data holder:
	ModelParams<1, AssetClass::Equity, double> params;
	params.Strike = 65.0;
	params.RiskFreeRate = 0.25;
	params.DividendRate = 0.0;
	params.Volatility = 0.3;
	params.Spot = 65.0;

	// Generate fixing dates:
	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods = 364;
	for (auto t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Create a binomial lattice:
	Lattice<LatticeType::Binomial, double, date> biLattice{ fixingDates };

	// Create a model:
	TrigeorgisModel<> trim{ params };

	// Name of the model:
	std::cout << decltype(trim)::name() << "\n";

	// Prepare delta times:
	double year{ 365.0 };
	auto fd = biLattice.fixingDates();
	std::vector<double> deltaTimes(fd.size() - 1);
	for (auto i = 0; i < deltaTimes.size(); ++i) {
		deltaTimes[i] = ((fd[i + 1] - fd[i]).days() / year);
	}

	// Build the lattice, i.e. populate it with asset prices from model:
	typedef ForwardInduction<LatticeType::Binomial, date, std::vector<double>,double> forward_induction;
	forward_induction fwd_induction;


	// Prepare payoffs:
	double K = params.Strike;
	auto callPayoff = [&K](double stock) {
		return std::max(stock - K, 0.0);
	};
	auto putPayoff = [&K](double stock) {
		return std::max(K - stock, 0.0);
	};

	// Prepare adjusters form american options:
	auto callAdjuster = [&callPayoff](double& optionValue, double stock) {
		optionValue = std::max(optionValue, callPayoff(stock));
	};

	auto putAdjuster = [&putPayoff](double& optionValue, double stock) {
		optionValue = std::max(optionValue, putPayoff(stock));
	};


	// spawn threads to do the forward induction:
	{
		scoped_thread t1(std::thread([&]() {
			fwd_induction(biLattice, trim, deltaTimes, params.Spot);
		}));
	}

	// create all the lattices:
	auto euroPutLattice = biLattice;
	auto euroCallLattice = biLattice;
	auto americanCallLattice = biLattice;
	auto americanPutLattice = biLattice;

	// Build the option price lattices, i.e. populate it with option prices from model:
	typedef BackwardInduction<LatticeType::Binomial, date, std::vector<double>> backward_induction;
	backward_induction bwd_induction;

	{
		scoped_thread t1(std::thread([&]() {
			bwd_induction( euroPutLattice, trim, deltaTimes, putPayoff);
		}));
		scoped_thread t2(std::thread([&]() {
			bwd_induction(euroCallLattice, trim, deltaTimes,callPayoff);
		}));
		scoped_thread t3(std::thread([&]() {
			bwd_induction(americanCallLattice, trim, deltaTimes, callPayoff, callAdjuster);
		}));
		scoped_thread t4(std::thread([&]() {
			bwd_induction(americanPutLattice, trim, deltaTimes, putPayoff, putAdjuster);
		}));
	}

	//// print the asset prices:
	//std::cout << "Asset price tree: \n";
	//auto first = biLattice.begin();
	//auto last = std::next(first, 5);
	//print(biLattice, first, last);
	//std::cout << "\n";

	//// print the option prices:
	//std::cout << "Euro call option price tree: \n";
	//first = euroCallLattice.begin();
	//last = std::next(first, 5);
	//print(euroCallLattice, first, last);

	//std::cout << "Euro put option price tree: \n";
	//first = euroPutLattice.begin();
	//last = std::next(first, 5);
	//print(euroPutLattice, first, last);

	//std::cout << "American call option price tree: \n";
	//first = americanCallLattice.begin();
	//last = std::next(first, 5);
	//print(americanCallLattice, first, last);

	//std::cout << "American put option price tree: \n";
	//first = americanPutLattice.begin();
	//last = std::next(first, 5);
	//print(americanPutLattice, first, last);

	std::cout << "European call price: " << euroCallLattice.apex() << "\n";
	std::cout << "European put price: " << euroPutLattice.apex() << "\n";
	std::cout << "American call price: " << americanCallLattice.apex() << "\n";
	std::cout << "American put price: " << americanPutLattice.apex() << "\n";
	std::cout << "time dimension: " << biLattice.timeDimension() << "\n";
}



void tmBinomialLatticeParallelPricingScoped() {

	std::cout << "\nPricing with scoped threads:\n";
	using lattice_structure::Lattice;
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_miscellaneous::scoped_thread;
	using lattice_model::TianModel;
	using lattice_algorithms::ForwardInduction;
	using lattice_algorithms::BackwardInduction;
	using lattice_utility::print;


	// Create option data holder:
	ModelParams<1, AssetClass::Equity, double> params;
	params.Strike = 65.0;
	params.RiskFreeRate = 0.25;
	params.DividendRate = 0.0;
	params.Volatility = 0.3;
	params.Spot = 65.0;

	// Generate fixing dates:
	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods = 364;
	for (auto t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Create a binomial lattice:
	Lattice<LatticeType::Binomial, double, date> biLattice{ fixingDates };

	// Create a model:
	TianModel<> tm{ params };

	// Name of the model:
	std::cout << decltype(tm)::name() << "\n";

	// Prepare delta times:
	double year{ 365.0 };
	auto fd = biLattice.fixingDates();
	std::vector<double> deltaTimes(fd.size() - 1);
	for (auto i = 0; i < deltaTimes.size(); ++i) {
		deltaTimes[i] = ((fd[i + 1] - fd[i]).days() / year);
	}

	// Build the lattice, i.e. populate it with asset prices from model:
	typedef ForwardInduction<LatticeType::Binomial, date, std::vector<double>, double> forward_induction;
	forward_induction fwd_induction;


	// Prepare payoffs:
	double K = params.Strike;
	auto callPayoff = [&K](double stock) {
		return std::max(stock - K, 0.0);
	};
	auto putPayoff = [&K](double stock) {
		return std::max(K - stock, 0.0);
	};

	// Prepare adjusters form american options:
	auto callAdjuster = [&callPayoff](double& optionValue, double stock) {
		optionValue = std::max(optionValue, callPayoff(stock));
	};

	auto putAdjuster = [&putPayoff](double& optionValue, double stock) {
		optionValue = std::max(optionValue, putPayoff(stock));
	};


	// spawn threads to do the forward induction:
	{
		scoped_thread t1(std::thread([&]() {
			fwd_induction(biLattice, tm, deltaTimes, params.Spot);
		}));
	}

	// create all the lattices:
	auto euroPutLattice = biLattice;
	auto euroCallLattice = biLattice;
	auto americanCallLattice = biLattice;
	auto americanPutLattice = biLattice;


	// Build the option price lattices, i.e. populate it with option prices from model:
	typedef BackwardInduction<LatticeType::Binomial, date, std::vector<double>> backward_induction;
	backward_induction bwd_induction;

	{
		scoped_thread t1(std::thread([&]() {
			bwd_induction(euroPutLattice, tm, deltaTimes, putPayoff);
		}));
		scoped_thread t2(std::thread([&]() {
			bwd_induction(euroCallLattice, tm, deltaTimes, callPayoff);
		}));
		scoped_thread t3(std::thread([&]() {
			bwd_induction( americanCallLattice, tm, deltaTimes, callPayoff, callAdjuster);
		}));
		scoped_thread t4(std::thread([&]() {
			bwd_induction( americanPutLattice, tm, deltaTimes, putPayoff, putAdjuster);
		}));
	}

	//// print the asset prices:
	//std::cout << "Asset price tree: \n";
	//auto first = biLattice.begin();
	//auto last = std::next(first, 5);
	//print(biLattice, first, last);
	//std::cout << "\n";

	//// print the option prices:
	//std::cout << "Euro call option price tree: \n";
	//first = euroCallLattice.begin();
	//last = std::next(first, 5);
	//print(euroCallLattice, first, last);

	//std::cout << "Euro put option price tree: \n";
	//first = euroPutLattice.begin();
	//last = std::next(first, 5);
	//print(euroPutLattice, first, last);

	//std::cout << "American call option price tree: \n";
	//first = americanCallLattice.begin();
	//last = std::next(first, 5);
	//print(americanCallLattice, first, last);

	//std::cout << "American put option price tree: \n";
	//first = americanPutLattice.begin();
	//last = std::next(first, 5);
	//print(americanPutLattice, first, last);

	std::cout << "European call price: " << euroCallLattice.apex() << "\n";
	std::cout << "European put price: " << euroPutLattice.apex() << "\n";
	std::cout << "American call price: " << americanCallLattice.apex() << "\n";
	std::cout << "American put price: " << americanPutLattice.apex() << "\n";
	std::cout << "time dimension: " << biLattice.timeDimension() << "\n";
}


void lrBinomialLatticeParallelPricingScoped() {

	std::cout << "\nPricing with scoped threads:\n";
	using lattice_structure::Lattice;
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_model_params::ModelParams;
	using lattice_miscellaneous::scoped_thread;
	using lattice_model::LeisenReimerModel;
	using lattice_model_components::leisen_reimer_inversion::PeizerPrattSecondInversion;
	using lattice_algorithms::ForwardInduction;
	using lattice_algorithms::BackwardInduction;
	using lattice_utility::print;


	// Create option data holder:
	ModelParams<1, AssetClass::Equity, double> params;
	params.Strike = 65.0;
	params.RiskFreeRate = 0.25;
	params.DividendRate = 0.0;
	params.Volatility = 0.3;
	params.Spot = 65.0;

	// Generate fixing dates:
	auto today = date(day_clock::local_day());
	std::set<date> fixingDates;
	fixingDates.emplace(today);
	std::size_t periods = 364;
	for (auto t = 1; t <= periods; ++t) {
		fixingDates.emplace(today + date_duration(t));
	}

	// Create a binomial lattice:
	Lattice<LatticeType::Binomial, double, date> biLattice{ fixingDates };

	// Prepare Inversion functor:
	std::size_t numberTimePoints{ periods + 1 };
	PeizerPrattSecondInversion<> ppi{ numberTimePoints };

	// Create a model:
	LeisenReimerModel<> lr{ params,periods,ppi };

	// Name of the model:
	std::cout << decltype(lr)::name() << "\n";

	// Prepare delta times:
	double year{ 365.0 };
	auto fd = biLattice.fixingDates();
	std::vector<double> deltaTimes(fd.size() - 1);
	for (auto i = 0; i < deltaTimes.size(); ++i) {
		deltaTimes[i] = ((fd[i + 1] - fd[i]).days() / year);
	}

	// Build the lattice, i.e. populate it with asset prices from model:
	typedef ForwardInduction<LatticeType::Binomial, date, std::vector<double>, double> forward_induction;
	forward_induction fwd_induction;


	// Prepare payoffs:
	double K = params.Strike;
	auto callPayoff = [&K](double stock) {
		return std::max(stock - K, 0.0);
	};
	auto putPayoff = [&K](double stock) {
		return std::max(K - stock, 0.0);
	};

	// Prepare adjusters form american options:
	auto callAdjuster = [&callPayoff](double& optionValue, double stock) {
		optionValue = std::max(optionValue, callPayoff(stock));
	};

	auto putAdjuster = [&putPayoff](double& optionValue, double stock) {
		optionValue = std::max(optionValue, putPayoff(stock));
	};


	// spawn threads to do the forward induction:
	{
		scoped_thread t1(std::thread([&]() {
			fwd_induction(biLattice, lr, deltaTimes, params.Spot);
		}));
	}

	// create all the lattices:
	auto euroPutLattice = biLattice;
	auto euroCallLattice = biLattice;
	auto americanCallLattice = biLattice;
	auto americanPutLattice = biLattice;

	// Build the option price lattices, i.e. populate it with option prices from model:
	typedef BackwardInduction<LatticeType::Binomial, date, std::vector<double>> backward_induction;
	backward_induction bwd_induction;

	{
		scoped_thread t1(std::thread([&]() {
			bwd_induction(euroPutLattice, lr, deltaTimes, putPayoff);
		}));
		scoped_thread t2(std::thread([&]() {
			bwd_induction( euroCallLattice, lr, deltaTimes, callPayoff);
		}));
		scoped_thread t3(std::thread([&]() {
			bwd_induction( americanCallLattice, lr, deltaTimes, callPayoff, callAdjuster);
		}));
		scoped_thread t4(std::thread([&]() {
			bwd_induction( americanPutLattice, lr, deltaTimes, putPayoff, putAdjuster);
		}));
	}

	//// print the asset prices:
	//std::cout << "Asset price tree: \n";
	//auto first = biLattice.begin();
	//auto last = std::next(first, 5);
	//print(biLattice, first, last);
	//std::cout << "\n";

	//// print the option prices:
	//std::cout << "Euro call option price tree: \n";
	//first = euroCallLattice.begin();
	//last = std::next(first, 5);
	//print(euroCallLattice, first, last);

	//std::cout << "Euro put option price tree: \n";
	//first = euroPutLattice.begin();
	//last = std::next(first, 5);
	//print(euroPutLattice, first, last);

	//std::cout << "American call option price tree: \n";
	//first = americanCallLattice.begin();
	//last = std::next(first, 5);
	//print(americanCallLattice, first, last);

	//std::cout << "American put option price tree: \n";
	//first = americanPutLattice.begin();
	//last = std::next(first, 5);
	//print(americanPutLattice, first, last);

	std::cout << "European call price: " << euroCallLattice.apex() << "\n";
	std::cout << "European put price: " << euroPutLattice.apex() << "\n";
	std::cout << "American call price: " << americanCallLattice.apex() << "\n";
	std::cout << "American put price: " << americanPutLattice.apex() << "\n";
	std::cout << "time dimension: " << biLattice.timeDimension() << "\n";
}

#endif //_LATTICE_EXAMPLES