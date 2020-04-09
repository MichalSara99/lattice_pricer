#pragma once
#if !defined(_LATTICE_ALGORITHMS_T)
#define _LATTICE_ALGORITHMS_T


#include"lattice_algorithms.h"
#include"lattice_utility.h"
#include<iostream>
#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>
#include<boost/math/special_functions/binomial.hpp>


using namespace boost::gregorian;

// ============================= Testing forward-induction algos =======================

void indexedLatticeBinomialForwardInduction() {

	const int N = 12;
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

	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value, double dt)->std::tuple<double,double> {
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

	std::cout << "Price: " << il.apex() << "\n\n";
}


void indexedLatticeBinomialForwardInductionNew() {

	const int N = 12;
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

	auto fwdgen = [](double value, double dt)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial, std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, fwdgen, dt,1.0);

	std::cout << "After forward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << il.apex() << "\n\n";
}

void indexedLatticeTrinomialForwardInduction() {
	const int N = 12;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	lattice_types::LeafForwardGenerator<double, double,double, double> fwdgen = [](double value,double dt)->std::tuple<double,double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value,d*value);
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

	std::cout << "Price: " << il.apex() << "\n\n";
}

void indexedLatticeTrinomialForwardInductionNew() {
	const int N = 12;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	auto fwdgen = [](double value, double dt)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial, std::size_t, double, double> forward_trinomial_induction;

	forward_trinomial_induction fwd_induction;
	fwd_induction(il, fwdgen, dt,1.0 );
	std::cout << "After forward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

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
	for (std::size_t t = 0; t < fixingDates.size()-1; ++t) {
		deltas[t] = ((fixingDates[t + 1] - fixingDates[t]).days() / daysInYear);
	}

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value,double dt)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	lattice_algorithms::forward_induction(la, fwdgen, 1.0, deltas);
	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << la.apex() << "\n\n";
}


void latticeBinomialForwardInductionNew() {

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

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	auto fwdgen = [](double value, double dt)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, fwdgen, deltas, 1.0);

	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

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

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	lattice_types::LeafForwardGenerator<double, double,double, double> fwdgen = [](double value,double dt)->std::tuple<double,double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	lattice_algorithms::forward_induction(la, fwdgen, 1.0, deltas);
	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << la.apex() << "\n\n";
}


void latticeTrinomialForwardInductionNew() {

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

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	auto fwdgen = [](double value, double dt)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>, double> forward_trinomial_induction;

	forward_trinomial_induction fwd_induction;
	fwd_induction(la, fwdgen, deltas, 1.0);
	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << la.apex() << "\n\n";
}

//
//// === discretely paid dividends ===
//
//
void indexedLatticeBinomialForwardInductionDividend() {

	const int N = 12;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);
	// dicrete dividends:
	std::map<std::size_t, double> dividends = { {3,0.01},{9,0.015} };

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

	lattice_algorithms::forward_induction<>(il, fwdgen, 1.0, dt, dividends);
	std::cout << "After forward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << il.apex() << "\n\n";
}

void indexedLatticeBinomialForwardInductionDividendNew() {

	const int N = 12;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);
	// dicrete dividends:
	std::map<std::size_t, double> dividends = { { 3,0.01 },{ 9,0.015 } };

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	auto fwdgen = [](double value, double dt)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, fwdgen, dt, 1.0, dividends);
	std::cout << "After forward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << il.apex() << "\n\n";
}

//

void indexedLatticeTrinomialForwardInductionDividends() {
	const int N = 12;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);
	std::map<std::size_t, double> dividends = { {3,0.01},{9,0.015} };

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	lattice_types::LeafForwardGenerator<double, double, double, double> fwdgen = [](double value, double dt)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	lattice_algorithms::forward_induction<>(il, fwdgen, 1.0, dt, dividends);
	std::cout << "After forward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << il.apex() << "\n\n";
}

void indexedLatticeTrinomialForwardInductionDividendsNew() {
	const int N = 12;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);
	std::map<std::size_t, double> dividends = { { 3,0.01 },{ 9,0.015 } };

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	auto fwdgen = [](double value, double dt)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double, double> forward_trinomial_induction;

	forward_trinomial_induction fwd_induction;
	fwd_induction(il, fwdgen, dt, 1.0, dividends);
	std::cout << "After forward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << il.apex() << "\n\n";
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

	std::map<date, double> dividends = { { today + date_duration(3),0.01},
										 { today + date_duration(9),0.015 } };

	// get fixing dates and compute deltas:
	auto fixingDates = la.fixingDates();
	double const daysInYear{ 365.0 };
	std::vector<double> deltas(fixingDates.size() - 1);
	for (std::size_t t = 0; t < fixingDates.size() - 1; ++t) {
		deltas[t] = ((fixingDates[t + 1] - fixingDates[t]).days() / daysInYear);
	}

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value, double dt)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	lattice_algorithms::forward_induction(la, fwdgen, 1.0, deltas, dividends);
	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << la.apex() << "\n\n";
}


void latticeBinomialForwardInductionDividendsNew() {

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

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	auto fwdgen = [](double value, double dt)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date,std::vector<double>, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, fwdgen, deltas, 1.0, dividends);

	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

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

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	lattice_types::LeafForwardGenerator<double, double, double, double> fwdgen = [](double value, double dt)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	lattice_algorithms::forward_induction(la, fwdgen, 1.0, deltas, dividends);
	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << la.apex() << "\n\n";
}


void latticeTrinomialForwardInductionDividendNew() {

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

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	auto fwdgen = [](double value, double dt)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>, double> forward_trinomial_induction;

	forward_trinomial_induction fwd_induction;
	fwd_induction(la, fwdgen, deltas, 1.0, dividends);

	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << la.apex() << "\n\n";
}

//
////
////// ============================= Testing backward-induction algos =======================
////
//
void indexedLatticeBinomialBackwardInduction() {

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

	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value,double dt)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	lattice_algorithms::forward_induction<>(il, fwdgen, 1.0,dt);
	std::cout << "After forward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	lattice_types::LeafBackwardGenerator<double, double, double,double> backgen = [](double upValue,double downValue,double dt) {
		double p = 0.5;
		return (p*upValue + (1.0 - p)*downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	lattice_algorithms::backward_induction<>(il, backgen, payoff,dt);
	std::cout << "After backward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << il.apex() << "\n\n";
}


void indexedLatticeBinomialBackwardInductionNew() {

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

	auto fwdgen = [](double value, double dt)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t,double, double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(il, fwdgen,dt, 1.0);

	std::cout << "After forward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	auto backgen = [](double upValue, double downValue, double dt) {
		double p = 0.5;
		return (p*upValue + (1.0 - p)*downValue);
	};

	auto payoff = [](double S) {
		return S;
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		std::size_t, double> backward_binomial_induction;

	backward_binomial_induction bwrd_induction;
	bwrd_induction(il, backgen, dt, payoff);

	std::cout << "After backward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

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

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	lattice_types::LeafForwardGenerator<double, double, double, double> fwdgen = [](double value,double dt)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
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

	lattice_types::LeafBackwardGenerator<double, double,double, double,double> backgen = [](double upValue, double midValue,double downValue,double dt) {
		double p = (1.0 / 3.0);
		return (p*upValue + p*midValue + p*downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	lattice_algorithms::backward_induction<>(il, backgen, payoff,dt);
	std::cout << "After backward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << il.apex() << "\n\n";

}

void indexedLatticeTrinomialBackwardInductionNew() {

	const int N = 3;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	auto fwdgen = [](double value, double dt)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double, double> forward_trinomial_induction;

	forward_trinomial_induction fwd_induction;
	fwd_induction(il, fwdgen, dt, 1.0);
	std::cout << "After forward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	auto backgen = [](double upValue, double midValue, double downValue, double dt) {
		double p = (1.0 / 3.0);
		return (p*upValue + p * midValue + p * downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial,
		std::size_t, double> backward_trinomial_induction;

	backward_trinomial_induction bwrd_induction;
	bwrd_induction(il, backgen, dt, payoff);

	std::cout << "After backward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << il.apex() << "\n\n";

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

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value,double dt)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	lattice_algorithms::forward_induction(la, fwdgen, 1.0, deltas);
	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}


	lattice_types::LeafBackwardGenerator<double, double, double,double> backgen = [](double upValue, double downValue,double dt) {
		double p = 0.5;
		return (p*upValue + (1.0 - p)*downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	lattice_algorithms::backward_induction<>(la, backgen, payoff, deltas);
	std::cout << "After backward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}


	std::cout << "Price: " << la.apex() << "\n\n";
}

void latticeBinomialBackwardInductionNew() {


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

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	auto fwdgen = [](double value, double dt)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>,double> forward_binomial_induction;

	forward_binomial_induction fwd_induction;
	fwd_induction(la, fwdgen, deltas, 1.0);

	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}


	auto backgen = [](double upValue, double downValue, double dt) {
		double p = 0.5;
		return (p*upValue + (1.0 - p)*downValue);
	};

	auto payoff = [](double S) {
		return S;
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Binomial,
		date, std::vector<double>> backward_binomial_induction;

	backward_binomial_induction bwrd_induction;
	bwrd_induction(la, backgen,deltas, payoff);

	std::cout << "After backward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}


	std::cout << "Price: " << la.apex() << "\n\n";
}

//
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

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	lattice_types::LeafForwardGenerator<double, double, double, double> fwdgen = [](double value,double dt)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	lattice_algorithms::forward_induction(la, fwdgen, 1.0, deltas);
	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	lattice_types::LeafBackwardGenerator<double, double, double, double,double> backgen = [](double upValue, double midValue, double downValue,double dt) {
		double p = (1.0 / 3.0);
		return (p*upValue + p * midValue + p * downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	lattice_algorithms::backward_induction<>(la, backgen, payoff, deltas);
	std::cout << "After backward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << la.apex() << "\n\n";
}

void latticeTrinomialBackwardInductionNew() {

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

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	auto fwdgen = [](double value, double dt)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	typedef lattice_algorithms::ForwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>, double> forward_trinomial_induction;

	forward_trinomial_induction fwd_induction;
	fwd_induction(la, fwdgen, deltas, 1.0);

	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	auto backgen = [](double upValue, double midValue, double downValue, double dt) {
		double p = (1.0 / 3.0);
		return (p*upValue + p * midValue + p * downValue);
	};

	auto payoff = [](double S) {
		return S;
	};

	typedef lattice_algorithms::BackwardInduction<lattice_types::LatticeType::Trinomial,
		date, std::vector<double>> backward_trinomial_induction;

	backward_trinomial_induction bwrd_induction;
	bwrd_induction(la, backgen, deltas, payoff);

	std::cout << "After backward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	std::cout << "Price: " << la.apex() << "\n\n";
}


//
//
//void mergeIndexedBinomial() {
//	const int N = 3;
//	const double maturity{ 1.0 };
//	double dt = maturity / ((double)N);
//
//	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> stockTree(N);
//	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";
//
//	lattice_utility::print(stockTree,stockTree.begin(), stockTree.end());
//
//	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value, double dt)->std::tuple<double, double> {
//		// Parameters to ensure recombining lattice (test purposes)
//		double u = 2.0;
//		double d = 1.0 / u;
//		return std::make_tuple(u*value, d*value);
//	};
//
//	lattice_algorithms::forward_induction<>(stockTree, fwdgen, 1.0, dt);
//	std::cout << "After forward induction step:\n";
//	lattice_utility::print(stockTree,stockTree.begin(), stockTree.end());
//
//	auto optionTree = stockTree;
//
//	lattice_types::LeafBackwardGenerator<double, double, double, double> backgen = [](double upValue, double downValue, double dt) {
//		double p = 0.5;
//		return (p*upValue + (1.0 - p)*downValue);
//	};
//
//	lattice_types::Payoff<double, double> payoff = [](double S) {
//		return S;
//	};
//
//	lattice_algorithms::backward_induction<>(optionTree, backgen, payoff, dt);
//	std::cout << "After backward induction step:\n";
//	lattice_utility::print(optionTree,optionTree.begin(), optionTree.end());
//
//	std::cout << "Price: " << optionTree.apex() << "\n\n";
//
//	std::cout << "Merging stockTree with optionTree:\n";
//	auto mergedTree = lattice_algorithms::merge(stockTree, optionTree,lattice_types::Launch::Parallel);
//	for (auto const &v : mergedTree.tree()) {
//		std::cout << "[";
//		for (auto const &e : v) {
//			std::cout << "(" << std::get<0>(e) << "," << std::get<1>(e) << "),";
//		}
//		std::cout << "]\n";
//	}
//
//
//}
//
//
//void mergeIndexedTrinomial() {
//
//	const int N = 3;
//	const double maturity{ 1.0 };
//	double dt = maturity / ((double)N);
//
//	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> stockTree(N);
//	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";
//
//	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());
//
//	lattice_types::LeafForwardGenerator<double, double, double, double> fwdgen = [](double value, double dt)->std::tuple<double, double, double> {
//		// Parameters to ensure recombining lattice (test purposes)
//		double u = 2.0;
//		double d = 1.0 / u;
//		double m = (u + d)*0.5;
//		return std::make_tuple(u*value, m*value, d*value);
//	};
//
//	lattice_algorithms::forward_induction<>(stockTree, fwdgen, 1.0, dt);
//	std::cout << "After forward induction step:\n";
//	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());
//
//	lattice_types::LeafBackwardGenerator<double, double, double, double, double> backgen = [](double upValue, double midValue, double downValue, double dt) {
//		double p = (1.0 / 3.0);
//		return (p*upValue + p * midValue + p * downValue);
//	};
//
//	lattice_types::Payoff<double, double> payoff = [](double S) {
//		return S;
//	};
//
//	auto optionTree = stockTree;
//
//	lattice_algorithms::backward_induction<>(optionTree, backgen, payoff, dt);
//	std::cout << "After backward induction step:\n";
//	lattice_utility::print(optionTree,optionTree.begin(), optionTree.end());
//
//	std::cout << "Price: " << optionTree.apex() << "\n\n";
//
//	std::cout << "Merged stockTree and optionTree:\n";
//	auto mergedTree = lattice_algorithms::merge(stockTree, optionTree);
//	for (auto const &v : mergedTree.tree()) {
//		std::cout << "[";
//		for (auto const &e : v) {
//			std::cout << "(" << std::get<0>(e) << "," << std::get<1>(e) << "),";
//		}
//		std::cout << "]\n";
//	}
//
//}
//
//
//void mergeBinomial() {
//
//	auto today = date(day_clock::local_day());
//	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date>
//		stockTree = { today,today + date_duration(2),today,today + date_duration(1) };
//
//	// get fixing dates and compute deltas:
//	auto fixingDates = stockTree.fixingDates();
//	double const daysInYear{ 365.0 };
//	std::vector<double> deltas(fixingDates.size() - 1);
//	for (std::size_t t = 0; t < fixingDates.size() - 1; ++t) {
//		deltas[t] = ((fixingDates[t + 1] - fixingDates[t]).days() / daysInYear);
//	}
//
//	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());
//
//	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value, double dt)->std::tuple<double, double> {
//		// Parameters to ensure recombining lattice (test purposes)
//		double u = 2.0;
//		double d = 1.0 / u;
//		return std::make_tuple(u*value, d*value);
//	};
//
//	lattice_algorithms::forward_induction(stockTree, fwdgen, 1.0, deltas);
//	std::cout << "After forward induction step:\n";
//	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());
//
//
//	lattice_types::LeafBackwardGenerator<double, double, double, double> backgen = [](double upValue, double downValue, double dt) {
//		double p = 0.5;
//		return (p*upValue + (1.0 - p)*downValue);
//	};
//
//	lattice_types::Payoff<double, double> payoff = [](double S) {
//		return S;
//	};
//
//	auto optionTree = stockTree;
//
//	lattice_algorithms::backward_induction<>(optionTree, backgen, payoff, deltas);
//	std::cout << "After backward induction step:\n";
//	lattice_utility::print(optionTree, optionTree.begin(), optionTree.end());
//
//
//	std::cout << "Price: " << optionTree.apex() << "\n\n";
//
//	auto mergedTree = lattice_algorithms::merge(stockTree, optionTree,lattice_types::Launch::Parallel);
//	std::cout << "Merged stockTree and optionTree:\n";
//
//	for (auto const &l : mergedTree.tree()) {
//		std::cout << "(" << l.first << "): [";
//		for (auto const &e : l.second) {
//			std::cout << "(" << std::get<0>(e) << "," << std::get<1>(e) << "),";
//		}
//		std::cout << "]\n";
//	}
//
//}
//
//
//void mergeTrinomial() {
//	auto today = date(day_clock::local_day());
//	lattice_structure::Lattice<lattice_types::LatticeType::Trinomial, double, date>
//		stockTree = { today,today + date_duration(2),today,today + date_duration(1) };
//
//	// get fixing dates and compute deltas:
//	auto fixingDates = stockTree.fixingDates();
//	double const daysInYear{ 365.0 };
//	std::vector<double> deltas(fixingDates.size() - 1);
//	for (std::size_t t = 0; t < fixingDates.size() - 1; ++t) {
//		deltas[t] = ((fixingDates[t + 1] - fixingDates[t]).days() / daysInYear);
//	}
//
//	lattice_utility::print(stockTree,stockTree.begin(), stockTree.end());
//
//	lattice_types::LeafForwardGenerator<double, double, double, double> fwdgen = [](double value, double dt)->std::tuple<double, double, double> {
//		// Parameters to ensure recombining lattice (test purposes)
//		double u = 2.0;
//		double d = 1.0 / u;
//		double m = (u + d)*0.5;
//		return std::make_tuple(u*value, m*value, d*value);
//	};
//
//	lattice_algorithms::forward_induction(stockTree, fwdgen, 1.0, deltas);
//	std::cout << "After forward induction step:\n";
//	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());
//
//	lattice_types::LeafBackwardGenerator<double, double, double, double, double> backgen = [](double upValue, double midValue, double downValue, double dt) {
//		double p = (1.0 / 3.0);
//		return (p*upValue + p * midValue + p * downValue);
//	};
//
//	lattice_types::Payoff<double, double> payoff = [](double S) {
//		return S;
//	};
//
//	auto optionTree = stockTree;
//
//	lattice_algorithms::backward_induction<>(optionTree, backgen, payoff, deltas);
//	std::cout << "After backward induction step:\n";
//	lattice_utility::print(optionTree, optionTree.begin(), optionTree.end());
//	std::cout << "Price: " << optionTree.apex() << "\n\n";
//
//	auto mergedTree = lattice_algorithms::merge(stockTree, optionTree);
//	std::cout << "Merged stockTree and optionTree:\n";
//
//	for (auto const &l : mergedTree.tree()) {
//		std::cout << "(" << l.first << "): [";
//		for (auto const &e : l.second) {
//			std::cout << "(" << std::get<0>(e) << "," << std::get<1>(e) << "),";
//		}
//		std::cout << "]\n";
//	}
//
//}


//
//void pascalTriangleIndexedTest() {
//
//	using lattice_structure::IndexedLattice;
//	using lattice_types::LatticeType;
//
//	std::size_t periods = 10;
//
//	IndexedLattice<LatticeType::Binomial, long> pascalLattice{ periods };
//
//	pascalLattice(0, 0) = 1;
//	for (auto t = pascalLattice.minIndex() + 1; t <= pascalLattice.maxIndex();++t) {
//
//		// Edges:
//		pascalLattice(t, 0) = pascalLattice(t - 1, 0);
//		pascalLattice(t, pascalLattice.nodesAt(t).size() - 1) = pascalLattice(t - 1, pascalLattice.nodesAt(t - 1).size() - 1);
//
//		// Inner nodes:
//		for (auto i = 1; i < pascalLattice.nodesAt(t).size() - 1; ++i) {
//			pascalLattice(t, i) = pascalLattice(t - 1, i - 1) + pascalLattice(t - 1, i);
//		}
//	}
//
//	// test against binomial coefficients in boost:
//	std::size_t counter{ 0 };
//	long latticeValue{ 0 };
//	double boostValue{};
//	for (auto t = pascalLattice.minIndex(); t <= pascalLattice.maxIndex(); ++t) {
//		for (auto k = 0; k < pascalLattice.nodesAt(t).size(); ++k) {
//			latticeValue = pascalLattice(t, k);
//			boostValue = boost::math::binomial_coefficient<double>(t, k);
//			if (static_cast<double>(latticeValue) != boostValue) {
//				counter++;
//			}
//		}
//	}
//	std::cout << "Number of differences: " << counter << "\n";
//
//}
//
//void pascalTriangleTest() {
//	using lattice_structure::Lattice;
//	using lattice_types::LatticeType;
//
//	const int periods = 10;
//
//	auto today = date(day_clock::local_day());
//	std::set<date> fDates;
//	fDates.emplace(today);
//
//	for (std::size_t p = 1; p <= periods; ++p) {
//		fDates.emplace(today + date_duration(p));
//	}
//
//	Lattice<LatticeType::Binomial, long, date> pascalLattice{ fDates };
//
//	pascalLattice(today, 0) = 1;
//
//	auto itr = pascalLattice.cbegin();
//	auto begin = std::next(pascalLattice.cbegin(),1);
//	auto end = pascalLattice.cend();
//
//	for (; begin != end; ++begin, ++itr) {
//
//		// Edges:
//		pascalLattice(begin->first, 0) = pascalLattice(itr->first, 0);
//		pascalLattice(begin->first, pascalLattice.nodesAt(begin->first).size() - 1) =
//			pascalLattice(itr->first, pascalLattice.nodesAt(itr->first).size() - 1);
//
//		// Inner nodes:
//		for (auto i = 1; i < pascalLattice.nodesAt(begin->first).size() - 1; ++i) {
//			pascalLattice(begin->first, i) = pascalLattice(itr->first, i - 1) + pascalLattice(itr->first, i);
//		}
//	}
//
//	// test against binomial coefficients in boost:
//	std::size_t counter{ 0 };
//	long latticeValue{ 0 };
//	double boostValue{};
//	std::size_t t{ 0 };
//	begin = pascalLattice.cbegin();
//	end = pascalLattice.cend();
//	for (; begin!= end; ++begin,++t) {
//		for (auto k = 0; k < pascalLattice.nodesAt(begin->first).size(); ++k) {
//			latticeValue = pascalLattice(begin->first, k);
//			boostValue = boost::math::binomial_coefficient<double>(t, k);
//			if (static_cast<double>(latticeValue) != boostValue) {
//				counter++;
//			}
//		}
//	}
//	std::cout << "Number of differences: " << counter << "\n";
//}



#endif ///_LATTICE_ALGORITHMS_T