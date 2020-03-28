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





#endif ///_LATTICE_GREEKS_T