#pragma once
#if !defined(_LATTICE_ALGORITHMS_T)
#define _LATTICE_ALGORITHMS_T


#include"lattice_algorithms.h"
#include<iostream>
#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>


using namespace boost::gregorian;

// ============================= Testing forward-induction algos =======================

void indexedLatticeBinomialForwardInduction() {

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(3);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value)->std::tuple<double,double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	lattice_algorithms::forward_induction<>(il, fwdgen, 1.0);
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

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(3);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	lattice_types::LeafForwardGenerator<double, double,double, double> fwdgen = [](double value)->std::tuple<double,double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value,d*value);
	};

	lattice_algorithms::forward_induction<>(il, fwdgen, 1.0);
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
		la = { today,today + date_duration(2),today,today + date_duration(1) };

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	lattice_algorithms::forward_induction(la, fwdgen, 1.0);
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
		la = { today,today + date_duration(2),today,today + date_duration(1) };

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	lattice_types::LeafForwardGenerator<double, double,double, double> fwdgen = [](double value)->std::tuple<double,double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	lattice_algorithms::forward_induction(la, fwdgen, 1.0);
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

// ============================= Testing backward-induction algos =======================

void indexedLatticeBinomialBackwardInduction() {

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> il(3);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	lattice_algorithms::forward_induction<>(il, fwdgen, 1.0);
	std::cout << "After forward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	lattice_types::LeafBackwardGenerator<double, double, double> backgen = [](double upValue,double downValue) {
		double p = 0.5;
		return (p*upValue + (1.0 - p)*downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	lattice_algorithms::backward_induction<>(il, backgen, payoff);
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

void indexedLatticeTrinomialBackwardInduction() {

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> il(3);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	lattice_types::LeafForwardGenerator<double, double, double, double> fwdgen = [](double value)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	lattice_algorithms::forward_induction<>(il, fwdgen, 1.0);
	std::cout << "After forward induction step:\n";
	for (auto const &v : il.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	lattice_types::LeafBackwardGenerator<double, double,double, double> backgen = [](double upValue, double midValue,double downValue) {
		double p = (1.0 / 3.0);
		return (p*upValue + p*midValue + p*downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	lattice_algorithms::backward_induction<>(il, backgen, payoff);
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



void latticeBinomialBackwardInduction() {

	auto today = date(day_clock::local_day());
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date>
		la = { today,today + date_duration(2),today,today + date_duration(1) };

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	lattice_algorithms::forward_induction(la, fwdgen, 1.0);
	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}


	lattice_types::LeafBackwardGenerator<double, double, double> backgen = [](double upValue, double downValue) {
		double p = 0.5;
		return (p*upValue + (1.0 - p)*downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	lattice_algorithms::backward_induction<>(la, backgen, payoff);
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

void latticeTrinomialBackwardInduction() {

	auto today = date(day_clock::local_day());
	lattice_structure::Lattice<lattice_types::LatticeType::Trinomial, double, date>
		la = { today,today + date_duration(2),today,today + date_duration(1) };

	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";

	lattice_types::LeafForwardGenerator<double, double, double, double> fwdgen = [](double value)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	lattice_algorithms::forward_induction(la, fwdgen, 1.0);
	std::cout << "After forward induction step:\n";
	for (auto const &l : la.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << e << ",";
		}
		std::cout << "]\n";
	}

	lattice_types::LeafBackwardGenerator<double, double, double, double> backgen = [](double upValue, double midValue, double downValue) {
		double p = (1.0 / 3.0);
		return (p*upValue + p * midValue + p * downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	lattice_algorithms::backward_induction<>(la, backgen, payoff);
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








#endif ///_LATTICE_ALGORITHMS_T