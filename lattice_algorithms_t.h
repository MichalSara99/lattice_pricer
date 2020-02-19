#pragma once
#if !defined(_LATTICE_ALGORITHMS_T)
#define _LATTICE_ALGORITHMS_T


#include"lattice_algorithms.h"
#include"lattice_utility.h"
#include<iostream>
#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>


using namespace boost::gregorian;

// ============================= Testing forward-induction algos =======================

void indexedLatticeBinomialForwardInduction() {

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

void indexedLatticeTrinomialForwardInduction() {
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



void latticeBinomialForwardInduction() {

	auto today = date(day_clock::local_day());
	lattice_structure::Lattice<lattice_types::LatticeType::Binomial, double, date>
		la = { today,today + date_duration(2),today,today + date_duration(1) };

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

void latticeTrinomialForwardInduction() {

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

//
//// ============================= Testing backward-induction algos =======================
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


void mergeIndexedBinomial() {
	const int N = 3;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double> stockTree(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Binomial, double>::TreeType).name() << "\n";

	lattice_utility::print(stockTree,stockTree.begin(), stockTree.end());

	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value, double dt)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	lattice_algorithms::forward_induction<>(stockTree, fwdgen, 1.0, dt);
	std::cout << "After forward induction step:\n";
	lattice_utility::print(stockTree,stockTree.begin(), stockTree.end());

	auto optionTree = stockTree;

	lattice_types::LeafBackwardGenerator<double, double, double, double> backgen = [](double upValue, double downValue, double dt) {
		double p = 0.5;
		return (p*upValue + (1.0 - p)*downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	lattice_algorithms::backward_induction<>(optionTree, backgen, payoff, dt);
	std::cout << "After backward induction step:\n";
	lattice_utility::print(optionTree,optionTree.begin(), optionTree.end());

	std::cout << "Price: " << optionTree.apex() << "\n\n";

	std::cout << "Merging stockTree with optionTree:\n";
	auto mergedTree = lattice_algorithms::merge(stockTree, optionTree,lattice_types::Launch::Parallel);
	for (auto const &v : mergedTree.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << "(" << std::get<0>(e) << "," << std::get<1>(e) << "),";
		}
		std::cout << "]\n";
	}


}


void mergeIndexedTrinomial() {

	const int N = 3;
	const double maturity{ 1.0 };
	double dt = maturity / ((double)N);

	lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double> stockTree(N);
	std::cout << "type of tree: " << typeid(lattice_structure::IndexedLattice<lattice_types::LatticeType::Trinomial, double>::TreeType).name() << "\n";

	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());

	lattice_types::LeafForwardGenerator<double, double, double, double> fwdgen = [](double value, double dt)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	lattice_algorithms::forward_induction<>(stockTree, fwdgen, 1.0, dt);
	std::cout << "After forward induction step:\n";
	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());

	lattice_types::LeafBackwardGenerator<double, double, double, double, double> backgen = [](double upValue, double midValue, double downValue, double dt) {
		double p = (1.0 / 3.0);
		return (p*upValue + p * midValue + p * downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	auto optionTree = stockTree;

	lattice_algorithms::backward_induction<>(optionTree, backgen, payoff, dt);
	std::cout << "After backward induction step:\n";
	lattice_utility::print(optionTree,optionTree.begin(), optionTree.end());

	std::cout << "Price: " << optionTree.apex() << "\n\n";

	std::cout << "Merged stockTree and optionTree:\n";
	auto mergedTree = lattice_algorithms::merge(stockTree, optionTree);
	for (auto const &v : mergedTree.tree()) {
		std::cout << "[";
		for (auto const &e : v) {
			std::cout << "(" << std::get<0>(e) << "," << std::get<1>(e) << "),";
		}
		std::cout << "]\n";
	}

}


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

	lattice_types::LeafForwardGenerator<double, double, double> fwdgen = [](double value, double dt)->std::tuple<double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		return std::make_tuple(u*value, d*value);
	};

	lattice_algorithms::forward_induction(stockTree, fwdgen, 1.0, deltas);
	std::cout << "After forward induction step:\n";
	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());


	lattice_types::LeafBackwardGenerator<double, double, double, double> backgen = [](double upValue, double downValue, double dt) {
		double p = 0.5;
		return (p*upValue + (1.0 - p)*downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	auto optionTree = stockTree;

	lattice_algorithms::backward_induction<>(optionTree, backgen, payoff, deltas);
	std::cout << "After backward induction step:\n";
	lattice_utility::print(optionTree, optionTree.begin(), optionTree.end());


	std::cout << "Price: " << optionTree.apex() << "\n\n";

	auto mergedTree = lattice_algorithms::merge(stockTree, optionTree,lattice_types::Launch::Parallel);
	std::cout << "Merged stockTree and optionTree:\n";

	for (auto const &l : mergedTree.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << "(" << std::get<0>(e) << "," << std::get<1>(e) << "),";
		}
		std::cout << "]\n";
	}

}


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

	lattice_types::LeafForwardGenerator<double, double, double, double> fwdgen = [](double value, double dt)->std::tuple<double, double, double> {
		// Parameters to ensure recombining lattice (test purposes)
		double u = 2.0;
		double d = 1.0 / u;
		double m = (u + d)*0.5;
		return std::make_tuple(u*value, m*value, d*value);
	};

	lattice_algorithms::forward_induction(stockTree, fwdgen, 1.0, deltas);
	std::cout << "After forward induction step:\n";
	lattice_utility::print(stockTree, stockTree.begin(), stockTree.end());

	lattice_types::LeafBackwardGenerator<double, double, double, double, double> backgen = [](double upValue, double midValue, double downValue, double dt) {
		double p = (1.0 / 3.0);
		return (p*upValue + p * midValue + p * downValue);
	};

	lattice_types::Payoff<double, double> payoff = [](double S) {
		return S;
	};

	auto optionTree = stockTree;

	lattice_algorithms::backward_induction<>(optionTree, backgen, payoff, deltas);
	std::cout << "After backward induction step:\n";
	lattice_utility::print(optionTree, optionTree.begin(), optionTree.end());
	std::cout << "Price: " << optionTree.apex() << "\n\n";

	auto mergedTree = lattice_algorithms::merge(stockTree, optionTree);
	std::cout << "Merged stockTree and optionTree:\n";

	for (auto const &l : mergedTree.tree()) {
		std::cout << "(" << l.first << "): [";
		for (auto const &e : l.second) {
			std::cout << "(" << std::get<0>(e) << "," << std::get<1>(e) << "),";
		}
		std::cout << "]\n";
	}

}



#endif ///_LATTICE_ALGORITHMS_T