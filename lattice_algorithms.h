#pragma once
#if !defined(_LATTICE_ALGORITHMS)
#define _LATTICE_ALGORITHMS

#include"lattice_structure.h"

namespace lattice_algorithms {

	using lattice_types::LatticeType;
	using lattice_types::LeafForwardGenerator;
	using lattice_types::LeafBackwardGenerator;
	using lattice_types::Payoff;
	using lattice_structure::IndexedLattice;
	using lattice_structure::Lattice;

	//==================== Forward Induction ============================

	template<typename Node>
	void forward_induction(IndexedLattice<LatticeType::Binomial,Node>& lattice,
		LeafForwardGenerator<Node, Node,Node> const &generator, Node apex) {
		lattice(0,0) = apex;
		std::tuple<Node,Node> tuple;
		for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex();++t){
			for (std::size_t l = 0; l < lattice.nodesAt(t-1).size(); ++l) {
				tuple = generator(lattice(t - 1, l));
				lattice(t, l) = std::get<0>(tuple);
				lattice(t, l + 1) = std::get<1>(tuple);
			}
		}
	}

	template<typename Node>
		void forward_induction(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafForwardGenerator<Node, Node,Node,Node> const &generator, Node apex) {
		lattice(0, 0) = apex;
		std::tuple<Node,Node,Node> tuple;
		for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex(); ++t) {
			for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
				tuple = generator(lattice(t - 1, l));
				lattice(t, l) = std::get<0>(tuple);
				lattice(t, l + 1) = std::get<1>(tuple);
				lattice(t, l + 2) = std::get<2>(tuple);
			}
		}
	}


	template<typename Node,typename TimeAxis>
	void forward_induction(Lattice<LatticeType::Binomial,Node,TimeAxis>& lattice,
		LeafForwardGenerator<Node, Node,Node> const &generator, Node apex) {
		auto fixingDates = lattice.fixingDates();
		lattice(fixingDates[0],0) = apex;
		std::tuple<Node,Node> tuple;
		for (std::size_t t = 1; t < fixingDates.size();++t){
			for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t-1]).size(); ++l) {
				tuple = generator(lattice(fixingDates[t - 1], l));
				lattice(fixingDates[t], l) = std::get<0>(tuple);
				lattice(fixingDates[t], l + 1) = std::get<1>(tuple);
			}
		}
	}

	template<typename Node,typename TimeAxis>
		void forward_induction(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node,Node,Node> const &generator, Node apex) {
		auto fixingDates = lattice.fixingDates();
		lattice(fixingDates[0], 0) = apex;
		std::tuple<Node,Node,Node> tuple;
		for (std::size_t t = 1; t < fixingDates.size(); ++t) {
			for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
				tuple = generator(lattice(fixingDates[t - 1], l));
				lattice(fixingDates[t], l) = std::get<0>(tuple);
				lattice(fixingDates[t], l + 1) = std::get<1>(tuple);
				lattice(fixingDates[t], l + 2) = std::get<2>(tuple);
			}
		}
	}

	//==================== Backward Induction ============================

	
	template<typename Node>
	void backward_induction(IndexedLattice<LatticeType::Binomial, Node>& lattice,
		LeafBackwardGenerator<Node, Node, Node>const &backwardGenerator, Payoff<Node,Node> const &payoff) {

		std::size_t lastIdx = lattice.maxIndex();
		for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
			lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
		}
		
		for (auto n = lastIdx - 1; n > 0; --n) {
			for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
				lattice(n,i) = backwardGenerator(lattice(n + 1,i), lattice(n + 1, i + 1));
			}
		}
		lattice(0,0) = backwardGenerator(lattice(1, 0), lattice(1,1));
	}

	template<typename Node, typename TimeAxis>
	void backward_induction(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
		LeafBackwardGenerator<Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff) {

		auto fixingDates = lattice.fixingDates();
		TimeAxis lastDate = fixingDates.back();
		for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
			lattice(lastDate, i) = payoff(lattice(lastDate, i));
		}

		for (auto n = fixingDates.size() - 2; n > 0; --n) {
			for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
				lattice(fixingDates[n], i) = backwardGenerator(lattice(fixingDates[n + 1], i), lattice(fixingDates[n + 1], i + 1));
			}
		}
		lattice(fixingDates[0], 0) = backwardGenerator(lattice(fixingDates[1], 0), lattice(fixingDates[1], 1));
	}

	template<typename Node>
	void backward_induction(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
		LeafBackwardGenerator<Node,Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff) {

		std::size_t lastIdx = lattice.maxIndex();
		for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
			lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
		}

		for (auto n = lastIdx - 1; n > 0; --n) {
			for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
				lattice(n, i) = backwardGenerator(lattice(n + 1, i), lattice(n + 1, i + 1), lattice(n + 1, i + 2));
			}
		}
		lattice(0, 0) = backwardGenerator(lattice(1, 0), lattice(1, 1), lattice(1, 2));
	}

	template<typename Node, typename TimeAxis>
	void backward_induction(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
		LeafBackwardGenerator<Node, Node,Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff) {

		auto fixingDates = lattice.fixingDates();
		TimeAxis lastDate = fixingDates.back();
		for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
			lattice(lastDate, i) = payoff(lattice(lastDate, i));
		}

		for (auto n = fixingDates.size() - 2; n > 0; --n) {
			for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
				lattice(fixingDates[n], i) = backwardGenerator(lattice(fixingDates[n + 1], i), 
																lattice(fixingDates[n + 1], i + 1),
																lattice(fixingDates[n + 1], i + 2));
			}
		}
		lattice(fixingDates[0], 0) = backwardGenerator(lattice(fixingDates[1], 0), 
														lattice(fixingDates[1], 1),
														lattice(fixingDates[1], 2));
	}



}




#endif //_LATTICE_ALGORITHMS