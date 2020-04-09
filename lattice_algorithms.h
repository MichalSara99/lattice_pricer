#pragma once
#if !defined(_LATTICE_ALGORITHMS)
#define _LATTICE_ALGORITHMS

#include<cassert>
#include<ppl.h>
#include"lattice_structure.h"
#include"lattice_forward_traversals.h"
#include"lattice_backward_traversals.h"
#include"lattice_macros.h"


namespace lattice_algorithms {

	using lattice_types::LatticeType;
	using lattice_types::LeafForwardGenerator;
	using lattice_types::LeafBackwardGenerator;
	using lattice_types::Payoff;
	using lattice_types::PayoffAdjuster;
	using lattice_structure::IndexedLattice;
	using lattice_structure::Lattice;
	using lattice_forward_traversals::ForwardTraversal;
	using lattice_backward_traversals::BackwardTraversal;
	

	// ==============================================================================
	// ==================== Forward Induction Algorithms ============================
	// ==============================================================================

	namespace {

		template<typename Node, typename DeltaTime>
		void _forward_induction_indexed_binomial_impl(IndexedLattice<LatticeType::Binomial, Node>& lattice,
			LeafForwardGenerator<Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,std::true_type) {
			assert(deltaTime.size() == lattice.maxIndex());
			lattice.at(0, 0) = apex;
			std::tuple<Node, Node> tuple;
			for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice.at(t - 1, l), deltaTime[t-1]);
					lattice.at(t, l) = std::get<0>(tuple);
					lattice.at(t, l + 1) = std::get<1>(tuple);
				}
			}
		}

		template<typename Node, typename DeltaTime>
		void _forward_induction_indexed_binomial_impl(IndexedLattice<LatticeType::Binomial, Node>& lattice,
			LeafForwardGenerator<Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime, std::false_type) {
			lattice.at(0, 0) = apex;
			std::tuple<Node, Node> tuple;
			for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice.at(t - 1, l), deltaTime);
					lattice.at(t, l) = std::get<0>(tuple);
					lattice.at(t, l + 1) = std::get<1>(tuple);
				}
			}
		}


		template<typename Node,typename DeltaTime>
		void _forward_induction_indexed_trinomial_impl(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,std::true_type) {
			assert(deltaTime.size() == lattice.maxIndex());
			lattice.at(0, 0) = apex;
			std::tuple<Node, Node, Node> tuple;
			for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice.at(t - 1, l),deltaTime[t-1]);
					lattice.at(t, l) = std::get<0>(tuple);
					lattice.at(t, l + 1) = std::get<1>(tuple);
					lattice.at(t, l + 2) = std::get<2>(tuple);
				}
			}
		}

		template<typename Node, typename DeltaTime>
		void _forward_induction_indexed_trinomial_impl(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex,DeltaTime deltaTime,std::false_type) {
			lattice.at(0, 0) = apex;
			std::tuple<Node, Node, Node> tuple;
			for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice.at(t - 1, l), deltaTime);
					lattice.at(t, l) = std::get<0>(tuple);
					lattice.at(t, l + 1) = std::get<1>(tuple);
					lattice.at(t, l + 2) = std::get<2>(tuple);
				}
			}
		}


		template<typename Node, typename TimeAxis,typename DeltaTime>
		void _forward_induction_lattice_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node, Node> const &generator, Node apex,DeltaTime deltaTime,std::true_type) {
			auto fixingDates = lattice.fixingDates();
			assert(deltaTime.size() == (fixingDates.size() - 1));
			lattice.at(fixingDates[0], 0) = apex;
			std::tuple<Node, Node> tuple;
			for (std::size_t t = 1; t < fixingDates.size(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice.at(fixingDates[t - 1], l), deltaTime[t - 1]);
					lattice.at(fixingDates[t], l) = std::get<0>(tuple);
					lattice.at(fixingDates[t], l + 1) = std::get<1>(tuple);
				}
			}
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _forward_induction_lattice_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime, std::false_type) {
			auto fixingDates = lattice.fixingDates();
			lattice.at(fixingDates[0], 0) = apex;
			std::tuple<Node, Node> tuple;
			for (std::size_t t = 1; t < fixingDates.size(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice.at(fixingDates[t - 1], l), deltaTime);
					lattice.at(fixingDates[t], l) = std::get<0>(tuple);
					lattice.at(fixingDates[t], l + 1) = std::get<1>(tuple);
				}
			}
		}

		template<typename Node, typename TimeAxis,typename DeltaTime>
		void _forward_induction_lattice_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex,DeltaTime deltaTime,std::true_type) {
			auto fixingDates = lattice.fixingDates();
			assert(deltaTime.size() == (fixingDates.size() - 1));
			lattice.at(fixingDates[0], 0) = apex;
			std::tuple<Node, Node, Node> tuple;
			for (std::size_t t = 1; t < fixingDates.size(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice.at(fixingDates[t - 1], l), deltaTime[t - 1]);
					lattice.at(fixingDates[t], l) = std::get<0>(tuple);
					lattice.at(fixingDates[t], l + 1) = std::get<1>(tuple);
					lattice.at(fixingDates[t], l + 2) = std::get<2>(tuple);
				}
			}
		}

		template<typename Node, typename TimeAxis,typename DeltaTime>
		void _forward_induction_lattice_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime, std::false_type) {
			auto fixingDates = lattice.fixingDates();
			lattice.at(fixingDates[0], 0) = apex;
			std::tuple<Node, Node, Node> tuple;
			for (std::size_t t = 1; t < fixingDates.size(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice.at(fixingDates[t - 1], l), deltaTime);
					lattice.at(fixingDates[t], l) = std::get<0>(tuple);
					lattice.at(fixingDates[t], l + 1) = std::get<1>(tuple);
					lattice.at(fixingDates[t], l + 2) = std::get<2>(tuple);
				}
			}
		}

		// =============================================
		// ===== overloads for discrete dividends ======
		// =============================================

		template<typename Node, typename DeltaTime>
		void _forward_induction_indexed_binomial_impl(IndexedLattice<LatticeType::Binomial, Node>& lattice,
			LeafForwardGenerator<Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,
			std::map<std::size_t, Node> const &dividendData, std::true_type) {
			assert(deltaTime.size() == lattice.maxIndex());
			lattice.at(0, 0) = apex;
			std::tuple<Node, Node> tuple;
			for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice.at(t - 1, l), deltaTime[t - 1]);
					lattice.at(t, l) = std::get<0>(tuple);
					lattice.at(t, l + 1) = std::get<1>(tuple);
				}
			}
			// Adjust the tree for discretely paying divdiends: 
			auto first = dividendData.begin();
			auto last = dividendData.end();
			if (first == last)return;
			std::size_t firstExIdx = first->first;
			if (firstExIdx > lattice.maxIndex())return;

			Node factor{ 1.0 };
			for (std::size_t t = firstExIdx; t <= lattice.maxIndex(); ++t) {
				auto nextExItr = dividendData.find(t);
				factor = 1.0;
				if (nextExItr != last) {
					factor *= (factor - nextExItr->second);
				}
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice.at(t - 1, l), deltaTime[t - 1]);
					lattice.at(t, l) = factor * std::get<0>(tuple);
					lattice.at(t, l + 1) = factor * std::get<1>(tuple);
				}
			}
		}

		template<typename Node, typename DeltaTime>
		void _forward_induction_indexed_binomial_impl(IndexedLattice<LatticeType::Binomial, Node>& lattice,
			LeafForwardGenerator<Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,
			std::map<std::size_t, Node> const &dividendData, std::false_type) {
			lattice.at(0, 0) = apex;
			std::tuple<Node, Node> tuple;
			for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice.at(t - 1, l), deltaTime);
					lattice.at(t, l) = std::get<0>(tuple);
					lattice.at(t, l + 1) = std::get<1>(tuple);
				}
			}
			// Adjust the tree for discretely paying divdiends: 
			auto first = dividendData.begin();
			auto last = dividendData.end();
			if (first == last)return;
			std::size_t firstExIdx = first->first;
			if (firstExIdx > lattice.maxIndex())return;

			Node factor{ 1.0 };
			for (std::size_t t = firstExIdx; t <= lattice.maxIndex(); ++t) {
				auto nextExItr = dividendData.find(t);
				factor = 1.0;
				if (nextExItr != last) {
					factor *= (factor - nextExItr->second);
				}
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice.at(t - 1, l), deltaTime);
					lattice.at(t, l) = factor * std::get<0>(tuple);
					lattice.at(t, l + 1) = factor * std::get<1>(tuple);
				}
			}
		}

		template<typename Node, typename DeltaTime>
		void _forward_induction_indexed_trinomial_impl(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,
			std::map<std::size_t, Node> const &dividendData,std::true_type) {
			assert(deltaTime.size() == lattice.maxIndex());
			lattice.at(0, 0) = apex;
			std::tuple<Node, Node, Node> tuple;
			for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice.at(t - 1, l), deltaTime[t - 1]);
					lattice.at(t, l) = std::get<0>(tuple);
					lattice.at(t, l + 1) = std::get<1>(tuple);
					lattice.at(t, l + 2) = std::get<2>(tuple);
				}
			}
			// Adjust the tree for discretely paying divdiends: 
			auto first = dividendData.begin();
			auto last = dividendData.end();
			if (first == last)return;
			std::size_t firstExIdx = first->first;
			if (firstExIdx > lattice.maxIndex())return;

			Node factor{ 1.0 };
			for (std::size_t t = firstExIdx; t <= lattice.maxIndex(); ++t) {
				auto nextExItr = dividendData.find(t);
				factor = 1.0;
				if (nextExItr != last) {
					factor *= (factor - nextExItr->second);
				}
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice.at(t - 1, l), deltaTime[t - 1]);
					lattice.at(t, l) = factor * std::get<0>(tuple);
					lattice.at(t, l + 1) = factor * std::get<1>(tuple);
					lattice.at(t, l + 2) = factor * std::get<2>(tuple);
				}
			}
		}


		template<typename Node, typename DeltaTime>
		void _forward_induction_indexed_trinomial_impl(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,
			std::map<std::size_t, Node> const &dividendData,std::false_type) {
			lattice.at(0, 0) = apex;
			std::tuple<Node, Node, Node> tuple;
			for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice.at(t - 1, l), deltaTime);
					lattice.at(t, l) = std::get<0>(tuple);
					lattice.at(t, l + 1) = std::get<1>(tuple);
					lattice.at(t, l + 2) = std::get<2>(tuple);
				}
			}
			// Adjust the tree for discretely paying divdiends: 
			auto first = dividendData.begin();
			auto last = dividendData.end();
			if (first == last)return;
			std::size_t firstExIdx = first->first;
			if (firstExIdx > lattice.maxIndex())return;

			Node factor{ 1.0 };
			for (std::size_t t = firstExIdx; t <= lattice.maxIndex(); ++t) {
				auto nextExItr = dividendData.find(t);
				factor = 1.0;
				if (nextExItr != last) {
					factor *= (factor - nextExItr->second);
				}
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice.at(t - 1, l), deltaTime);
					lattice.at(t, l) = factor * std::get<0>(tuple);
					lattice.at(t, l + 1) = factor * std::get<1>(tuple);
					lattice.at(t, l + 2) = factor * std::get<2>(tuple);
				}
			}
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _forward_induction_lattice_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,
			std::map<TimeAxis,Node> const &dividendData,std::true_type) {
			auto fixingDates = lattice.fixingDates();
			assert(deltaTime.size() == (fixingDates.size() - 1));
			lattice.at(fixingDates[0], 0) = apex;
			std::tuple<Node, Node> tuple;
			for (std::size_t t = 1; t < fixingDates.size(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice.at(fixingDates[t - 1], l), deltaTime[t - 1]);
					lattice.at(fixingDates[t], l) = std::get<0>(tuple);
					lattice.at(fixingDates[t], l + 1) = std::get<1>(tuple);
				}
			}

			// Adjust the tree for discretely paying divdiends: 
			auto first = dividendData.begin();
			auto last = dividendData.end();
			if (first == last)return;
			auto fd_itr = std::find(fixingDates.cbegin(), fixingDates.cend(), first->first);
			if (fd_itr == fixingDates.cend())return;

			Node factor{ 1.0 };
			std::size_t firstExIdx = std::distance(fixingDates.cbegin(), fd_itr);
			for (std::size_t t = firstExIdx; t < fixingDates.size(); ++t) {
				auto nextExItr = dividendData.find(fixingDates[t]);
				factor = 1.0;
				if (nextExItr != last) {
					factor *= (factor - nextExItr->second);
				}
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice.at(fixingDates[t - 1], l), deltaTime[t - 1]);
					lattice.at(fixingDates[t], l) = factor * std::get<0>(tuple);
					lattice.at(fixingDates[t], l + 1) = factor * std::get<1>(tuple);
				}
			}
		}
		
		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _forward_induction_lattice_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,
			std::map<TimeAxis, Node> const &dividendData,std::false_type) {
			auto fixingDates = lattice.fixingDates();
			lattice.at(fixingDates[0], 0) = apex;
			std::tuple<Node, Node> tuple;
			for (std::size_t t = 1; t < fixingDates.size(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice.at(fixingDates[t - 1], l), deltaTime);
					lattice.at(fixingDates[t], l) = std::get<0>(tuple);
					lattice.at(fixingDates[t], l + 1) = std::get<1>(tuple);
				}
			}

			// Adjust the tree for discretely paying divdiends: 
			auto first = dividendData.begin();
			auto last = dividendData.end();
			if (first == last)return;
			auto fd_itr = std::find(fixingDates.cbegin(), fixingDates.cend(), first->first);
			if (fd_itr == fixingDates.cend())return;

			Node factor{ 1.0 };
			std::size_t firstExIdx = std::distance(fixingDates.cbegin(), fd_itr);
			for (std::size_t t = firstExIdx; t < fixingDates.size(); ++t) {
				auto nextExItr = dividendData.find(fixingDates[t]);
				factor = 1.0;
				if (nextExItr != last) {
					factor *= (factor - nextExItr->second);
				}
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice.at(fixingDates[t - 1], l), deltaTime);
					lattice.at(fixingDates[t], l) = factor * std::get<0>(tuple);
					lattice.at(fixingDates[t], l + 1) = factor * std::get<1>(tuple);
				}
			}
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _forward_induction_lattice_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,
			std::map<TimeAxis,Node> const &dividendData,std::true_type) {
			auto fixingDates = lattice.fixingDates();
			assert(deltaTime.size() == (fixingDates.size() - 1));
			lattice.at(fixingDates[0], 0) = apex;
			std::tuple<Node, Node, Node> tuple;
			for (std::size_t t = 1; t < fixingDates.size(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice.at(fixingDates[t - 1], l), deltaTime[t - 1]);
					lattice.at(fixingDates[t], l) = std::get<0>(tuple);
					lattice.at(fixingDates[t], l + 1) = std::get<1>(tuple);
					lattice.at(fixingDates[t], l + 2) = std::get<2>(tuple);
				}
			}

			// Adjust the tree for discretely paying divdiends: 
			auto first = dividendData.begin();
			auto last = dividendData.end();
			if (first == last)return;
			auto fd_itr = std::find(fixingDates.cbegin(), fixingDates.cend(), first->first);
			if (fd_itr == fixingDates.cend())return;

			Node factor{ 1.0 };
			std::size_t firstExIdx = std::distance(fixingDates.cbegin(), fd_itr);
			for (std::size_t t = firstExIdx; t < fixingDates.size(); ++t) {
				auto nextExItr = dividendData.find(fixingDates[t]);
				factor = 1.0;
				if (nextExItr != last) {
					factor *= (factor - nextExItr->second);
				}
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice.at(fixingDates[t - 1], l), deltaTime[t - 1]);
					lattice.at(fixingDates[t], l) = factor * std::get<0>(tuple);
					lattice.at(fixingDates[t], l + 1) = factor * std::get<1>(tuple);
					lattice.at(fixingDates[t], l + 2) = factor * std::get<2>(tuple);
				}
			}
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _forward_induction_lattice_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,
			std::map<TimeAxis,Node> const &dividendData,std::false_type) {
			auto fixingDates = lattice.fixingDates();
			lattice.at(fixingDates[0], 0) = apex;
			std::tuple<Node, Node, Node> tuple;
			for (std::size_t t = 1; t < fixingDates.size(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice.at(fixingDates[t - 1], l), deltaTime);
					lattice.at(fixingDates[t], l) = std::get<0>(tuple);
					lattice.at(fixingDates[t], l + 1) = std::get<1>(tuple);
					lattice.at(fixingDates[t], l + 2) = std::get<2>(tuple);
				}
			}

			// Adjust the tree for discretely paying divdiends: 
			auto first = dividendData.begin();
			auto last = dividendData.end();
			if (first == last)return;
			auto fd_itr = std::find(fixingDates.cbegin(), fixingDates.cend(), first->first);
			if (fd_itr == fixingDates.cend())return;

			Node factor{ 1.0 };
			std::size_t firstExIdx = std::distance(fixingDates.cbegin(), fd_itr);
			for (std::size_t t = firstExIdx; t < fixingDates.size(); ++t) {
				auto nextExItr = dividendData.find(fixingDates[t]);
				factor = 1.0;
				if (nextExItr != last) {
					factor *= (factor - nextExItr->second);
				}
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice.at(fixingDates[t - 1], l), deltaTime);
					lattice.at(fixingDates[t], l) = factor * std::get<0>(tuple);
					lattice.at(fixingDates[t], l + 1) = factor * std::get<1>(tuple);
					lattice.at(fixingDates[t], l + 2) = factor * std::get<2>(tuple);
				}
			}
		}
	}



	template<typename Node,typename DeltaTime>
	void forward_induction(IndexedLattice<LatticeType::Binomial,Node>& lattice,
		LeafForwardGenerator<Node, Node,Node> const &generator, Node apex, DeltaTime deltaTime) {
		_forward_induction_indexed_binomial_impl(lattice, generator, apex, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node,typename DeltaTime>
		void forward_induction(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafForwardGenerator<Node, Node,Node,Node> const &generator, Node apex, DeltaTime deltaTime) {
		_forward_induction_indexed_trinomial_impl(lattice, generator, apex, deltaTime, std::is_compound<DeltaTime>());
	}


	template<typename Node,typename TimeAxis,typename DeltaTime>
	void forward_induction(Lattice<LatticeType::Binomial,Node,TimeAxis>& lattice,
		LeafForwardGenerator<Node, Node,Node> const &generator, Node apex, DeltaTime deltaTime) {
		_forward_induction_lattice_binomial_impl(lattice, generator, apex, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node,typename TimeAxis,typename DeltaTime>
		void forward_induction(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node,Node,Node> const &generator, Node apex, DeltaTime deltaTime) {
			_forward_induction_lattice_trinomial_impl(lattice, generator, apex, deltaTime, std::is_compound<DeltaTime>());
	}




	

	// ==============================================================================
	// ==================== Forward Induction Algorithms ============================
	// ==============================================================================

	template<LatticeType Type,
			typename TimeAxis,
			typename DeltaTime,
			typename Node>
			class ForwardInduction{};


	template<typename TimeAxis,
			typename DeltaTime,
			typename Node>
	class ForwardInduction<LatticeType::Binomial,TimeAxis,DeltaTime,Node> {
	public:
		template<typename LatticeObject,typename Generator>
		void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime,Node apex) {
			//LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
			ForwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime, Node>::
				traverse(lattice, std::forward<Generator>(generator), deltaTime, apex);
		}

		template<typename LatticeObject, typename Generator>
		void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
			std::map<TimeAxis,Node> const &dividendData) {
			//LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
			ForwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime, Node>::
				traverse(lattice, std::forward<Generator>(generator), deltaTime, apex, dividendData);
		}


	};

	template<typename TimeAxis,
		typename DeltaTime,
		typename Node>
		class ForwardInduction<LatticeType::Trinomial, TimeAxis, DeltaTime, Node> {
		public:
			template<typename LatticeObject, typename Generator>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex) {
				//LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				ForwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime, Node>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, apex);
			}
			template<typename LatticeObject, typename Generator>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis, Node> const &dividendData) {
				//LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				ForwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime, Node>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, apex, dividendData);
			}
	};





	// =============================================
	// ===== overloads for discrete dividends ======
	// =============================================

	template<typename Node,typename DeltaTime>
	void forward_induction(IndexedLattice<LatticeType::Binomial, Node>& lattice,
		LeafForwardGenerator<Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime, 
		std::map<std::size_t, Node> const &dividendData) {
		_forward_induction_indexed_binomial_impl(lattice, generator, apex, deltaTime,dividendData,std::is_compound<DeltaTime>());
	}

	template<typename Node, typename DeltaTime>
	void forward_induction(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
		LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,
		std::map<std::size_t, Node> const &dividendData) {
		_forward_induction_indexed_trinomial_impl(lattice, generator, apex, deltaTime, dividendData, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename TimeAxis, typename DeltaTime>
	void forward_induction(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
		LeafForwardGenerator<Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,
		std::map<TimeAxis,Node> const &dividendData) {
		_forward_induction_lattice_binomial_impl(lattice, generator, apex, deltaTime,dividendData,std::is_compound<DeltaTime>());
	}

	template<typename Node, typename TimeAxis, typename DeltaTime>
	void forward_induction(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
		LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,
		std::map<TimeAxis,Node> const &dividendData) {
		_forward_induction_lattice_trinomial_impl(lattice, generator, apex, deltaTime, dividendData, std::is_compound<DeltaTime>());
	}



	// ==============================================================================
	// =================== Backward Induction Algorithms ============================
	// ==============================================================================
	namespace {

		//==== Backward Induction with rewritable source lattice impls ==== 

		template<typename Node,typename DeltaTime>
		void _backward_induction_indexed_binomial_impl(IndexedLattice<LatticeType::Binomial, Node>& lattice,
			LeafBackwardGenerator<Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::true_type) {

			assert(deltaTime.size() == lattice.maxIndex());
			std::size_t lastIdx = lattice.maxIndex();
			for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
				lattice.at(lastIdx, i) = payoff(lattice.at(lastIdx, i));
			}

			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
					lattice.at(n, i) = backwardGenerator(lattice.at(n + 1, i), lattice.at(n + 1, i + 1), deltaTime[n]);
				}
			}
			lattice.at(0, 0) = backwardGenerator(lattice.at(1, 0), lattice.at(1, 1), deltaTime[0]);
		}

		template<typename Node,typename DeltaTime>
		void _backward_induction_indexed_binomial_impl(IndexedLattice<LatticeType::Binomial, Node>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::false_type) {

			std::size_t lastIdx = lattice.maxIndex();
			for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
				lattice.at(lastIdx, i) = payoff(lattice.at(lastIdx, i));
			}

			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
					lattice.at(n, i) = backwardGenerator(lattice.at(n + 1, i), lattice.at(n + 1, i + 1), deltaTime);
				}
			}
			lattice.at(0, 0) = backwardGenerator(lattice.at(1, 0), lattice.at(1, 1), deltaTime);
		}

		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_binomial_adj_impl(IndexedLattice<LatticeType::Binomial, Node>& lattice,
			LeafBackwardGenerator<Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&, Node>const &payoffAdjuster,DeltaTime deltaTime, std::true_type) {

			assert(deltaTime.size() == lattice.maxIndex());
			std::size_t lastIdx = lattice.maxIndex();
			for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
				lattice.at(lastIdx, i) = payoff(lattice.at(lastIdx, i));
			}
			Node value{};
			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
					value = backwardGenerator(lattice.at(n + 1, i), lattice.at(n + 1, i + 1), deltaTime[n]);
					payoffAdjuster(value, lattice.at(n, i));
					lattice.at(n, i) = value;
				}
			}
			value = backwardGenerator(lattice.at(1, 0), lattice.at(1, 1), deltaTime[0]);
			payoffAdjuster(value, lattice.at(0, 0));
			lattice.at(0, 0) = value;
		}

		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_binomial_adj_impl(IndexedLattice<LatticeType::Binomial, Node>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&,Node> const &payoffAdjuster,DeltaTime deltaTime, std::false_type) {

			std::size_t lastIdx = lattice.maxIndex();
			for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
				lattice.at(lastIdx, i) = payoff(lattice.at(lastIdx, i));
			}
			Node value{};
			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
					value = backwardGenerator(lattice.at(n + 1, i), lattice.at(n + 1, i + 1), deltaTime);
					payoffAdjuster(value, lattice.at(n, i));
					lattice.at(n, i) = value;
				}
			}
			value = backwardGenerator(lattice.at(1, 0), lattice.at(1, 1), deltaTime);
			payoffAdjuster(value, lattice.at(0, 0));
			lattice.at(0, 0) = value;
		}


		template<typename Node, typename TimeAxis,typename DeltaTime>
		void _backward_induction_lattice_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::true_type) {

			auto fixingDates = lattice.fixingDates();
			assert(deltaTime.size() == (fixingDates.size() - 1));
			TimeAxis lastDate = fixingDates.back();
			for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
				lattice.at(lastDate, i) = payoff(lattice.at(lastDate, i));
			}

			for (auto n = fixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
					lattice.at(fixingDates[n], i) = backwardGenerator(lattice.at(fixingDates[n + 1], i),
																	lattice.at(fixingDates[n + 1], i + 1),
																	deltaTime[n]);
				}
			}
			lattice.at(fixingDates[0], 0) = backwardGenerator(lattice.at(fixingDates[1], 0), lattice.at(fixingDates[1], 1), deltaTime[0]);
		}

		template<typename Node, typename TimeAxis,typename DeltaTime>
		void _backward_induction_lattice_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
			LeafBackwardGenerator<Node, Node, Node,Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::false_type) {

			auto fixingDates = lattice.fixingDates();
			TimeAxis lastDate = fixingDates.back();
			for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
				lattice.at(lastDate, i) = payoff(lattice.at(lastDate, i));
			}

			for (auto n = fixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
					lattice.at(fixingDates[n], i) = backwardGenerator(lattice.at(fixingDates[n + 1], i),
																	lattice.at(fixingDates[n + 1], i + 1),
																	deltaTime);
				}
			}
			lattice.at(fixingDates[0], 0) = backwardGenerator(lattice.at(fixingDates[1], 0), lattice.at(fixingDates[1], 1), deltaTime);
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_binomial_adj_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&,Node> const &payoffAdjuster,DeltaTime deltaTime, std::true_type) {

			auto fixingDates = lattice.fixingDates();
			assert(deltaTime.size() == (fixingDates.size() - 1));
			TimeAxis lastDate = fixingDates.back();
			for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
				lattice.at(lastDate, i) = payoff(lattice.at(lastDate, i));
			}
			Node value{};
			for (auto n = fixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
					value = backwardGenerator(lattice.at(fixingDates[n + 1], i),
											lattice.at(fixingDates[n + 1], i + 1),
											deltaTime[n]);
					payoffAdjuster(value, lattice.at(fixingDates[n], i));
					lattice.at(fixingDates[n], i) = value;
				}
			}
			value = backwardGenerator(lattice.at(fixingDates[1], 0), lattice.at(fixingDates[1], 1), deltaTime[0]);
			payoffAdjuster(value, lattice.at(fixingDates[0], 0));
			lattice.at(fixingDates[0], 0) = value;
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_binomial_adj_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&, Node> const &payoffAdjuster,DeltaTime deltaTime, std::false_type) {

			auto fixingDates = lattice.fixingDates();
			TimeAxis lastDate = fixingDates.back();
			for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
				lattice.at(lastDate, i) = payoff(lattice.at(lastDate, i));
			}
			Node value{};
			for (auto n = fixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
					value = backwardGenerator(lattice.at(fixingDates[n + 1], i),
											lattice.at(fixingDates[n + 1], i + 1),
											deltaTime);
					payoffAdjuster(value,lattice.at(fixingDates[n], i));
					lattice.at(fixingDates[n], i) = value;
				}
			}
			value = backwardGenerator(lattice.at(fixingDates[1], 0), lattice.at(fixingDates[1], 1), deltaTime);
			payoffAdjuster(value, lattice.at(fixingDates[0], 0));
			lattice(fixingDates[0], 0) = value;
		}


		template<typename Node,typename DeltaTime>
		void _backward_induction_indexed_trinomial_impl(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::true_type) {

			assert(deltaTime.size() == lattice.maxIndex());
			std::size_t lastIdx = lattice.maxIndex();
			for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
				lattice.at(lastIdx, i) = payoff(lattice.at(lastIdx, i));
			}

			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
					lattice.at(n, i) = backwardGenerator(lattice.at(n + 1, i), lattice.at(n + 1, i + 1), lattice.at(n + 1, i + 2), deltaTime[n]);
				}
			}
			lattice.at(0, 0) = backwardGenerator(lattice.at(1, 0), lattice.at(1, 1), lattice.at(1, 2), deltaTime[0]);
		}

		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_trinomial_impl(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::false_type) {

			std::size_t lastIdx = lattice.maxIndex();
			for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
				lattice.at(lastIdx, i) = payoff(lattice.at(lastIdx, i));
			}

			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
					lattice.at(n, i) = backwardGenerator(lattice.at(n + 1, i), lattice.at(n + 1, i + 1), lattice.at(n + 1, i + 2), deltaTime);
				}
			}
			lattice.at(0, 0) = backwardGenerator(lattice.at(1, 0), lattice.at(1, 1), lattice.at(1, 2), deltaTime);
		}

		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_trinomial_adj_impl(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&,Node> const &payoffAdjuster, DeltaTime deltaTime, std::true_type) {

			assert(deltaTime.size() == lattice.maxIndex());
			std::size_t lastIdx = lattice.maxIndex();
			for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
				lattice.at(lastIdx, i) = payoff(lattice.at(lastIdx, i));
			}
			Node value{};
			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
					value = backwardGenerator(lattice.at(n + 1, i), lattice.at(n + 1, i + 1), lattice.at(n + 1, i + 2), deltaTime[n]);
					payoffAdjuster(value, lattice.at(n, i));
					lattice.at(n, i) = value;
				}
			}
			value = backwardGenerator(lattice.at(1, 0), lattice.at(1, 1), lattice.at(1, 2), deltaTime[0]);
			payoffAdjuster(value, lattice.at(0, 0));
			lattice.at(0, 0) = value;
		}

		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_trinomial_adj_impl(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&,Node> const &payoffAdjuster,DeltaTime deltaTime, std::false_type) {

			std::size_t lastIdx = lattice.maxIndex();
			for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
				lattice.at(lastIdx, i) = payoff(lattice.at(lastIdx, i));
			}
			Node value{};
			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
					value = backwardGenerator(lattice.at(n + 1, i), lattice.at(n + 1, i + 1), lattice.at(n + 1, i + 2), deltaTime);
					payoffAdjuster(value, lattice.at(n, i));
					lattice.at(n, i) = value;
				}
			}
			value = backwardGenerator(lattice.at(1, 0), lattice.at(1, 1), lattice.at(1, 2), deltaTime);
			payoffAdjuster(value, lattice.at(0, 0));
			lattice.at(0, 0) = value;
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff, 
			DeltaTime deltaTime,std::true_type) {

			auto fixingDates = lattice.fixingDates();
			assert(deltaTime.size() == (fixingDates.size() - 1));
			TimeAxis lastDate = fixingDates.back();
			for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
				lattice.at(lastDate, i) = payoff(lattice.at(lastDate, i));
			}

			for (auto n = fixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
					lattice.at(fixingDates[n], i) = backwardGenerator(lattice.at(fixingDates[n + 1], i),
						lattice.at(fixingDates[n + 1], i + 1),
						lattice.at(fixingDates[n + 1], i + 2),
						deltaTime[n]);
				}
			}
			lattice.at(fixingDates[0], 0) = backwardGenerator(lattice.at(fixingDates[1], 0),
				lattice.at(fixingDates[1], 1),
				lattice.at(fixingDates[1], 2),
				deltaTime[0]);
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff, 
			DeltaTime deltaTime,std::false_type) {

			auto fixingDates = lattice.fixingDates();
			TimeAxis lastDate = fixingDates.back();
			for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
				lattice.at(lastDate, i) = payoff(lattice.at(lastDate, i));
			}

			for (auto n = fixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
					lattice.at(fixingDates[n], i) = backwardGenerator(lattice.at(fixingDates[n + 1], i),
						lattice.at(fixingDates[n + 1], i + 1),
						lattice.at(fixingDates[n + 1], i + 2),
						deltaTime);
				}
			}
			lattice.at(fixingDates[0], 0) = backwardGenerator(lattice.at(fixingDates[1], 0),
				lattice.at(fixingDates[1], 1),
				lattice.at(fixingDates[1], 2),
				deltaTime);
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_trinomial_adj_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&,Node> const &payoffAdjuster,DeltaTime deltaTime, std::true_type) {

			auto fixingDates = lattice.fixingDates();
			assert(deltaTime.size() == (fixingDates.size() - 1));
			TimeAxis lastDate = fixingDates.back();
			for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
				lattice.at(lastDate, i) = payoff(lattice.at(lastDate, i));
			}
			Node value{};
			for (auto n = fixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
					value = backwardGenerator(lattice.at(fixingDates[n + 1], i),
						lattice.at(fixingDates[n + 1], i + 1),
						lattice.at(fixingDates[n + 1], i + 2),
						deltaTime[n]);
					payoffAdjuster(value, lattice.at(fixingDates[n], i));
					lattice.at(fixingDates[n], i) = value;
				}
			}
			value = backwardGenerator(lattice.at(fixingDates[1], 0),
				lattice.at(fixingDates[1], 1),
				lattice.at(fixingDates[1], 2),
				deltaTime[0]);
			payoffAdjuster(value, lattice.at(fixingDates[0], 0));
			lattice.at(fixingDates[0], 0) = value;
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_trinomial_adj_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&,Node> const &payoffAdjuster,DeltaTime deltaTime, std::false_type) {

			auto fixingDates = lattice.fixingDates();
			TimeAxis lastDate = fixingDates.back();
			for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
				lattice.at(lastDate, i) = payoff(lattice.at(lastDate, i));
			}
			Node value{};
			for (auto n = fixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
					value = backwardGenerator(lattice.at(fixingDates[n + 1], i),
						lattice.at(fixingDates[n + 1], i + 1),
						lattice.at(fixingDates[n + 1], i + 2),
						deltaTime);
					payoffAdjuster(value, lattice.at(fixingDates[n], i));
					lattice.at(fixingDates[n], i) = value;
				}
			}
			value = backwardGenerator(lattice.at(fixingDates[1], 0),
				lattice.at(fixingDates[1], 1),
				lattice.at(fixingDates[1], 2),
				deltaTime);
			payoffAdjuster(value, lattice.at(fixingDates[0], 0));
			lattice.at(fixingDates[0], 0) = value;
		}


		//==== Backward Induction with source and modified lattice impls ====

		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_binomial_mod_impl(IndexedLattice<LatticeType::Binomial, Node>const &sourceLattice,
			IndexedLattice<LatticeType::Binomial, Node> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::true_type) {
			
			assert(deltaTime.size() == sourceLattice.maxIndex());
			assert(deltaTime.size() == modifiedLattice.maxIndex());

			std::size_t lastIdx = modifiedLattice.maxIndex();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastIdx).size(); ++i) {
				modifiedLattice.at(lastIdx, i) = payoff(modifiedLattice.at(lastIdx, i));
			}

			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(n).size(); ++i) {
					modifiedLattice.at(n, i) = backwardGenerator(modifiedLattice.at(n + 1, i), modifiedLattice.at(n + 1, i + 1), deltaTime[n]);
				}
			}
			modifiedLattice.at(0, 0) = backwardGenerator(modifiedLattice.at(1, 0), modifiedLattice.at(1, 1), deltaTime[0]);
		}

		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_binomial_mod_impl(IndexedLattice<LatticeType::Binomial, Node> const &sourceLattice,
			IndexedLattice<LatticeType::Binomial, Node> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::false_type) {

			assert(sourceLattice.maxIndex() == modifiedLattice.maxIndex());
			std::size_t lastIdx = modifiedLattice.maxIndex();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastIdx).size(); ++i) {
				modifiedLattice.at(lastIdx, i) = payoff(modifiedLattice.at(lastIdx, i));
			}

			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(n).size(); ++i) {
					modifiedLattice.at(n, i) = backwardGenerator(modifiedLattice.at(n + 1, i), modifiedLattice.at(n + 1, i + 1), deltaTime);
				}
			}
			modifiedLattice.at(0, 0) = backwardGenerator(modifiedLattice.at(1, 0), modifiedLattice.at(1, 1), deltaTime);
		}

		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_binomial_mod_adj_impl(IndexedLattice<LatticeType::Binomial, Node>const& sourceLattice,
			IndexedLattice<LatticeType::Binomial, Node> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&, Node>const &payoffAdjuster, DeltaTime deltaTime, std::true_type) {

			assert(deltaTime.size() == sourceLattice.maxIndex());
			assert(deltaTime.size() == modifiedLattice.maxIndex());

			std::size_t lastIdx = modifiedLattice.maxIndex();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastIdx).size(); ++i) {
				modifiedLattice.at(lastIdx, i) = payoff(modifiedLattice.at(lastIdx, i));
			}
			Node value{};
			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(n).size(); ++i) {
					value = backwardGenerator(modifiedLattice.at(n + 1, i), modifiedLattice.at(n + 1, i + 1), deltaTime[n]);
					payoffAdjuster(value, sourceLattice.at(n, i));
					modifiedLattice.at(n, i) = value;
				}
			}
			value = backwardGenerator(modifiedLattice.at(1, 0), modifiedLattice.at(1, 1), deltaTime[0]);
			payoffAdjuster(value, sourceLattice.at(0, 0));
			modifiedLattice.at(0, 0) = value;
		}

		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_binomial_mod_adj_impl(IndexedLattice<LatticeType::Binomial, Node>const& sourceLattice,
			IndexedLattice<LatticeType::Binomial, Node> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&, Node> const &payoffAdjuster, DeltaTime deltaTime, std::false_type) {

			assert(sourceLattice.maxIndex() == modifiedLattice.maxIndex());
			std::size_t lastIdx = modifiedLattice.maxIndex();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastIdx).size(); ++i) {
				modifiedLattice.at(lastIdx, i) = payoff(modifiedLattice.at(lastIdx, i));
			}
			Node value{};
			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(n).size(); ++i) {
					value = backwardGenerator(modifiedLattice.at(n + 1, i), modifiedLattice.at(n + 1, i + 1), deltaTime);
					payoffAdjuster(value, sourceLattice.at(n, i));
					modifiedLattice.at(n, i) = value;
				}
			}
			value = backwardGenerator(modifiedLattice.at(1, 0), modifiedLattice.at(1, 1), deltaTime);
			payoffAdjuster(value, sourceLattice.at(0, 0));
			modifiedLattice.at(0, 0) = value;
		}


		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_binomial_mod_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>const& sourceLattice,
			Lattice<LatticeType::Binomial, Node, TimeAxis> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::true_type) {

			auto sFixingDates = sourceLattice.fixingDates();
			auto mFixingDates = modifiedLattice.fixingDates();
			assert(sFixingDates == mFixingDates);
			assert(deltaTime.size() == (sFixingDates.size() - 1));


			TimeAxis lastDate = mFixingDates.back();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastDate).size(); ++i) {
				modifiedLattice.at(lastDate, i) = payoff(modifiedLattice.at(lastDate, i));
			}

			for (auto n = mFixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(mFixingDates[n]).size(); ++i) {
					modifiedLattice.at(mFixingDates[n], i) = backwardGenerator(modifiedLattice.at(mFixingDates[n + 1], i),
																			modifiedLattice.at(mFixingDates[n + 1], i + 1),
																			deltaTime[n]);
				}
			}
			modifiedLattice.at(mFixingDates[0], 0) = backwardGenerator(modifiedLattice.at(mFixingDates[1], 0),
																	modifiedLattice.at(mFixingDates[1], 1),
																	deltaTime[0]);
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_binomial_mod_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>const &sourceLattice,
			Lattice<LatticeType::Binomial, Node,TimeAxis> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::false_type) {

			auto sFixingDates = sourceLattice.fixingDates();
			auto mFixingDates = modifiedLattice.fixingDates();
			assert(sFixingDates == mFixingDates);

			TimeAxis lastDate = mFixingDates.back();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastDate).size(); ++i) {
				modifiedLattice.at(lastDate, i) = payoff(modifiedLattice.at(lastDate, i));
			}

			for (auto n = mFixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(mFixingDates[n]).size(); ++i) {
					modifiedLattice.at(mFixingDates[n], i) = backwardGenerator(modifiedLattice.at(mFixingDates[n + 1], i),
																			modifiedLattice.at(mFixingDates[n + 1], i + 1),
																			deltaTime);
				}
			}
			modifiedLattice.at(mFixingDates[0], 0) = backwardGenerator(modifiedLattice.at(mFixingDates[1], 0),
															modifiedLattice.at(mFixingDates[1], 1), deltaTime);
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_binomial_mod_adj_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>const &sourceLattice,
			Lattice<LatticeType::Binomial, Node,TimeAxis> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&, Node> const &payoffAdjuster, DeltaTime deltaTime, std::true_type) {

			auto sFixingDates = sourceLattice.fixingDates();
			auto mFixingDates = modifiedLattice.fixingDates();

			assert(sFixingDates == mFixingDates);
			assert(deltaTime.size() == (mFixingDates.size() - 1));

			TimeAxis lastDate = mFixingDates.back();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastDate).size(); ++i) {
				modifiedLattice.at(lastDate, i) = payoff(modifiedLattice.at(lastDate, i));
			}
			Node value{};
			for (auto n = mFixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(mFixingDates[n]).size(); ++i) {
					value = backwardGenerator(modifiedLattice.at(mFixingDates[n + 1], i),
											modifiedLattice.at(mFixingDates[n + 1], i + 1),
											deltaTime[n]);
					payoffAdjuster(value, sourceLattice.at(sFixingDates[n], i)   );
					modifiedLattice.at(mFixingDates[n], i) = value;
				}
			}
			value = backwardGenerator(modifiedLattice.at(mFixingDates[1], 0),
										modifiedLattice.at(mFixingDates[1], 1),
										deltaTime[0]);
			payoffAdjuster(value, sourceLattice.at(sFixingDates[0], 0));
			modifiedLattice.at(mFixingDates[0], 0) = value;
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_binomial_mod_adj_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>const &sourceLattice,
			Lattice<LatticeType::Binomial, Node, TimeAxis> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&, Node> const &payoffAdjuster, DeltaTime deltaTime, std::false_type) {

			auto sFixingDates = sourceLattice.fixingDates();
			auto mFixingDates = modifiedLattice.fixingDates();
			assert(sFixingDates == mFixingDates);

			TimeAxis lastDate = mFixingDates.back();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastDate).size(); ++i) {
				modifiedLattice.at(lastDate, i) = payoff(modifiedLattice.at(lastDate, i));
			}
			Node value{};
			for (auto n = mFixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(mFixingDates[n]).size(); ++i) {
					value = backwardGenerator(modifiedLattice.at(mFixingDates[n + 1], i),
											modifiedLattice.at(mFixingDates[n + 1], i + 1),
											deltaTime);
					payoffAdjuster(value, sourceLattice.at(mFixingDates[n], i));
					modifiedLattice.at(mFixingDates[n], i) = value;
				}
			}
			value = backwardGenerator(modifiedLattice.at(mFixingDates[1], 0), modifiedLattice.at(mFixingDates[1], 1), deltaTime);
			payoffAdjuster(value, sourceLattice.at(mFixingDates[0], 0));
			modifiedLattice.at(mFixingDates[0], 0) = value;
		}


		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_trinomial_mod_impl(IndexedLattice<LatticeType::Trinomial, Node>const &sourceLattice,
			IndexedLattice<LatticeType::Trinomial, Node> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::true_type) {

			assert(deltaTime.size() == sourceLattice.maxIndex());
			assert(deltaTime.size() == modifiedLattice.maxIndex());

			std::size_t lastIdx = modifiedLattice.maxIndex();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastIdx).size(); ++i) {
				modifiedLattice.at(lastIdx, i) = payoff(modifiedLattice.at(lastIdx, i));
			}

			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(n).size(); ++i) {
					modifiedLattice.at(n, i) = backwardGenerator(modifiedLattice.at(n + 1, i),
															modifiedLattice.at(n + 1, i + 1),
															modifiedLattice.at(n + 1, i + 2), deltaTime[n]);
				}
			}
			modifiedLattice.at(0, 0) = backwardGenerator(modifiedLattice.at(1, 0),
													modifiedLattice.at(1, 1),
													modifiedLattice.at(1, 2), deltaTime[0]);
		}

		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_trinomial_mod_impl(IndexedLattice<LatticeType::Trinomial, Node>const &sourceLattice,
			IndexedLattice<LatticeType::Trinomial, Node> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::false_type) {

			assert(sourceLattice.maxIndex() == modifiedLattice.maxIndex());
			std::size_t lastIdx = modifiedLattice.maxIndex();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastIdx).size(); ++i) {
				modifiedLattice.at(lastIdx, i) = payoff(modifiedLattice.at(lastIdx, i));
			}

			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(n).size(); ++i) {
					modifiedLattice.at(n, i) = backwardGenerator(modifiedLattice.at(n + 1, i),
																modifiedLattice.at(n + 1, i + 1),
																modifiedLattice.at(n + 1, i + 2), deltaTime);
				}
			}
			modifiedLattice.at(0, 0) = backwardGenerator(modifiedLattice.at(1, 0),
												modifiedLattice.at(1, 1),
												modifiedLattice.at(1, 2), deltaTime);
		}

		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_trinomial_mod_adj_impl(IndexedLattice<LatticeType::Trinomial, Node>const &sourceLattice,
			IndexedLattice<LatticeType::Trinomial, Node> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&, Node> const &payoffAdjuster, DeltaTime deltaTime, std::true_type) {

			assert(deltaTime.size() == sourceLattice.maxIndex());
			assert(deltaTime.size() == modifiedLattice.maxIndex());

			std::size_t lastIdx = modifiedLattice.maxIndex();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastIdx).size(); ++i) {
				modifiedLattice.at(lastIdx, i) = payoff(modifiedLattice.at(lastIdx, i));
			}
			Node value{};
			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(n).size(); ++i) {
					value = backwardGenerator(modifiedLattice.at(n + 1, i),
												modifiedLattice.at(n + 1, i + 1),
												modifiedLattice.at(n + 1, i + 2), deltaTime[n]);
					payoffAdjuster(value, sourceLattice.at(n, i));
					modifiedLattice.at(n, i) = value;
				}
			}
			value = backwardGenerator(modifiedLattice.at(1, 0), modifiedLattice.at(1, 1), modifiedLattice.at(1, 2), deltaTime[0]);
			payoffAdjuster(value, sourceLattice.at(0, 0));
			modifiedLattice.at(0, 0) = value;
		}

		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_trinomial_mod_adj_impl(IndexedLattice<LatticeType::Trinomial, Node>const &sourceLattice,
			IndexedLattice<LatticeType::Trinomial, Node> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&, Node> const &payoffAdjuster, DeltaTime deltaTime, std::false_type) {

			assert(sourceLattice.maxIndex() == modifiedLattice.maxIndex());

			std::size_t lastIdx = modifiedLattice.maxIndex();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastIdx).size(); ++i) {
				modifiedLattice.at(lastIdx, i) = payoff(modifiedLattice.at(lastIdx, i));
			}
			Node value{};
			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(n).size(); ++i) {
					value = backwardGenerator(modifiedLattice.at(n + 1, i),
											modifiedLattice.at(n + 1, i + 1),
											modifiedLattice.at(n + 1, i + 2), deltaTime);
					payoffAdjuster(value, sourceLattice.at(n, i));
					modifiedLattice.at(n, i) = value;
				}
			}
			value = backwardGenerator(modifiedLattice.at(1, 0), modifiedLattice.at(1, 1), modifiedLattice.at(1, 2), deltaTime);
			payoffAdjuster(value, sourceLattice.at(0, 0));
			modifiedLattice.at(0, 0) = value;
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_trinomial_mod_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>const &sourceLattice,
			Lattice<LatticeType::Trinomial, Node, TimeAxis> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::true_type) {

			auto sFixingDates = sourceLattice.fixingDates();
			auto mFixingDates = modifiedLattice.fixingDates();

			assert(sFixingDates == mFixingDates);
			assert(deltaTime.size() == (sFixingDates.size() - 1));

			TimeAxis lastDate = mFixingDates.back();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastDate).size(); ++i) {
				modifiedLattice.at(lastDate, i) = payoff(modifiedLattice.at(lastDate, i));
			}
			for (auto n = sFixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(mFixingDates[n]).size(); ++i) {
					modifiedLattice.at(mFixingDates[n], i) = backwardGenerator(modifiedLattice.at(mFixingDates[n + 1], i),
																			modifiedLattice.at(mFixingDates[n + 1], i + 1),
																			modifiedLattice.at(mFixingDates[n + 1], i + 2),
																			deltaTime[n]);
				}
			}
			modifiedLattice.at(mFixingDates[0], 0) = backwardGenerator(modifiedLattice.at(mFixingDates[1], 0),
																	modifiedLattice.at(mFixingDates[1], 1),
																	modifiedLattice.at(mFixingDates[1], 2),
																	deltaTime[0]);
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_trinomial_mod_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>const &sourceLattice,
			Lattice<LatticeType::Trinomial, Node, TimeAxis> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::false_type) {

			auto sFixingDates = sourceLattice.fixingDates();
			auto mFixingDates = modifiedLattice.fixingDates();
			assert(sFixingDates == mFixingDates);

			TimeAxis lastDate = mFixingDates.back();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastDate).size(); ++i) {
				modifiedLattice.at(lastDate, i) = payoff(modifiedLattice.at(lastDate, i));
			}

			for (auto n = sFixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(mFixingDates[n]).size(); ++i) {
					modifiedLattice.at(mFixingDates[n], i) = backwardGenerator(modifiedLattice.at(mFixingDates[n + 1], i),
																			modifiedLattice.at(mFixingDates[n + 1], i + 1),
																			modifiedLattice.at(mFixingDates[n + 1], i + 2),
																			deltaTime);
				}
			}
			modifiedLattice.at(mFixingDates[0], 0) = backwardGenerator(modifiedLattice.at(mFixingDates[1], 0),
																	modifiedLattice.at(mFixingDates[1], 1),
																	modifiedLattice.at(mFixingDates[1], 2),
																	deltaTime);
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_trinomial_mod_adj_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>const& sourceLattice,
			Lattice<LatticeType::Trinomial, Node, TimeAxis> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&, Node> const &payoffAdjuster, DeltaTime deltaTime, std::true_type) {

			auto sFixingDates = sourceLattice.fixingDates();
			auto mFixingDates = modifiedLattice.fixingDates();

			assert(sFixingDates == mFixingDates);
			assert(deltaTime.size() == (sFixingDates.size() - 1));

			TimeAxis lastDate = mFixingDates.back();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastDate).size(); ++i) {
				modifiedLattice.at(lastDate, i) = payoff(modifiedLattice.at(lastDate, i));
			}
			Node value{};
			for (auto n = mFixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(mFixingDates[n]).size(); ++i) {
					value = backwardGenerator(modifiedLattice.at(mFixingDates[n + 1], i),
												modifiedLattice.at(mFixingDates[n + 1], i + 1),
												modifiedLattice.at(mFixingDates[n + 1], i + 2),
												deltaTime[n]);
					payoffAdjuster(value, sourceLattice.at(mFixingDates[n], i));
					modifiedLattice.at(mFixingDates[n], i) = value;
				}
			}
			value = backwardGenerator(modifiedLattice.at(mFixingDates[1], 0),
										modifiedLattice.at(mFixingDates[1], 1),
										modifiedLattice.at(mFixingDates[1], 2),
										deltaTime[0]);
			payoffAdjuster(value, sourceLattice.at(mFixingDates[0], 0));
			modifiedLattice.at(mFixingDates[0], 0) = value;
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_trinomial_mod_adj_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>const& sourceLattice,
			Lattice<LatticeType::Trinomial, Node,TimeAxis> &modifiedLattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			PayoffAdjuster<Node&, Node> const &payoffAdjuster, DeltaTime deltaTime, std::false_type) {

			auto sFixingDates = sourceLattice.fixingDates();
			auto mFixingDates = modifiedLattice.fixingDates();
			assert(sFixingDates == mFixingDates);

			TimeAxis lastDate = mFixingDates.back();
			for (auto i = 0; i < modifiedLattice.nodesAt(lastDate).size(); ++i) {
				modifiedLattice.at(lastDate, i) = payoff(modifiedLattice.at(lastDate, i));
			}
			Node value{};
			for (auto n = mFixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < modifiedLattice.nodesAt(mFixingDates[n]).size(); ++i) {
					value = backwardGenerator(modifiedLattice.at(mFixingDates[n + 1], i),
											modifiedLattice.at(mFixingDates[n + 1], i + 1),
											modifiedLattice.at(mFixingDates[n + 1], i + 2),
											deltaTime);
					payoffAdjuster(value, sourceLattice.at(mFixingDates[n], i));
					modifiedLattice.at(mFixingDates[n], i) = value;
				}
			}
			value = backwardGenerator(modifiedLattice.at(mFixingDates[1], 0),
										modifiedLattice.at(mFixingDates[1], 1),
										modifiedLattice.at(mFixingDates[1], 2),
										deltaTime);
			payoffAdjuster(value, sourceLattice.at(mFixingDates[0], 0));
			modifiedLattice.at(mFixingDates[0], 0) = value;
		}

	}
	
	//==== Backward Induction with rewritable source lattice ====  

	template<typename Node,typename DeltaTime>
	void backward_induction(IndexedLattice<LatticeType::Binomial, Node>& lattice,
		LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node,Node> const &payoff,DeltaTime deltaTime) {
		_backward_induction_indexed_binomial_impl(lattice, backwardGenerator, payoff, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename DeltaTime>
	void backward_induction(IndexedLattice<LatticeType::Binomial, Node>& lattice,
		LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff, 
		PayoffAdjuster<Node&,Node> const &payoffAdjuster, DeltaTime deltaTime) {
		_backward_induction_indexed_binomial_adj_impl(lattice, backwardGenerator, payoff, payoffAdjuster, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename TimeAxis,typename DeltaTime>
	void backward_induction(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
		LeafBackwardGenerator<Node, Node, Node,Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,DeltaTime deltaTime) {
		_backward_induction_lattice_binomial_impl(lattice, backwardGenerator, payoff, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename TimeAxis, typename DeltaTime>
	void backward_induction(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
		LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
		PayoffAdjuster<Node&,Node> const &payoffAdjuster,DeltaTime deltaTime) {
		_backward_induction_lattice_binomial_adj_impl(lattice, backwardGenerator, payoff, payoffAdjuster, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node,typename DeltaTime>
	void backward_induction(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
		LeafBackwardGenerator<Node,Node, Node, Node,Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,DeltaTime deltaTime) {
		_backward_induction_indexed_trinomial_impl(lattice, backwardGenerator, payoff, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename DeltaTime>
	void backward_induction(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
		LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
		PayoffAdjuster<Node&,Node> const &payoffAdjuster,DeltaTime deltaTime) {
		_backward_induction_indexed_trinomial_adj_impl(lattice, backwardGenerator, payoff, payoffAdjuster, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename TimeAxis,typename DeltaTime>
	void backward_induction(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
		LeafBackwardGenerator<Node, Node,Node, Node,Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,DeltaTime deltaTime) {
		_backward_induction_lattice_trinomial_impl(lattice, backwardGenerator, payoff, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename TimeAxis, typename DeltaTime>
	void backward_induction(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
		LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
		PayoffAdjuster<Node&,Node> const &payoffAdjuster,DeltaTime deltaTime) {
		_backward_induction_lattice_trinomial_adj_impl(lattice, backwardGenerator, payoff, payoffAdjuster, deltaTime, std::is_compound<DeltaTime>());
	}

	//==== Backward Induction with source and modified lattice ====


	template<typename Node, typename DeltaTime>
	void backward_induction(IndexedLattice<LatticeType::Binomial, Node> const &sourceLattice,
		IndexedLattice<LatticeType::Binomial, Node> &modifiedLattice,
		LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff, DeltaTime deltaTime) {
		_backward_induction_indexed_binomial_mod_impl(sourceLattice,modifiedLattice, backwardGenerator, payoff, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename DeltaTime>
	void backward_induction(IndexedLattice<LatticeType::Binomial, Node>const &sourceLattice,
		IndexedLattice<LatticeType::Binomial,Node> &modifiedLattice,
		LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
		PayoffAdjuster<Node&, Node> const &payoffAdjuster, DeltaTime deltaTime) {
		_backward_induction_indexed_binomial_mod_adj_impl(sourceLattice,modifiedLattice, backwardGenerator, payoff, payoffAdjuster, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename TimeAxis, typename DeltaTime>
	void backward_induction(Lattice<LatticeType::Binomial, Node, TimeAxis>const &sourceLattice,
		Lattice<LatticeType::Binomial, Node, TimeAxis> &modifiedLattice,
		LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff, DeltaTime deltaTime) {
		_backward_induction_lattice_binomial_mod_impl(sourceLattice,modifiedLattice, backwardGenerator, payoff, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename TimeAxis, typename DeltaTime>
	void backward_induction(Lattice<LatticeType::Binomial, Node, TimeAxis> const &sourceLattice,
		Lattice<LatticeType::Binomial, Node,TimeAxis> &modifiedLattice,
		LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
		PayoffAdjuster<Node&, Node> const &payoffAdjuster, DeltaTime deltaTime) {
		_backward_induction_lattice_binomial_mod_adj_impl(sourceLattice,modifiedLattice, backwardGenerator, payoff, payoffAdjuster, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename DeltaTime>
	void backward_induction(IndexedLattice<LatticeType::Trinomial, Node>const &sourceLattice,
		IndexedLattice<LatticeType::Trinomial, Node> &modifiedLattice,
		LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff, DeltaTime deltaTime) {
		_backward_induction_indexed_trinomial_mod_impl(sourceLattice,modifiedLattice, backwardGenerator, payoff, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename DeltaTime>
	void backward_induction(IndexedLattice<LatticeType::Trinomial, Node>const &sourceLattice,
		IndexedLattice<LatticeType::Trinomial, Node> &modifiedLattice,
		LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
		PayoffAdjuster<Node&, Node> const &payoffAdjuster, DeltaTime deltaTime) {
		_backward_induction_indexed_trinomial_mod_adj_impl(sourceLattice,modifiedLattice, backwardGenerator, payoff, payoffAdjuster, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename TimeAxis, typename DeltaTime>
	void backward_induction(Lattice<LatticeType::Trinomial, Node, TimeAxis>const &sourceLattice,
		Lattice<LatticeType::Trinomial,Node, TimeAxis> &modifiedLattice,
		LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff, DeltaTime deltaTime) {
		_backward_induction_lattice_trinomial_mod_impl(sourceLattice,modifiedLattice, backwardGenerator, payoff, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename TimeAxis, typename DeltaTime>
	void backward_induction(Lattice<LatticeType::Trinomial, Node, TimeAxis>const &sourceLattice,
		Lattice<LatticeType::Trinomial, Node,TimeAxis> &modifiedLattice,
		LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
		PayoffAdjuster<Node&, Node> const &payoffAdjuster, DeltaTime deltaTime) {
		_backward_induction_lattice_trinomial_mod_adj_impl(sourceLattice,modifiedLattice, backwardGenerator, payoff, payoffAdjuster, deltaTime, std::is_compound<DeltaTime>());
	}



	// ==============================================================================
	// ==================== Backward Induction Algorithms ===========================
	// ==============================================================================

	
	template<LatticeType Type,
		typename TimeAxis,
		typename DeltaTime>
		class BackwardInduction {};


	template<typename TimeAxis,
		typename DeltaTime>
		class BackwardInduction<LatticeType::Binomial, TimeAxis, DeltaTime> {
		public:
			template<typename LatticeObject, typename Generator, typename Payoff>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff) {
				//LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				BackwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff));
			}

			template<typename LatticeObject, typename Generator, typename Payoff,typename PayoffAdjuster>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {
				//LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				BackwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff),std::forward<PayoffAdjuster>(payoffAdjuster));
			}

	};

	template<typename TimeAxis,
		typename DeltaTime>
		class BackwardInduction<LatticeType::Trinomial, TimeAxis, DeltaTime> {
		public:
			template<typename LatticeObject, typename Generator, typename Payoff>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff) {
				//LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				BackwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff));
			}

			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {
				//LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				BackwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff), std::forward<PayoffAdjuster>(payoffAdjuster));
			}
	};






	//==================== Merging algorithms ============================

	namespace {

		template<typename Node>
		IndexedLattice<LatticeType::Binomial, std::tuple<Node, Node>>
			_s_merge_binomial_indexed_impl(IndexedLattice<LatticeType::Binomial, Node>const &latticeA, IndexedLattice<LatticeType::Binomial, Node> const&latticeB) {
			assert(latticeA.timeDimension() == latticeB.timeDimension());
			IndexedLattice<LatticeType::Binomial, std::tuple<Node, Node>> result{ latticeA.maxIndex() };
			for (auto i = result.minIndex(); i <= result.maxIndex(); ++i) {
				for (auto j = 0; j < result.nodesAt(i).size(); ++j) {
					result(i, j) = std::make_tuple(latticeA(i, j), latticeB(i, j));
				}
			}
			return result;
		}

		template<typename Node>
		IndexedLattice<LatticeType::Binomial, std::tuple<Node, Node>>
			_p_merge_binomial_indexed_impl(IndexedLattice<LatticeType::Binomial, Node>const &latticeA, IndexedLattice<LatticeType::Binomial, Node> const&latticeB) {
			assert(latticeA.timeDimension() == latticeB.timeDimension());
			IndexedLattice<LatticeType::Binomial, std::tuple<Node, Node>> result{ latticeA.maxIndex() };
			concurrency::parallel_for(static_cast<std::size_t>(result.minIndex()), static_cast<std::size_t>(result.maxIndex() + 1), [&](std::size_t t) {
				for (auto j = 0; j < result.nodesAt(t).size(); ++j) {
					result(t, j) = std::make_tuple(latticeA(t, j), latticeB(t, j));
				}
			});
			return result;
		}

		template<typename Node>
		IndexedLattice<LatticeType::Trinomial, std::tuple<Node, Node>>
			_s_merge_trinomial_indexed_impl(IndexedLattice<LatticeType::Trinomial, Node>const &latticeA, IndexedLattice<LatticeType::Trinomial, Node> const&latticeB) {
			assert(latticeA.timeDimension() == latticeB.timeDimension());
			IndexedLattice<LatticeType::Trinomial, std::tuple<Node, Node>> result{ latticeA.maxIndex() };
			for (auto i = result.minIndex(); i <= result.maxIndex(); ++i) {
				for (auto j = 0; j < result.nodesAt(i).size(); ++j) {
					result(i, j) = std::make_tuple(latticeA(i, j), latticeB(i, j));
				}
			}
			return result;
		}

		template<typename Node>
		IndexedLattice<LatticeType::Trinomial, std::tuple<Node, Node>>
			_p_merge_trinomial_indexed_impl(IndexedLattice<LatticeType::Trinomial, Node>const &latticeA, IndexedLattice<LatticeType::Trinomial, Node> const&latticeB) {
			assert(latticeA.timeDimension() == latticeB.timeDimension());
			IndexedLattice<LatticeType::Trinomial, std::tuple<Node, Node>> result{ latticeA.maxIndex() };
			concurrency::parallel_for(static_cast<std::size_t>(result.minIndex()), static_cast<std::size_t>(result.maxIndex() + 1), [&](std::size_t t) {
				for (auto j = 0; j < result.nodesAt(t).size(); ++j) {
					result(t, j) = std::make_tuple(latticeA(t, j), latticeB(t, j));
				}
			});
			return result;
		}

		template<typename Node,typename TimeAxis>
		Lattice<LatticeType::Binomial,std::tuple<Node,Node>,TimeAxis>
			_s_merge_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis> const &latticeA, Lattice<LatticeType::Binomial, Node, TimeAxis> const &latticeB) {
			auto fdA = latticeA.fixingDates();
			assert(fdA == latticeB.fixingDates());
			Lattice<LatticeType::Binomial, std::tuple<Node, Node>, TimeAxis> result(std::set<TimeAxis>(fdA.begin(), fdA.end()));
			for (auto itr_d = fdA.begin(); itr_d != fdA.end(); ++itr_d) {
				for (auto i = 0; i < latticeA.nodesAt(*itr_d).size(); ++i) {
					result(*itr_d, i) = std::make_tuple(latticeA(*itr_d, i), latticeB(*itr_d, i));
				}
			}
			return result;
		}

		template<typename Node, typename TimeAxis>
		Lattice<LatticeType::Binomial, std::tuple<Node, Node>, TimeAxis>
			_p_merge_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis> const &latticeA, Lattice<LatticeType::Binomial, Node, TimeAxis> const &latticeB) {
			auto fdA = latticeA.fixingDates();
			assert(fdA == latticeB.fixingDates());
			Lattice<LatticeType::Binomial, std::tuple<Node, Node>, TimeAxis> result(std::set<TimeAxis>(fdA.begin(), fdA.end()));
			concurrency::parallel_for_each(fdA.cbegin(), fdA.cend(), [&](TimeAxis ta) {
				for (auto i = 0; i < latticeA.nodesAt(ta).size(); ++i) {
					result(ta, i) = std::make_tuple(latticeA(ta, i), latticeB(ta, i));
				}
			
			});
			return result;
		}

		template<typename Node, typename TimeAxis>
		Lattice<LatticeType::Trinomial, std::tuple<Node, Node>, TimeAxis>
			_s_merge_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>const &latticeA, Lattice<LatticeType::Trinomial, Node, TimeAxis> const &latticeB) {
			auto fdA = latticeA.fixingDates();
			assert(fdA == latticeB.fixingDates());
			Lattice<LatticeType::Trinomial, std::tuple<Node, Node>, TimeAxis> result(std::set<TimeAxis>(fdA.begin(), fdA.end()));
			for (auto itr_d = fdA.begin(); itr_d != fdA.end(); ++itr_d) {
				for (auto i = 0; i < latticeA.nodesAt(*itr_d).size(); ++i) {
					result(*itr_d, i) = std::make_tuple(latticeA(*itr_d, i), latticeB(*itr_d, i));
				}
			}
			return result;
		}

		template<typename Node, typename TimeAxis>
		Lattice<LatticeType::Trinomial, std::tuple<Node, Node>, TimeAxis>
			_p_merge_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>const &latticeA, Lattice<LatticeType::Trinomial, Node, TimeAxis> const &latticeB) {
			auto fdA = latticeA.fixingDates();
			assert(fdA == latticeB.fixingDates());
			Lattice<LatticeType::Trinomial, std::tuple<Node, Node>, TimeAxis> result(std::set<TimeAxis>(fdA.begin(), fdA.end()));
			concurrency::parallel_for_each(fdA.cbegin(), fdA.cend(), [&](TimeAxis ta) {
				for (auto i = 0; i < latticeA.nodesAt(ta).size(); ++i) {
					result(ta, i) = std::make_tuple(latticeA(ta, i), latticeB(ta, i));
				}
			});
			return result;
		}



	}

	template<typename Node>
	IndexedLattice<LatticeType::Binomial,std::tuple<Node,Node>> 
		merge(IndexedLattice<LatticeType::Binomial, Node> const &latticeA, IndexedLattice<LatticeType::Binomial, Node> const &latticeB,
			lattice_types::Launch launch = lattice_types::Launch::Sequential) {
		if (launch == lattice_types::Launch::Sequential) 
			return _s_merge_binomial_indexed_impl(latticeA, latticeB);
		else if(launch == lattice_types::Launch::Parallel)
			return _p_merge_binomial_indexed_impl(latticeA, latticeB);
	}

	template<typename Node>
	IndexedLattice<LatticeType::Trinomial, std::tuple<Node, Node>>
		merge(IndexedLattice<LatticeType::Trinomial, Node> const &latticeA, IndexedLattice<LatticeType::Trinomial, Node> const &latticeB,
			lattice_types::Launch launch = lattice_types::Launch::Sequential) {
		if (launch == lattice_types::Launch::Sequential)
			return _s_merge_trinomial_indexed_impl(latticeA, latticeB);
		else if (launch == lattice_types::Launch::Parallel)
			return _p_merge_trinomial_indexed_impl(latticeA, latticeB);
	}

	template<typename Node,typename TimeAxis>
	Lattice<LatticeType::Binomial, std::tuple<Node, Node>, TimeAxis>
		merge(Lattice<LatticeType::Binomial, Node, TimeAxis>const &latticeA, Lattice<LatticeType::Binomial, Node, TimeAxis>const &latticeB,
			lattice_types::Launch launch = lattice_types::Launch::Sequential) {
		if (launch == lattice_types::Launch::Sequential)
			return _s_merge_binomial_impl(latticeA, latticeB);
		else if (launch == lattice_types::Launch::Parallel)
			return _p_merge_binomial_impl(latticeA, latticeB);
	}

	template<typename Node, typename TimeAxis>
	Lattice<LatticeType::Trinomial, std::tuple<Node, Node>, TimeAxis>
		merge(Lattice<LatticeType::Trinomial, Node, TimeAxis>const &latticeA, Lattice<LatticeType::Trinomial, Node, TimeAxis>const &latticeB,
			lattice_types::Launch launch = lattice_types::Launch::Sequential) {
		if (launch == lattice_types::Launch::Sequential)
			return _s_merge_trinomial_impl(latticeA, latticeB);
		else if (launch == lattice_types::Launch::Parallel)
			return _p_merge_trinomial_impl(latticeA, latticeB);
	}



}




#endif //_LATTICE_ALGORITHMS