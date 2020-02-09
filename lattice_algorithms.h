#pragma once
#if !defined(_LATTICE_ALGORITHMS)
#define _LATTICE_ALGORITHMS

#include<cassert>
#include"lattice_structure.h"


namespace lattice_algorithms {

	using lattice_types::LatticeType;
	using lattice_types::LeafForwardGenerator;
	using lattice_types::LeafBackwardGenerator;
	using lattice_types::Payoff;
	using lattice_structure::IndexedLattice;
	using lattice_structure::Lattice;

	//==================== Forward Induction ============================

	namespace {

		template<typename Node, typename DeltaTime>
		void _forward_induction_indexed_binomial_impl(IndexedLattice<LatticeType::Binomial, Node>& lattice,
			LeafForwardGenerator<Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,std::true_type) {
			assert(deltaTime.size() == lattice.maxIndex());
			lattice(0, 0) = apex;
			std::tuple<Node, Node> tuple;
			for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice(t - 1, l), deltaTime[t-1]);
					lattice(t, l) = std::get<0>(tuple);
					lattice(t, l + 1) = std::get<1>(tuple);
				}
			}
		}

		template<typename Node, typename DeltaTime>
		void _forward_induction_indexed_binomial_impl(IndexedLattice<LatticeType::Binomial, Node>& lattice,
			LeafForwardGenerator<Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime, std::false_type) {
			lattice(0, 0) = apex;
			std::tuple<Node, Node> tuple;
			for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice(t - 1, l), deltaTime);
					lattice(t, l) = std::get<0>(tuple);
					lattice(t, l + 1) = std::get<1>(tuple);
				}
			}
		}


		template<typename Node,typename DeltaTime>
		void _forward_induction_indexed_trinomial_impl(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime,std::true_type) {
			assert(deltaTime.size() == lattice.maxIndex());
			lattice(0, 0) = apex;
			std::tuple<Node, Node, Node> tuple;
			for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice(t - 1, l),deltaTime[t-1]);
					lattice(t, l) = std::get<0>(tuple);
					lattice(t, l + 1) = std::get<1>(tuple);
					lattice(t, l + 2) = std::get<2>(tuple);
				}
			}
		}

		template<typename Node, typename DeltaTime>
		void _forward_induction_indexed_trinomial_impl(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex,DeltaTime deltaTime,std::false_type) {
			lattice(0, 0) = apex;
			std::tuple<Node, Node, Node> tuple;
			for (std::size_t t = lattice.minIndex() + 1; t <= lattice.maxIndex(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(t - 1).size(); ++l) {
					tuple = generator(lattice(t - 1, l), deltaTime);
					lattice(t, l) = std::get<0>(tuple);
					lattice(t, l + 1) = std::get<1>(tuple);
					lattice(t, l + 2) = std::get<2>(tuple);
				}
			}
		}


		template<typename Node, typename TimeAxis,typename DeltaTime>
		void _forward_induction_lattice_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node, Node> const &generator, Node apex,DeltaTime deltaTime,std::true_type) {
			auto fixingDates = lattice.fixingDates();
			assert(deltaTime.size() == (fixingDates.size() - 1));
			lattice(fixingDates[0], 0) = apex;
			std::tuple<Node, Node> tuple;
			for (std::size_t t = 1; t < fixingDates.size(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice(fixingDates[t - 1], l), deltaTime[t - 1]);
					lattice(fixingDates[t], l) = std::get<0>(tuple);
					lattice(fixingDates[t], l + 1) = std::get<1>(tuple);
				}
			}
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _forward_induction_lattice_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime, std::false_type) {
			auto fixingDates = lattice.fixingDates();
			lattice(fixingDates[0], 0) = apex;
			std::tuple<Node, Node> tuple;
			for (std::size_t t = 1; t < fixingDates.size(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice(fixingDates[t - 1], l), deltaTime);
					lattice(fixingDates[t], l) = std::get<0>(tuple);
					lattice(fixingDates[t], l + 1) = std::get<1>(tuple);
				}
			}
		}

		template<typename Node, typename TimeAxis,typename DeltaTime>
		void _forward_induction_lattice_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex,DeltaTime deltaTime,std::true_type) {
			auto fixingDates = lattice.fixingDates();
			assert(deltaTime.size() == (fixingDates.size() - 1));
			lattice(fixingDates[0], 0) = apex;
			std::tuple<Node, Node, Node> tuple;
			for (std::size_t t = 1; t < fixingDates.size(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice(fixingDates[t - 1], l), deltaTime[t - 1]);
					lattice(fixingDates[t], l) = std::get<0>(tuple);
					lattice(fixingDates[t], l + 1) = std::get<1>(tuple);
					lattice(fixingDates[t], l + 2) = std::get<2>(tuple);
				}
			}
		}

		template<typename Node, typename TimeAxis,typename DeltaTime>
		void _forward_induction_lattice_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafForwardGenerator<Node, Node, Node, Node> const &generator, Node apex, DeltaTime deltaTime, std::false_type) {
			auto fixingDates = lattice.fixingDates();
			lattice(fixingDates[0], 0) = apex;
			std::tuple<Node, Node, Node> tuple;
			for (std::size_t t = 1; t < fixingDates.size(); ++t) {
				for (std::size_t l = 0; l < lattice.nodesAt(fixingDates[t - 1]).size(); ++l) {
					tuple = generator(lattice(fixingDates[t - 1], l), deltaTime);
					lattice(fixingDates[t], l) = std::get<0>(tuple);
					lattice(fixingDates[t], l + 1) = std::get<1>(tuple);
					lattice(fixingDates[t], l + 2) = std::get<2>(tuple);
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

	//==================== Backward Induction ============================

	namespace {

		template<typename Node,typename DeltaTime>
		void _backward_induction_indexed_binomial_impl(IndexedLattice<LatticeType::Binomial, Node>& lattice,
			LeafBackwardGenerator<Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::true_type) {

			assert(deltaTime.size() == lattice.maxIndex());
			std::size_t lastIdx = lattice.maxIndex();
			for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
				lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
			}

			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
					lattice(n, i) = backwardGenerator(lattice(n + 1, i), lattice(n + 1, i + 1), deltaTime[n]);
				}
			}
			lattice(0, 0) = backwardGenerator(lattice(1, 0), lattice(1, 1), deltaTime[0]);
		}

		template<typename Node,typename DeltaTime>
		void _backward_induction_indexed_binomial_impl(IndexedLattice<LatticeType::Binomial, Node>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::false_type) {

			std::size_t lastIdx = lattice.maxIndex();
			for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
				lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
			}

			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
					lattice(n, i) = backwardGenerator(lattice(n + 1, i), lattice(n + 1, i + 1), deltaTime);
				}
			}
			lattice(0, 0) = backwardGenerator(lattice(1, 0), lattice(1, 1), deltaTime);
		}

		template<typename Node, typename TimeAxis,typename DeltaTime>
		void _backward_induction_lattice_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::true_type) {

			auto fixingDates = lattice.fixingDates();
			assert(deltaTime.size() == (fixingDates.size() - 1));
			TimeAxis lastDate = fixingDates.back();
			for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
				lattice(lastDate, i) = payoff(lattice(lastDate, i));
			}

			for (auto n = fixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
					lattice(fixingDates[n], i) = backwardGenerator(lattice(fixingDates[n + 1], i), 
																	lattice(fixingDates[n + 1], i + 1),
																	deltaTime[n]);
				}
			}
			lattice(fixingDates[0], 0) = backwardGenerator(lattice(fixingDates[1], 0), lattice(fixingDates[1], 1), deltaTime[0]);
		}

		template<typename Node, typename TimeAxis,typename DeltaTime>
		void _backward_induction_lattice_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
			LeafBackwardGenerator<Node, Node, Node,Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::false_type) {

			auto fixingDates = lattice.fixingDates();
			TimeAxis lastDate = fixingDates.back();
			for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
				lattice(lastDate, i) = payoff(lattice(lastDate, i));
			}

			for (auto n = fixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
					lattice(fixingDates[n], i) = backwardGenerator(lattice(fixingDates[n + 1], i), 
																	lattice(fixingDates[n + 1], i + 1),
																	deltaTime);
				}
			}
			lattice(fixingDates[0], 0) = backwardGenerator(lattice(fixingDates[1], 0), lattice(fixingDates[1], 1), deltaTime);
		}

		template<typename Node,typename DeltaTime>
		void _backward_induction_indexed_trinomial_impl(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::true_type) {

			assert(deltaTime.size() == lattice.maxIndex());
			std::size_t lastIdx = lattice.maxIndex();
			for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
				lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
			}

			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
					lattice(n, i) = backwardGenerator(lattice(n + 1, i), lattice(n + 1, i + 1), lattice(n + 1, i + 2), deltaTime[n]);
				}
			}
			lattice(0, 0) = backwardGenerator(lattice(1, 0), lattice(1, 1), lattice(1, 2), deltaTime[0]);
		}

		template<typename Node, typename DeltaTime>
		void _backward_induction_indexed_trinomial_impl(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,
			DeltaTime deltaTime, std::false_type) {

			std::size_t lastIdx = lattice.maxIndex();
			for (auto i = 0; i < lattice.nodesAt(lastIdx).size(); ++i) {
				lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
			}

			for (auto n = lastIdx - 1; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(n).size(); ++i) {
					lattice(n, i) = backwardGenerator(lattice(n + 1, i), lattice(n + 1, i + 1), lattice(n + 1, i + 2), deltaTime);
				}
			}
			lattice(0, 0) = backwardGenerator(lattice(1, 0), lattice(1, 1), lattice(1, 2), deltaTime);
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff, 
			DeltaTime deltaTime,std::true_type) {

			auto fixingDates = lattice.fixingDates();
			assert(deltaTime.size() == (fixingDates.size() - 1));
			TimeAxis lastDate = fixingDates.back();
			for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
				lattice(lastDate, i) = payoff(lattice(lastDate, i));
			}

			for (auto n = fixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
					lattice(fixingDates[n], i) = backwardGenerator(lattice(fixingDates[n + 1], i),
						lattice(fixingDates[n + 1], i + 1),
						lattice(fixingDates[n + 1], i + 2),
						deltaTime[n]);
				}
			}
			lattice(fixingDates[0], 0) = backwardGenerator(lattice(fixingDates[1], 0),
				lattice(fixingDates[1], 1),
				lattice(fixingDates[1], 2),
				deltaTime[0]);
		}

		template<typename Node, typename TimeAxis, typename DeltaTime>
		void _backward_induction_lattice_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
			LeafBackwardGenerator<Node, Node, Node, Node, Node>const &backwardGenerator, Payoff<Node, Node> const &payoff, 
			DeltaTime deltaTime,std::false_type) {

			auto fixingDates = lattice.fixingDates();
			TimeAxis lastDate = fixingDates.back();
			for (auto i = 0; i < lattice.nodesAt(lastDate).size(); ++i) {
				lattice(lastDate, i) = payoff(lattice(lastDate, i));
			}

			for (auto n = fixingDates.size() - 2; n > 0; --n) {
				for (auto i = 0; i < lattice.nodesAt(fixingDates[n]).size(); ++i) {
					lattice(fixingDates[n], i) = backwardGenerator(lattice(fixingDates[n + 1], i),
						lattice(fixingDates[n + 1], i + 1),
						lattice(fixingDates[n + 1], i + 2),
						deltaTime);
				}
			}
			lattice(fixingDates[0], 0) = backwardGenerator(lattice(fixingDates[1], 0),
				lattice(fixingDates[1], 1),
				lattice(fixingDates[1], 2),
				deltaTime);
		}

	}
	
	template<typename Node,typename DeltaTime>
	void backward_induction(IndexedLattice<LatticeType::Binomial, Node>& lattice,
		LeafBackwardGenerator<Node, Node, Node, Node>const &backwardGenerator, Payoff<Node,Node> const &payoff,DeltaTime deltaTime) {
		_backward_induction_indexed_binomial_impl(lattice, backwardGenerator, payoff, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename TimeAxis,typename DeltaTime>
	void backward_induction(Lattice<LatticeType::Binomial, Node, TimeAxis>& lattice,
		LeafBackwardGenerator<Node, Node, Node,Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,DeltaTime deltaTime) {
		_backward_induction_lattice_binomial_impl(lattice, backwardGenerator, payoff, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node,typename DeltaTime>
	void backward_induction(IndexedLattice<LatticeType::Trinomial, Node>& lattice,
		LeafBackwardGenerator<Node,Node, Node, Node,Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,DeltaTime deltaTime) {
		_backward_induction_indexed_trinomial_impl(lattice, backwardGenerator, payoff, deltaTime, std::is_compound<DeltaTime>());
	}

	template<typename Node, typename TimeAxis,typename DeltaTime>
	void backward_induction(Lattice<LatticeType::Trinomial, Node, TimeAxis>& lattice,
		LeafBackwardGenerator<Node, Node,Node, Node,Node>const &backwardGenerator, Payoff<Node, Node> const &payoff,DeltaTime deltaTime) {
		_backward_induction_lattice_trinomial_impl(lattice, backwardGenerator, payoff, deltaTime, std::is_compound<DeltaTime>());
	}



}




#endif //_LATTICE_ALGORITHMS