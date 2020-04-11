#pragma once
#if !defined(_LATTICE_ALGORITHMS)
#define _LATTICE_ALGORITHMS

#include<cassert>
#include<ppl.h>
#include"lattice_structure.h"
#include"lattice_forward_traversals.h"
#include"lattice_backward_traversals.h"
#include"lattice_macros.h"
#include"lattice_traits.h"


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
	using lattice_traits::MergeTraits;
	

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





	// ==============================================================================
	// ============================ Merging Algorithms ==============================
	// ==============================================================================



	template<typename TimeAxis>
	class MergeLattices {
	private:
		template<typename LatticeObject,
			typename ...LatticeObjects,
			typename Traits = MergeTraits<sizeof...(LatticeObjects), LatticeObject::Node_type>>
			auto construct_impl(LatticeObject latticeObject,std::true_type, LatticeObjects... latticeObjects)const {
			
			std::set<TimeAxis> dates;
			for (auto d : latticeObject.fixingDates()) {
				dates.emplace(d);
			}

			Lattice< LatticeObject::type(),
					typename Traits::NodeHolder,
					TimeAxis> result{ dates };

			return result;
		}

		template<typename LatticeObject,
			typename ...LatticeObjects,
			typename Traits = MergeTraits<sizeof...(LatticeObjects), LatticeObject::Node_type>>
			auto construct_impl(LatticeObject latticeObject, std::false_type, LatticeObjects... latticeObjects)const {

			IndexedLattice< LatticeObject::type(),
							typename Traits::NodeHolder> result{ latticeObject.maxIndex() };
			return result;
		}


	public:

		template<typename LatticeObject,
			typename ...LatticeObjects,
			typename Traits = MergeTraits<sizeof...(LatticeObjects), LatticeObject::Node_type>>
		auto operator()(LatticeObject latticeObject, LatticeObjects... latticeObjects)const {
			
			auto result = construct_impl(latticeObject,
					std::is_same<LatticeObject, Lattice< LatticeObject::type(),
				typename LatticeObject::Node_type,TimeAxis>>(), latticeObjects...);

			const std::size_t lastIdx = latticeObject.timeDimension() - 1;
			std::size_t nodesSize{ 0 };

			for (auto t = 0; t <= lastIdx; ++t) {
				nodesSize = latticeObject.nodesAtIdx(t).size();
				for (auto i = 0; i < nodesSize; ++i) {
					result(t, i) = std::make_tuple(latticeObject(t, i), latticeObjects(t, i)...);
				}
			}
			return result;
		}

	};




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