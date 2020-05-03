#pragma once
#if !defined(_LATTICE_ALGORITHMS)
#define _LATTICE_ALGORITHMS

#include<cassert>
#include<ppl.h>
#include"lattice_structure.h"
#include"lattice_forward_traversals.h"
#include"lattice_backward_traversals.h"
#include"lattice_types.h"
#include"lattice_macros.h"
#include"lattice_traits.h"
#include"lattice_merge.h"

namespace lattice_algorithms {

	using lattice_types::LatticeType;
	using lattice_types::LatticeClass;
	using lattice_types::LaunchPolicy;
	using lattice_structure::IndexedLattice;
	using lattice_structure::Lattice;
	using lattice_forward_traversals::ForwardTraversal;
	using lattice_backward_traversals::BackwardTraversal;
	using lattice_traits::MergeTraits;
	using lattice_merge::MergePolicy;
	

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
			LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
			ForwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime, Node>::
				traverse(lattice, std::forward<Generator>(generator), deltaTime, apex);
		}

		template<typename LatticeObject, typename Generator>
		void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
			std::map<TimeAxis,Node> const &dividendData) {
			LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
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
				LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				ForwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime, Node>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, apex);
			}

			template<typename LatticeObject, typename Generator>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis, Node> const &dividendData) {
				LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
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
				LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				BackwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff));
			}

			template<typename LatticeObject, typename Generator, typename Payoff,typename PayoffAdjuster>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {
				LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
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
				LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				BackwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff));
			}

			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {
				LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				BackwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff), std::forward<PayoffAdjuster>(payoffAdjuster));
			}
	};





	// ==============================================================================
	// ============================ Merging Algorithms ==============================
	// ==============================================================================

	template<typename TimeAxis>
	class MergeLattices {
	public:
		template<typename LatticeObject,
			typename ...LatticeObjects,
			typename Traits = MergeTraits<sizeof...(LatticeObjects), LatticeObject::Node_type>>
			auto operator()(LaunchPolicy launch,LatticeObject latticeObject, LatticeObjects... latticeObjects)const {

			auto result = MergePolicy<TimeAxis>::construct(latticeObject,
				std::is_same<LatticeObject, Lattice< LatticeObject::type(), typename LatticeObject::Node_type, TimeAxis>>(),
				latticeObjects...);

			if (launch == LaunchPolicy::Sequential) {
				MergePolicy<TimeAxis>::runSequential(result,latticeObject,latticeObjects...);
			}
			else {
				MergePolicy<TimeAxis>::runParallel(result, latticeObject, latticeObjects...);
			}
			return result;
		}

		template<typename LatticeObject>
			auto operator()(LaunchPolicy launch, LatticeObject latticeObject)const {
			return latticeObject;
		}

	};





	



}




#endif //_LATTICE_ALGORITHMS