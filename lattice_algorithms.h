#pragma once
#if !defined(_LATTICE_ALGORITHMS)
#define _LATTICE_ALGORITHMS

#include<cassert>
#include<ppl.h>
#include<type_traits>
#include"lattice_structure.h"
#include"lattice_multidimensional.h"
#include"lattice_forward_traversals.h"
#include"lattice_backward_traversals.h"
#include"lattice_calibrator_results.h"
#include"lattice_types.h"
#include"lattice_macros.h"
#include"lattice_traits.h"
#include"lattice_merge.h"

namespace lattice_algorithms {

	using lattice_types::LatticeType;
	using lattice_types::LatticeClass;
	using lattice_types::AssetClass;
	using lattice_types::LaunchPolicy;
	using lattice_types::BarrierType;
	using lattice_types::DiscountingStyle;
	using lattice_structure::IndexedLattice;
	using lattice_structure::Lattice;
	using lattice_structure::GeneralLattice;
	using lattice_calibrator_results::CalibratorTrinomialEquityResultsPtr;
	using lattice_forward_traversals::ForwardTraversal;
	using lattice_backward_traversals::BackwardTraversal;
	using lattice_backward_traversals::ImpliedBackwardTraversal;
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
			LASSERT(lattice.factors() == 1, "For one-dimensional lattice only");
			ForwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime, Node>::
				traverse(lattice, std::forward<Generator>(generator), deltaTime, apex);
		}

		template<typename LatticeObject, typename Generator>
		void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
			std::map<TimeAxis,Node> const &dividendData) {
			LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
			LASSERT(lattice.factors() == 1, "For one-dimensional lattice only");
			ForwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime, Node>::
				traverse(lattice, std::forward<Generator>(generator), deltaTime, apex, dividendData);
		}
	};


	template<typename TimeAxis,
		typename DeltaTime,
		typename Node>
		class ForwardInduction<LatticeType::TwoVariableBinomial, TimeAxis, DeltaTime, Node> {
		public:
			template<typename MultidimLatticeObject, typename Generator>
			void operator()(MultidimLatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, std::pair<Node, Node> const &apex) {
				LASSERT(lattice.getFactor(0).type() == generator.latticeType(), "Mismatch between lattice types");
				LASSERT(lattice.factors() > 1, "For multidimensional lattice only");
				ForwardTraversal<LatticeType::TwoVariableBinomial, TimeAxis, DeltaTime, Node>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, apex);
			}

			template<typename MultidimLatticeObject, typename Generator>
			void operator()(MultidimLatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, std::pair<Node, Node> const &apex,
				std::pair<std::map<TimeAxis, Node>, std::map<TimeAxis, Node>> const &dividendData) {
				LASSERT(lattice.getFactor(0).type() == generator.latticeType(), "Mismatch between lattice types");
				LASSERT(lattice.factors() > 1, "For multidimensional lattice only");
				ForwardTraversal<LatticeType::TwoVariableBinomial, TimeAxis, DeltaTime, Node>::
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
				LASSERT(lattice.factors() == 1, "For one-dimensional lattice only");
				ForwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime, Node>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, apex);
			}

			template<typename LatticeObject, typename Generator>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis, Node> const &dividendData) {
				LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				LASSERT(lattice.factors() == 1, "For one-dimensional lattice only");
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
				LASSERT(lattice.factors() == 1, "For one-dimensional lattice only");
				BackwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff));
			}

			template<typename LatticeObject, typename Generator, typename Payoff,typename PayoffAdjuster>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {
				LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				LASSERT(lattice.factors() == 1, "For one-dimensional lattice only");
				BackwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff),std::forward<PayoffAdjuster>(payoffAdjuster));
			}

			template<typename LatticeObject, typename Generator, typename Payoff>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff,
						BarrierType barrierType,typename LatticeObject::Node_type const &barrierValue,
						typename LatticeObject::Node_type const &rebateValue) {
				LASSERT(lattice.type() == generator.latticeType(), "Mismatch between lattice types");
				LASSERT(lattice.factors() == 1, "For one-dimensional lattice only");
				BackwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime>::
					traverseBarrier(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff),
						barrierType, barrierValue, rebateValue);
			}

			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, 
				PayoffAdjuster &&payoffAdjuster,BarrierType barrierType, typename LatticeObject::Node_type const &barrierValue,
				typename LatticeObject::Node_type const &rebateValue) {
				LASSERT(lattice.type() == generator.latticeType(), "Mismatch between lattice types");
				LASSERT(lattice.factors() == 1, "For one-dimensional lattice only");
				BackwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime>::
					traverseBarrier(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff),
						std::forward<PayoffAdjuster>(payoffAdjuster), barrierType, barrierValue, rebateValue);
			}

	};

	template<typename TimeAxis,
		typename DeltaTime>
		class BackwardInduction<LatticeType::TwoVariableBinomial, TimeAxis, DeltaTime> {
		public:
			template<typename LatticeObject,typename MultidimLatticeObject, typename Generator, typename Payoff>
			void operator()(LatticeObject &priceLattice, MultidimLatticeObject const &lattice, Generator &&generator,
				DeltaTime const &deltaTime, Payoff &&payoff) {
				LASSERT(priceLattice.type() == LatticeType::TwoVariableBinomial, "Mismatch between lattice types");
				LASSERT(generator.latticeType() == LatticeType::Binomial, "Mismatch between lattice types");
				LASSERT(lattice.factors() == 2, "For two-dimensional lattice only");
				BackwardTraversal<LatticeType::TwoVariableBinomial, TimeAxis, DeltaTime>::
					traverse(priceLattice, lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff));
			}

			template<typename LatticeObject, typename MultidimLatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			void operator()(LatticeObject &priceLattice, MultidimLatticeObject const &lattice, Generator &&generator,
				DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {
				LASSERT(priceLattice.type() == LatticeType::TwoVariableBinomial, "Mismatch between lattice types");
				LASSERT(generator.latticeType() == LatticeType::Binomial, "Mismatch between lattice types");
				LASSERT(lattice.factors() == 2, "For two-dimensional lattice only");
				BackwardTraversal<LatticeType::TwoVariableBinomial, TimeAxis, DeltaTime>::
					traverse(priceLattice,lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff), std::forward<PayoffAdjuster>(payoffAdjuster));
			}

	};

	template<typename TimeAxis,
		typename DeltaTime>
		class BackwardInduction<LatticeType::Trinomial, TimeAxis, DeltaTime> {
		public:
			template<typename LatticeObject, typename Generator, typename Payoff>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff) {
				LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				LASSERT(lattice.factors() == 1, "For one-dimensional lattice only");
				BackwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff));
			}

			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {
				LASSERT(lattice.type() == generator.latticeType(),"Mismatch between lattice types");
				LASSERT(lattice.factors() == 1, "For one-dimensional lattice only");
				BackwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff), std::forward<PayoffAdjuster>(payoffAdjuster));
			}

			template<typename LatticeObject, typename Generator, typename Payoff>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff,
				BarrierType barrierType, typename LatticeObject::Node_type const &barrierValue,
				typename LatticeObject::Node_type const &rebateValue) {
				LASSERT(lattice.type() == generator.latticeType(), "Mismatch between lattice types");
				LASSERT(lattice.factors() == 1, "For one-dimensional lattice only");
				BackwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime>::
					traverseBarrier(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff),
						barrierType, barrierValue, rebateValue);
			}

			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			void operator()(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff,
				PayoffAdjuster &&payoffAdjuster, BarrierType barrierType, typename LatticeObject::Node_type const &barrierValue,
				typename LatticeObject::Node_type const &rebateValue) {
				LASSERT(lattice.type() == generator.latticeType(), "Mismatch between lattice types");
				LASSERT(lattice.factors() == 1, "For one-dimensional lattice only");
				BackwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime>::
					traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff),
						std::forward<PayoffAdjuster>(payoffAdjuster), barrierType, barrierValue, rebateValue);
			}
	};


	// ======================================================================================
	// ==================== Implied Backward Induction Algorithms ===========================
	// ======================================================================================

	template<LatticeType Type,
		typename DeltaTime,
		typename RiskFreeRate>
	class ImpliedBackwardInduction {};


	template<typename DeltaTime,
			typename RiskFreeRate>
	class ImpliedBackwardInduction<LatticeType::Trinomial, DeltaTime, RiskFreeRate> {
	private:
		RiskFreeRate rate_;

	public:
		ImpliedBackwardInduction() = delete;
		explicit ImpliedBackwardInduction(RiskFreeRate rate)
			:rate_{rate}{}

		template<typename LatticeObject,typename Payoff>
		void operator()(LatticeObject &optionLattice, LatticeObject const &spotLattice,
			CalibratorTrinomialEquityResultsPtr<LatticeObject> const& calibrationResults,
			DeltaTime const &deltaTime, Payoff &&payoff, DiscountingStyle style = DiscountingStyle::Discrete) {

			ImpliedBackwardTraversal<LatticeType::Trinomial, DeltaTime, RiskFreeRate>::
				traverse(optionLattice, spotLattice, calibrationResults, deltaTime,
					rate_, std::forward<Payoff>(payoff), style);
		}

		template<typename LatticeObject, typename Payoff,typename PayoffAdjuster>
		void operator()(LatticeObject &optionLattice, LatticeObject const &spotLattice,
			CalibratorTrinomialEquityResultsPtr<LatticeObject> const& calibrationResults,
			DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster, DiscountingStyle style = DiscountingStyle::Discrete) {
			ImpliedBackwardTraversal<LatticeType::Trinomial, DeltaTime, RiskFreeRate>::
				traverse(optionLattice, spotLattice, calibrationResults, deltaTime,
					rate_, std::forward<Payoff>(payoff),std::forward<PayoffAdjuster>(payoffAdjuster), style);
		}

		template<typename LatticeObject,typename Payoff>
		void operator()(LatticeObject &optionLattice, LatticeObject const &spotLattice,
			CalibratorTrinomialEquityResultsPtr<LatticeObject> const& calibrationResults,
			DeltaTime const &deltaTime, Payoff &&payoff, BarrierType barrierType, 
			typename LatticeObject::Node_type const &barrier,typename LatticeObject::Node_type const &rebate,
			DiscountingStyle style = DiscountingStyle::Discrete) {
			ImpliedBackwardTraversal<LatticeType::Trinomial, DeltaTime, RiskFreeRate>::
				traverseBarrier(optionLattice, spotLattice, calibrationResults, deltaTime,
					rate_, std::forward<Payoff>(payoff),barrierType, barrier, rebate, style);
		}

		template<typename LatticeObject, typename Payoff,typename PayoffAdjuster>
		void operator()(LatticeObject &optionLattice, LatticeObject const &spotLattice,
			CalibratorTrinomialEquityResultsPtr<LatticeObject> const& calibrationResults,
			DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster,BarrierType barrierType,
			typename LatticeObject::Node_type const &barrier, typename LatticeObject::Node_type const &rebate,
			DiscountingStyle style = DiscountingStyle::Discrete) {
			ImpliedBackwardTraversal<LatticeType::Trinomial, DeltaTime, RiskFreeRate>::
				traverseBarrier(optionLattice, spotLattice, calibrationResults, deltaTime,
					rate_, std::forward<Payoff>(payoff),std::forward<PayoffAdjuster>(payoffAdjuster),
					barrierType, barrier, rebate, style);
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