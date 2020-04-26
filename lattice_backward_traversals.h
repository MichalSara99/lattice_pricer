#pragma once
#if !defined(_LATTICE_BACKWARD_TRAVERSALS)
#define _LATTICE_BACKWARD_TRAVERSALS

#include"lattice_types.h"
#include"lattice_macros.h"
#include"lattice_utility.h"

namespace lattice_backward_traversals {

	using lattice_types::LatticeType;
	using lattice_utility::DeltaTimeHolder;

	// ==============================================================================
	// ==================== Backward Traversal Algorithms ===========================
	// ==============================================================================

	template<LatticeType Type,
		typename TimeAxis,
		typename DeltaTime>
		struct BackwardTraversal {};

	//	LatticeType::Binomial specialization
	template<typename TimeAxis,
		typename DeltaTime>
		struct BackwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime> {
		private:
			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator, typename Payoff>
			static void _back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff);

			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void _back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff,
				PayoffAdjuster &&payoffAdjuster);


		public:
			template<typename LatticeObject, typename Generator, typename Payoff>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff) {
				_back_traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff));
			}

			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {
				_back_traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff), std::forward<PayoffAdjuster>(payoffAdjuster));
			}

	};




	//	LatticeType::Trinomial specialization
	template<typename TimeAxis,
		typename DeltaTime>
		struct BackwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime> {
		private:
			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator, typename Payoff>
			static void _back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff);

			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void _back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff,
				PayoffAdjuster &&payoffAdjuster);


		public:
			template<typename LatticeObject, typename Generator, typename Payoff>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff) {
				_back_traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff));
			}

			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {
				_back_traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff), std::forward<PayoffAdjuster>(payoffAdjuster));
			}

	};


}




template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime>::
_back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	const std::size_t lastIdx = lattice.timeDimension() - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();
	typename LatticeObject::Node_type dt{};

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}
	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = lattice.nodesAtIdx(n).size();
		dt = DT::deltaTime(n, deltaTime);
		for (auto i = 0; i < nodesSize; ++i) {
			lattice(n, i) = generator(lattice(n, i),lattice(n + 1, i), lattice(n + 1, i + 1), dt);
		}
	}
	dt = DT::deltaTime(0, deltaTime);
	lattice(0, 0) = generator(lattice(0, 0),lattice(1, 0), lattice(1, 1), dt);
}



template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime>::
_back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	const std::size_t lastIdx = lattice.timeDimension() - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();
	typename LatticeObject::Node_type dt{};

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = lattice.nodesAtIdx(n).size();
		dt = DT::deltaTime(n, deltaTime);
		for (auto i = 0; i < nodesSize; ++i) {
			auto value = generator(lattice(n, i),lattice(n + 1, i), lattice(n + 1, i + 1), dt);
			payoffAdjuster(value, lattice(n, i));
			lattice(n, i) = value;
		}
	}
	dt = DT::deltaTime(0, deltaTime);
	auto lastValue = generator(lattice(0, 0),lattice(1, 0), lattice(1, 1), dt);
	payoffAdjuster(lastValue, lattice(0, 0));
	lattice(0, 0) = lastValue;
}


template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
_back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	const std::size_t lastIdx = lattice.timeDimension() - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();
	typename LatticeObject::Node_type dt{};

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}
	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		dt = DT::deltaTime(n, deltaTime);
		nodesSize = lattice.nodesAtIdx(n).size();
		for (auto i = 0; i < nodesSize; ++i) {
			lattice(n, i) = generator(lattice(n, i), lattice(n + 1, i), lattice(n + 1, i + 1), lattice(n + 1, i + 2), dt);
		}
	}
	dt = DT::deltaTime(0, deltaTime);
	lattice(0, 0) = generator(lattice(0, 0),lattice(1, 0), lattice(1, 1), lattice(1, 2), dt);

}


template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
_back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	const std::size_t lastIdx = lattice.timeDimension() - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();
	typename LatticeObject::Node_type dt{};

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		dt = DT::deltaTime(n, deltaTime);
		nodesSize = lattice.nodesAtIdx(n).size();
		for (auto i = 0; i < nodesSize; ++i) {
			auto value  = generator(lattice(n, i),lattice(n + 1, i), lattice(n + 1, i + 1), lattice(n + 1, i + 2), dt);
			payoffAdjuster(value, lattice(n, i));
			lattice(n, i) = value;
		}
	}
	dt = DT::deltaTime(0, deltaTime);
	auto lastValue = generator(lattice(0, 0),lattice(1, 0), lattice(1, 1), lattice(1, 2), dt);
	payoffAdjuster(lastValue, lattice(0, 0));
	lattice(0, 0) = lastValue;
}



#endif ///_LATTICE_BACKWARD_TRAVERSALS