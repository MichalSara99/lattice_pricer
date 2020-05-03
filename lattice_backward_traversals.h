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
			template<typename LatticeObject, typename Generator, typename Payoff>
			static void _backTraverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff);

			//	This one is for payoffadjusted
			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void _backTraverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff,
				PayoffAdjuster &&payoffAdjuster);


		public:
			template<typename LatticeObject, typename Generator, typename Payoff>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff) {
				_backTraverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff));
			}

			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {
				_backTraverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff), std::forward<PayoffAdjuster>(payoffAdjuster));
			}

	};




	//	LatticeType::Trinomial specialization
	template<typename TimeAxis,
		typename DeltaTime>
		struct BackwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime> {
		private:

			template<typename LatticeObject, typename Generator>
			static void _backTraverseNormal(std::size_t timeIdx,LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime);

			template<typename LatticeObject, typename Generator>
			static void _backTraverseReverting(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime);

			template<typename LatticeObject, typename Generator, typename Payoff>
			static void _backTraverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff);


			template<typename LatticeObject, typename Generator,typename PayoffAdjuster>
			static void _backTraverseNormal(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, PayoffAdjuster &&payoffAdjuster);

			template<typename LatticeObject, typename Generator, typename PayoffAdjuster>
			static void _backTraverseReverting(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, PayoffAdjuster &&payoffAdjuster);

			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void _backTraverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff,
				PayoffAdjuster &&payoffAdjuster);


		public:
			template<typename LatticeObject, typename Generator, typename Payoff>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff) {
				_backTraverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff));
			}

			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {
				_backTraverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff), std::forward<PayoffAdjuster>(payoffAdjuster));
			}

	};


}




template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime>::
_backTraverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff) {

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
_backTraverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	const std::size_t lastIdx = lattice.timeDimension() - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();
	typename LatticeObject::Node_type dt{};
	typename LatticeObject::Node_type value{};

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = lattice.nodesAtIdx(n).size();
		dt = DT::deltaTime(n, deltaTime);
		for (auto i = 0; i < nodesSize; ++i) {
			value = generator(lattice(n, i),lattice(n + 1, i), lattice(n + 1, i + 1), dt);
			payoffAdjuster(value, lattice(n, i));
			lattice(n, i) = value;
		}
	}
	dt = DT::deltaTime(0, deltaTime);
	value = generator(lattice(0, 0),lattice(1, 0), lattice(1, 1), dt);
	payoffAdjuster(value, lattice(0, 0));
	lattice(0, 0) = value;
}





template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
_backTraverseNormal(std::size_t timeIdx,LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	typename LatticeObject::Node_type dt{};

	std::size_t nodesSize{ 0 };
	for (auto n = timeIdx - 1; n > 0; --n) {
		dt = DT::deltaTime(n, deltaTime);
		nodesSize = lattice.nodesAtIdx(n).size();
		for (auto i = 0; i < nodesSize; ++i) {
			lattice(n, i) = generator(lattice(n, i), lattice(n + 1, i), lattice(n + 1, i + 1), lattice(n + 1, i + 2), dt);
		}
	}
	dt = DT::deltaTime(0, deltaTime);
	lattice(0, 0) = generator(lattice(0, 0), lattice(1, 0), lattice(1, 1), lattice(1, 2), dt);

}

template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
_backTraverseReverting(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	typename LatticeObject::Node_type dt{};
	const std::size_t lastIdx = lattice.timeDimension() - 1;

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n >= timeIdx; --n) {
		dt = DT::deltaTime(n, deltaTime);
		nodesSize = lattice.nodesAtIdx(n).size();
		lattice(n, 0) = generator(lattice(n, 0), lattice(n + 1, 0), lattice(n + 1, 1), lattice(n + 1, 2), dt);
		for (auto i = 1; i < nodesSize - 1; ++i) {
			lattice(n, i) = generator(lattice(n, i), lattice(n + 1, i - 1), lattice(n + 1, i), lattice(n + 1, i + 1), dt);
		}
		lattice(n, nodesSize - 1) = generator(lattice(n, nodesSize - 1), lattice(n + 1, nodesSize - 3), lattice(n + 1, nodesSize - 2), lattice(n + 1, nodesSize - 1), dt);
	}
}



template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
_backTraverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff) {

	const std::size_t firstRevertIdx = lattice.firstRevertingIdx();
	const std::size_t treeSize = lattice.timeDimension();

	const std::size_t lastIdx = treeSize - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}

	if (firstRevertIdx == 0) {
		// This trinomial tree does not have reverting property:
		_backTraverseNormal(lastIdx, lattice, std::forward<Generator>(generator), deltaTime);
	}
	else {
		// This trinomial tree does have reverting property:
		_backTraverseReverting(firstRevertIdx - 1, lattice, std::forward<Generator>(generator), deltaTime);
		_backTraverseNormal(firstRevertIdx - 1, lattice, std::forward<Generator>(generator), deltaTime);

	}
}




template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
_backTraverseNormal(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, PayoffAdjuster &&payoffAdjuster) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	typename LatticeObject::Node_type dt{};
	typename LatticeObject::Node_type value{};

	std::size_t nodesSize{ 0 };
	for (auto n = timeIdx - 1; n > 0; --n) {
		dt = DT::deltaTime(n, deltaTime);
		nodesSize = lattice.nodesAtIdx(n).size();
		for (auto i = 0; i < nodesSize; ++i) {
			value = generator(lattice(n, i), lattice(n + 1, i), lattice(n + 1, i + 1), lattice(n + 1, i + 2), dt);
			payoffAdjuster(value, lattice(n, i));
			lattice(n, i) = value;
		}
	}
	dt = DT::deltaTime(0, deltaTime);
	value = generator(lattice(0, 0), lattice(1, 0), lattice(1, 1), lattice(1, 2), dt);
	payoffAdjuster(value, lattice(0, 0));
	lattice(0, 0) = value;
}


template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
_backTraverseReverting(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, PayoffAdjuster &&payoffAdjuster) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	typename LatticeObject::Node_type dt{};
	typename LatticeObject::Node_type value{};
	const std::size_t lastIdx = lattice.timeDimension() - 1;

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n >= timeIdx; --n) {
		dt = DT::deltaTime(n, deltaTime);
		nodesSize = lattice.nodesAtIdx(n).size();

		value = generator(lattice(n, 0), lattice(n + 1, 0), lattice(n + 1, 1), lattice(n + 1, 2), dt);
		payoffAdjuster(value, lattice(n, 0));
		lattice(n, 0) = value;

		for (auto i = 1; i < nodesSize - 1; ++i) {
			value = generator(lattice(n, i), lattice(n + 1, i), lattice(n + 1, i + 1), lattice(n + 1, i + 2), dt);
			payoffAdjuster(value, lattice(n, i));
			lattice(n, i) = value;
		}

		value = generator(lattice(n, nodesSize - 1), lattice(n + 1, nodesSize - 3), lattice(n + 1, nodesSize - 2), lattice(n + 1, nodesSize - 1), dt);
		payoffAdjuster(value, lattice(n, nodesSize - 1));
		lattice(n, nodesSize - 1) = value;
	}

}


template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
_backTraverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {

	const std::size_t firstRevertIdx = lattice.firstRevertingIdx();
	const std::size_t treeSize = lattice.timeDimension();

	const std::size_t lastIdx = treeSize - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}

	if (firstRevertIdx == 0) {
		// This trinomial tree does not have reverting property:
		_backTraverseNormal(lastIdx, lattice, std::forward<Generator>(generator), deltaTime, std::forward<PayoffAdjuster>(payoffAdjuster));
	}
	else {
		// This trinomial tree does have reverting property:
		_backTraverseReverting(firstRevertIdx - 1, lattice, std::forward<Generator>(generator), deltaTime,std::forward<PayoffAdjuster>(payoffAdjuster));
		_backTraverseNormal(firstRevertIdx - 1, lattice, std::forward<Generator>(generator), deltaTime, std::forward<PayoffAdjuster>(payoffAdjuster));

	}
}



#endif ///_LATTICE_BACKWARD_TRAVERSALS