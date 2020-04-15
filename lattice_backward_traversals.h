#pragma once
#if !defined(_LATTICE_BACKWARD_TRAVERSALS)
#define _LATTICE_BACKWARD_TRAVERSALS

#include"lattice_types.h"
#include"lattice_macros.h"

namespace lattice_backward_traversals {

	using lattice_types::LatticeType;

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
			static void _back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, std::true_type);
			//	This one is for non-compound DeltaTime object
			template<typename LatticeObject, typename Generator, typename Payoff>
			static void _back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, std::false_type);
			
			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void _back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff,
				PayoffAdjuster &&payoffAdjuster,std::true_type);
			//	This one is for non-compound DeltaTime object
			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void _back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff,
				PayoffAdjuster &&payoffAdjuster, std::false_type);


		public:
			template<typename LatticeObject, typename Generator, typename Payoff>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff) {
				_back_traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff), std::is_compound<DeltaTime>());
			}

			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {
				_back_traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff), std::forward<PayoffAdjuster>(payoffAdjuster), std::is_compound<DeltaTime>());
			}

	};




	//	LatticeType::Trinomial specialization
	template<typename TimeAxis,
		typename DeltaTime>
		struct BackwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime> {
		private:
			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator, typename Payoff>
			static void _back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, std::true_type);
			//	This one is for non-compound DeltaTime object
			template<typename LatticeObject, typename Generator, typename Payoff>
			static void _back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, std::false_type);

			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void _back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff,
				PayoffAdjuster &&payoffAdjuster, std::true_type);
			//	This one is for non-compound DeltaTime object
			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void _back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff,
				PayoffAdjuster &&payoffAdjuster, std::false_type);

		public:
			template<typename LatticeObject, typename Generator, typename Payoff>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff) {
				_back_traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff), std::is_compound<DeltaTime>());
			}

			template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster) {
				_back_traverse(lattice, std::forward<Generator>(generator), deltaTime, std::forward<Payoff>(payoff), std::forward<PayoffAdjuster>(payoffAdjuster), std::is_compound<DeltaTime>());
			}

	};


}




template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime>::
_back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, std::true_type) {

	LASSERT((deltaTime.size() + 1) == lattice.timeDimension(), "deltaTime is always one element shorter then time dimension");
	const std::size_t lastIdx = lattice.timeDimension() - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}
	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = lattice.nodesAtIdx(n).size();
		for (auto i = 0; i < nodesSize; ++i) {
			lattice(n, i) = generator(lattice(n, i),lattice(n + 1, i), lattice(n + 1, i + 1), deltaTime[n]);
		}
	}
	lattice(0, 0) = generator(lattice(0, 0),lattice(1, 0), lattice(1, 1), deltaTime[0]);
}


template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime>::
_back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, std::false_type) {

	const std::size_t lastIdx = lattice.timeDimension() - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}
	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = lattice.nodesAtIdx(n).size();
		for (auto i = 0; i < nodesSize; ++i) {
			lattice(n, i) = generator(lattice(n, i),lattice(n + 1, i), lattice(n + 1, i + 1), deltaTime);
		}
	}
	lattice(0, 0) = generator(lattice(0, 0),lattice(1, 0), lattice(1, 1), deltaTime);

}

template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime>::
_back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster, std::true_type) {

	LASSERT((deltaTime.size() + 1) == lattice.timeDimension(), "deltaTime is always one element shorter then time dimension");
	const std::size_t lastIdx = lattice.timeDimension() - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = lattice.nodesAtIdx(n).size();
		for (auto i = 0; i < nodesSize; ++i) {
			auto value = generator(lattice(n, i),lattice(n + 1, i), lattice(n + 1, i + 1), deltaTime[n]);
			payoffAdjuster(value, lattice(n, i));
			lattice(n, i) = value;
		}
	}
	auto lastValue = generator(lattice(0, 0),lattice(1, 0), lattice(1, 1), deltaTime[0]);
	payoffAdjuster(lastValue, lattice(0, 0));
	lattice(0, 0) = lastValue;
}


template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff,typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime>::
_back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster, std::false_type) {

	const std::size_t lastIdx = lattice.timeDimension() - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = lattice.nodesAtIdx(n).size();
		for (auto i = 0; i < nodesSize; ++i) {
			auto value = generator(lattice(n, i),lattice(n + 1, i), lattice(n + 1, i + 1), deltaTime);
			payoffAdjuster(value, lattice(n, i));
			lattice(n, i) = value;
		}
	}
	auto lastValue = generator(lattice(0, 0),lattice(1, 0), lattice(1, 1), deltaTime);
	payoffAdjuster(lastValue, lattice(0, 0));
	lattice(0, 0) = lastValue;
}

template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
_back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, std::true_type) {

	LASSERT((deltaTime.size() + 1) == lattice.timeDimension(), "deltaTime is always one element shorter then time dimension");
	const std::size_t lastIdx = lattice.timeDimension() - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}
	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = lattice.nodesAtIdx(n).size();
		for (auto i = 0; i < nodesSize; ++i) {
			lattice(n, i) = generator(lattice(n, i), lattice(n + 1, i), lattice(n + 1, i + 1), lattice(n + 1, i + 2), deltaTime[n]);
		}
	}
	lattice(0, 0) = generator(lattice(0, 0),lattice(1, 0), lattice(1, 1), lattice(1, 2), deltaTime[0]);

}


template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
_back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, std::false_type) {
	const std::size_t lastIdx = lattice.timeDimension() - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}
	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = lattice.nodesAtIdx(n).size();
		for (auto i = 0; i < nodesSize; ++i) {
			lattice(n, i) = generator(lattice(n, i),lattice(n + 1, i), lattice(n + 1, i + 1), lattice(n + 1, i + 2), deltaTime);
		}
	}
	lattice(0, 0) = generator(lattice(0, 0),lattice(1, 0), lattice(1, 1), lattice(1, 2), deltaTime);

}

template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
_back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster, std::true_type) {

	LASSERT((deltaTime.size() + 1) == lattice.timeDimension(), "deltaTime is always one element shorter then time dimension");
	const std::size_t lastIdx = lattice.timeDimension() - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = lattice.nodesAtIdx(n).size();
		for (auto i = 0; i < nodesSize; ++i) {
			auto value  = generator(lattice(n, i),lattice(n + 1, i), lattice(n + 1, i + 1), lattice(n + 1, i + 2), deltaTime[n]);
			payoffAdjuster(value, lattice(n, i));
			lattice(n, i) = value;
		}
	}
	auto lastValue = generator(lattice(0, 0),lattice(1, 0), lattice(1, 1), lattice(1, 2), deltaTime[0]);
	payoffAdjuster(lastValue, lattice(0, 0));
	lattice(0, 0) = lastValue;
}


template<typename TimeAxis,
	typename DeltaTime>
	template<typename LatticeObject, typename Generator, typename Payoff, typename PayoffAdjuster>
void lattice_backward_traversals::BackwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime>::
_back_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Payoff &&payoff, PayoffAdjuster &&payoffAdjuster, std::false_type) {

	const std::size_t lastIdx = lattice.timeDimension() - 1;
	const std::size_t lastNodesSize = lattice.nodesAtIdx(lastIdx).size();

	for (auto i = 0; i < lastNodesSize; ++i) {
		lattice(lastIdx, i) = payoff(lattice(lastIdx, i));
	}

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = lattice.nodesAtIdx(n).size();
		for (auto i = 0; i < nodesSize; ++i) {
			auto value = generator(lattice(n, i),lattice(n + 1, i), lattice(n + 1, i + 1), lattice(n + 1, i + 2), deltaTime);
			payoffAdjuster(value, lattice(n, i));
			lattice(n, i) = value;
		}
	}
	auto lastValue = generator(lattice(0, 0),lattice(1, 0), lattice(1, 1), lattice(1, 2), deltaTime);
	payoffAdjuster(lastValue, lattice(0, 0));
	lattice(0, 0) = lastValue;
}


#endif ///_LATTICE_BACKWARD_TRAVERSALS