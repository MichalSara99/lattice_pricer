#pragma once
#if !defined(_LATTICE_TRAVERSALS)
#define _LATTICE_TRAVERSALS

#include"lattice_types.h"
#include"lattice_macros.h"
#include"lattice_utility.h"

namespace lattice_forward_traversals {

	using lattice_types::LatticeType;
	using lattice_utility::DeltaTimeHolder;

	// ==============================================================================
	// ==================== Forward Traversal Algorithms ============================
	// ==============================================================================


	template<LatticeType Type,
		typename TimeAxis,
		typename DeltaTime,
		typename Node>
		struct ForwardTraversal {};

	//	LatticeType::Binomial specialization
	template<typename TimeAxis,
		typename DeltaTime,
		typename Node>
		struct ForwardTraversal<LatticeType::Binomial, TimeAxis, DeltaTime, Node> {
		private:
			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverse(LatticeObject &lattice, Generator &&generator,DeltaTime const &deltaTime, Node apex);

			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis,Node>const &dividendData);


		public:
			template<typename LatticeObject, typename Generator>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex) {
				_traverse(lattice, std::forward<Generator>(generator), deltaTime, apex);
			}

			template<typename LatticeObject, typename Generator>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis,Node> const &dividendData) {
				_traverse(lattice, std::forward<Generator>(generator), deltaTime, apex, dividendData);
			}
	};

	//	LatticeType::TwoVariableBinomial specialization
	template<typename TimeAxis,
		typename DeltaTime,
		typename Node>
		struct ForwardTraversal<LatticeType::TwoVariableBinomial, TimeAxis, DeltaTime, Node> {
		private:
			//	This one is for compound DeltaTime object
			template<typename MultidimLatticeObject, typename Generator>
			static void _traverse(MultidimLatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, std::pair<Node,Node> const &apex);

			//	This one is for compound DeltaTime object
			template<typename MultidimLatticeObject, typename Generator>
			static void _traverse(MultidimLatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, std::pair<Node, Node> const &apex,
				std::pair<std::map<TimeAxis, Node>, std::map<TimeAxis, Node>> const &dividendData);


		public:
			template<typename MultidimLatticeObject, typename Generator>
			static void traverse(MultidimLatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, std::pair<Node, Node> const &apex) {
				_traverse(lattice, std::forward<Generator>(generator), deltaTime, apex);
			}

			template<typename MultidimLatticeObject, typename Generator>
			static void traverse(MultidimLatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, std::pair<Node, Node> const &apex,
				std::pair<std::map<TimeAxis, Node>, std::map<TimeAxis, Node>> const &dividendData) {
				_traverse(lattice, std::forward<Generator>(generator), deltaTime, apex, dividendData);
			}
	};


	//	LatticeType::Trinomial specialization
	template<typename TimeAxis,
		typename DeltaTime,
		typename Node>
		struct ForwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime, Node> {
		private:
			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverseNormal(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex);

			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverseReverting(std::size_t timeIdx,LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex);


			template<typename LatticeObject, typename Generator>
			static void _traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex);


			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverseNormal(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis, Node>const &dividendData);

			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverseReverting(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis, Node>const &dividendData);


			template<typename LatticeObject, typename Generator>
			static void _traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis, Node>const &dividendData);


		public:
			template<typename LatticeObject, typename Generator>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex) {
				_traverse(lattice, std::forward<Generator>(generator), deltaTime, apex);
			}

			template<typename LatticeObject, typename Generator>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis, Node> const &dividendData) {
				_traverse(lattice, std::forward<Generator>(generator), deltaTime, apex, dividendData);
			}
	};



}


template<typename TimeAxis,
	typename DeltaTime,
	typename Node>
template<typename MultidimLatticeObject, typename Generator>
void lattice_forward_traversals::ForwardTraversal<lattice_types::LatticeType::TwoVariableBinomial,TimeAxis,DeltaTime,Node>::
_traverse(MultidimLatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, std::pair<Node, Node> const &apex) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	Node dt{};
	std::size_t nodesSize{};
	std::size_t const treeSize = lattice.getFactor(0).timeDimension();
	auto gens = generator.forwardGenerator();
	auto generator0 = gens.first;
	auto generator1 = gens.second;

	lattice(0, 0, 0) = apex.first;
	lattice(1, 0, 0) = apex.second;
	std::tuple<Node, Node> tuple;

	for (std::size_t t = 1; t < treeSize; ++t) {
		dt = DT::deltaTime(t - 1, deltaTime);
		nodesSize = lattice.getFactor(0).nodesAtIdx(t - 1).size();
		for (std::size_t l = 0; l < nodesSize; ++l) {
			tuple = generator0(lattice(0, t - 1, l), dt, l, t, false);
			lattice(0, t, l) = std::get<0>(tuple);
			lattice(0, t, l + 1) = std::get<1>(tuple);
			tuple = generator1(lattice(1, t - 1, l), dt, l, t, false);
			lattice(1, t, l) = std::get<0>(tuple);
			lattice(1, t, l + 1) = std::get<1>(tuple);
		}
	}
}

template<typename TimeAxis,
	typename DeltaTime,
	typename Node>
	template<typename MultidimLatticeObject, typename Generator>
void lattice_forward_traversals::ForwardTraversal<lattice_types::LatticeType::TwoVariableBinomial, TimeAxis, DeltaTime, Node>::
_traverse(MultidimLatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, std::pair<Node, Node> const &apex,
	std::pair<std::map<TimeAxis, Node>, std::map<TimeAxis, Node>> const &dividendData) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	Node dt{};
	std::size_t nodesSize{};
	std::size_t const treeSize = lattice.getFactor(0).timeDimension();
	auto gens = generator.forwardGenerator();
	auto generator0 = gens.first;
	auto generator1 = gens.second;

	lattice(0, 0, 0) = apex.first;
	lattice(1, 0, 0) = apex.second;
	std::tuple<Node, Node> tuple;

	for (std::size_t t = 1; t < treeSize; ++t) {
		dt = DT::deltaTime(t - 1, deltaTime);
		nodesSize = lattice.getFactor(0).nodesAtIdx(t - 1).size();
		for (std::size_t l = 0; l < nodesSize; ++l) {
			tuple = generator0(lattice(0, t - 1, l), dt, l, t, false);
			lattice(0, t, l) = std::get<0>(tuple);
			lattice(0, t, l + 1) = std::get<1>(tuple);
			tuple = generator1(lattice(1, t - 1, l), dt, l, t, false);
			lattice(1, t, l) = std::get<0>(tuple);
			lattice(1, t, l + 1) = std::get<1>(tuple);
		}
	}
	// Adjust the tree for discretely paying divdiends:
	auto divData = dividendData.first;
	auto curr = divData.begin();
	auto last = divData.end();
	if (curr == last)
		return;
	std::size_t firstExIdx = lattice.getFactor(0).indexOf(curr->first);
	if (firstExIdx >= lattice.getFactor(0).timeDimension())
		return;

	Node factor{ 1.0 };
	for (std::size_t t = firstExIdx; t < treeSize; ++t) {
		auto nextExItr = divData.find(lattice.getFactor(0).timeAt(t));
		factor = 1.0;
		if (nextExItr != last) {
			factor *= (factor - nextExItr->second);
		}
		dt = DT::deltaTime(t - 1, deltaTime);
		nodesSize = lattice.getFactor(0).nodesAtIdx(t - 1).size();
		for (std::size_t l = 0; l < nodesSize; ++l) {
			tuple = generator0(lattice(0, t - 1, l), dt, l, t, false);
			lattice(0,t, l) = factor * std::get<0>(tuple);
			lattice(0,t, l + 1) = factor * std::get<1>(tuple);
		}
	}

	// Adjust the tree for discretely paying divdiends:
	auto divData = dividendData.second;
	auto curr = divData.begin();
	auto last = divData.end();
	if (curr == last)
		return;
	std::size_t firstExIdx = lattice.getFactor(1).indexOf(curr->first);
	if (firstExIdx >= lattice.getFactor(1).timeDimension())
		return;

	Node factor{ 1.0 };
	for (std::size_t t = firstExIdx; t < treeSize; ++t) {
		auto nextExItr = divData.find(lattice.getFactor(1).timeAt(t));
		factor = 1.0;
		if (nextExItr != last) {
			factor *= (factor - nextExItr->second);
		}
		dt = DT::deltaTime(t - 1, deltaTime);
		nodesSize = lattice.getFactor(1).nodesAtIdx(t - 1).size();
			for (std::size_t l = 0; l < nodesSize; ++l) {
				tuple = generator1(lattice(1, t - 1, l), dt, l, t, false);
				lattice(1, t, l) = factor * std::get<0>(tuple);
				lattice(1, t, l + 1) = factor * std::get<1>(tuple);
			}
	}
}


template<typename TimeAxis,
	typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_forward_traversals::ForwardTraversal<lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime, Node>::
_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	Node dt{};

	lattice(0, 0) = apex;
	std::tuple<Node, Node> tuple;
	for (std::size_t t = 1; t < lattice.timeDimension(); ++t) {
		dt = DT::deltaTime(t - 1, deltaTime);
		for (std::size_t l = 0; l < lattice.nodesAtIdx(t - 1).size(); ++l) {
			tuple = generator(lattice(t - 1, l), dt,l,t, false);
			lattice(t, l) = std::get<0>(tuple);
			lattice(t, l + 1) = std::get<1>(tuple);
		}
	}
}





template<typename TimeAxis,
	typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_forward_traversals::ForwardTraversal<lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime, Node>::
_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
	std::map<TimeAxis,Node>const &dividendData) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	Node dt{};

	lattice(0, 0) = apex;
	std::tuple<Node, Node> tuple;
	for (std::size_t t = 1; t < lattice.timeDimension(); ++t) {
		dt = DT::deltaTime(t - 1, deltaTime);
		for (std::size_t l = 0; l < lattice.nodesAtIdx(t - 1).size(); ++l) {
			tuple = generator(lattice(t - 1, l), dt,l,t, false);
			lattice(t, l) = std::get<0>(tuple);
			lattice(t, l + 1) = std::get<1>(tuple);
		}
	}
	// Adjust the tree for discretely paying divdiends: 
	auto curr = dividendData.begin();
	auto last = dividendData.end();
	if (curr == last)
		return;
	std::size_t firstExIdx = lattice.indexOf(curr->first);
	if (firstExIdx >= lattice.timeDimension())
		return;

	Node factor{ 1.0 };
	for (std::size_t t = firstExIdx; t < lattice.timeDimension(); ++t) {
		auto nextExItr = dividendData.find(lattice.timeAt(t));
		factor = 1.0;
		if (nextExItr != last) {
			factor *= (factor - nextExItr->second);
		}
		dt = DT::deltaTime(t - 1, deltaTime);
		for (std::size_t l = 0; l < lattice.nodesAtIdx(t - 1).size(); ++l) {
			tuple = generator(lattice(t - 1, l), dt,l,t, false);
			lattice(t, l) = factor * std::get<0>(tuple);
			lattice(t, l + 1) = factor * std::get<1>(tuple);
		}
	}
}






template<typename TimeAxis,
	typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_forward_traversals::ForwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime, Node>::
_traverseNormal(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	Node dt{};

	lattice(0, 0) = apex;
	std::size_t nodesSize{ 0 };
	std::tuple<Node, Node, Node> tuple;
	for (std::size_t t = 1; t < timeIdx; ++t) {
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		dt = DT::deltaTime(t - 1, deltaTime);
		for (std::size_t l = 0; l < nodesSize; ++l) {
			tuple = generator(lattice(t - 1, l), dt,l,t, false);
			lattice(t, l) = std::get<0>(tuple);
			lattice(t, l + 1) = std::get<1>(tuple);
			lattice(t, l + 2) = std::get<2>(tuple);
		}
	}

}



template<typename TimeAxis,
	typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_forward_traversals::ForwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime, Node>::
_traverseReverting(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	Node dt{};

	std::size_t const treeSize = lattice.timeDimension();
	std::size_t nodesSize{ 0 };
	std::tuple<Node, Node, Node> tuple;
	for (std::size_t t = timeIdx; t < treeSize; ++t) {
		dt = DT::deltaTime(t - 1, deltaTime);
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		//for lowest node:
		tuple = generator(lattice(t - 1, 0), dt, 0, t,true);
		lattice(t, 0) = std::get<0>(tuple); //low
		lattice(t, 1) = std::get<1>(tuple); //mid
		lattice(t, 2) = std::get<2>(tuple); //high

		for (std::size_t l = 1; l < nodesSize - 1; ++l) {
			tuple = generator(lattice(t - 1, l), dt, l, t, false);
			lattice(t, l - 1) = std::get<0>(tuple);
			lattice(t, l) = std::get<1>(tuple);
			lattice(t, l + 1) = std::get<2>(tuple);
		}
		//for highest node:
		tuple = generator(lattice(t - 1, nodesSize - 1), dt, nodesSize - 1, t, true);
		lattice(t, nodesSize - 1 - 2) = std::get<0>(tuple);
		lattice(t, nodesSize - 1 - 1) = std::get<1>(tuple);
		lattice(t, nodesSize - 1) = std::get<2>(tuple);
	}
}



template<typename TimeAxis,
	typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_forward_traversals::ForwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime, Node>::
_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex) {

	const std::size_t firstRevertIdx = lattice.firstRevertingIdx();
	const std::size_t treeSize = lattice.timeDimension();

	if (firstRevertIdx == 0) {
		// This trinomial tree does not have reverting property:
		_traverseNormal(treeSize, lattice, std::forward<Generator>(generator), deltaTime, apex);
	}
	else {
		// This trinomial tree does have reverting property:
		_traverseNormal(firstRevertIdx, lattice, std::forward<Generator>(generator), deltaTime, apex);
		_traverseReverting(firstRevertIdx, lattice, std::forward<Generator>(generator), deltaTime, apex);
	}
}





template<typename TimeAxis,
	typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_forward_traversals::ForwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime, Node>::
_traverseNormal(std::size_t timeIdx,LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
	std::map<TimeAxis, Node>const &dividendData) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	Node dt{};
	std::size_t nodesSize{ 0 };

	lattice(0, 0) = apex;
	std::tuple<Node, Node, Node> tuple;
	for (std::size_t t = 1; t < timeIdx; ++t) {
		dt = DT::deltaTime(t - 1, deltaTime);
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		for (std::size_t l = 0; l < nodesSize; ++l) {
			tuple = generator(lattice(t - 1, l), dt,l,t, false);
			lattice(t, l) = std::get<0>(tuple);
			lattice(t, l + 1) = std::get<1>(tuple);
			lattice(t, l + 2) = std::get<2>(tuple);
		}
	}
	// Adjust the tree for discretely paying divdiends: 
	auto curr = dividendData.begin();
	auto last = dividendData.end();
	if (curr == last)
		return;
	std::size_t firstExIdx = lattice.indexOf(curr->first);
	if (firstExIdx >= timeIdx)
		return;

	Node factor{ 1.0 };
	for (std::size_t t = firstExIdx; t < timeIdx; ++t) {
		auto nextExItr = dividendData.find(lattice.timeAt(t));
		factor = 1.0;
		if (nextExItr != last) {
			factor *= (factor - nextExItr->second);
		}
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		dt = DT::deltaTime(t - 1, deltaTime);
		for (std::size_t l = 0; l < nodesSize; ++l) {
			tuple = generator(lattice(t - 1, l), dt,l,t, false);
			lattice(t, l) = factor * std::get<0>(tuple);
			lattice(t, l + 1) = factor * std::get<1>(tuple);
			lattice(t, l + 2) = factor * std::get<2>(tuple);
		}
	}
}




template<typename TimeAxis,
	typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_forward_traversals::ForwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime, Node>::
_traverseReverting(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
	std::map<TimeAxis, Node>const &dividendData) {

	typedef DeltaTimeHolder<DeltaTime> DT;
	Node dt{};
	std::size_t const treeSize = lattice.timeDimension();
	std::size_t nodesSize{ 0 };

	lattice(0, 0) = apex;
	std::tuple<Node, Node, Node> tuple;
	for (std::size_t t = timeIdx; t < treeSize; ++t) {
		dt = DT::deltaTime(t - 1, deltaTime);
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		//for lowest node:
		tuple = generator(lattice(t - 1, 0), dt, 0, t, true);
		lattice(t, 0) = std::get<0>(tuple); //low
		lattice(t, 1) = std::get<1>(tuple); //mid
		lattice(t, 2) = std::get<2>(tuple); //high

		for (std::size_t l = 1; l < nodesSize - 1; ++l) {
			tuple = generator(lattice(t - 1, l), dt, l, t, false);
			lattice(t, l - 1) = std::get<0>(tuple);
			lattice(t, l) = std::get<1>(tuple);
			lattice(t, l + 1) = std::get<2>(tuple);
		}
		//for highest node:
		tuple = generator(lattice(t - 1, nodesSize - 1), dt, nodesSize - 1, t, true);
		lattice(t, nodesSize - 1 - 2) = std::get<0>(tuple);
		lattice(t, nodesSize - 1 - 1) = std::get<1>(tuple);
		lattice(t, nodesSize - 1) = std::get<2>(tuple);
	}

	// Adjust the tree for discretely paying divdiends: 
	auto curr = dividendData.begin();
	auto last = dividendData.end();
	if (curr == last)
		return;
	std::size_t firstExIdx = lattice.indexOf(curr->first);
	if ((firstExIdx >= treeSize) || (firstExIdx < timeIdx))
		return;

	Node factor{ 1.0 };
	for (std::size_t t = firstExIdx; t < treeSize; ++t) {
		auto nextExItr = dividendData.find(lattice.timeAt(t));
		factor = 1.0;
		if (nextExItr != last) {
			factor *= (factor - nextExItr->second);
		}
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		dt = DT::deltaTime(t - 1, deltaTime);
		//for lowest node:
		tuple = generator(lattice(t - 1, 0), dt, 0, t, true);
		lattice(t, 0) = factor * std::get<0>(tuple); //low
		lattice(t, 1) = factor * std::get<1>(tuple); //mid
		lattice(t, 2) = factor * std::get<2>(tuple); //high

		for (std::size_t l = 1; l < nodesSize - 1; ++l) {
			tuple = generator(lattice(t - 1, l), dt, l, t, false);
			lattice(t, l - 1) = factor * std::get<0>(tuple);
			lattice(t, l) = factor * std::get<1>(tuple);
			lattice(t, l + 1) = factor * std::get<2>(tuple);
		}
		//for highest node:
		tuple = generator(lattice(t - 1, nodesSize - 1), dt, nodesSize - 1, t, true);
		lattice(t, nodesSize - 1 - 2) = factor * std::get<0>(tuple);
		lattice(t, nodesSize - 1 - 1) = factor * std::get<1>(tuple);
		lattice(t, nodesSize - 1) = factor * std::get<2>(tuple);

	}
}


template<typename TimeAxis,
	typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_forward_traversals::ForwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime, Node>::
_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex, std::map<TimeAxis, Node>const &dividendData) {

	const std::size_t firstRevertIdx = lattice.firstRevertingIdx();
	const std::size_t treeSize = lattice.timeDimension();

	if (firstRevertIdx == 0) {
		// This trinomial tree does not have reverting property:
		_traverseNormal(treeSize, lattice, std::forward<Generator>(generator), deltaTime, apex, dividendData);
	}
	else {
		// This trinomial tree does have reverting property:
		_traverseNormal(firstRevertIdx, lattice, std::forward<Generator>(generator), deltaTime, apex, dividendData);
		_traverseReverting(firstRevertIdx, lattice, std::forward<Generator>(generator), deltaTime, apex, dividendData);
	}
}



#endif ///_LATTICE_FORWARD_TRAVERSALS