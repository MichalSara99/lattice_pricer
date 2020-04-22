#pragma once
#if !defined(_LATTICE_TRAVERSALS)
#define _LATTICE_TRAVERSALS

#include"lattice_types.h"
#include"lattice_macros.h"

namespace lattice_forward_traversals {

	using lattice_types::LatticeType;


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
			static void _traverse(LatticeObject &lattice, Generator &&generator,DeltaTime const &deltaTime, Node apex, std::true_type);
			//	This one is for non-compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex, std::false_type);

			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis,Node>const &dividendData,std::true_type);
			//	This one is for non-compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverse(LatticeObject &lattice, Generator &&generator,DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis, Node>const &dividendData, std::false_type);

		public:
			template<typename LatticeObject, typename Generator>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex) {
				_traverse(lattice, std::forward<Generator>(generator), deltaTime, apex, std::is_compound<DeltaTime>());
			}

			template<typename LatticeObject, typename Generator>
			static void traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis,Node> const &dividendData) {
				_traverse(lattice, std::forward<Generator>(generator), deltaTime, apex, dividendData,std::is_compound<DeltaTime>());
			}

	};




	//	LatticeType::Binomial specialization
	template<typename TimeAxis,
		typename DeltaTime,
		typename Node>
		struct ForwardTraversal<LatticeType::Trinomial, TimeAxis, DeltaTime, Node> {
		private:
			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverseNormal(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex, std::true_type);
			//	This one is for non-compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverseNormal(std::size_t timeIdx,LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex, std::false_type);

			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverseReverting(std::size_t timeIdx,LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex, std::true_type);
			//	This one is for non-compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverseReverting(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex, std::false_type);

			template<typename LatticeObject, typename Generator>
			static void _traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex);



			
			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverseNormal(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis, Node>const &dividendData, std::true_type);
			//	This one is for non-compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverseNormal(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis, Node>const &dividendData, std::false_type);

			//	This one is for compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverseReverting(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis, Node>const &dividendData, std::true_type);
			//	This one is for non-compound DeltaTime object
			template<typename LatticeObject, typename Generator>
			static void _traverseReverting(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
				std::map<TimeAxis, Node>const &dividendData, std::false_type);

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
	template<typename LatticeObject, typename Generator>
void lattice_forward_traversals::ForwardTraversal<lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime, Node>::
_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex, std::true_type) {

	LASSERT((deltaTime.size() + 1) == lattice.timeDimension(), "deltaTime is always one element shorter then time dimension");

	lattice(0, 0) = apex;
	std::tuple<Node, Node> tuple;
	for (std::size_t t = 1; t < lattice.timeDimension(); ++t) {
		for (std::size_t l = 0; l < lattice.nodesAtIdx(t - 1).size(); ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime[t - 1],l,t);
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
_traverse(LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex, std::false_type) {

	lattice(0, 0) = apex;
	std::tuple<Node, Node> tuple;
	for (std::size_t t = 1; t < lattice.timeDimension(); ++t) {
		for (std::size_t l = 0; l < lattice.nodesAtIdx(t - 1).size(); ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime,l,t);
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
	std::map<TimeAxis,Node>const &dividendData, std::true_type) {

	LASSERT((deltaTime.size() + 1) == lattice.timeDimension(), "deltaTime is always one element shorter then time dimension");
	
	lattice(0, 0) = apex;
	std::tuple<Node, Node> tuple;
	for (std::size_t t = 1; t < lattice.timeDimension(); ++t) {
		for (std::size_t l = 0; l < lattice.nodesAtIdx(t - 1).size(); ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime[t - 1],l,t);
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
		for (std::size_t l = 0; l < lattice.nodesAtIdx(t - 1).size(); ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime[t - 1],l,t);
			lattice(t, l) = factor * std::get<0>(tuple);
			lattice(t, l + 1) = factor * std::get<1>(tuple);
		}
	}
}


template<typename TimeAxis,
	typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_forward_traversals::ForwardTraversal<lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime, Node>::
_traverse(LatticeObject &lattice, Generator &&generator,  DeltaTime const &deltaTime, Node apex,
	std::map<TimeAxis, Node>const &dividendData, std::false_type) {

	lattice(0, 0) = apex;
	std::tuple<Node, Node> tuple;
	for (std::size_t t = 1; t < lattice.timeDimension(); ++t) {
		for (std::size_t l = 0; l < lattice.nodesAtIdx(t - 1).size(); ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime,l,t);
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
	TimeAxis time;
	Node factor{ 1.0 };
	for (std::size_t t = firstExIdx; t < lattice.timeDimension(); ++t) {
		time = lattice.timeAt(t);
		auto nextExItr = dividendData.find(time);
		factor = 1.0;
		if (nextExItr != last) {
			factor *= (factor - nextExItr->second);
		}
		for (std::size_t l = 0; l < lattice.nodesAtIdx(t - 1).size(); ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime,l,t);
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
_traverseNormal(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex, std::true_type) {

	LASSERT((deltaTime.size() + 1) == lattice.timeDimension(), "deltaTime is always one element shorter then time dimension");

	lattice(0, 0) = apex;
	std::size_t nodesSize{ 0 };
	std::tuple<Node, Node, Node> tuple;
	for (std::size_t t = 1; t < timeIdx; ++t) {
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		for (std::size_t l = 0; l < nodesSize; ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime[t - 1],l,t);
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
_traverseNormal(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex, std::false_type) {

	lattice(0, 0) = apex;
	std::size_t nodesSize{ 0 };
	std::tuple<Node, Node, Node> tuple;
	for (std::size_t t = 1; t < timeIdx; ++t) {
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		for (std::size_t l = 0; l < nodesSize; ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime,l,t);
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
_traverseReverting(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex, std::true_type) {

	LASSERT((deltaTime.size() + 1) == lattice.timeDimension(), "deltaTime is always one element shorter then time dimension");

	std::size_t const treeSize = lattice.timeDimension();
	std::size_t nodesSize{ 0 };
	std::tuple<Node, Node, Node> tuple;
	for (std::size_t t = timeIdx; t < treeSize; ++t) {
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		//for lowest node:
		tuple = generator(lattice(t - 1, 0), deltaTime[t - 1], 0, t,true);
		lattice(t, 0) = std::get<0>(tuple); //low
		lattice(t, 1) = std::get<1>(tuple); //mid
		lattice(t, 2) = std::get<2>(tuple); //high

		for (std::size_t l = 1; l < nodesSize - 1; ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime[t - 1], l, t);
			lattice(t, l - 1) = std::get<0>(tuple);
			lattice(t, l) = std::get<1>(tuple);
			lattice(t, l + 1) = std::get<2>(tuple);
		}
		//for highest node:
		tuple = generator(lattice(t - 1, nodesSize - 1), deltaTime[t - 1], nodesSize - 1, t, true);
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
_traverseReverting(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex, std::false_type) {

	std::size_t const treeSize = lattice.timeDimension();
	std::size_t nodesSize{ 0 };
	std::tuple<Node, Node, Node> tuple;
	for (std::size_t t = timeIdx; t < treeSize; ++t) {
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		//for lowest node:
		tuple = generator(lattice(t - 1, 0), deltaTime, 0, t, true);
		lattice(t, 0) = std::get<0>(tuple); //low
		lattice(t, 1) = std::get<1>(tuple); //mid
		lattice(t, 2) = std::get<2>(tuple); //high

		for (std::size_t l = 1; l < nodesSize - 1; ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime, l, t);
			lattice(t, l - 1) = std::get<0>(tuple);
			lattice(t, l) = std::get<1>(tuple);
			lattice(t, l + 1) = std::get<2>(tuple);
		}
		//for highest node:
		tuple = generator(lattice(t - 1, nodesSize - 1), deltaTime, nodesSize - 1, t, true);
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
		_traverseNormal(treeSize, lattice, std::forward<Generator>(generator), deltaTime, apex, std::is_compound<DeltaTime>());
	}
	else {
		// This trinomial tree does have reverting property:
		_traverseNormal(firstRevertIdx, lattice, std::forward<Generator>(generator), deltaTime, apex, std::is_compound<DeltaTime>());
		_traverseReverting(firstRevertIdx, lattice, std::forward<Generator>(generator), deltaTime, apex, std::is_compound<DeltaTime>());
	}
}





template<typename TimeAxis,
	typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_forward_traversals::ForwardTraversal<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime, Node>::
_traverseNormal(std::size_t timeIdx,LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
	std::map<TimeAxis, Node>const &dividendData, std::true_type) {

	LASSERT((deltaTime.size() + 1) == lattice.timeDimension(), "deltaTime is always one element shorter then time dimension");
	std::size_t nodesSize{ 0 };

	lattice(0, 0) = apex;
	std::tuple<Node, Node, Node> tuple;
	for (std::size_t t = 1; t < timeIdx; ++t) {
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		for (std::size_t l = 0; l < nodesSize; ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime[t - 1],l,t);
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
		for (std::size_t l = 0; l < nodesSize; ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime[t - 1],l,t);
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
_traverseNormal(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
	std::map<TimeAxis, Node>const &dividendData, std::false_type) {

	std::size_t nodesSize{ 0 };

	lattice(0, 0) = apex;
	std::tuple<Node, Node, Node> tuple;
	for (std::size_t t = 1; t < timeIdx; ++t) {
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		for (std::size_t l = 0; l < nodesSize; ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime,l,t);
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
	TimeAxis time;
	Node factor{ 1.0 };
	for (std::size_t t = firstExIdx; t < timeIdx; ++t) {
		time = lattice.timeAt(t);
		auto nextExItr = dividendData.find(time);
		factor = 1.0;
		if (nextExItr != last) {
			factor *= (factor - nextExItr->second);
		}
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		for (std::size_t l = 0; l < nodesSize; ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime,l,t);
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
	std::map<TimeAxis, Node>const &dividendData, std::true_type) {

	LASSERT((deltaTime.size() + 1) == lattice.timeDimension(), "deltaTime is always one element shorter then time dimension");
	std::size_t const treeSize = lattice.timeDimension();
	std::size_t nodesSize{ 0 };

	lattice(0, 0) = apex;
	std::tuple<Node, Node, Node> tuple;
	for (std::size_t t = timeIdx; t < treeSize; ++t) {
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		//for lowest node:
		tuple = generator(lattice(t - 1, 0), deltaTime[t - 1], 0, t, true);
		lattice(t, 0) = std::get<0>(tuple); //low
		lattice(t, 1) = std::get<1>(tuple); //mid
		lattice(t, 2) = std::get<2>(tuple); //high

		for (std::size_t l = 1; l < nodesSize - 1; ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime[t - 1], l, t);
			lattice(t, l - 1) = std::get<0>(tuple);
			lattice(t, l) = std::get<1>(tuple);
			lattice(t, l + 1) = std::get<2>(tuple);
		}
		//for highest node:
		tuple = generator(lattice(t - 1, nodesSize - 1), deltaTime[t - 1], nodesSize - 1, t, true);
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
		//for lowest node:
		tuple = generator(lattice(t - 1, 0), deltaTime[t - 1], 0, t, true);
		lattice(t, 0) = factor * std::get<0>(tuple); //low
		lattice(t, 1) = factor * std::get<1>(tuple); //mid
		lattice(t, 2) = factor * std::get<2>(tuple); //high

		for (std::size_t l = 1; l < nodesSize - 1; ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime[t - 1], l, t);
			lattice(t, l - 1) = factor * std::get<0>(tuple);
			lattice(t, l) = factor * std::get<1>(tuple);
			lattice(t, l + 1) = factor * std::get<2>(tuple);
		}
		//for highest node:
		tuple = generator(lattice(t - 1, nodesSize - 1), deltaTime[t - 1], nodesSize - 1, t, true);
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
_traverseReverting(std::size_t timeIdx, LatticeObject &lattice, Generator &&generator, DeltaTime const &deltaTime, Node apex,
	std::map<TimeAxis, Node>const &dividendData, std::false_type) {

	std::size_t const treeSize = lattice.timeDimension();
	std::size_t nodesSize{ 0 };

	lattice(0, 0) = apex;
	std::tuple<Node, Node, Node> tuple;
	for (std::size_t t = timeIdx; t < treeSize; ++t) {
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		//for lowest node:
		tuple = generator(lattice(t - 1, 0), deltaTime, 0, t, true);
		lattice(t, 0) = std::get<0>(tuple); //low
		lattice(t, 1) = std::get<1>(tuple); //mid
		lattice(t, 2) = std::get<2>(tuple); //high

		for (std::size_t l = 1; l < nodesSize - 1; ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime, l, t);
			lattice(t, l - 1) = std::get<0>(tuple);
			lattice(t, l) = std::get<1>(tuple);
			lattice(t, l + 1) = std::get<2>(tuple);
		}
		//for highest node:
		tuple = generator(lattice(t - 1, nodesSize - 1), deltaTime, nodesSize - 1, t, true);
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
	if ((firstExIdx >= timeIdx) || (firstExIdx < timeIdx))
		return;
	TimeAxis time;
	Node factor{ 1.0 };
	for (std::size_t t = firstExIdx; t < treeSize; ++t) {
		time = lattice.timeAt(t);
		auto nextExItr = dividendData.find(time);
		factor = 1.0;
		if (nextExItr != last) {
			factor *= (factor - nextExItr->second);
		}
		nodesSize = lattice.nodesAtIdx(t - 1).size();
		//for lowest node:
		tuple = generator(lattice(t - 1, 0), deltaTime, 0, t, true);
		lattice(t, 0) = factor * std::get<0>(tuple); //low
		lattice(t, 1) = factor * std::get<1>(tuple); //mid
		lattice(t, 2) = factor * std::get<2>(tuple); //high

		for (std::size_t l = 1; l < nodesSize - 1; ++l) {
			tuple = generator(lattice(t - 1, l), deltaTime, l, t);
			lattice(t, l - 1) = factor * std::get<0>(tuple);
			lattice(t, l) = factor * std::get<1>(tuple);
			lattice(t, l + 1) = factor * std::get<2>(tuple);
		}
		//for highest node:
		tuple = generator(lattice(t - 1, nodesSize - 1), deltaTime, nodesSize - 1, t, true);
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
		_traverseNormal(treeSize, lattice, std::forward<Generator>(generator), deltaTime, apex, dividendData, std::is_compound<DeltaTime>());
	}
	else {
		// This trinomial tree does have reverting property:
		_traverseNormal(firstRevertIdx, lattice, std::forward<Generator>(generator), deltaTime, apex, dividendData, std::is_compound<DeltaTime>());
		_traverseReverting(firstRevertIdx, lattice, std::forward<Generator>(generator), deltaTime, apex, dividendData, std::is_compound<DeltaTime>());
	}
}



#endif ///_LATTICE_FORWARD_TRAVERSALS