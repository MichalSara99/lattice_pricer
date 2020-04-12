#pragma once
#if !defined(_LATTICE_MERGE)
#define _LATTICE_MERGE

#include<set>
#include<ppl.h>
#include"lattice_structure.h"
#include"lattice_traits.h"

namespace lattice_merge {


	using lattice_structure::IndexedLattice;
	using lattice_structure::Lattice;
	using lattice_traits::MergeTraits;


	template<typename TimeAxis>
	struct MergePolicy {
		template<typename LatticeObject,
			typename ...LatticeObjects,
			typename Traits = MergeTraits<sizeof...(LatticeObjects), LatticeObject::Node_type>>
			static auto construct(LatticeObject latticeObject, std::true_type, LatticeObjects... latticeObjects);
		
		template<typename LatticeObject,
			typename ...LatticeObjects,
			typename Traits = MergeTraits<sizeof...(LatticeObjects), LatticeObject::Node_type>>
			static auto construct(LatticeObject latticeObject, std::false_type, LatticeObjects... latticeObjects);

		template<typename ResultLattice,
			typename LatticeObject,
			typename ...LatticeObjects,
			typename Traits = MergeTraits<sizeof...(LatticeObjects), LatticeObject::Node_type>>
			static void runParallel(ResultLattice& resultLattice,LatticeObject latticeObject, LatticeObjects... latticeObjects);

		template<typename ResultLattice,
			typename LatticeObject,
			typename ...LatticeObjects,
			typename Traits = MergeTraits<sizeof...(LatticeObjects), LatticeObject::Node_type>>
			static void runSequential(ResultLattice& resultLattice,LatticeObject latticeObject, LatticeObjects... latticeObjects);

	};


}


template<typename TimeAxis>
template<typename LatticeObject,
	typename ...LatticeObjects,
	typename Traits = lattice_traits::MergeTraits<sizeof...(LatticeObjects), LatticeObject::Node_type>>
	auto lattice_merge::MergePolicy<TimeAxis>::construct(LatticeObject latticeObject, std::true_type, LatticeObjects... latticeObjects) {

	std::set<TimeAxis> dates;
	for (auto d : latticeObject.fixingDates()) {
		dates.emplace(d);
	}

	lattice_structure::Lattice< LatticeObject::type(),
		typename Traits::NodeHolder,
		TimeAxis> result{ dates };

	return result;
}

template<typename TimeAxis>
template<typename LatticeObject,
	typename ...LatticeObjects,
	typename Traits = lattice_traits::MergeTraits<sizeof...(LatticeObjects), LatticeObject::Node_type>>
	auto lattice_merge::MergePolicy<TimeAxis>::construct(LatticeObject latticeObject, std::false_type, LatticeObjects... latticeObjects) {

	lattice_structure::IndexedLattice< LatticeObject::type(),
		typename Traits::NodeHolder> result{ latticeObject.maxIndex() };
	return result;
}

template<typename TimeAxis>
template<typename ResultLattice,
	typename LatticeObject,
	typename ...LatticeObjects,
	typename Traits = lattice_traits::MergeTraits<sizeof...(LatticeObjects), LatticeObject::Node_type>>
	void lattice_merge::MergePolicy<TimeAxis>::runParallel(ResultLattice& resultLattice,LatticeObject latticeObject, LatticeObjects... latticeObjects) {
	
	const std::size_t lastIdx = latticeObject.timeDimension() - 1;
	std::size_t nodesSize{ 0 };

	concurrency::parallel_for(static_cast<std::size_t>(0), static_cast<std::size_t>(lastIdx + 1),
		[&](std::size_t t) {
		nodesSize = latticeObject.nodesAtIdx(t).size();
		for (auto i = 0; i < nodesSize; ++i) {
			resultLattice(t, i) = std::make_tuple(latticeObject(t, i), latticeObjects(t, i)...);
		}
	});
}

template<typename TimeAxis>
template<typename ResultLattice,
	typename LatticeObject,
	typename ...LatticeObjects,
	typename Traits = lattice_traits::MergeTraits<sizeof...(LatticeObjects), LatticeObject::Node_type>>
	void lattice_merge::MergePolicy<TimeAxis>::runSequential(ResultLattice& resultLattice, LatticeObject latticeObject, LatticeObjects... latticeObjects) {

	const std::size_t lastIdx = latticeObject.timeDimension() - 1;
	std::size_t nodesSize{ 0 };

	for (auto t = 0; t <= lastIdx; ++t) {
		nodesSize = latticeObject.nodesAtIdx(t).size();
		for (auto i = 0; i < nodesSize; ++i) {
			resultLattice(t, i) = std::make_tuple(latticeObject(t, i), latticeObjects(t, i)...);
		}
	}
}


#endif ///_LATTICE_MERGE