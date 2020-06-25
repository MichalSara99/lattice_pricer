#pragma once
#if !defined(_LATTICE_BUILDER)
#define _LATTICE_BUILDER

#include"lattice_types.h"
#include"lattice_structure.h"
#include"lattice_miscellaneous.h"

namespace lattice_builder {

	using lattice_types::LatticeType;
	using lattice_structure::IndexedLattice;
	using lattice_structure::MeanRevertingIndexedLattice;
	using lattice_structure::MeanRevertingLattice;
	using lattice_miscellaneous::MeanRevertingParams;
	using lattice_structure::Lattice;

	// ==============================================================================
	// =========================== LatticeBuilder ===================================
	// ==============================================================================

	template<std::size_t FactorCount>
	class LatticeBuilder {};


	// ==============================================================================
	// ============= LatticeBuilder Partial Spec for one-factor models ==============
	// ==============================================================================

	template<>
	class LatticeBuilder<1> {
	public:
		template<LatticeType Type,typename Node>
		static IndexedLattice<Type, Node> 
			createIndexedBasedLattice(std::size_t numberPeriods) {
			return IndexedLattice<Type, Node>(numberPeriods);
		}

		template<typename Node,typename DeltaTime>
		static MeanRevertingIndexedLattice<Node>
			createIndexedBasedMRLattice(std::size_t numberPeriods, MeanRevertingParams<Node> const &params, DeltaTime const &deltaTime) {
			return MeanRevertingIndexedLattice<Node>(numberPeriods, params, deltaTime);
		}

		template<typename Node>
		static IndexedLattice<LatticeType::Binomial, Node> 
			createIndexedBasedBinomialLattice(std::size_t numberPeriods) {
			return IndexedLattice<LatticeType::Binomial, Node>{numberPeriods};
		}

		template<typename Node>
		static IndexedLattice<LatticeType::Trinomial, Node> 
			createIndexedBasedTrinomialLattice(std::size_t numberPeriods) {
			return IndexedLattice<LatticeType::Trinomial, Node>{numberPeriods};
		}

		template<LatticeType Type,typename Node,typename TimeAxis>
		static Lattice<Type,Node,TimeAxis>
			createLattice(std::set<TimeAxis> const &fixingDatesSet) {
			return Lattice<Type, Node, TimeAxis>{fixingDatesSet};
		}

		template<typename Node,typename TimeAxis,typename DeltaTime>
		static MeanRevertingLattice<Node,TimeAxis>
			createMRLattice(std::set<TimeAxis> const &fixingDatesSet, MeanRevertingParams<Node> const &params, DeltaTime const &deltaTime) {
			return MeanRevertingLattice<Node, TimeAxis>(fixingDatesSet, params, deltaTime);
		}


		template<typename Node,typename TimeAxis>
		static Lattice<LatticeType::Binomial,Node,TimeAxis>
			createBinomialLattice(std::set<TimeAxis> const &fixingDatesSet) {
			return Lattice<LatticeType::Binomial, Node, TimeAxis>{fixingDatesSet};
		}

		template<typename Node,typename TimeAxis>
		static Lattice<LatticeType::Trinomial,Node,TimeAxis>
			createTrinomialLattice(std::set<TimeAxis> const &fixingDatesSet) {
			return Lattice<LatticeType::Trinomial, Node, TimeAxis>{std::move(fixingDatesSet)};
		}

		template<LatticeType Type,typename Node,typename TimeAxis,
				template<typename T,typename Alloc> typename Container,
				typename T = TimeAxis,
				typename Alloc = std::allocator<T>>
		static Lattice<Type,Node,TimeAxis>
			createLattice(Container<T,Alloc> const &fixindDatesContainer) {
			std::set<TimeAxis> fixings{ fixindDatesContainer.cbegin(),fixindDatesContainer.cend() };
			return Lattice<Type, Node, TimeAxis>{std::move(fixings)};
		}

		template<typename Node, typename TimeAxis, typename DeltaTime,
			template<typename T, typename Alloc> typename Container,
			typename T = TimeAxis,
			typename Alloc = std::allocator<T>>
			static MeanRevertingLattice<Node, TimeAxis>
			createMRLattice(Container<T,Alloc> const &fixingDatesContainer, MeanRevertingParams<Node> const &params, DeltaTime const &deltaTime) {
			std::set<TimeAxis> fixings(fixingDatesContainer.cbegin(), fixingDatesContainer.cend());
			return MeanRevertingLattice<Node, TimeAxis>(std::move(fixings), params, deltaTime);
		}

		template<typename Node, typename TimeAxis,
				template<typename T,typename Alloc> typename Container,
				typename T = TimeAxis,
				typename Alloc = std::allocator<T>>
		static Lattice<LatticeType::Binomial, Node, TimeAxis>
			createBinomialLattice(Container<T,Alloc> const &fixindDatesContainer) {
			std::set<TimeAxis> fixings{ fixindDatesContainer.cbegin(),fixindDatesContainer.cend() };
			return Lattice<LatticeType::Binomial, Node, TimeAxis>{std::move(fixings)};
		}

		template<typename Node, typename TimeAxis,
				template<typename T,typename Alloc> typename Container,
				typename T = TimeAxis,
				typename Alloc = std::allocator<T>>
		static Lattice<LatticeType::Trinomial, Node, TimeAxis>
			createTrinomialLattice(Container<T,Alloc> const &fixingDatesContainer) {
			std::set<TimeAxis> fixings{ fixingDatesContainer.cbegin(),fixingDatesContainer.cend() };
			return Lattice<LatticeType::Trinomial, Node, TimeAxis>{std::move(fixings)};
		}

	};




	// ==============================================================================
	// ============= LatticeBuilder Partial Spec for two-factor models ==============
	// ==============================================================================

	template<>
	class LatticeBuilder<2> {


	};






}













#endif //_LATTICE_BUILDER