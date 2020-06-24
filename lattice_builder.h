#pragma once
#if !defined(_LATTICE_BUILDER)
#define _LATTICE_BUILDER

#include"lattice_types.h"
#include"lattice_structure.h"

namespace lattice_builder {

	using lattice_types::LatticeType;
	using lattice_structure::IndexedLattice;
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
			return Lattice<Type, Node, TimeAxis>{std::move(fixingDatesSet)};
		}

		template<typename Node,typename TimeAxis>
		static Lattice<LatticeType::Binomial,Node,TimeAxis>
			createBinomialLattice(std::set<TimeAxis> const &fixingDatesSet) {
			return Lattice<LatticeType::Binomial, Node, TimeAxis>{std::move(fixingDatesSet)};
		}

		template<typename Node,typename TimeAxis>
		static Lattice<LatticeType::Trinomial,Node,TimeAxis>
			createTrinomialLattice(std::set<TimeAxis> const &fixingDatesSet) {
			return Lattice<LatticeType::Trinomial, Node, TimeAxis>{std::move(fixingDatesSet)};
		}

		template<LatticeType Type,typename Node,typename TimeAxis>
		static Lattice<Type,Node,TimeAxis>
			createLattice(std::vector<TimeAxis> const &fixindDatesVec) {
			std::set<TimeAxis> fixings{ fixindDatesVec.cbegin(),fixindDatesVec.cend() };
			return Lattice<Type, Node, TimeAxis>{std::move(fixings)};
		}

		template<typename Node, typename TimeAxis>
		static Lattice<LatticeType::Binomial, Node, TimeAxis>
			createBinomialLattice(std::vector<TimeAxis> const &fixindDatesVec) {
			std::set<TimeAxis> fixings{ fixindDatesVec.cbegin(),fixindDatesVec.cend() };
			return Lattice<LatticeType::Binomial, Node, TimeAxis>{std::move(fixings)};
		}

		template<typename Node, typename TimeAxis>
		static Lattice<LatticeType::Trinomial, Node, TimeAxis>
			createTrinomialLattice(std::vector<TimeAxis> const &fixingDatesVec) {
			std::set<TimeAxis> fixings{ fixindDatesVec.cbegin(),fixindDatesVec.cend() };
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