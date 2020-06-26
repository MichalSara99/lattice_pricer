#pragma once
#if !defined(_LATTICE_MULTIDIMENSIONAL)
#define _LATTICE_MULTIDIMENSIONAL

#include"lattice_types.h"
#include"lattice_macros.h"
#include"lattice_structure.h"
#include"lattice_miscellaneous.h"
#include<type_traits>
#include<array>

namespace lattice_multidimensional {

	using lattice_types::LatticeType;
	using lattice_structure::IndexedLattice;
	using lattice_structure::MeanRevertingLattice;
	using lattice_structure::Lattice;
	using lattice_structure::MeanRevertingIndexedLattice;
	using lattice_miscellaneous::MeanRevertingParams;

	// ==============================================================================
	// ============================= MultidimGeneralLattice =========================
	// ==============================================================================

	template<std::size_t LatticeDimension,
			LatticeType Type,
			typename Node,
			typename TimeAxis,
			typename NodeContainerType,
			template<LatticeType,typename...> typename LatticeObject,
			typename  = typename std::enable_if<(LatticeDimension > 1)>::type>
	class MultidimGeneralLattice {
	public:
		virtual ~MultidimGeneralLattice(){}

		std::size_t constexpr factors()const { return LatticeDimension; }
	
	};

	// ==============================================================================
	// ====================== MultidimMeanRevertingGeneralLattice ===================
	// ==============================================================================

	template<std::size_t LatticeDimension,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType,
		template<typename...> typename LatticeObject,
		typename = typename std::enable_if<(LatticeDimension > 1)>::type>
		class MultidimMeanRevertingGeneralLattice {
		public:
			virtual ~MultidimMeanRevertingGeneralLattice() {}

			std::size_t constexpr factors()const { return LatticeDimension; }

	};


	// ==============================================================================
	// ============================= MultidimIndexedLattice =========================
	// ==============================================================================


	template<std::size_t Dimension,
			LatticeType Type,
			typename Node>
	class MultidimIndexedLattice :public MultidimGeneralLattice<Dimension,Type, Node, std::size_t, std::vector<Node>, IndexedLattice> {
	protected:
		std::vector<IndexedLattice<Type, Node>> multiTree_;

	public:
		explicit MultidimIndexedLattice(std::size_t numberPeriods) {
			for (std::size_t t = 0; t < Dimension; ++t)
				multiTree_.emplace_back(std::move(IndexedLattice<Type, Node>{ numberPeriods }));
		}

		IndexedLattice<Type, Node> getFactor(std::size_t factorIdx) {
			LASSERT(factorIdx < Dimension, "factor index is out of range");
			return multiTree_[factorIdx];
		}

		Node const& operator()(std::size_t factorIdx, std::size_t timeIdx, std::size_t leafIdx)const {
			LASSERT(factorIdx < Dimension, "factor index is out of range");
			return multiTree_[factorIdx](timeIdx, leafIdx);
		}

		Node &operator()(std::size_t factorIdx, std::size_t timeIdx, std::size_t leafIdx) {
			LASSERT(factorIdx < Dimension, "factor index is out of range");
			return multiTree_[factorIdx](timeIdx, leafIdx);
		}

		Node const &at(std::size_t factorIdx, std::size_t timeIdx, std::size_t leafIdx)const {
			LASSERT(factorIdx < Dimension, "factor index is out of range");
			return multiTree_[factorIdx](timeIdx, leafIdx);
		}

		Node apex(std::size_t factorIdx)const {
			LASSERT(factorIdx < Dimension, "factor index is out of range");
			return multiTree_[factorIdx].apex();
		}
	};


	// ==============================================================================
	// =================================== MultidimLattice ==========================
	// ==============================================================================
	 
	template<std::size_t Dimension,
		LatticeType Type,
		typename Node,
		typename TimeAxis>
	class MultidimLattice :public MultidimGeneralLattice<Dimension,Type, Node, TimeAxis, std::vector<Node>,Lattice> {
	protected:
		std::vector<Lattice<Type, Node, TimeAxis>> multiTree_;

	public:
		explicit MultidimLattice(std::set<TimeAxis> const &fixingDatesSet) {
			for (std::size_t t = 0; t < Dimension; ++t)
				multiTree_.emplace_back(std::move(Lattice<Type, Node, TimeAxis>{ fixingDatesSet }));
		}

		explicit MultidimLattice(std::initializer_list<TimeAxis> const &fixingDates) {
			for (std::size_t t = 0; t < Dimension; ++t)
				multiTree_.emplace_back(std::move(Lattice<Type, Node, TimeAxis>{ fixingDates }));
		}

		Lattice<Type, Node, TimeAxis> getFactor(std::size_t factorIdx) {
			LASSERT(factorIdx < Dimension, "factor index is out of range");
			return multiTree_[factorIdx];
		}

		Node const& operator()(std::size_t factorIdx, std::size_t timeIdx, std::size_t leafIdx)const {
			LASSERT(factorIdx < Dimension, "factor index is out of range");
			return multiTree_[factorIdx](timeIdx, leafIdx);
		}

		Node &operator()(std::size_t factorIdx, std::size_t timeIdx, std::size_t leafIdx) {
			LASSERT(factorIdx < Dimension, "factor index is out of range");
			return multiTree_[factorIdx](timeIdx, leafIdx);
		}

		Node const &at(std::size_t factorIdx, std::size_t timeIdx, std::size_t leafIdx)const {
			LASSERT(factorIdx < Dimension, "factor index is out of range");
			return multiTree_[factorIdx](timeIdx, leafIdx);
		}

		Node apex(std::size_t factorIdx)const {
			LASSERT(factorIdx < Dimension, "factor index is out of range");
			return multiTree_[factorIdx].apex();
		}


	};


	// ==============================================================================
	// ===================== MultidimMeanRevertingIndexedLattice ====================
	// ==============================================================================


	template<std::size_t Dimension,
		typename Node>
		class MultidimMeanRevertingIndexedLattice :public MultidimMeanRevertingGeneralLattice<Dimension,Node, std::size_t,
												std::vector<Node>, MeanRevertingIndexedLattice> {
		protected:
			std::vector<MeanRevertingIndexedLattice<Node>> multiTree_;

		public:
			template<typename DeltaTime>
			explicit MultidimMeanRevertingIndexedLattice(std::size_t numberPeriods, MeanRevertingParams<Node> const &params, DeltaTime const &deltaTime) {
				for (std::size_t t = 0; t < Dimension; ++t)
					multiTree_.emplace_back(std::move(MeanRevertingIndexedLattice<Node>{ numberPeriods,params,deltaTime }));
			}

			MeanRevertingIndexedLattice<Node> getFactor(std::size_t factorIdx) {
				LASSERT(factorIdx < Dimension, "factor index is out of range");
				return multiTree_[factorIdx];
			}

			Node const& operator()(std::size_t factorIdx, std::size_t timeIdx, std::size_t leafIdx)const {
				LASSERT(factorIdx < Dimension, "factor index is out of range");
				return multiTree_[factorIdx](timeIdx, leafIdx);
			}

			Node &operator()(std::size_t factorIdx, std::size_t timeIdx, std::size_t leafIdx) {
				LASSERT(factorIdx < Dimension, "factor index is out of range");
				return multiTree_[factorIdx](timeIdx, leafIdx);
			}

			Node const &at(std::size_t factorIdx, std::size_t timeIdx, std::size_t leafIdx)const {
				LASSERT(factorIdx < Dimension, "factor index is out of range");
				return multiTree_[factorIdx](timeIdx, leafIdx);
			}

			Node apex(std::size_t factorIdx)const {
				LASSERT(factorIdx < Dimension, "factor index is out of range");
				return multiTree_[factorIdx].apex();
			}
	};


	// ==============================================================================
	// ===================== MultidimMeanRevertingLattice ===========================
	// ==============================================================================


	template<std::size_t Dimension,
		typename Node,
		typename TimeAxis>
		class MultidimMeanRevertingLattice :public MultidimMeanRevertingGeneralLattice<Dimension, Node, TimeAxis,
		std::vector<Node>, MeanRevertingLattice> {
		protected:
			std::vector<MeanRevertingLattice<Node,TimeAxis>> multiTree_;

		public:
			template<typename DeltaTime>
			explicit MultidimMeanRevertingLattice(std::set<TimeAxis> const &fixingDatesSet, MeanRevertingParams<Node> const &params, DeltaTime const &deltaTime) {
				for (std::size_t t = 0; t < Dimension; ++t)
					multiTree_.emplace_back(std::move(MeanRevertingLattice<Node,TimeAxis>{ numberPeriods, params, deltaTime }));
			}

			MeanRevertingLattice<Node, TimeAxis> getFactor(std::size_t factorIdx) {
				LASSERT(factorIdx < Dimension, "factor index is out of range");
				return multiTree_[factorIdx];
			}

			Node const& operator()(std::size_t factorIdx, std::size_t timeIdx, std::size_t leafIdx)const {
				LASSERT(factorIdx < Dimension, "factor index is out of range");
				return multiTree_[factorIdx](timeIdx, leafIdx);
			}

			Node &operator()(std::size_t factorIdx, std::size_t timeIdx, std::size_t leafIdx) {
				LASSERT(factorIdx < Dimension, "factor index is out of range");
				return multiTree_[factorIdx](timeIdx, leafIdx);
			}

			Node const &at(std::size_t factorIdx, std::size_t timeIdx, std::size_t leafIdx)const {
				LASSERT(factorIdx < Dimension, "factor index is out of range");
				return multiTree_[factorIdx](timeIdx, leafIdx);
			}

			Node apex(std::size_t factorIdx)const {
				LASSERT(factorIdx < Dimension, "factor index is out of range");
				return multiTree_[factorIdx].apex();
			}
	};


}






#endif ///_LATTICE_MULTIDIMENSIONAL