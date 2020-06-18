#pragma once
#if !defined(_LATTICE_CALIBRATOR_RESULTS)
#define _LATTICE_CALIBRATOR_RESULTS

#include<vector>
#include<tuple>
#include"lattice_types.h"
#include"lattice_macros.h"

namespace lattice_calibrator_results {


	using lattice_types::AssetClass;
	using lattice_types::LatticeType;


	// ==============================================================================
	// ============================ CalibratorResults ===============================
	// ==============================================================================

	template<AssetClass AClass,
			typename LatticeObject,
			typename ...Others>
	struct CalibratorResults {};


	// ==============================================================================
	// ======== CalibratorResults Partial Specialization for InterestRate ===========
	// ==============================================================================

	template<typename LatticeObject>
	struct CalibratorResults<AssetClass::InterestRate, LatticeObject> {
		typedef typename LatticeObject::Node_type Node;

		std::vector<std::tuple<Node, Node, Node, std::size_t>> 
			thetaOptimizers;

		LatticeObject 
			arrowDebreuLattice;

		explicit CalibratorResults(LatticeObject const &arrowDebreu,
				std::vector<std::tuple<Node, Node, Node, std::size_t>> const &theta)
			:thetaOptimizers{ theta }, arrowDebreuLattice{ arrowDebreu }{}

	};


	// ==============================================================================
	// ========== CalibratorResults Partial Specialization for Equity ===============
	// ==============================================================================


	template<typename LatticeObject,typename ProbabilityHolder>
	struct CalibratorResults<AssetClass::Equity, LatticeObject, ProbabilityHolder> {

		LatticeObject resultLattice;

		std::vector<std::vector<ProbabilityHolder>> impliedProbabilities;

		explicit CalibratorResults(LatticeObject const &results)
			:resultLattice{ results }{}

		explicit CalibratorResults(LatticeObject const &results,
			std::vector<std::vector<ProbabilityHolder>> const &impliedProbs)
			:resultLattice{ results },impliedProbabilities(impliedProbs) {}

	};

}





#endif ///_LATTICE_CALIBRATOR_RESULTS