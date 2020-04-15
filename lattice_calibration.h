#pragma once
#if !defined(_LATTICE_CALIBRATION)
#define _LATTICE_CALIBRATION

#include"lattice_types.h"


namespace lattice_calibration {

	using lattice_types::AssetClass;
	using lattice_types::LatticeType;

	// ==============================================================================
	// =============================== Calibrator ===================================
	// ==============================================================================

	template<LatticeType Type,
			AssetClass AClass,
			typename TimeAxis,
			typename DeltaTime,
			typename Node>
		class Calibrator {};



	template<typename TimeAxis,
		typename DeltaTime,
		typename Node>
	class Calibrator<LatticeType::Binomial,AssetClass::InterestRate,
						TimeAxis,DeltaTime,Node> {
	
	
	};


	template<typename TimeAxis,
		typename DeltaTime,
		typename Node>
	class Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		TimeAxis, DeltaTime, Node> {


	};





}





#endif ///_LATTICE_CALIBRATION