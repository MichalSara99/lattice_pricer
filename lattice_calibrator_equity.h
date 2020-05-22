#pragma once
#if !defined(_LATTICE_CALIBRATOR_EQUITY)
#define _LATTICE_CALIBRATOR_EQUITY


#include"lattice_types.h"
#include"lattice_calibrator_results.h"
#include<memory>

namespace lattice_calibrator_equity {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_calibrator_results::CalibratorResults;


	// ==============================================================================
	// ====================== CalibratorEquity ======================================
	// ==============================================================================


	template<LatticeType Type,
		typename TimeAxis,
		typename DeltaTime,
		typename RiskFreeRate,
		typename OptionSurface>
	struct CalibratorEquity {};


	template<typename TimeAxis,
			typename DeltaTime,
			typename RiskFreeRate,
			typename OptionSurface>
	struct CalibratorEquity<LatticeType::Binomial, TimeAxis, DeltaTime, RiskFreeRate, OptionSurface> {


	};



	template<typename TimeAxis,
			typename DeltaTime,
			typename RiskFreeRate,
			typename OptionSurface>
	struct CalibratorEquity<LatticeType::Trinomial, TimeAxis, DeltaTime, RiskFreeRate, OptionSurface> {
	private:
		template<typename LatticeObject>
		static std::shared_ptr<CalibratorResults<AssetClass::Equity, LatticeObject>> const
			_statePriceLattice_impl(LatticeObject &statePriceLattice, DeltaTime const &deltaTime,
				RiskFreeRate const &riskFreeRate, OptionSurface const &optionSurface);
	
	public:
		template<typename LatticeObject>
		static std::shared_ptr<CalibratorResults<AssetClass::Equity, LatticeObject>> const
			statePriceLattice(LatticeObject &statePriceLattice, DeltaTime const &deltaTime,
				RiskFreeRate const &riskFreeRate, OptionSurface const &optionSurface) {
			_statePriceLattice_impl(statePriceLattice, deltaTime, riskFreeRate, optionSurface);
		}



	};



}





#endif ///_LATTICE_CALIBRATOR_EQUITY