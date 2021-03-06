#pragma once
#if! defined(_LATTICE_CALIBRATOR)
#define _LATTICE_CALIBRATOR

#include"lattice_types.h"
#include"lattice_macros.h"
#include"lattice_calibrator_results.h"
#include"lattice_calibrator_ir.h"
#include"lattice_calibrator_equity.h"
#include<memory>

namespace lattice_calibrator {

	using lattice_types::AssetClass;
	using lattice_types::LatticeType;
	using lattice_calibrator_results::CalibratorIRResultsPtr;
	using lattice_calibrator_results::CalibratorBinomialEquityResultsPtr;
	using lattice_calibrator_results::CalibratorTrinomialEquityResultsPtr;
	using lattice_calibrator_ir::CalibratorIR;
	using lattice_calibrator_equity::CalibratorEquity;

	// ==============================================================================
	// =============================== Calibrator ===================================
	// ==============================================================================

	template<LatticeType Type,
			AssetClass AClass,
			typename TimeAxis,
			typename DeltaTime,
			typename ...Others>
	class Calibrator {};


	// ==============================================================================
	// ============= Calibrator Partial Specialization for InterestRate =============
	// ==============================================================================

	template<typename TimeAxis,
		typename DeltaTime,
		typename DiscountCurve>
	class Calibrator<LatticeType::Binomial,AssetClass::InterestRate,
						TimeAxis,DeltaTime, DiscountCurve> {
	private:
		DiscountCurve discountCurve_;

	public:
		explicit Calibrator(DiscountCurve const &discountCurve)
			:discountCurve_{ discountCurve }{}


		template<typename LatticeObject,typename Generator>
		CalibratorIRResultsPtr<LatticeObject> const
			operator()(LatticeObject &rateLattice, 
				Generator &&generator,
				DeltaTime const &deltaTime) const {
			LASSERT(rateLattice.type() == generator.latticeType(), "Mismatch between lattice types");
			LASSERT(generator.assetClass() == AssetClass::InterestRate, "Mismatch between asset classes");
			return CalibratorIR<LatticeType::Binomial, TimeAxis, DeltaTime, DiscountCurve>::
				calibrate(rateLattice, std::forward<Generator>(generator), deltaTime, this->discountCurve_);
		}
	


	};


	template<typename TimeAxis,
		typename DeltaTime,
		typename DiscountCurve>
	class Calibrator<LatticeType::Trinomial, AssetClass::InterestRate,
		TimeAxis, DeltaTime, DiscountCurve> {
	private:
		DiscountCurve discountCurve_;

	public:
		explicit Calibrator(DiscountCurve const &discountCurve)
			:discountCurve_{ discountCurve } {}


		template<typename LatticeObject, typename Generator>
		CalibratorIRResultsPtr<LatticeObject> const
			operator()(LatticeObject &rateLattice, 
				Generator &&generator,
				DeltaTime const &deltaTime) const {
			LASSERT(rateLattice.type() == generator.latticeType(), "Mismatch between lattice types");
			LASSERT(generator.assetClass() == AssetClass::InterestRate, "Mismatch between asset classes");
			return CalibratorIR<LatticeType::Trinomial, TimeAxis, DeltaTime, DiscountCurve>::
				calibrate(rateLattice, std::forward<Generator>(generator), deltaTime, this->discountCurve_);
		}

	};



	// ==============================================================================
	// ================= Calibrator Partial Specialization for Equity ===============
	// ==============================================================================

	template<typename TimeAxis,
			typename DeltaTime,
			typename RiskFreeRate,
			typename OptionData>
	class Calibrator<LatticeType::Binomial, AssetClass::Equity,
		TimeAxis, DeltaTime, RiskFreeRate, OptionData> {
	private:
		OptionData optionData_;
		RiskFreeRate rate_;

	public:
		explicit Calibrator(OptionData const &optionData,RiskFreeRate const &rate)
			:rate_{rate},optionData_{optionData}{}

		template<typename LatticeObject>
		CalibratorBinomialEquityResultsPtr<LatticeObject> const
			operator()(LatticeObject &stockPriceLattice, DeltaTime const &deltaTime,
				RiskFreeRate const &riskFreeRate, typename LatticeObject::Node_type const &apexPrice) {

			return lattice_calibrator_equity::CalibratorEquity<LatticeType::Binomial, TimeAxis,
				DeltaTime, RiskFreeRate, OptionData>::
				implyTree(stockPriceLattice, deltaTime, riskFreeRate, apexPrice, this->optionData_);
		}


	};


	template<typename TimeAxis,
		typename DeltaTime,
		typename RiskFreeRate,
		typename OptionData>
	class Calibrator<LatticeType::Trinomial, AssetClass::Equity,
		TimeAxis, DeltaTime, RiskFreeRate, OptionData> {
	private:
		OptionData optionData_;
		RiskFreeRate rate_;

	public:
		explicit Calibrator(OptionData const &optionData,RiskFreeRate const &rate)
			:rate_{ rate }, optionData_{ optionData } {}

		template<typename LatticeObject>
		CalibratorTrinomialEquityResultsPtr<LatticeObject> const
			operator()(LatticeObject &stockPriceLattice, DeltaTime const &deltaTime,
				typename LatticeObject::Node_type const &apexPrice,bool areCallPricesLiquid = true)const {
			
			return lattice_calibrator_equity::CalibratorEquity<LatticeType::Trinomial, TimeAxis,
				DeltaTime, RiskFreeRate, OptionData>::
				statePriceLattice(stockPriceLattice, deltaTime, apexPrice, this->optionData_, areCallPricesLiquid);

		}

		template<typename LatticeObject>
		CalibratorTrinomialEquityResultsPtr<LatticeObject> const
			operator()(LatticeObject const &statePriceLattice, LatticeObject const &stockPriceLattice,
				DeltaTime const &deltaTime, RiskFreeRate const &riskFreeRate, bool areCallPricesLiquid = true) {
			return lattice_calibrator::CalibratorEquity<LatticeType::Trinomial, TimeAxis,
				DeltaTime, RiskFreeRate, OptionData>::
				impliedProbability(statePriceLattice, stockPriceLattice, deltaTime, riskFreeRate, areCallPricesLiquid);

		}




	};




}


#endif ///_LATTICE_CALIBRATOR