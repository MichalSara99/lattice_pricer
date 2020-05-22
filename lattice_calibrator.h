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
	using lattice_calibrator_results::CalibratorResults;
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
		std::shared_ptr<CalibratorResults<AssetClass::InterestRate, LatticeObject>> const
			operator()(LatticeObject &rateLattice, Generator &&generator,
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
		std::shared_ptr<CalibratorResults<AssetClass::InterestRate, LatticeObject>> const
			operator()(LatticeObject &rateLattice, Generator &&generator,
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
			typename OptionSurface>
	class Calibrator<LatticeType::Binomial, AssetClass::Equity,
		TimeAxis, DeltaTime, RiskFreeRate,OptionSurface> {




	};


	template<typename TimeAxis,
		typename DeltaTime,
		typename RiskFreeRate,
		typename OptionSurface>
	class Calibrator<LatticeType::Trinomial, AssetClass::Equity,
		TimeAxis, DeltaTime, RiskFreeRate, OptionSurface> {
	private:
		OptionSurface optionSurface_;
		RiskFreeRate rate_;

	public:
		explicit Calibrator(OptionSurface const &optionSurface,RiskFreeRate const &rate)
			:rate_{ rate }, optionSurface_{ optionSurface } {}

		template<typename LatticeObject>
		std::shared_ptr<CalibratorResults<AssetClass::Equity,LatticeObject>> const
			operator()(LatticeObject &statePriceLattice, DeltaTime const &deltaTime)const {
			
			return 


		}




	};




}


#endif ///_LATTICE_CALIBRATOR