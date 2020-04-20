#pragma once
#if! defined(_LATTICE_CALIBRATOR)
#define _LATTICE_CALIBRATOR

#include"lattice_types.h"
#include"lattice_macros.h"
#include"lattice_calibrator_results.h"
#include"lattice_calibrator_ir.h"
#include<memory>

namespace lattice_calibrator {

	using lattice_types::AssetClass;
	using lattice_types::LatticeType;
	using lattice_calibrator_results::CalibratorResults;
	using lattice_calibrator_ir::CalibratorIR;

	// ==============================================================================
	// =============================== Calibrator ===================================
	// ==============================================================================

	template<LatticeType Type,
			AssetClass AClass,
			typename TimeAxis,
			typename DeltaTime,
			typename DiscountCurve>
	class Calibrator {};



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



}


#endif ///_LATTICE_CALIBRATOR