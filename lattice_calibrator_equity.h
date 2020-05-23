#pragma once
#if !defined(_LATTICE_CALIBRATOR_EQUITY)
#define _LATTICE_CALIBRATOR_EQUITY

#include"lattice_macros.h"
#include"lattice_types.h"
#include"lattice_calibrator_results.h"
#include<memory>

namespace lattice_calibrator_equity {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_utility::DeltaTimeHolder;
	using lattice_calibrator_results::CalibratorResults;

	// ==============================================================================
	// ====================== VolSurfaceHolder ======================================
	// ==============================================================================


	template<typename VolSurface>
	struct VolSurfaceHolder {
	private:
		static inline typename VolSurface::value_type const _maxVolatility_impl(std::size_t maturityIdx,std::size_t maturitySize,
			VolSurface const &volSurface, std::true_type) {

			using vol = typename VolSurface::value_type;
			std::size_t const volSize = volSurface.size();
			vol result = volSurface.at(maturityIdx);
			vol max = result;
			for (std::size_t i = maturityIdx + maturitySize; i < volSize; i += maturitySize) {
				max = volSurface.at(i);
				if (max > result) {
					result = max;
				}
			}
			return result;
		}

		static inline VolSurface const _maxVolatility_impl(std::size_t maturityIdx,std::size_t maturitySize,
			VolSurface const &volSurface, std::false_type) {
			return volSurface;
		}


	public:
		static inline auto const maxVolatility(std::size_t maturityIdx,std::size_t maturitySize, VolSurface const &volSurface) {
			return _maxVolatility_impl(maturityIdx, maturitySize,volSurface,std::is_compound<VolSurface>());
		}

	};



	// ==============================================================================
	// ====================== CalibratorEquity ======================================
	// ==============================================================================


	template<LatticeType Type,
		typename TimeAxis,
		typename DeltaTime,
		typename RiskFreeRate,
		typename OptionVolatilitySurface>
	struct CalibratorEquity {};


	template<typename TimeAxis,
			typename DeltaTime,
			typename RiskFreeRate,
			typename OptionSurface>
	struct CalibratorEquity<LatticeType::Binomial, TimeAxis, DeltaTime, RiskFreeRate, OptionVolatilitySurface> {


	};



	template<typename TimeAxis,
			typename DeltaTime,
			typename RiskFreeRate,
			typename OptionVolatilitySurface>
	struct CalibratorEquity<LatticeType::Trinomial, TimeAxis, DeltaTime, RiskFreeRate, OptionVolatilitySurface> {
	private:
		template<typename LatticeObject>
		static std::shared_ptr<CalibratorResults<AssetClass::Equity, LatticeObject>> const
			_statePriceLattice_impl(LatticeObject &stockPriceLattice, DeltaTime const &deltaTime,
				typename LatticeObject::Node_type const &apexPrice,
				OptionVolatilitySurface const &optionVolSurface);
	
	public:
		template<typename LatticeObject>
		static std::shared_ptr<CalibratorResults<AssetClass::Equity, LatticeObject>> const
			statePriceLattice(LatticeObject &stockPriceLattice, DeltaTime const &deltaTime,
				typename LatticeObject::Node_type const &apexPrice,
				OptionVolatilitySurface const &optionVolSurface) {
			return _statePriceLattice_impl(stockPriceLattice, deltaTime, apexPrice, optionVolSurface);
		}

	};
}





template<typename TimeAxis,
		typename DeltaTime,
		typename RiskFreeRate,
		typename OptionVolatilitySurface>
template<typename LatticeObject>
std::shared_ptr<lattice_calibrator_results::CalibratorResults<lattice_types::AssetClass::Equity,LatticeObject>> const
lattice_calibrator_equity::CalibratorEquity<lattice_types::LatticeType::Trinomial,TimeAxis,DeltaTime,RiskFreeRate, OptionVolatilitySurface>::
_statePriceLattice_impl(LatticeObject &stockPriceLattice, DeltaTime const &deltaTime, typename LatticeObject::Node_type const &apexPrice,
	OptionVolatilitySurface const &optionVolSurface) {

	// unpack the strikes:
	auto strikes = std::get<0>(optionVolSurface);
	// unpack the option maturities:
	auto maturities = std::get<1>(optionVolSurface);
	// unpack the option price surface:
	auto optionSurface = std::get<2>(optionVolSurface);
	// unpack the implied volatility surface:
	auto volSurface = std::get<3>(optionVolSurface);

	// Since strikes, maturities, and optionSurface must be containers: 
	std::size_t const strikeSize = strikes.size();
	std::size_t const maturitySize = maturities.size();
	std::size_t const optionSurfaceSize = optionSurface.size();
	// we must check that the dimension matches:
	LASSERT((strikeSize*maturitySize) == optionSurfaceSize,
		"option surface must have size strikes x maturities");
	// Get the trees time dimension:
	std::size_t const treeSize = stockPriceLattice.timeDimension();
	// Tree dimension must match the size of maturities:
	LASSERT((treeSize == maturitySize),
		"Time dimension of the stockPriceLattice must match the size of maturities");

	// typedef the vol surface holder:
	typedef VolSurfaceHolder<decltype(volSurface)> VS;
	// typedef deltaTime holder:
	typedef DeltaTimeHolder<DeltaTime> DT;
	// typedef the node type:
	typedef typename LatticeObject::Node_type Node;


	// create statePriceLattice from empty stockPriceLattice:
	LatticeObject statePricelattice(stockPriceLattice);
	std::size_t nodesSize{ 0 };
	std::size_t k{ 0 };
	Node dt{};
	Node vol{};
	Node fact{};

	// First populate stockPriceLattice:
	stockPriceLattice(0, 0) = apexPrice;
	for (std::size_t t = 1; t < treeSize; ++t) {
		dt = DT::deltaTime(t - 1, deltaTime);
		vol = VS::maxVolatility(t, maturitySize, volSurface);
		fact = vol * std::sqrt(3.0*dt);
		nodesSize = stockPriceLattice.nodesAtIdx(t).size();
		for (std::size_t l = 0; l < nodesSize; ++l) {
			stockPriceLattice(t, l) = apexPrice * std::exp((l - t)*fact);
		}
	}

	// populate statePriceLattice here:

}








#endif ///_LATTICE_CALIBRATOR_EQUITY