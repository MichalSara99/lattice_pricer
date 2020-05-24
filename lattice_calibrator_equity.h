#pragma once
#if !defined(_LATTICE_CALIBRATOR_EQUITY)
#define _LATTICE_CALIBRATOR_EQUITY

#include"lattice_macros.h"
#include"lattice_types.h"
#include"lattice_utility.h"
#include"lattice_calibrator_results.h"
#include<memory>

namespace lattice_calibrator_equity {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_utility::DeltaTimeHolder;
	using lattice_utility::LinearInterpolator;
	using lattice_calibrator_results::CalibratorResults;

	// ==============================================================================
	// ====================== OptionPriceLocator ====================================
	// ==============================================================================

	template<typename OptionPriceContainer>
	struct OptionPriceLocator {
	public:
		static inline std::pair<OptionPriceContainer, OptionPriceContainer>
			const callPutPrices(std::size_t maturityIdx, std::size_t maturitySize,
				OptionPriceContainer const &optionPriceContainer) {

			std::size_t const containerSize = optionPriceContainer.size();
			OptionPriceContainer calls(containerSize / maturitySize);
			OptionPriceContainer puts(containerSize / maturitySize);

			for (std::size_t i = maturityIdx; i < containerSize; i += maturitySize) {
				calls[i] = optionPriceContainer.at(i).first;
				puts[i] = optionPriceContainer.at(i).second;
			}

			return std::make_pair(calls, puts);
		}


	};




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
		typename OptionData>
	struct CalibratorEquity {};


	template<typename TimeAxis,
			typename DeltaTime,
			typename RiskFreeRate,
			typename OptionData>
	struct CalibratorEquity<LatticeType::Binomial, TimeAxis, DeltaTime, RiskFreeRate, OptionData> {


	};



	template<typename TimeAxis,
			typename DeltaTime,
			typename RiskFreeRate,
			typename OptionData>
	struct CalibratorEquity<LatticeType::Trinomial, TimeAxis, DeltaTime, RiskFreeRate, OptionData> {
	private:
		template<typename LatticeObject>
		static std::shared_ptr<CalibratorResults<AssetClass::Equity, LatticeObject>> const
			_statePriceLattice_impl(LatticeObject &stockPriceLattice, DeltaTime const &deltaTime,
				typename LatticeObject::Node_type const &apexPrice,
				OptionData const &optionData);
	
	public:
		template<typename LatticeObject>
		static std::shared_ptr<CalibratorResults<AssetClass::Equity, LatticeObject>> const
			statePriceLattice(LatticeObject &stockPriceLattice, DeltaTime const &deltaTime,
				typename LatticeObject::Node_type const &apexPrice,
				OptionData const &optionData) {
			return _statePriceLattice_impl(stockPriceLattice, deltaTime, apexPrice, optionData);
		}

	};
}





template<typename TimeAxis,
		typename DeltaTime,
		typename RiskFreeRate,
		typename OptionData>
template<typename LatticeObject>
std::shared_ptr<lattice_calibrator_results::CalibratorResults<lattice_types::AssetClass::Equity,LatticeObject>> const
lattice_calibrator_equity::CalibratorEquity<lattice_types::LatticeType::Trinomial,TimeAxis,DeltaTime,RiskFreeRate, OptionData>::
_statePriceLattice_impl(LatticeObject &stockPriceLattice, DeltaTime const &deltaTime, typename LatticeObject::Node_type const &apexPrice,
	OptionData const &optionData) {

	// unpack the strikes:
	auto strikes = std::get<0>(optionData);
	// unpack the option maturities:
	auto maturities = std::get<1>(optionData);
	// unpack the option price surface:
	auto optionSurface = std::get<2>(optionData);
	// unpack the implied volatility surface:
	auto volSurface = std::get<3>(optionData);

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

	// typedef the LinearInterpolator:
	typedef LinearInterpolator<decltype(strikes)> LERP;
	// typedef the OptionPriceLocator:
	typedef OptionPriceLocator<decltype(optionSurface)> OPL;
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
	Node lerpPrice{};

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

	// prepare pair prices of calls and puts:
	std::pair<decltype(optionSurface), decltype(optionSurface)> pairPrices;
	// populate statePriceLattice here:
	statePricelattice(0, 0) = 1.0;
	// initialize LinearInterpolator here:
	LERP lerp;
	for (std::size_t t = 1; t < treeSize; ++t) {
		nodesSize = statePricelattice.nodesAtIdx(t).size();
		pairPrices = OPL::callPutPrices(t, maturitySize, optionSurface);




	}



}








#endif ///_LATTICE_CALIBRATOR_EQUITY