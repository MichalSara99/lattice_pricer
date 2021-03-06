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
	using lattice_utility::RiskFreeRateHolder;
	using lattice_utility::LinearInterpolator;
	using lattice_utility::probFloorCapper;
	using lattice_calibrator_results::CalibratorBinomialEquityResultsPtr;
	using lattice_calibrator_results::CalibratorBinomialEquityResultsT;
	using lattice_calibrator_results::CalibratorTrinomialEquityResultsPtr;
	using lattice_calibrator_results::CalibratorTrinomialEquityResultsT;


	// ==============================================================================
	// ====================== OptionPriceLocator ====================================
	// ==============================================================================

	template<typename Container>
	struct OptionPriceLocator {
	public:
		template<typename OptionPriceContainer>
		static inline std::pair<Container, Container>
			const callPutPrices(std::size_t maturityIdx, std::size_t maturitySize,
				OptionPriceContainer const &optionPriceContainer) {

			std::size_t const containerSize = optionPriceContainer.size();
			Container calls;
			calls.reserve(containerSize / maturitySize);
			Container puts;
			puts.reserve(containerSize / maturitySize);

			for (std::size_t i = maturityIdx; i < containerSize; i += maturitySize) {
				calls.push_back(optionPriceContainer.at(i).first);
				puts.push_back(optionPriceContainer.at(i).second);
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
		static auto const _maxVolatility_impl(std::size_t maturityIdx,std::size_t maturitySize,
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

		static auto const _maxVolatility_impl(std::size_t maturityIdx,std::size_t maturitySize,
			VolSurface const &volSurface, std::false_type) {
			return volSurface;
		}


	public:
		static auto const maxVolatility(std::size_t maturityIdx,std::size_t maturitySize, VolSurface const &volSurface) {
			return _maxVolatility_impl(maturityIdx, maturitySize,volSurface,std::is_compound<VolSurface>());
		}

	};




	// ==============================================================================
	// ====================== ImpliedProbability ====================================
	// ==============================================================================

	template<typename T>
	struct ImpliedProbability {
		static std::tuple<std::function<T(T,T,T)>,std::function<T(T,T,T,T,T)>,
			std::function<T(T,T,T,T,T,T,T)>,std::function<T(T,T,T,T,T,T)>> 
			const probabilityFunctions() {

			std::function<T(T, T, T)> upProb3 = [](T infl, T q1, T q2)->T {
				return (infl * q1 / q2);
			};

			std::function<T(T, T, T, T, T)> upProb5 = [](T infl, T q1,T pm, T q2,T q3)->T {
				return ((infl*q1 - pm * q2) / q3);
			};

			std::function<T(T, T, T, T, T, T, T)> upProb7 = [](T infl, T q1, T pd, T q2, T pm,T q3,T q4)->T {
				return ((infl*q1 - pd * q2 - pm * q3) / q4);
			};

			std::function<T(T, T, T, T, T, T)> midProb6 = [](T infl, T s1, T s2, T pu, T s3, T s4)->T {
				return ((infl*s1 - s2 - pu * (s3 - s2)) / (s4 - s2));
			};

			return std::make_tuple(upProb3, upProb5, upProb7, midProb6);
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
	private:
		template<typename LatticeObject>
		static CalibratorBinomialEquityResultsPtr<LatticeObject> const
			_implyTree_impl(LatticeObject &stockPriceLattice,
				DeltaTime const &deltaTime, 
				RiskFreeRate const &riskFreeRate,
				typename LatticeObject::Node_type const &apexPrice, 
				OptionData const &optionData);

	public:
		template<typename LatticeObject>
		static CalibratorBinomialEquityResultsPtr<LatticeObject> const
			implyTree(LatticeObject &stockPriceLattice, 
				DeltaTime const &deltaTime, 
				RiskFreeRate const &riskFreeRate,
				typename LatticeObject::Node_type const &apexPrice, 
				OptionData const &optionData) {
			return _implyTree_impl(stockPriceLattice, deltaTime, riskFreeRate, apexPrice, optionData);
		}

	};



	template<typename TimeAxis,
			typename DeltaTime,
			typename RiskFreeRate,
			typename OptionData>
	struct CalibratorEquity<LatticeType::Trinomial, TimeAxis, DeltaTime, RiskFreeRate, OptionData> {
	private:
		// ============= state price algos ===============
		// for liquid call option prices
		template<typename LatticeObject>
		static CalibratorTrinomialEquityResultsPtr<LatticeObject> const
			_statePriceLatticeCallKernel_impl(LatticeObject &stockPriceLattice, 
				DeltaTime const &deltaTime,
				typename LatticeObject::Node_type const &apexPrice,
				OptionData const &optionData);
		// for liquid put option prices
		template<typename LatticeObject>
		static CalibratorTrinomialEquityResultsPtr<LatticeObject> const
			_statePriceLatticePutKernel_impl(LatticeObject &stockPriceLattice,
				DeltaTime const &deltaTime,
				typename LatticeObject::Node_type const &apexPrice,
				OptionData const &optionData);

		// ============= implied probabilities algos ===============

		template<typename LatticeObject>
		static CalibratorTrinomialEquityResultsPtr<LatticeObject> const
			_impliedProbabilityCallKernel_impl(LatticeObject const &statePriceLattice, 
				LatticeObject const &stockPriceLattice,
				DeltaTime const &deltaTime, 
				RiskFreeRate const &riskFreeRate);

		template<typename LatticeObject>
		static CalibratorTrinomialEquityResultsPtr<LatticeObject> const
			_impliedProbabilityPutKernel_impl(LatticeObject const &statePriceLattice, 
				LatticeObject const &stockPriceLattice,
				DeltaTime const &deltaTime, 
				RiskFreeRate const &riskFreeRate);

	
	public:
		template<typename LatticeObject>
		static CalibratorTrinomialEquityResultsPtr<LatticeObject> const
			statePriceLattice(LatticeObject &stockPriceLattice,
				DeltaTime const &deltaTime,
				typename LatticeObject::Node_type const &apexPrice,
				OptionData const &optionData,
				bool areCallPricesLiquid = true) {
			if(areCallPricesLiquid)
				return _statePriceLatticeCallKernel_impl(stockPriceLattice, deltaTime, apexPrice, optionData);
			else
				return _statePriceLatticePutKernel_impl(stockPriceLattice, deltaTime, apexPrice, optionData);
		}

		template<typename LatticeObject>
		static CalibratorTrinomialEquityResultsPtr<LatticeObject> const
			impliedProbability(LatticeObject const &statePriceLattice,
				LatticeObject const &stockPriceLattice,
				DeltaTime const &deltaTime,
				RiskFreeRate const &riskFreeRate,
				bool areCallPricesLiquid = true) {
			LASSERT(statePriceLattice.timeDimension() == stockPriceLattice.timeDimension(),
				"statePriceLattice must be of the same time dimension as stockPriceLattice");
			if(areCallPricesLiquid)
				return _impliedProbabilityCallKernel_impl(statePriceLattice, stockPriceLattice, deltaTime, riskFreeRate);
			else
				return _impliedProbabilityPutKernel_impl(statePriceLattice, stockPriceLattice, deltaTime, riskFreeRate);
		}


	};
}





template<typename TimeAxis,
		typename DeltaTime,
		typename RiskFreeRate,
		typename OptionData>
template<typename LatticeObject>
lattice_calibrator_equity::CalibratorTrinomialEquityResultsPtr<LatticeObject> const
lattice_calibrator_equity::CalibratorEquity<lattice_types::LatticeType::Trinomial,TimeAxis,DeltaTime,RiskFreeRate, OptionData>::
_statePriceLatticeCallKernel_impl(LatticeObject &stockPriceLattice, DeltaTime const &deltaTime, typename LatticeObject::Node_type const &apexPrice,
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
	typedef OptionPriceLocator<decltype(strikes)> OPL;
	// typedef the vol surface holder:
	typedef VolSurfaceHolder<decltype(volSurface)> VS;
	// typedef deltaTime holder:
	typedef DeltaTimeHolder<DeltaTime> DT;
	// typedef the node type:
	typedef typename LatticeObject::Node_type Node;
	// typedef probability triplet:
	typedef std::tuple<Node, Node, Node> Triplet;

	// create statePriceLattice from empty stockPriceLattice:
	LatticeObject statePriceLattice(stockPriceLattice);
	std::size_t nodesSize{ 0 };
	Node dt{};
	Node vol{};
	Node fact{};
	Node lerpPrice{};
	Node priceUp{}, priceDown{};
	Node sum{};

	// First populate stockPriceLattice:
	stockPriceLattice(0, 0) = apexPrice;
	for (long long t = 1; t < treeSize; ++t) {
		dt = DT::deltaTime(t - 1, deltaTime);
		vol = VS::maxVolatility(t, maturitySize, volSurface);
		fact = vol * std::sqrt(3.0*dt);
		nodesSize = stockPriceLattice.nodesAtIdx(t).size();
		for (long long l = 0; l < nodesSize; ++l) {
			stockPriceLattice(t, l) = apexPrice * std::exp((l - t)*fact);
		}
	}

	// prepare pair prices of calls and puts:
	std::pair<decltype(strikes), decltype(strikes)> pairPrices;
	// populate statePriceLattice here:
	statePriceLattice(0, 0) = 1.0;
	// initialize LinearInterpolator here:
	LERP lerp;
	for (std::size_t t = 1; t < treeSize; ++t) {
		pairPrices = OPL::callPutPrices(t, maturitySize, optionSurface);

		// Set linear interpolation for put prices first
		lerp.setPoints(strikes, pairPrices.second);
		// Using put prices first as we traverse from down up
		// clean-up cached value
		sum = 0.0;
		for (std::size_t l = 0; l < t; ++l) {
			priceDown = stockPriceLattice(t, l);
			priceUp = stockPriceLattice(t, l + 1);
			lerpPrice = lerp.getValue(priceUp);
			statePriceLattice(t, l) = (lerpPrice  -  sum)/ (priceUp - priceDown);

			if (l != (t - 1)) {
				sum = 0.0;
				for (std::size_t k = 0; k < l + 1; ++k) {
					sum += (stockPriceLattice(t, l + 2) - stockPriceLattice(t, k)) * statePriceLattice(t, k);
				}
			}
		}
		// Set linear interpolation for call prices
		lerp.setPoints(strikes, pairPrices.first);
		// Using call prices as we traverse from up down to the center of tree
		nodesSize = statePriceLattice.nodesAtIdx(t).size();
		// clean-up cached value
		sum = 0.0;
		for (std::size_t l = nodesSize - 1; l >= t; --l) {
			priceDown = stockPriceLattice(t, l - 1);
			priceUp = stockPriceLattice(t, l);
			lerpPrice = lerp.getValue(priceDown);
			statePriceLattice(t, l) = (lerpPrice - sum) / (priceUp - priceDown);

			if (l != t) {
				sum = 0.0;
				for (std::size_t k = nodesSize - 1; k > l - 1; --k) {
					sum += (stockPriceLattice(t, k) - stockPriceLattice(t, l - 2)) * statePriceLattice(t, k);
				}
			}
		}
	}

	return lattice_calibrator_equity::CalibratorTrinomialEquityResultsPtr<LatticeObject>{new lattice_calibrator_equity::
		CalibratorTrinomialEquityResultsT<LatticeObject>(statePriceLattice)};
}


template<typename TimeAxis,
	typename DeltaTime,
	typename RiskFreeRate,
	typename OptionData>
	template<typename LatticeObject>
lattice_calibrator_equity::CalibratorTrinomialEquityResultsPtr<LatticeObject> const
lattice_calibrator_equity::CalibratorEquity<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime, RiskFreeRate, OptionData>::
_statePriceLatticePutKernel_impl(LatticeObject &stockPriceLattice, DeltaTime const &deltaTime, typename LatticeObject::Node_type const &apexPrice,
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
	typedef OptionPriceLocator<decltype(strikes)> OPL;
	// typedef the vol surface holder:
	typedef VolSurfaceHolder<decltype(volSurface)> VS;
	// typedef deltaTime holder:
	typedef DeltaTimeHolder<DeltaTime> DT;
	// typedef the node type:
	typedef typename LatticeObject::Node_type Node;
	// typedef probability triplet:
	typedef std::tuple<Node,Node,Node> Triplet;

	// create statePriceLattice from empty stockPriceLattice:
	LatticeObject statePriceLattice(stockPriceLattice);
	std::size_t nodesSize{ 0 };
	Node dt{};
	Node vol{};
	Node fact{};
	Node lerpPrice{};
	Node priceUp{}, priceDown{};
	Node sum{};

	// First populate stockPriceLattice:
	stockPriceLattice(0, 0) = apexPrice;
	for (long long t = 1; t < treeSize; ++t) {
		dt = DT::deltaTime(t - 1, deltaTime);
		vol = VS::maxVolatility(t, maturitySize, volSurface);
		fact = vol * std::sqrt(3.0*dt);
		nodesSize = stockPriceLattice.nodesAtIdx(t).size();
		for (long long l = 0; l < nodesSize; ++l) {
			stockPriceLattice(t, l) = apexPrice * std::exp((l - t)*fact);
		}
	}

	// prepare pair prices of calls and puts:
	std::pair<decltype(strikes), decltype(strikes)> pairPrices;
	// populate statePriceLattice here:
	statePriceLattice(0, 0) = 1.0;
	// initialize LinearInterpolator here:
	LERP lerp;
	for (std::size_t t = 1; t < treeSize; ++t) {
		pairPrices = OPL::callPutPrices(t, maturitySize, optionSurface);

		// Set linear interpolation for put prices first
		lerp.setPoints(strikes, pairPrices.second);
		// Using put prices first as we traverse from down up
		// clean-up cached value
		sum = 0.0;
		for (std::size_t l = 0; l <= t; ++l) {
			priceDown = stockPriceLattice(t, l);
			priceUp = stockPriceLattice(t, l + 1);
			lerpPrice = lerp.getValue(priceUp);
			statePriceLattice(t, l) = (lerpPrice - sum) / (priceUp - priceDown);

			if (l != t) {
				sum = 0.0;
				for (std::size_t k = 0; k < l + 1; ++k) {
					sum += (stockPriceLattice(t, l + 2) - stockPriceLattice(t, k)) * statePriceLattice(t, k);
				}
			}
		}
		// Set linear interpolation for call prices
		lerp.setPoints(strikes, pairPrices.first);
		// Using call prices as we traverse from up down to the center of tree
		nodesSize = statePriceLattice.nodesAtIdx(t).size();
		// clean-up cached value
		sum = 0.0;
		for (std::size_t l = nodesSize - 1; l > t; --l) {
			priceDown = stockPriceLattice(t, l - 1);
			priceUp = stockPriceLattice(t, l);
			lerpPrice = lerp.getValue(priceDown);
			statePriceLattice(t, l) = (lerpPrice - sum) / (priceUp - priceDown);

			if (l != (t + 1)) {
				sum = 0.0;
				for (std::size_t k = nodesSize - 1; k > l - 1; --k) {
					sum += (stockPriceLattice(t, k) - stockPriceLattice(t, l - 2)) * statePriceLattice(t, k);
				}
			}
		}
	}

	return lattice_calibrator_equity::CalibratorTrinomialEquityResultsPtr<LatticeObject>{new lattice_calibrator_equity::
		CalibratorTrinomialEquityResultsT<LatticeObject>(statePriceLattice)};
}




template<typename TimeAxis,
	typename DeltaTime,
	typename RiskFreeRate,
	typename OptionData>
template<typename LatticeObject>
lattice_calibrator_equity::CalibratorTrinomialEquityResultsPtr<LatticeObject> const
lattice_calibrator_equity::CalibratorEquity<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime, RiskFreeRate, OptionData>::
_impliedProbabilityCallKernel_impl(LatticeObject const &statePriceLattice, LatticeObject const &stockPriceLattice,
	DeltaTime const &deltaTime, RiskFreeRate const &riskFreeRate) {

	// typedef deltaTime holder:
	typedef DeltaTimeHolder<DeltaTime> DT;
	// typedef the node type:
	typedef typename LatticeObject::Node_type Node;
	// typedef riskFreeRate holder:
	typedef RiskFreeRateHolder<RiskFreeRate> RFR;
	// typedef probability triplet
	typedef std::tuple<Node, Node, Node> Triplet;
	
	// unpack the implied probability functions:
	auto probTpl = ImpliedProbability<Node>::probabilityFunctions();
	// unpack them into up_3,up_5,up_7 and mid_6:
	auto up_3 = std::get<0>(probTpl);
	auto up_5 = std::get<1>(probTpl);
	auto up_7 = std::get<2>(probTpl);
	auto mid_6 = std::get<3>(probTpl);

	// get timeDimension of the lattice:
	std::size_t const treeSize = statePriceLattice.timeDimension();
	std::size_t nodesSize{ 0 };
	// prepare container for probabilities:
	std::vector<std::vector<Triplet>> impliedProbs(treeSize - 1);
	// Prepare temporary values:
	Node up{}, mid{}, down{}, p_m{}, p_u{}, p_d{};
	Node dt{}, rfr{}, infl{};


	// going through time
	for (long long t = 0; t < treeSize - 1; ++t) {
		dt = DT::deltaTime(t, deltaTime);
		rfr = RFR::rate(t, riskFreeRate);
		infl = std::exp(rfr*dt);
		nodesSize = statePriceLattice.nodesAtIdx(t).size();
		std::vector<Triplet> probs(nodesSize);
		// first going from bottom of the branch to the node before the center of tree:
		for (long long l = 0; l < t; ++l) {
			if (l == 0) {
				down = probFloorCapper(up_3(infl, statePriceLattice(t + 1, l), statePriceLattice(t, l)));
			}
			else if (l == 1) {
				p_m = std::get<1>(probs.at(l - 1));
				down = probFloorCapper(up_5(infl, statePriceLattice(t + 1, l), p_m,
					statePriceLattice(t, l-1), statePriceLattice(t, l)));
			}
			else {
				p_u = std::get<2>(probs.at(l - 2));
				p_m = std::get<1>(probs.at(l - 1));
				down = probFloorCapper(up_7(infl, statePriceLattice(t + 1, l), p_u, statePriceLattice(t, l - 2),
					p_m, statePriceLattice(t, l - 1), statePriceLattice(t, l)));
			}
			mid = probFloorCapper(mid_6(infl, stockPriceLattice(t, l), stockPriceLattice(t + 1, l + 2), down, stockPriceLattice(t + 1, l),
				stockPriceLattice(t + 1, l + 1)));
			up = 1.0 -  mid -  down;
			probs[l] = std::make_tuple(down, mid, up);
		}

		// first going from top of the branch to the center of tree:
		for (long long l = nodesSize - 1; l >= t; l--) {
			if (l == (nodesSize - 1)) {
				up = probFloorCapper(up_3(infl, statePriceLattice(t + 1, l + 2), statePriceLattice(t, l)));
			}
			else if (l == (nodesSize - 2)) {
				p_m = std::get<1>(probs.at(l + 1));
				up = probFloorCapper(up_5(infl, statePriceLattice(t + 1, l + 2), p_m,
					statePriceLattice(t, l + 1), statePriceLattice(t, l)));
			}
			else {
				p_d = std::get<0>(probs.at(l + 2));
				p_m = std::get<1>(probs.at(l + 1));
				up = probFloorCapper(up_7(infl, statePriceLattice(t + 1, l+2), p_d, statePriceLattice(t, l + 2),
					p_m, statePriceLattice(t, l + 1), statePriceLattice(t, l)));
			}
			mid = probFloorCapper(mid_6(infl, stockPriceLattice(t, l), stockPriceLattice(t + 1, l), up, stockPriceLattice(t + 1, l+2),
				stockPriceLattice(t + 1, l + 1)));
			down = 1.0 - mid - up;
			probs[l] = std::make_tuple(down, mid, up);

		}
		impliedProbs[t] = std::move(probs);
	}


	return lattice_calibrator_equity::CalibratorTrinomialEquityResultsPtr<LatticeObject>{new lattice_calibrator_equity::
		CalibratorTrinomialEquityResultsT<LatticeObject>(statePriceLattice,impliedProbs)};

}


template<typename TimeAxis,
	typename DeltaTime,
	typename RiskFreeRate,
	typename OptionData>
	template<typename LatticeObject>
lattice_calibrator_equity::CalibratorTrinomialEquityResultsPtr<LatticeObject> const
lattice_calibrator_equity::CalibratorEquity<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime, RiskFreeRate, OptionData>::
_impliedProbabilityPutKernel_impl(LatticeObject const &statePriceLattice, LatticeObject const &stockPriceLattice,
	DeltaTime const &deltaTime, RiskFreeRate const &riskFreeRate) {

	// typedef deltaTime holder:
	typedef DeltaTimeHolder<DeltaTime> DT;
	// typedef the node type:
	typedef typename LatticeObject::Node_type Node;
	// typedef riskFreeRate holder:
	typedef RiskFreeRateHolder<RiskFreeRate> RFR;
	// typedef probability triplet
	typedef std::tuple<Node, Node, Node> Triplet;

	// unpack the implied probability functions:
	auto probTpl = ImpliedProbability<Node>::probabilityFunctions();
	// unpack them into up_3,up_5,up_7 and mid_6:
	auto up_3 = std::get<0>(probTpl);
	auto up_5 = std::get<1>(probTpl);
	auto up_7 = std::get<2>(probTpl);
	auto mid_6 = std::get<3>(probTpl);

	// get timeDimension of the lattice:
	std::size_t const treeSize = statePriceLattice.timeDimension();
	std::size_t nodesSize{ 0 };
	// prepare container for probabilities:
	std::vector<std::vector<Triplet>> impliedProbs(treeSize - 1);
	// Prepare temporary values:
	Node up{}, mid{}, down{}, p_m{}, p_u{}, p_d{};
	Node dt{}, rfr{}, infl{};


	// going through time
	for (long long t = 0; t < treeSize - 1; ++t) {
		dt = DT::deltaTime(t, deltaTime);
		rfr = RFR::rate(t, riskFreeRate);
		infl = std::exp(rfr*dt);
		nodesSize = statePriceLattice.nodesAtIdx(t).size();
		std::vector<Triplet> probs(nodesSize);
		// first going from bottom of the branch to the node before the center of tree:
		for (long long l = 0; l <= t; ++l) {
			if (l == 0) {
				down = probFloorCapper(up_3(infl, statePriceLattice(t + 1, l), statePriceLattice(t, l)));
			}
			else if (l == 1) {
				p_m = std::get<1>(probs.at(l - 1));
				down = probFloorCapper(up_5(infl, statePriceLattice(t + 1, l), p_m,
					statePriceLattice(t, l - 1), statePriceLattice(t, l)));
			}
			else {
				p_u = std::get<2>(probs.at(l - 2));
				p_m = std::get<1>(probs.at(l - 1));
				down = probFloorCapper(up_7(infl, statePriceLattice(t + 1, l), p_u, statePriceLattice(t, l - 2),
					p_m, statePriceLattice(t, l - 1), statePriceLattice(t, l)));
			}
			mid = probFloorCapper(mid_6(infl, stockPriceLattice(t, l), stockPriceLattice(t + 1, l + 2), down, stockPriceLattice(t + 1, l),
				stockPriceLattice(t + 1, l + 1)));
			up = 1.0 - mid - down;
			probs[l] = std::make_tuple(down, mid, up);
		}

		// first going from top of the branch to the center of tree:
		for (long long l = nodesSize - 1; l > t; l--) {
			if (l == (nodesSize - 1)) {
				up = probFloorCapper(up_3(infl, statePriceLattice(t + 1, l + 2), statePriceLattice(t, l)));
			}
			else if (l == (nodesSize - 2)) {
				p_m = std::get<1>(probs.at(l + 1));
				up = probFloorCapper(up_5(infl, statePriceLattice(t + 1, l + 2), p_m,
					statePriceLattice(t, l + 1), statePriceLattice(t, l)));
			}
			else {
				p_d = std::get<0>(probs.at(l + 2));
				p_m = std::get<1>(probs.at(l + 1));
				up = probFloorCapper(up_7(infl, statePriceLattice(t + 1, l + 2), p_d, statePriceLattice(t, l + 2),
					p_m, statePriceLattice(t, l + 1), statePriceLattice(t, l)));
			}
			mid = probFloorCapper(mid_6(infl, stockPriceLattice(t, l), stockPriceLattice(t + 1, l), up, stockPriceLattice(t + 1, l + 2),
				stockPriceLattice(t + 1, l + 1)));
			down = 1.0 - mid - up;
			probs[l] = std::make_tuple(down, mid, up);

		}
		impliedProbs[t] = std::move(probs);
	}


	return lattice_calibrator_equity::CalibratorTrinomialEquityResultsPtr<LatticeObject>{new lattice_calibrator_equity::
		CalibratorTrinomialEquityResultsT<LatticeObject>(statePriceLattice, impliedProbs)};

}


template<typename TimeAxis,
	typename DeltaTime,
	typename RiskFreeRate,
	typename OptionData>
template<typename LatticeObject>
static lattice_calibrator_equity::CalibratorBinomialEquityResultsPtr<LatticeObject> const
lattice_calibrator_equity::CalibratorEquity<lattice_types::LatticeType::Binomial,TimeAxis,DeltaTime, RiskFreeRate,OptionData>::
_implyTree_impl(LatticeObject &stockPriceLattice, DeltaTime const &deltaTime, RiskFreeRate const &riskFreeRate,
	typename LatticeObject::Node_type const &apexPrice, OptionData const &optionData) {

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
	LASSERT((strikeSize * maturitySize) == optionSurfaceSize,
		"option surface must have size strikes x maturities");
	// Get the trees time dimension:
	std::size_t const treeSize = stockPriceLattice.timeDimension();
	// Tree dimension must match the size of maturities:
	LASSERT((treeSize == maturitySize),
		"Time dimension of the stockPriceLattice must match the size of maturities");

	// typedef the LinearInterpolator:
	typedef LinearInterpolator<decltype(strikes)> LERP;
	// typedef the OptionPriceLocator:
	typedef OptionPriceLocator<decltype(strikes)> OPL;
	// typedef the vol surface holder:
	typedef VolSurfaceHolder<decltype(volSurface)> VS;
	// typedef deltaTime holder:
	typedef DeltaTimeHolder<DeltaTime> DT;
	// typedef RiskFreeRate holder:
	typedef RiskFreeRateHolder<RiskFreeRate> RFR;
	// typedef the node type:
	typedef typename LatticeObject::Node_type Node;
	// typedef probability holder:
	typedef std::pair<Node, Node> Pair;

	// prepare container for probabilities:
	std::vector<std::vector<Pair>> impliedProbs(treeSize - 1);
	// create statePriceLattice from empty stockPriceLattice:
	LatticeObject statePriceLattice(stockPriceLattice);
	std::size_t nodesSize{ 0 };
	long startIdx{ 0 };
	Node dt{};
	Node rfr{};
	Node infl{};
	Node stock{};
	Node call{}, put{};
	Node sum{};
	Node p{}, q{};

	// First populate stockPriceLattice apex:
	stockPriceLattice(0, 0) = apexPrice;
	// Second populate statePriceLattice apex:
	statePriceLattice(0, 0) = 1.0;
	// prepare pair prices of calls and puts:
	std::pair<decltype(strikes), decltype(strikes)> pairPrices;
	// initialize LinearInterpolator here:
	LERP lerp;
	// first populate stockPriceLattice:
	for (std::size_t t = 1; t < treeSize; ++t) {
		dt = DT::deltaTime(t, deltaTime);
		rfr = RFR::rate(t, riskFreeRate);
		infl = std::exp(rfr*dt);
		pairPrices = OPL::callPutPrices(t, maturitySize, optionSurface);
		// Set linear interpolation for call prices first
		lerp.setPoints(strikes, pairPrices.first,true);
		// Using call prices first as we traverse from top to the center of tree:
		nodesSize = stockPriceLattice.nodesAtIdx(t).size();
		// first with even number of nodes:
		if ((nodesSize % 2) == 0) {
			for (std::size_t j = (nodesSize / 2); j < nodesSize; ++j) {
				// clear the cached value:
				sum = 0.0;
				for (std::size_t k = j; k < t; ++k) {
					sum += (statePriceLattice(t - 1, k)*
						(infl*stockPriceLattice(t - 1, k) - stockPriceLattice(t - 1, j - 1)));
				}
				stock = stockPriceLattice(t - 1, j - 1);
				call = lerp.getValue(stock);
				stockPriceLattice(t, j) = ((stock*(infl*call + statePriceLattice(t - 1, j - 1)*stock - sum)) /
					(statePriceLattice(t - 1, j - 1)*infl*stock - infl * call + sum));
				if (j == (nodesSize / 2)) {
					stockPriceLattice(t, j - 1) = (stock*stock / stockPriceLattice(t, j));
				}
			}
		}
		// then with odd number of nodes:
		else {
			startIdx = (nodesSize - 1) / 2;
			stockPriceLattice(t, startIdx) = apexPrice;
			for (std::size_t j = startIdx + 1; j < nodesSize; ++j) {
				// clear the cached value:
				sum = 0.0;
				for (std::size_t k = j; k < t; ++k) {
					sum += (statePriceLattice(t - 1, k)*
						(infl*stockPriceLattice(t - 1, k) - stockPriceLattice(t - 1, j - 1)));
				}
				stock = stockPriceLattice(t - 1, j - 1);
				call = lerp.getValue(stock);
				stockPriceLattice(t, j) = ((stockPriceLattice(t, j - 1)*(infl*call - sum) -
					statePriceLattice(t - 1, j - 1)*stock*(infl*stock - stockPriceLattice(t, j - 1))) /
					((infl*call - sum) - statePriceLattice(t - 1, j - 1)*(infl*stock - stockPriceLattice(t, j - 1))));
			}
		}
		// Using put prices next as we traverse from center of tree down:
		if (t > 1) {
			// Set linear interpolation for put prices first
			lerp.setPoints(strikes, pairPrices.second,true);

			if ((nodesSize % 2) != 0) {
				startIdx = ((nodesSize - 1) / 2) - 1;
			}
			else {
				startIdx = (nodesSize / 2) - 2;
			}
			for (long j = startIdx; j >= 0; --j) {
				// clear the cached value:
				sum = 0.0;
				for (std::size_t k = 0; k < j; ++k) {
					sum += (statePriceLattice(t - 1, k)*
						(stockPriceLattice(t - 1, j) - infl * stockPriceLattice(t - 1, k)));
				}

				stock = stockPriceLattice(t - 1, j);
				put = lerp.getValue(stock);

				stockPriceLattice(t, j) = ((stockPriceLattice(t, j + 1)*(infl*put - sum) +
					statePriceLattice(t - 1, j)*stock*(infl*stock - stockPriceLattice(t, j + 1))) /
					((infl*put - sum) + statePriceLattice(t - 1, j)*(infl*stock - stockPriceLattice(t, j + 1))));
			}
		}
		// next, compute probabilities:
		nodesSize = stockPriceLattice.nodesAtIdx(t - 1).size();
		std::vector<Pair> probs(nodesSize);
		for (std::size_t j = 0; j < nodesSize; ++j) {
			p = ((infl*stockPriceLattice(t - 1, j) - stockPriceLattice(t, j)) /
				(stockPriceLattice(t, j + 1) - stockPriceLattice(t, j)));
			q = 1.0 - p;
			probs[j] = std::make_pair(q, p);
		}
		impliedProbs[t - 1] = std::move(probs);

		// finally compute state price lattice:
		statePriceLattice(t, 0) = (statePriceLattice(t - 1, 0)*impliedProbs[t - 1][0].first) / infl;
		statePriceLattice(t, t) = (statePriceLattice(t - 1, t - 1)*impliedProbs[t - 1][t - 1].second) / infl;
		nodesSize = stockPriceLattice.nodesAtIdx(t).size();
		for (std::size_t j = 1; j < nodesSize - 1; ++j) {
			statePriceLattice(t, j) = ((statePriceLattice(t - 1, j) * impliedProbs[t - 1][j].first +
				statePriceLattice(t - 1, j - 1) * impliedProbs[t - 1][j - 1].second) / infl);
		}
	}

	return lattice_calibrator_equity::CalibratorBinomialEquityResultsPtr<LatticeObject>{new lattice_calibrator_equity::
		CalibratorBinomialEquityResultsT<LatticeObject>(statePriceLattice, impliedProbs)};

}


#endif ///_LATTICE_CALIBRATOR_EQUITY