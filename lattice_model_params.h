#pragma once
#if !defined(_LATTICE_MODEL_PARAMS)
#define _LATTICE_MODEL_PARAMS

#include"lattice_types.h"

namespace lattice_model_params {

	using lattice_types::AssetClass;

	template<std::size_t FactorCount,AssetClass Asset,typename T>
	struct ModelParams {
	};


	template<typename T>
	struct ModelParams<1,AssetClass::InterestRate,T> {
		T ReversionSpeed;
		T Volatility;
	};


	template<typename T>
	struct ModelParams<1,AssetClass::Equity, T> {
		T RiskFreeRate;
		T Volatility;
		T DividendRate;
		T Spot;
		T Strike;
	};

	template<typename T>
	struct ModelParams<2, AssetClass::Equity, T> {
		T RiskFreeRate;
		T DividendRate1;
		T DividendRate2;
		T Volatility1;
		T Volatility2;
		T Spot1;
		T Spot2;
		T Strike;
		T Correlation;
	};

}

#endif ///_LATTICE_MODEL_PARAMS