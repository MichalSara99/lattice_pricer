#pragma once
#if !defined(_LATTICE_MODEL)
#define _LATTICE_MODEL

#include"lattice_miscellaneous.h"

namespace lattice_model {

	using lattice_miscellaneous::OptionData;


	// ===================== Cox-Rubinstein-Ross model ===========================

	template<typename T=double>
	class CoxRubinsteinRossModel {
	private:
		T upFactor_;
		T downFactor_;
		T prob_;
		T discount_;
		OptionData<T> option_;

	public:
		CoxRubinsteinRossModel(OptionData<T>const &option)
			:option_{ option }, prob_{0.5} {
		}

		// Forward generator
		std::tuple<T, T> operator()(T value, T dt) {
			T s = option_.Volatility;
			T r = option_.RiskFreeRate;
			T r1 = (r - 0.5*s*s)*dt;
			T r2 = s * std::sqrt(dt);
			upFactor_ = std::exp(r1 + r2);
			downFactor_ = std::exp(r1 - r2);
			return std::make_tuple(upFactor_*value, downFactor_*value);
		}

		// Backward generator
		T operator()(T upValue, T downValue, T dt) {
			T r = option_.RiskFreeRate;
			discount_ = std::exp(-1.0*r*dt);
			return (discount_ * (prob_*upValue + (1.0 - prob_)*downValue));
		}
	};







}




#endif //_LATTICE_MODEL