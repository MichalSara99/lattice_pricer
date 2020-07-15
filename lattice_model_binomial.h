#pragma once
#if !defined(_LATTICE_MODEL_BINOMIAL)
#define  _LATTICE_MODEL_BINOMIAL

#include<tuple>
#include"lattice_model_interface.h"
#include"lattice_types.h"
#include"lattice_utility.h"
#include"lattice_model_params.h"

namespace lattice_model {

	using lattice_model_params::ModelParams;
	using lattice_types::DiscountingStyle;
	using lattice_utility::DiscountingFactor;

	// =============================================================================================
	// ===================== Cox-Rubinstein-Ross model (binomial lattice) ==========================
	// =============================================================================================

	template<typename T = double>
	class CoxRubinsteinRossModel :public BinomialModel<1, T> {
	private:
		T prob_;
		ModelParams<1,AssetClass::Equity,T> params_;

	public:
		CoxRubinsteinRossModel(ModelParams<1,AssetClass::Equity,T>const &params)
			:params_{ params }, prob_{ 0.5 } {
		}

		// Forward generator
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting) const override {
			T const q = params_.DividendRate;
			T const sig = params_.Volatility;
			T const r = params_.RiskFreeRate;
			T const r1 = (r - q - 0.5*sig*sig)*dt;
			T const r2 = sig * std::sqrt(dt);
			T const up = std::exp(r1 + r2);
			T const down = std::exp(r1 - r2);
			return std::make_tuple(down*value,up*value);
		}

		// Backward generator
		T operator()(T currValue, T downValue, T upValue, T dt) override {
			T const r = params_.RiskFreeRate;
			T const disc = std::exp(-1.0*r*dt);
			return (disc * (prob_*upValue + (1.0 - prob_)*downValue));
		}

		static std::string const name() {
			return std::string{ "Cox-Rubinstein-Ross model" };
		}

		static constexpr AssetClass assetClass() { return AssetClass::Equity; }
	};


	// =============================================================================================
	// ================ Two-factor Cox-Rubinstein-Ross model (binomial lattice) ====================
	// =============================================================================================

	template<typename T=double>
	class CoxRubinsteinRossModel2Factor :public BinomialModel<2, T> {
	private:
		T rho_;
		ModelParams<2, AssetClass::Equity, T> params_;

	public:
		explicit CoxRubinsteinRossModel2Factor(ModelParams<2, AssetClass::Equity, T> const& params)
			:params_(params) {}

		// Returns risk-neutral probability:
		std::tuple<T,T,T,T> nodeRiskNeutralProb(T dt)const {
			T const r = params_.RiskFreeRate;
			T const q1 = params_.DividendRate1;
			T const q2 = params_.DividendRate2;
			T const sig1 = params_.Volatility1;
			T const sig2 = params_.Volatility2;
			T const nu1 = r - q1 - 0.5*sig1*sig1;
			T const nu2 = r - q2 - 0.5*sig2*sig2;
			T const dx1 = sig1 * std::sqrt(dt);
			T const dx2 = sig2 * std::sqrt(dt);
			T const mix = dx1*dx2;
			T const corr = params_.Correlation * sig1*sig2;

			T const pdd = (mix + (-1.0*dx2*nu1 - dx1 * nu2 + corr)*dt) / (4.0*mix);
			T const pdu = (mix + (-1.0*dx2*nu1 + dx1 * nu2 - corr)*dt) / (4.0*mix);
			T const pud = (mix + (dx2*nu1 - dx1 * nu2 - corr)*dt) / (4.0*mix);
			T const puu = (mix + (dx2*nu1 + dx1 * nu2 + corr)*dt) / (4.0*mix);
			return std::make_tuple(pdd, pdu, pud, puu);
		}
		
		// Forward generators:
		std::pair<LeafForwardGenerator<T, T, T>,
			LeafForwardGenerator<T, T, T>> forwardGenerator() const override {
			ModelParams<1, AssetClass::Equity, T> params1;
			params1.RiskFreeRate = params_.RiskFreeRate;
			params1.DividendRate = params_.DividendRate1;
			params1.Volatility = params_.Volatility1;
			params1.Spot = params_.Spot1;
			params1.Strike = params_.Strike;
			ModelParams<1, AssetClass::Equity, T> params2;
			params2.RiskFreeRate = params_.RiskFreeRate;
			params2.DividendRate = params_.DividendRate2;
			params2.Volatility = params_.Volatility2;
			params2.Spot = params_.Spot2;
			params2.Strike = params_.Strike;
			CoxRubinsteinRossModel<T> factor1{ params1 };
			CoxRubinsteinRossModel<T> factor2{ params2 };
			LeafForwardGenerator<T, T, T> first = 
				[=](T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting)->std::tuple<T,T> {
				return factor1(value, dt, leafIdx, timeIdx, isMeanReverting);
			};
			LeafForwardGenerator<T, T, T> second = 
				[=](T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting)->std::tuple<T,T> {
				return factor2(value, dt, leafIdx, timeIdx, isMeanReverting);
			};
			return std::make_pair(first, second);
		}

		// Backward generator:
		T operator()(T currValue, T downDownValue, T downUpValue, T upDownValue, T upUpValue, T dt) {
			// taking risk-free rate for discounting from first factor data:
			T const r = params_.RiskFreeRate;
			T const disc = std::exp(-1.0*r*dt);
			std::tuple<T, T, T, T> const prob = nodeRiskNeutralProb(dt);
			T const value = (std::get<0>(prob) *downDownValue +
				std::get<1>(prob) *downUpValue +
				std::get<2>(prob) * upDownValue +
				std::get<3>(prob) * upUpValue);
			return (disc * value);
		}

		static std::string const name() {
			return std::string{ "Cox-Rubinstein-Ross 2-factor model" };
		}

		static constexpr AssetClass assetClass() { return AssetClass::Equity; }

	};




	// =======================================================================================================
	// ===================== Modified Cox-Rubinstein-Ross model (binomial lattice) ===========================
	// =======================================================================================================

	template<typename T = double>
	class ModifiedCoxRubinsteinRossModel :public BinomialModel<1, T> {
	private:
		ModelParams<1, AssetClass::Equity, T> params_;
		std::size_t n_;

	public:
		ModifiedCoxRubinsteinRossModel(ModelParams<1, AssetClass::Equity, T> const &params, std::size_t numberPeriods)
			:params_{ params }, n_{ numberPeriods } {
		}

		// Forward generator:
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting)const override {
			T const s = params_.Spot;
			T const sig = params_.Volatility;
			T const k = params_.Strike;
			T const K_n = (std::log(k / s) / static_cast<T>(n_));
			T const V_n = sig * std::sqrt(dt);
			T const up = std::exp(K_n + V_n);
			T const down = std::exp(K_n - V_n);
			return std::make_tuple(down*value, up*value);
		}

		// Backward generator:
		T operator()(T currValue, T downValue, T upValue, T dt) override {
			T const q = params_.DividendRate;
			T const s = params_.Spot;
			T const sig = params_.Volatility;
			T const r = params_.RiskFreeRate;
			T const k = params_.Strike;
			T const K_n = (std::log(k / s) / static_cast<T>(n_));
			T const V_n = sig * std::sqrt(dt);
			T const up = std::exp(K_n + V_n);
			T const down = std::exp(K_n - V_n);
			T const disc = std::exp(-1.0*r *dt);
			T const prob = (std::exp((r - q)*dt) - down) / (up - down);
			return (disc * (prob*upValue + (1.0 - prob)*downValue));
		}

		static std::string const name() {
			return std::string{ "Modified Cox-Rubinstein-Ross model" };
		}
		static constexpr AssetClass assetClass() { return AssetClass::Equity; }
	};

	// =====================================================================================
	// ===================== Jarrow-Rudd model (binomial lattice) ==========================
	// =====================================================================================

	template<typename T = double>
	class JarrowRuddModel :public BinomialModel<1, T> {
	private:
		T prob_;
		ModelParams<1, AssetClass::Equity, T> params_;

	public:
		JarrowRuddModel(ModelParams<1, AssetClass::Equity, T> const &params)
			:params_{ params }, prob_{ 0.5 } {
		}

		// Forward generator
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting) const override {
			T const q = params_.DividendRate;
			T const sig = params_.Volatility;
			T const r = params_.RiskFreeRate;
			T const d = (r - q - 0.5*sig*sig)*dt;
			T const x1 = d + sig * std::sqrt(dt);
			T const x2 = d - sig * std::sqrt(dt);
			T const up = std::exp(x1);
			T const down = std::exp(x2);
			return std::make_tuple(down*value,up*value);
		}

		// Backward generator
		T operator()(T currValue, T downValue, T upValue, T dt) override {
			T const r = params_.RiskFreeRate;
			T const disc = std::exp(-1.0*r *dt);
			return (disc * (prob_*upValue + (1.0 - prob_)*downValue));
		}

		static std::string const name() {
			return std::string{ "Jarrow-Rudd model" };
		}
		static constexpr AssetClass assetClass() { return AssetClass::Equity; }
	};

	// =====================================================================================
	// ===================== Trigeorgis model (binomial lattice) ===========================
	// =====================================================================================

	template<typename T = double>
	class TrigeorgisModel :public BinomialModel<1, T> {
	private:
		T gamma_;
		ModelParams<1, AssetClass::Equity, T> params_;

	public:
		TrigeorgisModel(ModelParams<1, AssetClass::Equity, T> const &params)
			:params_{ params },
			gamma_{ params.RiskFreeRate - params.DividendRate - 0.5*params.Volatility*params.Volatility } {
		}

		// Forward generator:
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting) const override {
			T const sig = params_.Volatility;
			T const x = std::sqrt(sig*sig*dt + gamma_ * gamma_*dt*dt);
			T const up = std::exp(x);
			T const down = (1.0 / up);
			return std::make_tuple(down*value,up*value);
		}

		// Backward generator:
		T operator()(T currValue, T downValue, T upValue, T dt)override {
			T const sig = params_.Volatility;
			T const r = params_.RiskFreeRate;
			T const x = std::sqrt(sig*sig*dt + gamma_ * gamma_*dt*dt);
			T const disc = std::exp(-1.0*r *dt);
			T const prob = 0.5*(1.0 + (gamma_*(dt / x)));
			return (disc * (prob*upValue + (1.0 - prob)*downValue));
		}

		static std::string const name() {
			return std::string{ "Trigeorgis model" };
		}

		static constexpr AssetClass assetClass() { return AssetClass::Equity; }
	};

	// ===============================================================================
	// ===================== Tian model (binomial lattice) ===========================
	// ===============================================================================

	template<typename T = double>
	class TianModel :public BinomialModel<1, T> {
	private:
		ModelParams<1, AssetClass::Equity, T> params_;

	public:
		TianModel(ModelParams<1, AssetClass::Equity, T> const &params)
			:params_{ params } {
		}

		// Forward generator:
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting)const  override {
			T const q = params_.DividendRate;
			T const sig = params_.Volatility;
			T const r = params_.RiskFreeRate;
			T const v = std::exp(sig*sig*dt);
			T const x = std::sqrt(v *v + 2.0*v - 3.0);
			T const up = (0.5*std::exp((r - q)*dt)*v*(v + 1.0 + x));
			T const down = (0.5*std::exp((r - q)*dt)*v*(v + 1.0 - x));
			return std::make_tuple(down*value, up*value);
		}

		// Backward generator:
		T operator()(T currValue, T downValue, T upValue, T dt)override {
			T const q = params_.DividendRate;
			T const sig = params_.Volatility;
			T const r = params_.RiskFreeRate;
			T const v = std::exp(sig*sig*dt);
			T const x = std::sqrt(v *v + 2.0*v - 3.0);
			T const up = (0.5*std::exp((r - q)*dt)*v*(v + 1.0 + x));
			T const down = (0.5*std::exp((r - q)*dt)*v*(v + 1.0 - x));
			T const disc = std::exp(-1.0*r*dt);
			T const prob = ((std::exp((r - q)*dt) - down) / (up - down));
			return (disc * (prob*upValue + (1.0 - prob)*downValue));
		}

		static std::string const name() {
			return std::string{ "Tian model" };
		}

		static constexpr AssetClass assetClass() { return AssetClass::Equity; }
	};

	// ========================================================================================
	// ===================== Leisen-Reimer model (binomial lattice) ===========================
	// ========================================================================================

	template<typename T = double>
	class LeisenReimerModel :public BinomialModel<1, T> {
	private:
		ModelParams<1, AssetClass::Equity, T> params_;
		std::size_t numberPeriods_;
		std::function<T(T)> inversion_;

	public:
		LeisenReimerModel(ModelParams<1, AssetClass::Equity, T> const& params,
			std::size_t numberPeriods,
			std::function<T(T)> const &inversionFormula)
			:params_{ params },
			numberPeriods_{ numberPeriods },
			inversion_{ inversionFormula } {
		}

		// Forward generator:
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting) const override {
			T const sig = params_.Volatility;
			T const r = params_.RiskFreeRate;
			T const q = params_.DividendRate;
			T const s = params_.Spot;
			T const k = params_.Strike;
			T const d1 = ((std::log(s / k) + (r - q + 0.5 * sig * sig)) * (numberPeriods_ * dt)) /
				(sig * std::sqrt(numberPeriods_ * dt));
			T const d2 = d1 - sig * std::sqrt(numberPeriods_ * dt);
			T const e = std::exp((r - q) * dt);
			T const h1 = inversion_(d1);
			T const h2 = inversion_(d2);
			T const p = h2;
			T const up = e * (h1 / h2);
			T const down = (e - p * up) / (1.0 - p);
			return std::make_tuple(down * value, up * value);
		}

		// Backward generator:
		T operator()(T currValue, T downValue, T upValue, T dt)override {
			T const sig = params_.Volatility;
			T const r = params_.RiskFreeRate;
			T const q = params_.DividendRate;
			T const s = params_.Spot;
			T const k = params_.Strike;
			T const d1 = ((std::log(s / k) + (r - q + 0.5 * sig * sig)) * (numberPeriods_ * dt)) /
				(sig * std::sqrt(numberPeriods_ * dt));
			T const d2 = d1 - sig * std::sqrt(numberPeriods_ * dt);
			T const h2 = inversion_(d2);
			T const prob = h2;
			T const disc = std::exp(-1.0 * r * dt);
			return (disc * (prob * upValue + (1.0 - prob) * downValue));
		}

		static std::string const name() {
			return std::string{ "Leisen-Reimer model" };
		}

		static constexpr AssetClass assetClass() { return AssetClass::Equity; }

	};


	// =============================================================================================
	// ======================== Black-Derman-Toy model (binomial lattice) ==========================
	// ==	MODEL:																				  ==
	// ==	d ln(r(t)) = theta(t)*dt + sigma*dw(t)												  ==
	// =============================================================================================

	template<typename T = double>
	class BlackDermanToyModel :public BinomialModel<1, T> {
	private:
		// typedef discounting factor:
		typedef DiscountingFactor<T> DCF;

		std::function<T(T, T)> dcf_;
		DiscountingStyle ds_;
		T prob_;
		std::vector<T> theta_;
		ModelParams<1, AssetClass::InterestRate, T> params_;

		std::tuple<T, T> _branching(T theta, T sig, T sqrtdt, std::size_t leafIdx)const {
			T const up = std::exp(sig*(leafIdx + 1)*sqrtdt);
			T const down = std::exp(sig*leafIdx*sqrtdt);
			return std::make_tuple(down*theta, up*theta);
		}

	public:
		BlackDermanToyModel(ModelParams<1, AssetClass::InterestRate, T> const &params, DiscountingStyle style = DiscountingStyle::Continuous)
			:params_{ params }, prob_{ 0.5 }, ds_{style} {
			dcf_ = DCF::function(style);
		}

		DiscountingStyle discountingStyle()const { return ds_; }

		// Returns risk-neutral probability:
		T nodeRiskNeutralProb()const {
			return this->prob_;
		}

		// Forward generator
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting) const override {
			LASSERT(!theta_.empty(), "Populate theta via setTheta() member function!");
			T const sig = params_.Volatility;
			T const theta = theta_.at(timeIdx - 1);
			T const sqrtdt = std::sqrt(dt);
			return _branching(theta, sig, sqrtdt, leafIdx);
		}

		// Backward generator
		T operator()(T currValue, T downValue, T upValue, T dt) override {
			T const  disc = dcf_(currValue, dt);
			return (disc * (prob_*upValue + (1.0 - prob_)*downValue));
		}

		void setTheta(std::vector<T> const &theta) { theta_ = theta; }

		static std::string const name() {
			return std::string{ "Black-Derman-Toy model" };
		}

		static constexpr AssetClass assetClass() { return AssetClass::InterestRate; }

		// Calibration minimizer:
		auto calibrationMinimizer()const {
			// analytical theta via continuous/discrete discounting is not available here 
		}


		// Calibration objective function:
		auto calibrationObjective()const {
			T const sig = params_.Volatility;

			return [=](T theta, T dt, T marketDiscount, std::vector<T> const &prevRateStates,
				std::vector<T> const &arrowDebreuStates, std::function<T(T, T)>const &dcf)->T {

				T const sqrtdt = std::sqrt(dt);
				T sum{ 0.0 };
				T rate{ 0.0 };
				std::size_t const statesSize = arrowDebreuStates.size();
				for (std::size_t i = 0; i < statesSize; ++i) {
					rate = theta * std::exp(sig*i*sqrtdt);
					sum += arrowDebreuStates.at(i) * dcf(rate, dt);
				}
				return ((sum - marketDiscount)*(sum - marketDiscount));
			};
		}

		// Calibration forward function:
		auto calibrationForwardGenerator()const {
			T const sig = params_.Volatility;

			return [=](T theta, T value, T dt, std::size_t leafIdx)->std::tuple<T, T> {
				T const sqrtdt = std::sqrt(dt);
				return _branching(theta, sig, sqrtdt, leafIdx);

			};

		}

	};

	// =============================================================================================
	// ================================ Ho-Lee model (binomial lattice) ============================
	// ==	MODEL:																				  ==
	// ==	dr(t) = theta(t)*dt + sigma*dw(t)													  ==
	// =============================================================================================

	template<typename T = double>
	class HoLeeModel :public BinomialModel<1, T> {
	private:
		// typedef discounting factor:
		typedef DiscountingFactor<T> DCF;

		std::function<T(T, T)> dcf_;
		DiscountingStyle ds_;
		T prob_;
		std::vector<T> theta_;
		ModelParams<1, AssetClass::InterestRate, T> params_;


		std::tuple<T, T> _branching(T mean, T sigdt)const {
			T const down = mean - sigdt;
			T const up = mean + sigdt;
			return std::make_tuple(down, up);
		}

	public:
		HoLeeModel(ModelParams<1, AssetClass::InterestRate, T> const &params, DiscountingStyle style = DiscountingStyle::Continuous)
			:params_{ params }, prob_{ 0.5 }, ds_{style} {
			dcf_ = DCF::function(style);
		}


		DiscountingStyle discountingStyle()const { return ds_; }

		// Returns risk-neutral probability:
		T nodeRiskNeutralProb()const {
			return this->prob_;
		}


		// Forward generator
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting ) const override {
			LASSERT(!theta_.empty(), "Populate theta via setTheta() member function!");
			T const sig = params_.Volatility;
			T const sqrtdt = std::sqrt(dt);
			T const sigdt = sig * sqrtdt;
			T const mean = value + theta_.at(timeIdx - 1)*dt;
			return _branching(mean, sigdt);
		}


		// Backward generator
		T operator()(T currValue, T downValue, T upValue, T dt) override {
			T const disc = dcf_(currValue, dt);
			return (disc * (prob_*upValue + (1.0 - prob_)*downValue));
		}

		void setTheta(std::vector<T> const &theta) { theta_ = theta; }

		static std::string const name() {
			return std::string{ "Ho-Lee model" };
		}


		static constexpr AssetClass assetClass() { return AssetClass::InterestRate; }


		// Calibration minimizer:
		auto calibrationMinimizer()const {
			T const sig = params_.Volatility;

			// analytical theta via continuous discounting:
			return [=](T theta, T dt, T marketDiscount, std::vector<T> const &prevRateStates,
				std::vector<T> const &arrowDebreuStates)->T {

				const std::size_t prevNodeSize = prevRateStates.size();
				std::vector<T> rateStates(prevNodeSize + 1);
				T const sqrtdt = std::sqrt(dt);
				T sum{ 0.0 };

				for (std::size_t l = 0; l < prevNodeSize; ++l) {
					rateStates[l] = prevRateStates.at(l) - sig * sqrtdt;     //down l
					rateStates[l + 1] = prevRateStates.at(l) + sig * sqrtdt; //up l + 1  
				}

				std::size_t const statesSize = arrowDebreuStates.size();
				for (std::size_t i = 0; i < statesSize; ++i) {
					sum += arrowDebreuStates.at(i) * std::exp(-1.0*rateStates.at(i) * dt);
				}
				sum = (sum / marketDiscount);
				return std::log(sum) / (dt*dt);
			};
		}


		// Calibration objective function:
		auto calibrationObjective()const {

			T const sig = params_.Volatility;

			return [=](T theta, T dt, T marketDiscount, std::vector<T> const &prevRateStates,
				std::vector<T> const &arrowDebreuStates, std::function<T(T, T)>const &dcf)->T {

				const std::size_t prevNodeSize = prevRateStates.size();
				std::vector<T> rateStates(prevNodeSize + 1);
				T const sqrtdt = std::sqrt(dt);
				T sum{ 0.0 };

				for (std::size_t l = 0; l < prevNodeSize; ++l) {
					rateStates[l] = prevRateStates.at(l) + (theta * dt) - sig * sqrtdt;     //down l
					rateStates[l + 1] = prevRateStates.at(l) + (theta * dt) + sig * sqrtdt; //up l + 1  
				}

				std::size_t const statesSize = arrowDebreuStates.size();
				for (std::size_t i = 0; i < statesSize; ++i) {
					sum += arrowDebreuStates.at(i) * dcf(rateStates.at(i), dt);
				}
				return ((sum - marketDiscount)*(sum - marketDiscount));
			};
		}

		// Calibration forward function:
		auto calibrationForwardGenerator()const {
			T const sig = params_.Volatility;

			return [=](T theta, T value, T dt, std::size_t leafIdx)->std::tuple<T, T> {
				T const sqrtdt = std::sqrt(dt);
				T const sigdt = sig * sqrtdt;
				T const mean = value + theta * dt;
				return _branching(mean, sigdt);
			};

		}

	};



}






#endif  ///_LATTICE_MODEL_BINOMIAL