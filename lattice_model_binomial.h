#pragma once
#if !defined(_LATTICE_MODEL_BINOMIAL)
#define  _LATTICE_MODEL_BINOMIAL

#include"lattice_model_interface.h"


namespace lattice_model {



	// =============================================================================================
	// ===================== Cox-Rubinstein-Ross model (binomial lattice) ==========================
	// =============================================================================================

	template<typename T = double>
	class CoxRubinsteinRossModel :public BinomialModel<1, T> {
	private:
		T prob_;
		OptionData<T> option_;

	public:
		CoxRubinsteinRossModel(OptionData<T>const &optionData)
			:option_{ optionData }, prob_{ 0.5 } {
		}

		// Forward generator
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false) override {
			T const q = option_.DividentRate;
			T const sig = option_.Volatility;
			T const r = option_.RiskFreeRate;
			T const r1 = (r - q - 0.5*sig*sig)*dt;
			T const r2 = sig * std::sqrt(dt);
			T const up = std::exp(r1 + r2);
			T const down = std::exp(r1 - r2);
			return std::make_tuple(up*value, down*value);
		}

		// Backward generator
		T operator()(T currValue, T upValue, T downValue, T dt) override {
			T const r = option_.RiskFreeRate;
			T const disc = std::exp(-1.0*r*dt);
			return (disc * (prob_*upValue + (1.0 - prob_)*downValue));
		}

		static std::string const name() {
			return std::string{ "Cox-Rubinstein-Ross model" };
		}

		static constexpr AssetClass assetClass() { return AssetClass::Equity; }
	};


	// =======================================================================================================
	// ===================== Modified Cox-Rubinstein-Ross model (binomial lattice) ===========================
	// =======================================================================================================

	template<typename T = double>
	class ModifiedCoxRubinsteinRossModel :public BinomialModel<1, T> {
	private:
		OptionData<T> option_;
		std::size_t n_;

	public:
		ModifiedCoxRubinsteinRossModel(OptionData<T> const &optionData, std::size_t numberPeriods)
			:option_{ optionData }, n_{ numberPeriods } {
		}

		// Forward generator:
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false)override {
			T const s = option_.Underlying;
			T const sig = option_.Volatility;
			T const k = option_.Strike;
			T const K_n = (std::log(k / s) / static_cast<T>(n_));
			T const V_n = sig * std::sqrt(dt);
			T const up = std::exp(K_n + V_n);
			T const down = std::exp(K_n - V_n);
			return std::make_tuple(up*value, down*value);
		}

		// Backward generator:
		T operator()(T currValue, T upValue, T downValue, T dt) override {
			T const q = option_.DividentRate;
			T const s = option_.Underlying;
			T const sig = option_.Volatility;
			T const r = option_.RiskFreeRate;
			T const k = option_.Strike;
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
		OptionData<T> option_;

	public:
		JarrowRuddModel(OptionData<T>const &optionData)
			:option_{ optionData }, prob_{ 0.5 } {
		}

		// Forward generator
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false)override {
			T const q = option_.DividentRate;
			T const sig = option_.Volatility;
			T const r = option_.RiskFreeRate;
			T const d = (r - q - 0.5*sig*sig)*dt;
			T const x1 = d + sig * std::sqrt(dt);
			T const x2 = d - sig * std::sqrt(dt);
			T const up = std::exp(x1);
			T const down = std::exp(x2);
			return std::make_tuple(up*value, down*value);
		}

		// Backward generator
		T operator()(T currValue, T upValue, T downValue, T dt) override {
			T const r = option_.RiskFreeRate;
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
		OptionData<T> option_;

	public:
		TrigeorgisModel(OptionData<T> const &optionData)
			:option_{ optionData },
			gamma_{ optionData.RiskFreeRate - optionData.DividentRate - 0.5*optionData.Volatility*optionData.Volatility } {
		}

		// Forward generator:
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false) override {
			T const sig = option_.Volatility;
			T const x = std::sqrt(sig*sig*dt + gamma_ * gamma_*dt*dt);
			T const up = std::exp(x);
			T const down = (1.0 / up);
			return std::make_tuple(up*value, down*value);
		}

		// Backward generator:
		T operator()(T currValue, T upValue, T downValue, T dt)override {
			T const sig = option_.Volatility;
			T const r = option_.RiskFreeRate;
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
		OptionData<T> option_;

	public:
		TianModel(OptionData<T> const &optionData)
			:option_{ optionData } {
		}

		// Forward generator:
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false) override {
			T const q = option_.DividentRate;
			T const sig = option_.Volatility;
			T const r = option_.RiskFreeRate;
			T const v = std::exp(sig*sig*dt);
			T const x = std::sqrt(v *v + 2.0*v - 3.0);
			T const up = (0.5*std::exp((r - q)*dt)*v*(v + 1.0 + x));
			T const down = (0.5*std::exp((r - q)*dt)*v*(v + 1.0 - x));
			return std::make_tuple(up*value, down*value);
		}

		// Backward generator:
		T operator()(T currValue, T upValue, T downValue, T dt)override {
			T const q = option_.DividentRate;
			T const sig = option_.Volatility;
			T const r = option_.RiskFreeRate;
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
		OptionData<T> option_;
		std::size_t numberPeriods_;
		std::function<T(T)> inversion_;

	public:
		LeisenReimerModel(OptionData<T> const& optionData,
			std::size_t numberPeriods,
			std::function<T(T)> const &inversionFormula)
			:option_{ optionData },
			numberPeriods_{ numberPeriods },
			inversion_{ inversionFormula } {
		}

		// Forward generator:
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false) override {
			T const sig = option_.Volatility;
			T const r = option_.RiskFreeRate;
			T const q = option_.DividentRate;
			T const s = option_.Underlying;
			T const k = option_.Strike;
			T const d1 = ((std::log(s / k) + (r - q + 0.5 * sig * sig)) * (numberPeriods_ * dt)) /
				(sig * std::sqrt(numberPeriods_ * dt));
			T const d2 = d1 - sig * std::sqrt(numberPeriods_ * dt);
			T const e = std::exp((r - q) * dt);
			T const h1 = inversion_(d1);
			T const h2 = inversion_(d2);
			T const p = h2;
			T const up = e * (h1 / h2);
			T const down = (e - p * up) / (1.0 - p);
			return std::make_tuple(up * value, down * value);
		}

		// Backward generator:
		T operator()(T currValue, T upValue, T downValue, T dt)override {
			T const sig = option_.Volatility;
			T const r = option_.RiskFreeRate;
			T const q = option_.DividentRate;
			T const s = option_.Underlying;
			T const k = option_.Strike;
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
		T prob_;
		std::vector<T> theta_;
		OptionData<T> option_;

		std::tuple<T, T> _branching(T theta, T sig, T sqrtdt, std::size_t leafIdx)const {
			T const up = std::exp(sig*(leafIdx + 1)*sqrtdt);
			T const down = std::exp(sig*leafIdx*sqrtdt);
			return std::make_tuple(down*theta, up*theta);
		}

	public:
		BlackDermanToyModel(OptionData<T>const &optionData)
			:option_{ optionData }, prob_{ 0.5 } {
		}

		// Returns risk-neutral probability:
		T nodeRiskNeutralProb()const {
			return this->prob_;
		}

		// Forward generator
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false) override {
			LASSERT(!theta_.empty(), "Populate theta via setTheta() member function!");
			T const sig = option_.Volatility;
			T const theta = theta_.at(timeIdx - 1);
			T const sqrtdt = std::sqrt(dt);
			return _branching(theta, sig, sqrtdt, leafIdx);
		}

		// Backward generator
		T operator()(T currValue, T upValue, T downValue, T dt) override {
			T const disc = std::exp(-1.0*currValue*dt);
			return (disc * (prob_*upValue + (1.0 - prob_)*downValue));
		}

		void setTheta(std::vector<T> const &theta) { theta_ = theta; }

		static std::string const name() {
			return std::string{ "Black-Derman-Toy model" };
		}

		OptionData<T> const &optionData()const { return option_; }

		static constexpr AssetClass assetClass() { return AssetClass::InterestRate; }

		// Calibration minimizer:
		auto calibrationMinimizer()const {
			// analytical theta via continuous/discrete discounting is not available here 
		}


		// Calibration objective function:
		auto calibrationObjective()const {
			T const sig = option_.Volatility;

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
			T const sig = option_.Volatility;

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
		T prob_;
		std::vector<T> theta_;
		OptionData<T> option_;


		std::tuple<T, T> _branching(T mean, T sigdt)const {
			T const down = mean - sigdt;
			T const up = mean + sigdt;
			return std::make_tuple(down, up);
		}

	public:
		HoLeeModel(OptionData<T>const &optionData)
			:option_{ optionData }, prob_{ 0.5 } {
		}


		// Returns risk-neutral probability:
		T nodeRiskNeutralProb()const {
			return this->prob_;
		}


		// Forward generator
		std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false) override {
			LASSERT(!theta_.empty(), "Populate theta via setTheta() member function!");
			T const sig = option_.Volatility;
			T const sqrtdt = std::sqrt(dt);
			T const sigdt = sig * sqrtdt;
			T const mean = value + theta_.at(timeIdx - 1)*dt;
			return _branching(mean, sigdt);
		}


		// Backward generator
		T operator()(T currValue, T upValue, T downValue, T dt) override {
			T const disc = std::exp(-1.0*currValue*dt);
			return (disc * (prob_*upValue + (1.0 - prob_)*downValue));
		}

		void setTheta(std::vector<T> const &theta) { theta_ = theta; }

		static std::string const name() {
			return std::string{ "Ho-Lee model" };
		}

		OptionData<T> const &optionData()const { return option_; }

		static constexpr AssetClass assetClass() { return AssetClass::InterestRate; }


		// Calibration minimizer:
		auto calibrationMinimizer()const {
			T const sig = option_.Volatility;

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

			T const sig = option_.Volatility;

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
			T const sig = option_.Volatility;

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