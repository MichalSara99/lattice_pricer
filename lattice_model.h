#pragma once
#if !defined(_LATTICE_MODEL)
#define _LATTICE_MODEL

#include"lattice_miscellaneous.h"

namespace lattice_model {

	using lattice_miscellaneous::OptionData;


	template<typename T=double>
	class BinomialModel {
	public:
		// Forward generator:
		virtual std::tuple<T, T> operator()(T value, T dt) = 0;

		// Backward generator:
		virtual T operator()(T upValue, T downValue, T dt) = 0;

	};

	// ===================== Cox-Rubinstein-Ross model (binomial lattice) ==========================
	template<typename T=double>
	class CoxRubinsteinRossModel:public BinomialModel<T> {
	private:
		T prob_;
		OptionData<T> option_;

	public:
		CoxRubinsteinRossModel(OptionData<T>const &optionData)
			:option_{ optionData }, prob_{0.5} {
		}

		// Forward generator
		std::tuple<T, T> operator()(T value, T dt) override {
			T s = option_.Volatility;
			T r = option_.RiskFreeRate;
			T r1 = (r - 0.5*s*s)*dt;
			T r2 = s * std::sqrt(dt);
			T up = std::exp(r1 + r2);
			T down = std::exp(r1 - r2);
			return std::make_tuple(up*value, down*value);
		}

		// Backward generator
		T operator()(T upValue, T downValue, T dt) override {
			T r = option_.RiskFreeRate;
			T disc = std::exp(-1.0*r*dt);
			return (disc * (prob_*upValue + (1.0 - prob_)*downValue));
		}

		static std::string const name() {
			return std::string{ "Cox-Rubinstein-Ross model" };
		}
	};


	// ===================== Modified Cox-Rubinstein-Ross model (binomial lattice) ===========================
	template<typename T=double>
	class ModifiedCoxRubinsteinRossModel:public BinomialModel<T> {
	private:
		OptionData<T> option_;
		std::size_t n_;

	public:
		ModifiedCoxRubinsteinRossModel(OptionData<T> const &optionData, std::size_t numberPeriods)
			:option_{ optionData },n_{numberPeriods} {
		}

		// Forward generator:
		std::tuple<T, T> operator()(T value, T dt)override {
			T s = option_.Underlying;
			T sig = option_.Volatility;
			T k = option_.Strike;
			T K_n = (std::log(k / s) / static_cast<T>(n_));
			T V_n = sig * std::sqrt(dt);
			T up = std::exp(K_n + V_n);
			T down = std::exp(K_n - V_n);
			return std::make_tuple(up*value, down*value);
		}

		// Backward generator:
		T operator()(T upValue, T downValue, T dt) override {
			T s = option_.Underlying;
			T sig = option_.Volatility;
			T r = option_.RiskFreeRate;
			T k = option_.Strike;
			T K_n = (std::log(k / s) / static_cast<T>(n_));
			T V_n = sig * std::sqrt(dt);
			T up = std::exp(K_n + V_n);
			T down = std::exp(K_n - V_n);
			T disc = std::exp(-1.0*r*dt);
			T prob = (std::exp(r*dt) - down) / (up - down);
			return (disc * (prob*upValue + (1.0 - prob)*downValue));
		}

		static std::string const name() {
			return std::string{ "Modified Cox-Rubinstein-Ross model" };
		}
	};

	// ===================== Jarrow-Rudd model (binomial lattice) ==========================
	template<typename T = double>
	class JarrowRuddModel:public BinomialModel<T> {
	private:
		T prob_;
		OptionData<T> option_;

	public:
		JarrowRuddModel(OptionData<T>const &optionData)
			:option_{ optionData }, prob_{ 0.5 } {
		}

		// Forward generator
		std::tuple<T, T> operator()(T value, T dt)override {
			T sig = option_.Volatility;
			T r = option_.RiskFreeRate;
			T d = (r - 0.5*sig*sig)*dt;
			T x1 = d + sig * std::sqrt(dt);
			T x2 = d - sig * std::sqrt(dt);
			T up = std::exp(x1);
			T down = std::exp(x2);
			return std::make_tuple(up*value, down*value);
		}

		// Backward generator
		T operator()(T upValue, T downValue, T dt) override {
			T r = option_.RiskFreeRate;
			T disc = std::exp(-1.0*r*dt);
			return (disc * (prob_*upValue + (1.0 - prob_)*downValue));
		}

		static std::string const name(){
			return std::string{ "Jarrow-Rudd model" };
		}
	};

	// ===================== Trigeorgis model (binomial lattice) ===========================
	template<typename T=double>
	class TrigeorgisModel:public BinomialModel<T> {
	private:
		T gamma_;
		OptionData<T> option_;

	public:
		TrigeorgisModel(OptionData<T> const &optionData)
			:option_{ optionData },
			gamma_{ optionData.RiskFreeRate - 0.5*optionData.Volatility*optionData.Volatility } {
		}

		// Forward generator:
		std::tuple<T, T> operator()(T value, T dt) override {
			T sig = option_.Volatility;
			T x = std::sqrt(sig*sig*dt + gamma_ * gamma_*dt*dt);
			T up = std::exp(x);
			T down = (1.0 / up);
			return std::make_tuple(up*value, down*value);
		}

		// Backward generator:
		T operator()(T upValue, T downValue, T dt)override {
			T sig = option_.Volatility;
			T r = option_.RiskFreeRate;
			T x = std::sqrt(sig*sig*dt + gamma_ * gamma_*dt*dt);
			T disc = std::exp(-1.0*r*dt);
			T prob = 0.5*(1.0 + (gamma_*(dt / x)));
			return (disc * (prob*upValue + (1.0 - prob)*downValue));
		}

		static std::string const name()  {
			return std::string{ "Trigeorgis model" };
		}
	};

	// ===================== Tian model (binomial lattice) ===========================
	template<typename T = double>
	class TianModel:public BinomialModel<T> {
	private:
		OptionData<T> option_;

	public:
		TianModel(OptionData<T> const &optionData)
			:option_{ optionData } {
		}

		// Forward generator:
		std::tuple<T, T> operator()(T value, T dt) override {
			T sig = option_.Volatility;
			T r = option_.RiskFreeRate;
			T v = std::exp(sig*sig*dt);
			T x = std::sqrt(v *v + 2.0*v - 3.0);
			T up = (0.5*std::exp(r*dt)*v*(v + 1.0 + x));
			T down = (0.5*std::exp(r*dt)*v*(v + 1.0 - x));
			return std::make_tuple(up*value, down*value);
		}

		// Backward generator:
		T operator()(T upValue, T downValue, T dt)override {
			T sig = option_.Volatility;
			T r = option_.RiskFreeRate;
			T v = std::exp(sig*sig*dt);
			T x = std::sqrt(v *v + 2.0*v - 3.0);
			T up = (0.5*std::exp(r*dt)*v*(v + 1.0 + x));
			T down = (0.5*std::exp(r*dt)*v*(v + 1.0 - x));
			T disc = std::exp(-1.0*r*dt);
			T prob = ((std::exp(r*dt) - down) / (up - down));
			return (disc * (prob*upValue + (1.0 - prob)*downValue));
		}

		static std::string const name(){
			return std::string{ "Tian model" };
		}
	};

	// ===================== Leisen-Reimer model (binomial lattice) ===========================
	template<typename T =double>
	class LeisenReimerModel :public BinomialModel<T> {
	private:
		OptionData<T> option_;
		std::size_t numberPeriods_;
		std::function<T(T)> inversion_;

	public:
		LeisenReimerModel(OptionData<T> const& optionData,
			std::size_t numberPeriods,
			std::function<T(T)> const &inversionFormula)
			:option_{ optionData }, 
			numberPeriods_{numberPeriods},
			inversion_{ inversionFormula }{
		}

		// Forward generator:
		std::tuple<T, T> operator()(T value, T dt) override {
			T sig = option_.Volatility;
			T r = option_.RiskFreeRate;
			T q = option_.DividentRate;
			T s = option_.Underlying;
			T k = option_.Strike;
			T d1 = ((std::log(s / k) + (r - q + 0.5 * sig * sig)) * (numberPeriods_ * dt)) /
				(sig * std::sqrt(numberPeriods_ * dt));
			T d2 = d1 - sig * std::sqrt(numberPeriods_ * dt);
			T e = std::exp((r - q) * dt);
			T h1 = inversion_(d1);
			T h2 = inversion_(d2);
			T p = h2;
			T up = e * (h1 / h2);
			T down = (e - p * up) / (1.0 - p);
			return std::make_tuple(up * value, down * value);
		}

		// Backward generator:
		T operator()(T upValue, T downValue, T dt)override {
			T sig = option_.Volatility;
			T r = option_.RiskFreeRate;
			T q = option_.DividentRate;
			T s = option_.Underlying;
			T k = option_.Strike;
			T d1 = ((std::log(s / k) + (r - q + 0.5 * sig * sig)) * (numberPeriods_ * dt)) /
				(sig * std::sqrt(numberPeriods_ * dt));
			T d2 = d1 - sig * std::sqrt(numberPeriods_ * dt);
			T h2 = inversion_(d2);
			T prob = h2;
			T disc = std::exp(-1.0 * r * dt);
			return (disc * (prob * upValue + (1.0 - prob) * downValue));
		}

		static std::string const name() {
			return std::string{ "Leisen-Reimer model" };
		}

	};









}




#endif //_LATTICE_MODEL