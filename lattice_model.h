#pragma once
#if !defined(_LATTICE_MODEL)
#define _LATTICE_MODEL

#include"lattice_types.h"
#include"lattice_miscellaneous.h"

namespace lattice_model {

	using lattice_miscellaneous::OptionData;
	using lattice_types::LeafForwardGenerator;
	using lattice_types::LeafBackwardGenerator;
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;


	template<std::size_t FactorCount,typename T>
	class BinomialModel{};

	template<std::size_t FactorCount,typename T>
	class TrinomialModel{};

	template<typename T>
	class BinomialModel<1,T> {
	public:
		// Forward generator:
		virtual std::tuple<T, T> operator()(T value, T dt) = 0;

		// Backward generator:
		virtual T operator()(T currValue,T upValue, T downValue, T dt) = 0;

		// Factor count:
		enum { FactorCount = 1 };

		// LatticeType:
		LatticeType latticeType()const { return LatticeType::Binomial; }


	};

	template<typename T>
	class BinomialModel<2, T> {
	public:
		// Forward generators:
		virtual std::pair<LeafForwardGenerator<T, T, T>,
							LeafForwardGenerator<T, T, T>> forwardGenerator()const = 0;

		// Forward generator 2:
		virtual std::pair<LeafBackwardGenerator<T, T ,T, T, T>,
							LeafBackwardGenerator<T, T, T, T, T>> backwardGenerator()const = 0;
		// Factor count:
		enum { FactorCount = 2 };

		// LatticeType:
		LatticeType latticeType()const { return LatticeType::Binomial; }

	};

	template<typename T>
	class TrinomialModel<1,T> {
	public:
		// Forward generator:
		virtual std::tuple<T, T, T> operator()(T value, T dt) = 0;

		// Backward generator:
		virtual T operator()(T currValue,T upValue,T midValue, T downValue, T dt) = 0;

		// Factor count:
		enum { FactorCount = 1 };

		// LatticeType:
		LatticeType latticeType()const { return LatticeType::Trinomial; }

	};

	template<typename T>
	class TrinomialModel<2, T> {
	public:

		// Forward generators:
		virtual std::pair<LeafForwardGenerator<T, T, T, T>,
			LeafForwardGenerator<T, T, T,T>> forwardGenerator()const = 0;

		// Forward generator 2:
		virtual std::pair<LeafBackwardGenerator<T,T, T, T, T,T>,
			LeafBackwardGenerator<T,T, T, T, T,T>> backwardGenerator()const = 0;

		// Factor count:
		enum { FactorCount = 2 };

		// LatticeType:
		LatticeType latticeType()const { return LatticeType::Trinomial; }

	};

	// =============================================================================================
	// ===================== Cox-Rubinstein-Ross model (binomial lattice) ==========================
	// =============================================================================================


	template<typename T = double>
	class CoxRubinsteinRossModel:public BinomialModel<1,T> {
	private:
		T prob_;
		OptionData<T> option_;

	public:
		CoxRubinsteinRossModel(OptionData<T>const &optionData)
			:option_{ optionData }, prob_{0.5} {
		}

		// Forward generator
		std::tuple<T, T> operator()(T value, T dt) override {
			T q = option_.DividentRate;
			T sig = option_.Volatility;
			T r = option_.RiskFreeRate;
			T r1 = (r - q - 0.5*sig*sig)*dt;
			T r2 = sig * std::sqrt(dt);
			T up = std::exp(r1 + r2);
			T down = std::exp(r1 - r2);
			return std::make_tuple(up*value, down*value);
		}

		// Backward generator
		T operator()(T currValue, T upValue, T downValue, T dt) override {
			T r = option_.RiskFreeRate;
			T disc = std::exp(-1.0*r*dt);
			return (disc * (prob_*upValue + (1.0 - prob_)*downValue));
		}

		static std::string const name() {
			return std::string{ "Cox-Rubinstein-Ross model" };
		}

		static AssetClass const assetClass() { return AssetClass::Equity; }
	};


	// =======================================================================================================
	// ===================== Modified Cox-Rubinstein-Ross model (binomial lattice) ===========================
	// =======================================================================================================
	template<typename T=double>
	class ModifiedCoxRubinsteinRossModel:public BinomialModel<1,T> {
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
		T operator()(T currValue, T upValue, T downValue, T dt) override {
			T q = option_.DividentRate;
			T s = option_.Underlying;
			T sig = option_.Volatility;
			T r = option_.RiskFreeRate;
			T k = option_.Strike;
			T K_n = (std::log(k / s) / static_cast<T>(n_));
			T V_n = sig * std::sqrt(dt);
			T up = std::exp(K_n + V_n);
			T down = std::exp(K_n - V_n);
			T disc = std::exp(-1.0*r *dt);
			T prob = (std::exp((r-q)*dt) - down) / (up - down);
			return (disc * (prob*upValue + (1.0 - prob)*downValue));
		}

		static std::string const name() {
			return std::string{ "Modified Cox-Rubinstein-Ross model" };
		}
		static AssetClass const assetClass() { return AssetClass::Equity; }
	};

	// =====================================================================================
	// ===================== Jarrow-Rudd model (binomial lattice) ==========================
	// =====================================================================================
	template<typename T = double>
	class JarrowRuddModel:public BinomialModel<1,T> {
	private:
		T prob_;
		OptionData<T> option_;

	public:
		JarrowRuddModel(OptionData<T>const &optionData)
			:option_{ optionData }, prob_{ 0.5 } {
		}

		// Forward generator
		std::tuple<T, T> operator()(T value, T dt)override {
			T q = option_.DividentRate;
			T sig = option_.Volatility;
			T r = option_.RiskFreeRate;
			T d = (r - q - 0.5*sig*sig)*dt;
			T x1 = d + sig * std::sqrt(dt);
			T x2 = d - sig * std::sqrt(dt);
			T up = std::exp(x1);
			T down = std::exp(x2);
			return std::make_tuple(up*value, down*value);
		}

		// Backward generator
		T operator()(T currValue, T upValue, T downValue, T dt) override {
			T r = option_.RiskFreeRate;
			T disc = std::exp(-1.0*r *dt);
			return (disc * (prob_*upValue + (1.0 - prob_)*downValue));
		}

		static std::string const name(){
			return std::string{ "Jarrow-Rudd model" };
		}
		static AssetClass const assetClass() { return AssetClass::Equity; }
	};

	// =====================================================================================
	// ===================== Trigeorgis model (binomial lattice) ===========================
	// =====================================================================================
	template<typename T=double>
	class TrigeorgisModel:public BinomialModel<1, T> {
	private:
		T gamma_;
		OptionData<T> option_;

	public:
		TrigeorgisModel(OptionData<T> const &optionData)
			:option_{ optionData },
			gamma_{ optionData.RiskFreeRate - optionData.DividentRate - 0.5*optionData.Volatility*optionData.Volatility } {
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
		T operator()(T currValue, T upValue, T downValue, T dt)override {
			T sig = option_.Volatility;
			T r = option_.RiskFreeRate;
			T x = std::sqrt(sig*sig*dt + gamma_ * gamma_*dt*dt);
			T disc = std::exp(-1.0*r *dt);
			T prob = 0.5*(1.0 + (gamma_*(dt / x)));
			return (disc * (prob*upValue + (1.0 - prob)*downValue));
		}

		static std::string const name()  {
			return std::string{ "Trigeorgis model" };
		}

		static AssetClass const assetClass() { return AssetClass::Equity; }
	};

	// ===============================================================================
	// ===================== Tian model (binomial lattice) ===========================
	// ===============================================================================
	template<typename T = double>
	class TianModel:public BinomialModel<1, T> {
	private:
		OptionData<T> option_;

	public:
		TianModel(OptionData<T> const &optionData)
			:option_{ optionData } {
		}

		// Forward generator:
		std::tuple<T, T> operator()(T value, T dt) override {
			T q = option_.DividentRate;
			T sig = option_.Volatility;
			T r = option_.RiskFreeRate;
			T v = std::exp(sig*sig*dt);
			T x = std::sqrt(v *v + 2.0*v - 3.0);
			T up = (0.5*std::exp((r- q)*dt)*v*(v + 1.0 + x));
			T down = (0.5*std::exp((r-q)*dt)*v*(v + 1.0 - x));
			return std::make_tuple(up*value, down*value);
		}

		// Backward generator:
		T operator()(T currValue,T upValue, T downValue, T dt)override {
			T q = option_.DividentRate;
			T sig = option_.Volatility;
			T r = option_.RiskFreeRate;
			T v = std::exp(sig*sig*dt);
			T x = std::sqrt(v *v + 2.0*v - 3.0);
			T up = (0.5*std::exp((r - q)*dt)*v*(v + 1.0 + x));
			T down = (0.5*std::exp((r - q)*dt)*v*(v + 1.0 - x));
			T disc = std::exp(-1.0*r*dt);
			T prob = ((std::exp((r - q)*dt) - down) / (up - down));
			return (disc * (prob*upValue + (1.0 - prob)*downValue));
		}

		static std::string const name(){
			return std::string{ "Tian model" };
		}

		static AssetClass const assetClass() { return AssetClass::Equity; }
	};

	// ========================================================================================
	// ===================== Leisen-Reimer model (binomial lattice) ===========================
	// ========================================================================================
	template<typename T =double>
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
		T operator()(T currValue, T upValue, T downValue, T dt)override {
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

		static AssetClass const assetClass() { return AssetClass::Equity; }

	};


	// =============================================================================================
	// ======================== Black-Derman-Toy model (binomial lattice) ==========================
	// =============================================================================================


	template<typename T = double>
	class BlackDermanToyModel :public BinomialModel<1, T> {
	private:
		T prob_;
		std::vector<T> theta_;
		OptionData<T> option_;
		std::size_t thetaCounter_;

	public:
		BlackDermanToyModel(OptionData<T>const &optionData,std::vector<T> const &theta)
			:option_{ optionData }, prob_{ 0.5 },
			theta_{ theta }, thetaCounter_{ 0 } {
		}

		// Forward generator
		std::tuple<T, T> operator()(T value, T dt) override {
			T sig = option_.Volatility;
			T sqrtdt = std::sqrt(dt);
			T up = std::exp(sig*sqrtdt);
			T down = 1.0 / up;
			T theta = theta_.at(thetaCounter_++);
			return std::make_tuple(up*theta*value, down*theta*value);
		}

		// Backward generator
		T operator()(T currValue, T upValue, T downValue, T dt) override {
			T disc = std::exp(-1.0*currValue*dt);
			return (disc * (prob_*upValue + (1.0 - prob_)*downValue));
		}

		static std::string const name() {
			return std::string{ "Black-Derman-Toy model" };
		}

		static AssetClass const assetClass() { return AssetClass::InterestRate; }
	};


	// ================================================================================
	// ===================== Boyle model (trinomial lattice) ==========================
	// ================================================================================
	template<typename T = double>
	class BoyleModel :public TrinomialModel<1, T> {
	private:
		OptionData<T> option_;

	public:
		BoyleModel(OptionData<T>const &optionData)
			:option_{ optionData } {
		}

		// Forward generator
		std::tuple<T, T, T> operator()(T value, T dt) override {
			T sig = option_.Volatility;
			T expon = sig * std::sqrt(2.0*dt);
			T up = std::exp(expon);
			T mid = 1.0;
			T down = 1.0 / up;
			return std::make_tuple(up*value, mid*value, down*value);
		}

		// Backward generator
		T operator()(T currValue, T upValue,T midValue,T downValue, T dt) override {
			T q = option_.DividentRate;
			T r = option_.RiskFreeRate;
			T sig = option_.Volatility;
			T mu = r - q;
			T e_sig = std::exp(sig * std::sqrt(0.5*dt));
			T e_mu = std::exp(mu * 0.5*dt);
			T p_u = std::pow(((e_mu - (1.0/e_sig))/(e_sig - (1.0/e_sig))), 2.0);
			T p_d = std::pow(((e_sig - e_mu) / (e_sig - (1.0 / e_sig))), 2.0);
			T p_m = 1.0 - (p_u + p_d);
			T disc = std::exp(-1.0*r*dt);
			return (disc * (p_u * upValue + p_m * midValue + p_d * downValue));
		}

		static std::string const name() {
			return std::string{ "Boyle model" };
		}

		static AssetClass const assetClass() { return AssetClass::Equity; }
	};

}




#endif //_LATTICE_MODEL