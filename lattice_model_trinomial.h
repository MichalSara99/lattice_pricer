#pragma once
#if !defined(_LATTICE_MODEL_TRINOMIAL)
#define _LATTICE_MODEL_TRINOMIAL


#include"lattice_model_interface.h"
#include"lattice_types.h"
#include"lattice_utility.h"


namespace lattice_model {

	using lattice_types::DiscountingStyle;
	using lattice_utility::DiscountingFactor;


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
		std::tuple<T, T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false) override {
			T const sig = option_.Volatility;
			T const expon = sig * std::sqrt(2.0*dt);
			T const up = std::exp(expon);
			T const mid = 1.0;
			T const down = 1.0 / up;
			return std::make_tuple(up*value, mid*value, down*value);
		}

		// Backward generator
		T operator()(T currValue, T upValue, T midValue, T downValue, T dt) override {
			T const q = option_.DividentRate;
			T const r = option_.RiskFreeRate;
			T const sig = option_.Volatility;
			T const mu = r - q;
			T const e_sig = std::exp(sig * std::sqrt(0.5*dt));
			T const e_mu = std::exp(mu * 0.5*dt);
			T const p_u = std::pow(((e_mu - (1.0 / e_sig)) / (e_sig - (1.0 / e_sig))), 2.0);
			T const p_d = std::pow(((e_sig - e_mu) / (e_sig - (1.0 / e_sig))), 2.0);
			T const p_m = 1.0 - (p_u + p_d);
			T const disc = std::exp(-1.0*r*dt);
			return (disc * (p_u * upValue + p_m * midValue + p_d * downValue));
		}

		static std::string const name() {
			return std::string{ "Boyle model" };
		}

		static constexpr AssetClass assetClass() { return AssetClass::Equity; }
	};



	// =============================================================================================
	// ============================ Hull-White model (trinomial lattice) ===========================
	// ==	MODEL:																				  ==
	// ==	dr(t) = (theta(t) - a*r(t))*dt + sigma * dw(t)										  ==
	// =============================================================================================

	template<typename T = double>
	class HullWhiteModel :public TrinomialModel<1, T> {
	private:
		// typedef discounting factor:
		typedef DiscountingFactor<T> DCF;

		std::function<T(T, T)> dcf_;
		DiscountingStyle ds_;
		std::vector<T> theta_;
		OptionData<T> option_;

		std::tuple<T, T, T> _riskNeutralProb(std::size_t revertBranchesSize, std::size_t nodesSize, std::size_t leafIdx, T dt)const {
			T const sig = option_.Volatility;
			T const a = option_.ReversionSpeed;
			T const sqrtdt = std::sqrt(3.0*dt);
			T const dr = sig * sqrtdt;
			T const denom = dr * dr;
			std::size_t const decrement{ (nodesSize - 1) / 2 };
			long const i = (leafIdx - decrement);
			T const eps = sign(i)*std::floor(std::abs<T>(i) / revertBranchesSize);
			T const nu = (-1.0*a*dt*i + eps)*dr;
			T const numer = (nu * nu + sig * sig*dt);
			return std::make_tuple(0.5*((numer / denom) - (nu / dr)),
				(1.0 - (numer / denom)),
				0.5*((numer / denom) + (nu / dr)));
		}

		std::tuple<T, T, T> _branching(T mean, T dr, std::size_t leafIdx, bool isMeanReverting = false)const {
			if ((leafIdx == 0) && (isMeanReverting == true)) {
				// we are at the bottom of the tree -> going up with branching:
				return std::make_tuple(mean/*low*/, mean + 1.0*dr/*mid*/, mean + 2.0*dr/*high*/);
			}
			if (isMeanReverting) {
				// we are at the top of the tree -> going down with branching:
				return std::make_tuple(mean - 2.0*dr /*low*/, mean - 1.0*dr/*mid*/, mean/*high*/);
			}
			// Normal branching here:
			return std::make_tuple(mean - 1.0*dr /*low*/, mean/*mid*/, mean + 1.0*dr/*high*/);
		}

	public:
		HullWhiteModel(OptionData<T>const &optionData, DiscountingStyle style = DiscountingStyle::Continuous)
			:option_{ optionData }, ds_{style} 
		{
			dcf_ = DCF::function(style);
		}

		DiscountingStyle discountingStyle()const { return ds_; }

		// Returns tuple of probabilities:
		std::tuple<T, T, T> nodeRiskNeutralProb(std::size_t revertBranchesSize, std::size_t nodesSize, std::size_t leafIdx, T dt)const {
			return _riskNeutralProb(revertBranchesSize, nodesSize, leafIdx, dt);
		}

		// Forward generator
		std::tuple<T, T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false) override {
			LASSERT(!theta_.empty(), "Populate theta via setTheta() member function!");
			T const sig = option_.Volatility;
			T const a = option_.ReversionSpeed;
			T const sqrtdt = std::sqrt(3.0*dt);
			T const dr = sig * sqrtdt;
			T const mean = value + theta_.at(timeIdx - 1);
			return _branching(mean, dr, leafIdx, isMeanReverting);
		}

		// Backward generator
		T operator()(T currValue, T upValue, T midValue, T downValue, T dt) override {
			T const disc = dcf_(currValue,dt);
			return (disc * (0.5*upValue + (1.0 - 0.5)*downValue));
		}

		void setTheta(std::vector<T> const &theta) { theta_ = theta; }

		static std::string const name() {
			return std::string{ "Hull-White model" };
		}

		OptionData<T> const &optionData()const { return option_; }

		static constexpr AssetClass assetClass() { return AssetClass::InterestRate; }


		// Calibration objective function:
		std::function<T(T, T, T, std::vector<T> const&, std::vector<T> const&, std::function<T(T, T)>const &)>
			calibrationObjective(BranchingStyle branchingStyle)const {
			T const sig = option_.Volatility;
			T const a = option_.ReversionSpeed;

			if (branchingStyle == BranchingStyle::Normal) {
				// Calibration normal:
				return [=](T theta, T dt, T marketDiscount, std::vector<T> const &prevRateStates,
					std::vector<T> const &arrowDebreuStates, std::function<T(T, T)>const &dcf)->T {

					T const sqrtdt = std::sqrt(3.0*dt);
					T const dr = sig * sqrtdt;
					T mean{};
					const std::size_t prevNodeSize = prevRateStates.size();
					std::vector<T> rateStates(prevNodeSize + 2);

					for (std::size_t l = 0; l < prevNodeSize; ++l) {
						mean = prevRateStates.at(l) + theta;
						rateStates[l] = mean - dr;  // down
						rateStates[l + 1] = mean;  //mid
						rateStates[l + 2] = mean + dr; //high
					}

					T sum{};
					std::size_t const statesSize = arrowDebreuStates.size();
					for (std::size_t i = 0; i < statesSize; ++i) {
						sum += arrowDebreuStates.at(i) * dcf(rateStates.at(i), dt);
					}
					return ((sum - marketDiscount)*(sum - marketDiscount));
				};
			}
			else {

				// Calibration reverting:
				return [=](T theta, T dt, T marketDiscount, std::vector<T> const &prevRateStates,
					std::vector<T> const &arrowDebreuStates, std::function<T(T, T)>const &dcf)->T {

					T const sqrtdt = std::sqrt(3.0*dt);
					T const dr = sig * sqrtdt;
					T mean{};
					const std::size_t prevNodeSize = prevRateStates.size();
					std::vector<T> rateStates(prevNodeSize);

					// Branching upward here:
					mean = prevRateStates.at(0) + theta;
					rateStates[0] = mean;  // down
					rateStates[1] = mean + dr;  //mid
					rateStates[2] = mean + 2.0 * dr; //high
													 // normal branching here:
					for (std::size_t l = 1; l < prevNodeSize - 1; ++l) {
						mean = prevRateStates.at(l) + theta;
						rateStates[l - 1] = mean - dr;  // down
						rateStates[l] = mean;  //mid
						rateStates[l + 1] = mean + dr; //high
					}
					// Branching downward here:
					mean = prevRateStates.at(prevNodeSize - 1) + theta;
					rateStates[prevNodeSize - 1 - 2] = mean - 2.0*dr;  // down
					rateStates[prevNodeSize - 1 - 1] = mean - dr;  //mid
					rateStates[prevNodeSize - 1] = mean; //high

					T sum{};
					std::size_t const statesSize = arrowDebreuStates.size();
					for (std::size_t i = 0; i < statesSize; ++i) {
						sum += arrowDebreuStates.at(i) * dcf(rateStates.at(i), dt);
					}
					return ((sum - marketDiscount)*(sum - marketDiscount));
				};
			}
		}

		// Calibration forward function:
		auto calibrationForwardGenerator()const {
			T const sig = option_.Volatility;
			T const a = option_.ReversionSpeed;

			return [=](T theta, T value, T dt, std::size_t leafIdx, bool isMeanReverting = false)->std::tuple<T, T, T> {
				T const sqrtdt = std::sqrt(3.0*dt);
				T const dr = sig * sqrtdt;
				T const mean = value + theta;
				return _branching(mean, dr, leafIdx, isMeanReverting);
			};
		}

	};


	// =============================================================================================
	// ======================= Black-Karasinski model (trinomial lattice) ==========================
	// ==	MODEL:																				  ==
	// ==	d ln(r(t)) = (theta(t) - a*ln(r(t)))*dt + sigma * dw(t)							      ==
	// =============================================================================================

	template<typename T = double>
	class BlackKarasinskiModel :public TrinomialModel<1, T> {
	private:
		// typedef discounting factor:
		typedef DiscountingFactor<T> DCF;

		std::function<T(T, T)> dcf_;
		DiscountingStyle ds_;
		std::vector<T> theta_;
		OptionData<T> option_;

		std::tuple<T, T, T> _riskNeutralProb(std::size_t revertBranchesSize, std::size_t nodesSize, std::size_t leafIdx, T dt)const {
			T const sig = option_.Volatility;
			T const a = option_.ReversionSpeed;
			T const sqrtdt = std::sqrt(3.0*dt);
			T const dr = sig * sqrtdt;
			T const denom = dr * dr;
			std::size_t const decrement{ (nodesSize - 1) / 2 };
			long const i = (leafIdx - decrement);
			T const eps = sign(i)*std::floor(std::abs<T>(i) / revertBranchesSize);
			T const nu = (-1.0*a*dt*i + eps)*dr;
			T const numer = (nu * nu + sig * sig*dt);
			return std::make_tuple(0.5*((numer / denom) - (nu / dr)),
				(1.0 - (numer / denom)),
				0.5*((numer / denom) + (nu / dr)));
		}

		std::tuple<T, T, T> _branching(T mean, T dr, std::size_t leafIdx, bool isMeanReverting = false)const {
			if ((leafIdx == 0) && (isMeanReverting == true)) {
				// we are at the bottom of the tree -> going up with branching:
				return std::make_tuple(std::exp(mean)/*low*/, std::exp(mean + 1.0*dr)/*mid*/, std::exp(mean + 2.0*dr)/*high*/);
			}
			if (isMeanReverting) {
				// we are at the top of the tree -> going down with branching:
				return std::make_tuple(std::exp(mean - 2.0*dr) /*low*/, std::exp(mean - 1.0*dr)/*mid*/, std::exp(mean)/*high*/);
			}
			// Normal branching here:
			return std::make_tuple(std::exp(mean - 1.0*dr) /*low*/, std::exp(mean)/*mid*/, std::exp(mean + 1.0*dr)/*high*/);
		}

	public:
		BlackKarasinskiModel(OptionData<T>const &optionData, DiscountingStyle style = DiscountingStyle::Continuous)
			:option_{ optionData }, ds_{style} {
			dcf_ = DCF::function(style);
		}

		DiscountingStyle discountingStyle()const { return ds_; }

		// Returns tuple of probabilities:
		std::tuple<T, T, T> nodeRiskNeutralProb(std::size_t revertBranchesSize, std::size_t nodesSize, std::size_t leafIdx, T dt)const {
			return _riskNeutralProb(revertBranchesSize, nodesSize, leafIdx, dt);
		}

		// Forward generator
		std::tuple<T, T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting = false) override {
			LASSERT(!theta_.empty(), "Populate theta via setTheta() member function!");
			T const sig = option_.Volatility;
			T const a = option_.ReversionSpeed;
			T const sqrtdt = std::sqrt(3.0*dt);
			T const dr = sig * sqrtdt;
			T const mean = std::log(value) + theta_.at(timeIdx - 1);
			return _branching(mean, dr, leafIdx, isMeanReverting);
		}

		// Backward generator
		T operator()(T currValue, T upValue, T midValue, T downValue, T dt) override {
			T const disc = dcf_(currValue, dt);
			return (disc * (0.5*upValue + (1.0 - 0.5)*downValue));
		}

		void setTheta(std::vector<T> const &theta) { theta_ = theta; }

		static std::string const name() {
			return std::string{ "Black-Karasinski model" };
		}

		OptionData<T> const &optionData()const { return option_; }

		static constexpr AssetClass assetClass() { return AssetClass::InterestRate; }


		// Calibration objective function:
		std::function<T(T, T, T, std::vector<T> const&, std::vector<T> const&, std::function<T(T, T)>const &)>
			calibrationObjective(BranchingStyle branchingStyle)const {
			T const sig = option_.Volatility;
			T const a = option_.ReversionSpeed;

			if (branchingStyle == BranchingStyle::Normal) {
				// Calibration normal:
				return [=](T theta, T dt, T marketDiscount, std::vector<T> const &prevRateStates,
					std::vector<T> const &arrowDebreuStates, std::function<T(T, T)>const &dcf)->T {

					T const sqrtdt = std::sqrt(3.0*dt);
					T const dr = sig * sqrtdt;
					T mean{};
					const std::size_t prevNodeSize = prevRateStates.size();
					std::vector<T> rateStates(prevNodeSize + 2);

					for (std::size_t l = 0; l < prevNodeSize; ++l) {
						mean = std::log(prevRateStates.at(l)) + theta;
						rateStates[l] = std::exp(mean - dr);  // down
						rateStates[l + 1] = std::exp(mean);  //mid
						rateStates[l + 2] = std::exp(mean + dr); //high
					}

					T sum{};
					std::size_t const statesSize = arrowDebreuStates.size();
					for (std::size_t i = 0; i < statesSize; ++i) {
						sum += arrowDebreuStates.at(i) * dcf(rateStates.at(i), dt);
					}
					return ((sum - marketDiscount)*(sum - marketDiscount));
				};
			}
			else {

				// Calibration reverting:
				return [=](T theta, T dt, T marketDiscount, std::vector<T> const &prevRateStates,
					std::vector<T> const &arrowDebreuStates, std::function<T(T, T)>const &dcf)->T {

					T const sqrtdt = std::sqrt(3.0*dt);
					T const dr = sig * sqrtdt;
					T mean{};
					const std::size_t prevNodeSize = prevRateStates.size();
					std::vector<T> rateStates(prevNodeSize);

					// Branching upward here:
					mean = std::log(prevRateStates.at(0)) + theta;
					rateStates[0] = std::exp(mean);  // down
					rateStates[1] = std::exp(mean + dr);  //mid
					rateStates[2] = std::exp(mean + 2.0 * dr); //high
													 // normal branching here:
					for (std::size_t l = 1; l < prevNodeSize - 1; ++l) {
						mean = std::log(prevRateStates.at(l)) + theta;
						rateStates[l - 1] = std::exp(mean - dr);  // down
						rateStates[l] = std::exp(mean);  //mid
						rateStates[l + 1] = std::exp(mean + dr); //high
					}
					// Branching downward here:
					mean = std::log(prevRateStates.at(prevNodeSize - 1)) + theta;
					rateStates[prevNodeSize - 1 - 2] = std::exp(mean - 2.0*dr);  // down
					rateStates[prevNodeSize - 1 - 1] = std::exp(mean - dr);  //mid
					rateStates[prevNodeSize - 1] = std::exp(mean); //high

					T sum{};
					std::size_t const statesSize = arrowDebreuStates.size();
					for (std::size_t i = 0; i < statesSize; ++i) {
						sum += arrowDebreuStates.at(i) * dcf(rateStates.at(i), dt);
					}
					return ((sum - marketDiscount)*(sum - marketDiscount));
				};
			}
		}

		// Calibration forward function:
		auto calibrationForwardGenerator()const {
			T const sig = option_.Volatility;
			T const a = option_.ReversionSpeed;

			return [=](T theta, T value, T dt, std::size_t leafIdx, bool isMeanReverting = false)->std::tuple<T, T, T> {
				T const sqrtdt = std::sqrt(3.0*dt);
				T const dr = sig * sqrtdt;
				T const mean = std::log(value) + theta;
				return _branching(mean, dr, leafIdx, isMeanReverting);
			};
		}

	};




}



#endif ///_LATTICE_MODEL_TRINOMIAL