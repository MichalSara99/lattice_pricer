#pragma once
#if !defined(_LATTICE_PRODUCT)
#define _LATTICE_PRODUCT

#include<set>
#include<limits>
#include<memory>
#include<map>

#include"lattice_types.h"
#include"lattice_model_params.h"

namespace lattice_product {

	using lattice_model_params::ModelParams;
	using lattice_types::AssetClass;

	template<typename T>
	class Product {
	protected:
		std::string name_;

	public:
		virtual ~Product(){}


		inline void setName(std::string const &name) { name_ = name; }
		inline std::string const &name()const { return name_; }
	};

	// ===========================================================================
	// ================================= Option ==================================
	// ===========================================================================
	

	template<typename T>
	class Option:public Product<T>{
	protected:
		T strike_;
		T spot_;
		T rate_;
		T div_;
		std::size_t periods_;
		T vol_;

	public:
		explicit Option() {}

		virtual ~Option(){}

		ModelParams<1,AssetClass::Equity,T> modelParams()const{
			return { rate(),volatility(),dividend(),spot(),strike(),T{} };
		}

		inline void setStrike(T strike) { strike_ = strike; }
		inline T strike()const { return strike_; }

		inline void setSpot(T spot) { spot_ = spot; }
		inline T spot()const { return spot_; }

		inline void setRate(T rate) { rate_ = rate; }
		inline T rate()const { return rate_; }

		inline void setDividend(T div) { div_ = div; }
		inline T dividend()const { return div_; }

		inline void setPeriods(std::size_t periods) { periods_ = periods; }
		inline std::size_t periods()const { return periods_; }

		inline void setVolatility(T vol) { vol_ = vol; }
		inline T volatility()const { return vol_; }

	};

	// ===========================================================================
	// =========================== SpreadOption ==================================
	// ===========================================================================

	template<typename T>
	class SpreadOption :public Product<T> {
	private:
		T strike_;
		T spot1_;
		T spot2_;
		T rate_;
		T div1_;
		T div2_;
		T vol1_;
		T vol2_;
		T corr_;
		std::size_t periods_;

	public:
		explicit SpreadOption() {}

		virtual ~SpreadOption() {}

		ModelParams<2, AssetClass::Equity, T> modelParams()const {
			return { rate(),dividend1(),dividend2(),
				volatility1(),volatility2(),
				spot1(),spot2(),strike(),correlation() };
		}

		inline void setPeriods(std::size_t periods) { periods_ = periods; }
		inline std::size_t periods()const { return periods_; }

		inline void setStrike(T strike) { strike_ = strike; }
		inline T strike()const { return strike_; }

		inline void setSpot1(T spot) { spot1_ = spot; }
		inline T spot1()const { return spot1_; }

		inline void setSpot2(T spot) { spot2_ = spot; }
		inline T spot2()const { return spot2_; }

		inline void setRate(T rate) { rate_ = rate; }
		inline T rate()const { return rate_; }

		inline void setDividend1(T div) { div1_ = div; }
		inline T dividend1()const { return div1_; }

		inline void setDividend2(T div) { div2_ = div; }
		inline T dividend2()const { return div2_; }

		inline void setVolatility1(T vol) { vol1_ = vol; }
		inline T volatility1()const { return vol1_; }

		inline void setVolatility2(T vol) { vol2_ = vol; }
		inline T volatility2()const { return vol2_; }

		inline void setCorrelation(T value) { corr_ = value; }
		inline T correlation()const { return corr_; }

	};


	// ===========================================================================
	// ============================ BarrierOption ================================
	// ===========================================================================


	template<typename T>
	class BarrierOption :public Option<T> {
	protected:
		T barrier_;

	public:
		explicit BarrierOption() {}

		virtual ~BarrierOption() {}

		ModelParams<1, AssetClass::Equity, T> modelParams()const {
			return { this->rate(),this->volatility(),this->dividend(),this->spot(),this->strike(),barrier() };
		}

		inline void setBarrier(T barrier) { barrier_ = barrier; }
		inline T barrier()const { return barrier_; }

	};

	// ===========================================================================
	// =============================== PureDiscountBond ==========================
	// ===========================================================================


	template<typename T>
	class PureDiscountBond :public Product<T> {
	protected:
		T nominal_;
		std::size_t periods_;
	public:

		explicit PureDiscountBond() {}
		virtual ~PureDiscountBond() {}

		inline void setNominal(T value) { nominal_ = value; }
		inline T nominal()const { return nominal_; }

		inline void setPeriods(std::size_t periods) { periods_ = periods; }
		inline std::size_t periods()const { return periods_; }

	};


	// ===========================================================================
	// ================================= CouponBond ==============================
	// ===========================================================================

	template<typename T,typename TimeAxis>
	class CouponBond :public PureDiscountBond<T> {
	private:
		T lastCoupon_;
		std::map<TimeAxis, T> coupons_;

	public:

		explicit CouponBond() {}
		virtual ~CouponBond() {}

		inline void setLastCoupon(T value) { lastCoupon_ = value; }
		inline T lastCoupon()const { return lastCoupon_; }

		inline void addCoupon(TimeAxis time, T value) { coupons_.emplace(time, value); }
		inline T coupon(TimeAxis time) { 
			auto it = coupons_.find(time);
			return ((it != coupons_.end()) ? (*it)->first : std::numeric_limits<T>::signaling_NaN());
		}

		inline std::map<TimeAxis, T> const &coupons()const { return coupons_; }

	};

	// ===========================================================================
	// =========================== OptionOnPureDiscountBond  =====================
	// ===========================================================================


	template<typename T>
	class OptionOnPureDiscountBond :public PureDiscountBond<T> {
	private:
		T strike_;

	public:
		explicit OptionOnPureDiscountBond() {}
		virtual ~OptionOnPureDiscountBond() {}

		inline void setStrike(T value) { strike_ = value; }
		inline T strike()const { return strike_; }

	};


	// ===========================================================================
	// =========================== OptionOnCouponBond ============================
	// ===========================================================================
	
	template<typename T, typename TimeAxis>
	class OptionOnCouponBond :public CouponBond<T,TimeAxis> {
	private:
		T lastCoupon_;
		std::set<std::pair<TimeAxis, T>> coupons_;

	public:

		explicit OptionOnCouponBond() {}
		virtual ~OptionOnCouponBond() {}

		inline void setStrike(T value) { strike_ = value; }
		inline T strike()const { return strike_; }

		inline void setLastCoupon(T value) { lastCoupon_ = value; }
		inline T lastCoupon()const { return lastCoupon_; }

		inline void addCoupon(TimeAxis time, T value) { coupons_.emplace(std::make_pair(time, value)); }
		inline T coupon(TimeAxis time) {
			auto it = coupons_.find(time);
			return ((it != coupons_.end()) ? (*it)->first : std::numeric_limits<T>::signaling_NaN());
		}

		inline std::set<std::pair<TimeAxis, T>> const &coupons()const { return coupons_; }

	};











}


#endif ///_LATTICE_PRODUCT