#pragma once
#if !defined(_LATTICE_PRODUCT)
#define _LATTICE_PRODUCT

#include<set>
#include<limits>
#include<memory>

#include"lattice_types.h"

namespace lattice_product {


	template<typename T>
	class Product {
	protected:
		std::string name_;

	public:
		virtual ~Product(){}


		inline void setName(std::string const &name) { return name_ = name; }
		inline std::string const &name()const { return name_; }
	};

	// ===========================================================================
	// ====================== Option | OptionBuilder =============================
	// ===========================================================================
	
	// forward class declaration
	template<typename T>
	class OptionBuilder;

	template<typename T>
	class Option:public Product<T>{
	protected:
		T strike_;
		T spot_;
		T rate_;
		T div_;
		T maturity_;
		T vol_;

		explicit Option() {}

	public:
		static std::unique_ptr<OptionBuilder<T>> build() {
			return std::make_unique<OptionBuilder<T>>();
		}

		inline void setStrike(T strike) { strike_ = strike; }
		inline T strike()const { return strike_; }

		inline void setSpot(T spot) { spot_ = spot; }
		inline T spot()const { return spot_; }

		inline void setRate(T rate) { rate_ = rate; }
		inline T rate()const { return rate_; }

		inline void setDividend(T div) { div_ = div; }
		inline T dividend()const { return div_; }

		inline void setMaturity(T maturity) { mat_ = maturity; }
		inline T maturity()const { return mat_; }

		inline void setVolatility(T vol) { vol_ = vol; }
		inline T volatility()const { return vol_; }
	};

	template<typename T>
	class OptionBuilder {
	private:
		typedef OptionBuilder<T> &self;
		Option<T> option_;

	public:
		~OptionBuilder(){}

		operator Option<T>() {return option_;}

		self withName(std::string const &name) { option_.setName(name); }
		self withStrike(T value) { option_.setStrike(value); }
		self withSpot(T value) { option_.setSpot(value); }
		self withRate(T value) { option_.setRate(value); }
		self withDividend(T value) { option_.setDividend(value); }
		self withMaturity(T value) { option_.setMaturity(value); }
		self withVolatility(T value) { option_.setVolatility(value); }

	};


	// ===========================================================================
	// ============== BarrierOption | BarrierOptionBuilder =======================
	// ===========================================================================

	template<typename T>
	class BarrierOptionBuilder;

	template<typename T>
	class BarrierOption :public Option<T> {
	protected:
		T barrier_;

		explicit BarrierOption() {}

	public:
		static std::unique_ptr<BarrierOptionBuilder<T>> build() {
			return std::make_unique<BarrierOptionBuilder<T>>();
		}

		inline void setBarrier(T barrier) { barrier_ = barrier; }
		inline T barrier()const { return barrier_; }

	};

	template<typename T>
	class BarrierOptionBuilder{
	private:
		typedef BarrierOptionBuilder<T> &self;
		BarrierOption<T> option_;

	public:
		~BarrierOptionBuilder() {}

		operator BarrierOption<T>() { return option_; }

		self withName(std::string const &name) { option_.setName(name); }
		self withStrike(T value) { option_.setStrike(value); }
		self withBarrier(T value) { option_.setBarrier(value); }
		self withSpot(T value) { option_.setSpot(value); }
		self withRate(T value) { option_.setRate(value); }
		self withDividend(T value) { option_.setDividend(value); }
		self withMaturity(T value) { option_.setMaturity(value); }
		self withVolatility(T value) { option_.setVolatility(value); }

	};

	// ===========================================================================
	// ======================= PDBond | PDBondBuilder ============================
	// ===========================================================================

	template<typename T>
	class PDBondBuilder;

	template<typename T>
	class PDBond :public Product<T> {
	protected:
		T nominal_;
		explicit PDBond() {}

	public:
		static std::unique_ptr<PDBond<T>> build() {
			return std::make_unique<PDBond<T>>();
		}

		inline void setNominal(T value) { nominal_ = value; }
		inline T nominal()const { return nominal_; }

	};

	template<typename T>
	class PDBondBuilder {
	private:
		typedef PDBondBuilder<T> &self;
		PDBond<T> bond_;

	public:
		~PDBondBuilder(){}

		operator PDBond<T>() { return bond_; }

		self withName(std::string const &name) { bond_.setName(name); }
		self withNominal(T value) { bond_.setNominal(value); }

	};


	// ===========================================================================
	// ======================= CouponBond | CouponBondBuilder ====================
	// ===========================================================================

	template<typename T,typename TimeAxis>
	class CouponBondBuilder;

	template<typename T,typename TimeAxis>
	class CouponBond :public PDBond<T> {
	private:
		T lastCoupon_;
		std::set<std::pair<TimeAxis, T>> coupons_;

		explicit CouponBond() {}

	public:
		static std::unique_ptr<CouponBondBuilder<T,TimeAxis>> build() {
			return std::make_unique<CouponBondBuilder<T, TimeAxis>>();
		}

		inline void setLastCoupon(T value) { lastCoupon_ = value; }
		inline T lastCoupon()const { return lastCoupon_; }

		inline void addCoupon(TimeAxis time, T value) { coupons_.emplace(std::make_pair(time, value)); }
		inline T coupon(TimeAxis time) { 
			auto it = coupons_.find(time);
			return ((it != coupons_.end()) ? (*it)->first : std::numeric_limits<T>::signaling_NaN());
		}

		inline std::set<std::pair<TimeAxis, T>> const &coupons()const { return coupons_; }

	};

	template<typename T,typename TimeAxis>
	class CouponBondBuilder {
	private:
		typedef CouponBondBuilder<T, TimeAxis> &self;
		CouponBond<T, TimeAxis> bond_;

	public:
		operator CouponBond<T, TimeAxis>() { return bond_; }

		self withName(std::string const &name) { bond_.setName(name); }
		self withNominal(T value) { bond_.setNominal(value); }
		self withLastCoupon(T value) { bond_.setLastCoupon(value); }
		self withCouponPair(TimeAxis time, T value) { bond_.addCoupon(time, value); }

	};

	// ===========================================================================
	// =================== OptionOnPDBond | OptionOnPDBondBuilder ================
	// ===========================================================================

	template<typename T>
	class OptionOnPDBondBuilder;

	template<typename T>
	class OptionOnPDBond :public PDBond<T> {
	private:
		T strike_;
		explicit OptionOnPDBond() {}

	public:
		static std::unique_ptr<OptionOnPDBondBuilder<T>> build() {
			return std::make_unique<OptionOnPDBondBuilder<T>>();
		}

		inline void setStrike(T value) { strike_ = value; }
		inline T strike()const { return strike_; }

	};

	template<typename T>
	class OptionOnPDBondBuilder {
	private:
		typedef OptionOnPDBondBuilder<T> &self;
		OptionOnPDBond<T> option_;

	public:
		operator OptionOnPDBond<T>() { return option_; }

		self withName(std::string const &name) { option_.setName(name); }
		self withNominal(T value) { option_.setNominal(value); }
		self withStrike(T value) { option_.setStrike(value); }

	};

	// ===========================================================================
	// =========================== OptionOnCouponBond ============================
	// ===========================================================================

	template<typename T, typename TimeAxis>
	class OptionOnCouponBond :public CouponBond<T,TimeAxis> {
	private:
		T strike_;
	public:
		explicit OptionOnCouponBond() {}

		inline void setStrike(T value) { strike_ = value; }
		inline T strike()const { return strike_; }

	};











}


#endif ///_LATTICE_PRODUCT