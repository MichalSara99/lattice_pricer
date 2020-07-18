#pragma once
#if !defined(_LATTICE_PRODUCT_BUILDER)
#define _LATTICE_PRODUCT_BUILDER

#include"lattice_product.h"

namespace lattice_product_builder {

	using lattice_product::Option;
	using lattice_product::SpreadOption;
	using lattice_product::BarrierOption;
	using lattice_product::PureDiscountBond;
	using lattice_product::CouponBond;
	using lattice_product::OptionOnPureDiscountBond;
	using lattice_product::OptionOnCouponBond;


	// ===========================================================================
	// =========================== OptionBuilder =================================
	// ===========================================================================

	template<typename T>
	class OptionBuilder {
	private:
		typedef OptionBuilder<T> &self;
		Option<T> option_;

	public:

		~OptionBuilder() {}

		Option<T> build() { return std::move(option_); }

		self withName(std::string const &name) { option_.setName(name); return *this; }
		self withStrike(T value) { option_.setStrike(value); return *this;}
		self withSpot(T value) { option_.setSpot(value); return *this;}
		self withRate(T value) { option_.setRate(value); return *this;}
		self withDividend(T value) { option_.setDividend(value); return *this;}
		self withPeriods(std::size_t periods) { option_.setPeriods(periods); return *this;}
		self withVolatility(T value) { option_.setVolatility(value); return *this;}

	};

	// ===========================================================================
	// ======================= SpreadOptionBuilder ===============================
	// ===========================================================================

	template<typename T>
	class SpreadOptionBuilder {
	private:
		typedef SpreadOptionBuilder<T> &self;
		SpreadOption<T> option_;

	public:

		~SpreadOptionBuilder() {}

		SpreadOption<T> build() { return std::move(option_); }

		self withName(std::string const &name) { option_.setName(name); return *this; }
		self withStrike(T value) { option_.setStrike(value); return *this; }
		self withSpot1(T value) { option_.setSpot1(value); return *this; }
		self withSpot2(T value) { option_.setSpot2(value); return *this; }
		self withRate(T value) { option_.setRate(value); return *this; }
		self withDividend1(T value) { option_.setDividend1(value); return *this; }
		self withDividend2(T value) { option_.setDividend2(value); return *this; }
		self withPeriods(std::size_t periods) { option_.setPeriods(periods); return *this; }
		self withVolatility1(T value) { option_.setVolatility1(value); return *this; }
		self withVolatility2(T value) { option_.setVolatility2(value); return *this; }
		self withCorrelation(T value) { option_.setCorrelation(value); return *this; }

	};

	// ===========================================================================
	// ======================== BarrierOptionBuilder =============================
	// ===========================================================================

	template<typename T>
	class BarrierOptionBuilder {
	private:
		typedef BarrierOptionBuilder<T> &self;
		BarrierOption<T> option_;

	public:
		~BarrierOptionBuilder() {}

		BarrierOption<T> build() { return std::move(option_); }

		self withName(std::string const &name) { option_.setName(name); return *this;}
		self withStrike(T value) { option_.setStrike(value); return *this;}
		self withBarrier(T value) { option_.setBarrier(value); return *this;}
		self withSpot(T value) { option_.setSpot(value); return *this;}
		self withRate(T value) { option_.setRate(value); return *this;}
		self withDividend(T value) { option_.setDividend(value); return *this;}
		self withPeriods(std::size_t periods) { option_.setPeriods(periods); return *this;}
		self withVolatility(T value) { option_.setVolatility(value); return *this;}

	};

	// ===========================================================================
	// ============================ PureDiscountBondBuilder ======================
	// ===========================================================================

	template<typename T>
	class PureDiscountBondBuilder {
	private:
		typedef PureDiscountBondBuilder<T> &self;
		PureDiscountBond<T> bond_;

	public:
		~PureDiscountBondBuilder() {}

		PureDiscountBond<T> build() { return std::move(bond_); }

		self withName(std::string const &name) { bond_.setName(name); return *this;}
		self withNominal(T value) { bond_.setNominal(value); return *this;}
		self withPeriods(std::size_t periods) { bond_.setPeriods(periods); return *this; }


	};

	// ===========================================================================
	// =============================== CouponBondBuilder =========================
	// ===========================================================================

	template<typename T, typename TimeAxis>
	class CouponBondBuilder {
	private:
		typedef CouponBondBuilder<T, TimeAxis> &self;
		CouponBond<T, TimeAxis> bond_;

	public:
		~CouponBondBuilder() {}

		CouponBond<T, TimeAxis> build() { return std::move(bond_); }

		self withName(std::string const &name) { bond_.setName(name); return *this;}
		self withPeriods(std::size_t periods) { bond_.setPeriods(periods); return *this; }
		self withNominal(T value) { bond_.setNominal(value); return *this;}
		self withLastCoupon(T value) { bond_.setLastCoupon(value); return *this;}
		self withCouponPair(TimeAxis time, T value) { bond_.addCoupon(time, value); return *this;}
		self withCouponData(std::map<TimeAxis, T> const &data) { bond_.setCouponData(data); return *this; }
	};


	// ===========================================================================
	// ========================= OptionOnPureDiscountBondBuilder =================
	// ===========================================================================

	template<typename T>
	class OptionOnPureDiscountBondBuilder {
	private:
		typedef OptionOnPureDiscountBondBuilder<T> &self;
		OptionOnPureDiscountBond<T> option_;

	public:
		~OptionOnPureDiscountBondBuilder() {}

		OptionOnPureDiscountBond<T> build() { return std::move(option_); }

		self withName(std::string const &name) { option_.setName(name); return *this;}
		self withStrike(T value) { option_.setStrike(value); return *this;}
		self withPeriods(std::size_t periods) { option_.setPeriods(periods); return *this; }

	};

	// ===========================================================================
	// ========================= OptionOnCouponBondBuilder =======================
	// ===========================================================================

	template<typename T>
	class OptionOnCouponBondBuilder {
	private:
		typedef OptionOnCouponBondBuilder<T> &self;
		OptionOnCouponBond<T> option_;

	public:
		~OptionOnCouponBondBuilder() {}

		OptionOnCouponBond<T> build() { return std::move(option_); }

		self withName(std::string const &name) { option_.setName(name); return *this; }
		self withStrike(T value) { option_.setStrike(value); return *this; }
		self withPeriods(std::size_t periods) { option_.setPeriods(periods); return *this; }

	};

}









#endif ///_LATTICE_PRODUCT_BUILDER