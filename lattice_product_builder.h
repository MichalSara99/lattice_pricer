#pragma once
#if !defined(_LATTICE_PRODUCT_BUILDER)
#define _LATTICE_PRODUCT_BUILDER

#include"lattice_product.h"

namespace lattice_product_builder {

	using lattice_product::Option;
	using lattice_product::BarrierOption;
	using lattice_product::PureDiscountBond;
	using lattice_product::CouponBond;
	using lattice_product::OptionOnPureDiscountBond;


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
		self withMaturity(T value) { option_.setMaturity(value); return *this;}
		self withVolatility(T value) { option_.setVolatility(value); return *this;}

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
		self withMaturity(T value) { option_.setMaturity(value); return *this;}
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
		self withNominal(T value) { bond_.setNominal(value); return *this;}
		self withLastCoupon(T value) { bond_.setLastCoupon(value); return *this;}
		self withCouponPair(TimeAxis time, T value) { bond_.addCoupon(time, value); return *this;}

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
		self withNominal(T value) { option_.setNominal(value); return *this;}
		self withStrike(T value) { option_.setStrike(value); return *this;}

	};

}









#endif ///_LATTICE_PRODUCT_BUILDER