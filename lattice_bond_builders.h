#pragma once
#if !defined(_LATTICE_BOND_BUILDERS)
#define _LATTICE_BOND_BUILDERS

#include"lattice_types.h"
#include"lattice_utility.h"
#include"lattice_macros.h"

namespace lattice_bond_builders {

	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_types::LaunchPolicy;
	using lattice_utility::DiscountingFactor;
	using lattice_utility::DeltaTimeHolder;
	using lattice_utility::DiscountingStyle;

	// ==============================================================================
	// =============================== BondBuilder  =================================
	// ==============================================================================

	template<LatticeType Type,
		typename DeltaTime,
		typename Node>
		class BondBuilder {};

	template<typename DeltaTime,
		typename Node>
	class BondBuilder<LatticeType::Binomial, DeltaTime, Node>{
	private:
		template<typename LatticeObject, typename Generator>
		static void _buildPureDiscountBondTree(std::size_t lastIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice,
			Generator const &generator, Node nominal, DeltaTime const &deltaTime, lattice_types::DiscountingStyle style);
		
		template<typename LatticeObject, typename Generator>
		static void _buildCouponBondTree(std::size_t lastIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
			Node nominal, Node couponRate, std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
			DeltaTime const &deltaTime, lattice_types::DiscountingStyle style);

		template<typename LatticeObject,typename Generator>
		static void _buildBondTree(LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
			Node nominal, Node couponRate, Node lastCoupon,std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
			DeltaTime const &deltaTime, DiscountingStyle style);

	public:
		template<typename LatticeObject,typename Generator>
		void operator()(LatticeObject &bondLattice, LatticeObject const &calibartedRateLattice, Generator const &generator,
			Node nominal, Node couponRate,Node lastCoupon, std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
			DeltaTime const &deltaTime, DiscountingStyle style = DiscountingStyle::Continuous) {
			LASSERT(bondLattice.timeDimension() > 0, "Passed bondLattice must have at least 1 period.");
			LASSERT(bondLattice.timeDimension() <= calibartedRateLattice.timeDimension(),
				"bondLattice's dimension must be equal or less then calibratedRateLatice's dimension.");
			LASSERT(bondLattice.type() == generator.latticeType(),
				"Mismatch between lattice types");
			LASSERT(Generator::assetClass() == AssetClass::InterestRate,
				"The passed model must be interest-rate model");
			_buildBondTree(bondLattice, calibartedRateLattice, generator,
				nominal, couponRate,lastCoupon,couponDates, deltaTime, style);
		}

	};


	template<typename DeltaTime,
			typename Node>
	class BondBuilder<LatticeType::Trinomial, DeltaTime, Node> {
	private:
		template<typename LatticeObject, typename Generator>
		static void _buildNormalPureDiscountBondTree(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice,
			Generator const &generator, Node nominal, DeltaTime const &deltaTime, lattice_types::DiscountingStyle style);

		template<typename LatticeObject, typename Generator>
		static void _buildNormalCouponBondTree(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
			Node nominal, Node couponRate, std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
			DeltaTime const &deltaTime, DiscountingStyle style);

		template<typename LatticeObject,typename Generator>
		static void _buildNormal(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
			Node nominal, Node couponRate, std::set<typename LatticeObject::TimeAxis_type> const &couponDates, 
			DeltaTime const &deltaTime, DiscountingStyle style);

		template<typename LatticeObject, typename Generator>
		static void _buildRevertingPureDiscountBondTree(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice,
			Generator const &generator, Node nominal, DeltaTime const &deltaTime, DiscountingStyle style);

		template<typename LatticeObject, typename Generator>
		static void _buildRevertingCouponBondTree(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
			Node nominal, Node couponRate, std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
			DeltaTime const &deltaTime, DiscountingStyle style);


		template<typename LatticeObject,typename Generator>
		static void _buildReverting(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
			Node nominal, Node couponRate,std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
			DeltaTime const &deltaTime, DiscountingStyle style);


		template<typename LatticeObject, typename Generator>
		static void _buildBondTree(LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
			Node nominal, Node couponRate,Node lastCoupon, std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
			DeltaTime const &deltaTime, DiscountingStyle style);

	public:
		template<typename LatticeObject, typename Generator>
		void operator()(LatticeObject &bondLattice, LatticeObject const &calibartedRateLattice, Generator const &generator,
			Node nominal, Node couponRate,Node lastCoupon, std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
			DeltaTime const &deltaTime, DiscountingStyle style = DiscountingStyle::Continuous) {
			LASSERT(bondLattice.timeDimension() > 0, "Passed bondLattice must have at least 1 period.");
			LASSERT(bondLattice.timeDimension() <= calibartedRateLattice.timeDimension(),
				"bondLattice's dimension must be equal or less then calibratedRateLatice's dimension.");
			LASSERT(bondLattice.type() == generator.latticeType(),
				"Mismatch between lattice types");
			LASSERT(Generator::assetClass() == AssetClass::InterestRate,
				"The passed model must be interest-rate model");
			_buildBondTree(bondLattice, calibartedRateLattice, generator, nominal, couponRate,lastCoupon,couponDates, deltaTime, style);
		}

	};



	// ==============================================================================
	// ======================== OptionOnBondBuilder  ================================
	// ==============================================================================


	template<LatticeType Type,
		typename DeltaTime,
		typename Node>
	class OptionOnBondBuilder {};


	template<typename DeltaTime,
			typename Node>
	class OptionOnBondBuilder<LatticeType::Binomial, DeltaTime, Node> {

	};


	template<typename DeltaTime,
		typename Node>
	class OptionOnBondBuilder<LatticeType::Trinomial, DeltaTime, Node> {

	};



}


template<typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, DeltaTime, Node>::
_buildPureDiscountBondTree(std::size_t lastIdx,LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice,
	Generator const &generator,Node nominal, DeltaTime const &deltaTime, lattice_types::DiscountingStyle style) {

	typedef DiscountingFactor<Node> DCF;
	typedef DeltaTimeHolder<DeltaTime> DT;
	auto dcf = DCF::function(style);

	Node dt{};
	Node prob = generator.nodeRiskNeutralProb();

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = bondLattice.nodesAtIdx(n).size();
		dt = DT::deltaTime(n, deltaTime);

		for (auto i = 0; i < nodesSize; ++i) {
			bondLattice(n, i) =
				prob * (bondLattice(n + 1, i) + bondLattice(n + 1, i + 1))*
				dcf(calibratedRateLattice(n, i), dt);
		}
	}

	dt = DT::deltaTime(0, deltaTime);
	bondLattice(0, 0) = prob * (bondLattice(1, 0) + bondLattice(1, 1))*
		dcf(calibratedRateLattice(0, 0), dt);
}


template<typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial, DeltaTime, Node>::
_buildCouponBondTree(std::size_t lastIdx,LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
	Node nominal, Node couponRate, std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
	DeltaTime const &deltaTime, lattice_types::DiscountingStyle style) {

	typedef DiscountingFactor<Node> DCF;
	typedef DeltaTimeHolder<DeltaTime> DT;
	using TimeAxis = typename LatticeObject::TimeAxis_type;
	auto dcf = DCF::function(style);

	Node dt{};
	Node coupon{};
	TimeAxis date{};
	Node prob = generator.nodeRiskNeutralProb();


	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = bondLattice.nodesAtIdx(n).size();
		dt = DT::deltaTime(n, deltaTime);
		coupon = 0.0;
		date = bondLattice.timeAt(n);
		if (couponDates.find(date) != couponDates.end()) {
			coupon = nominal * couponRate * dt;
		}
		for (auto i = 0; i < nodesSize; ++i) {
			bondLattice(n, i) =
				prob * (bondLattice(n + 1, i) + bondLattice(n + 1, i + 1) + 2.0 * coupon)*
				dcf(calibratedRateLattice(n, i), dt);
		}
	}

	dt = DT::deltaTime(0, deltaTime);
	coupon = 0.0;
	date = bondLattice.timeAt(0);
	if (couponDates.find(date) != couponDates.end()) {
		coupon = nominal * couponRate * dt;
	}
	bondLattice(0, 0) = prob * (bondLattice(1, 0) + bondLattice(1, 1) + 2.0*coupon)*
		dcf(calibratedRateLattice(0, 0), dt);
}


template<typename DeltaTime,
	typename Node>
template<typename LatticeObject, typename Generator>
void lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial,DeltaTime,Node>::
_buildBondTree(LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
	Node nominal, Node couponRate, Node lastCoupon, std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
	DeltaTime const &deltaTime, lattice_types::DiscountingStyle style) {

	// declare for further use in last payoff:
	const std::size_t lastIdx = bondLattice.timeDimension() - 1;
	const std::size_t lastNodesSize = bondLattice.nodesAtIdx(lastIdx).size();

	for (auto i = 0; i < lastNodesSize; ++i) {
		bondLattice(lastIdx, i) = nominal + lastCoupon;
	}

	// If couponDates is an empty set -> ignore possible non-zero couponRate variable
	// and  compute pure discount bond tree
	if (couponDates.empty()) {
		_buildPureDiscountBondTree(lastIdx, bondLattice, calibratedRateLattice,
			generator, nominal, deltaTime, style);
	}
	else {
		_buildCouponBondTree(lastIdx, bondLattice, calibratedRateLattice,
			generator, nominal, couponRate, couponDates, deltaTime, style);
	}
}


template<typename DeltaTime,
		typename Node>
template<typename LatticeObject,typename Generator>
void lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial,DeltaTime,Node>::
_buildNormalPureDiscountBondTree(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice,
	Generator const &generator, Node nominal, DeltaTime const &deltaTime, lattice_types::DiscountingStyle style) {
	
	typedef DiscountingFactor<Node> DCF;
	typedef DeltaTimeHolder<DeltaTime> DT;
	auto dcf = DCF::function(style);

	std::size_t const revertBranchesSize = timeIdx;

	std::tuple<Node, Node, Node> prob;
	Node dt{};

	std::size_t nodesSize{ 0 };
	for (auto n = timeIdx - 1; n > 0; --n) {
		nodesSize = bondLattice.nodesAtIdx(n).size();
		dt = DT::deltaTime(n, deltaTime);
		for (auto l = 0; l < nodesSize; ++l) {
			prob = generator.nodeRiskNeutralProb(revertBranchesSize, nodesSize, l, dt);
			bondLattice(n, l) = (std::get<0>(prob)*bondLattice(n + 1, l) +
				std::get<1>(prob)*bondLattice(n + 1, l + 1) +
				std::get<2>(prob)*bondLattice(n + 1, l + 2))*
				dcf(calibratedRateLattice(n, l), dt);

		}
	}

	dt = DT::deltaTime(0, deltaTime);
	prob = generator.nodeRiskNeutralProb(revertBranchesSize, 1, 0, dt);
	bondLattice(0, 0) = (std::get<0>(prob)*bondLattice(1, 0) +
		std::get<1>(prob)*bondLattice(1, 1) +
		std::get<2>(prob)*bondLattice(1, 2))*
		dcf(calibratedRateLattice(0, 0), dt);
}


template<typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial,DeltaTime,Node>::
_buildNormalCouponBondTree(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
	Node nominal, Node couponRate, std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
	DeltaTime const &deltaTime, DiscountingStyle style) {

	typedef DiscountingFactor<Node> DCF;
	typedef DeltaTimeHolder<DeltaTime> DT;
	using TimeAxis = typename LatticeObject::TimeAxis_type;
	auto dcf = DCF::function(style);

	std::size_t const revertBranchesSize = timeIdx;

	Node dt{};
	Node coupon{};
	TimeAxis date{};
	std::tuple<Node, Node, Node> prob;


	std::size_t nodesSize{ 0 };
	for (auto n = timeIdx - 1; n > 0; --n) {
		nodesSize = bondLattice.nodesAtIdx(n).size();
		dt = DT::deltaTime(n, deltaTime);
		coupon = 0.0;
		date = bondLattice.timeAt(n);
		if (couponDates.find(date) != couponDates.end()) {
			coupon = nominal * couponRate * dt;
		}
		for (auto l = 0; l < nodesSize; ++l) {
			prob = generator.nodeRiskNeutralProb(revertBranchesSize, nodesSize, l, dt);
			bondLattice(n, l) = (std::get<0>(prob)*(bondLattice(n + 1, l) + coupon) +
				std::get<1>(prob)*(bondLattice(n + 1, l + 1) + coupon) +
				std::get<2>(prob)*(bondLattice(n + 1, l + 2) + coupon))*
				dcf(calibratedRateLattice(n, l), dt);

		}
	}

	dt = DT::deltaTime(0, deltaTime);
	coupon = 0.0;
	date = bondLattice.timeAt(0);
	if (couponDates.find(date) != couponDates.end()) {
		coupon = nominal * couponRate * dt;
	}
	prob = generator.nodeRiskNeutralProb(revertBranchesSize, 1, 0, dt);
	bondLattice(0, 0) = (std::get<0>(prob)*(bondLattice(1, 0) + coupon) +
		std::get<1>(prob)*(bondLattice(1, 1) + coupon) +
		std::get<2>(prob)*(bondLattice(1, 2) + coupon))*
		dcf(calibratedRateLattice(0, 0), dt);

}


template<typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, DeltaTime, Node>::
_buildNormal(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
	Node nominal, Node couponRate, std::set<typename LatticeObject::TimeAxis_type> const &couponDates, 
	DeltaTime const &deltaTime, lattice_types::DiscountingStyle style) {


	// If couponDates is an empty set -> ignore possible non-zero couponRate variable
	// and  compute pure discount bond tree
	if (couponDates.empty()) {
		_buildNormalPureDiscountBondTree(timeIdx, bondLattice, calibratedRateLattice,
			generator, nominal, deltaTime, style);
	}
	else {
		_buildNormalCouponBondTree(timeIdx, bondLattice, calibratedRateLattice,
			generator, nominal, couponRate, couponDates, deltaTime, style);
	}
}


template<typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void  lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, DeltaTime, Node>::
_buildRevertingPureDiscountBondTree(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice,
	Generator const &generator, Node nominal, DeltaTime const &deltaTime, DiscountingStyle style) {

	typedef DiscountingFactor<Node> DCF;
	typedef DeltaTimeHolder<DeltaTime> DT;
	auto dcf = DCF::function(style);

	std::size_t const lastIdx = bondLattice.timeDimension() - 1;
	std::size_t const revertBranchesSize = timeIdx;

	std::tuple<Node, Node, Node> prob;
	Node dt{};

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n >= timeIdx; --n) {
		nodesSize = bondLattice.nodesAtIdx(n).size();
		dt = DT::deltaTime(n, deltaTime);
		prob = generator.nodeRiskNeutralProb(revertBranchesSize, nodesSize, 0, dt);
		bondLattice(n, 0) = (std::get<0>(prob)*bondLattice(n + 1, 0) +
			std::get<1>(prob)*bondLattice(n + 1, 1)+
			std::get<2>(prob)*bondLattice(n + 1, 2))*
			dcf(calibratedRateLattice(n, 0), dt);

		for (auto l = 1; l < nodesSize - 1; ++l) {
			prob = generator.nodeRiskNeutralProb(revertBranchesSize, nodesSize, l, dt);
			bondLattice(n, l) = (std::get<0>(prob)*bondLattice(n + 1, l - 1) +
				std::get<1>(prob)*bondLattice(n + 1, l) +
				std::get<2>(prob)*bondLattice(n + 1, l + 1) )*
				dcf(calibratedRateLattice(n, l), dt);

		}

		prob = generator.nodeRiskNeutralProb(revertBranchesSize, nodesSize, nodesSize - 1, dt);
		bondLattice(n, nodesSize - 1) = (std::get<0>(prob)*bondLattice(n + 1, nodesSize - 3) +
			std::get<1>(prob)*bondLattice(n + 1, nodesSize - 2) +
			std::get<2>(prob)*bondLattice(n + 1, nodesSize - 1))*
			dcf(calibratedRateLattice(n, nodesSize - 1), dt);
	}
}

template<typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, DeltaTime, Node>::
_buildRevertingCouponBondTree(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
	Node nominal, Node couponRate, std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
	DeltaTime const &deltaTime, DiscountingStyle style) {

	typedef DiscountingFactor<Node> DCF;
	using TimeAxis = typename LatticeObject::TimeAxis_type;
	typedef DeltaTimeHolder<DeltaTime> DT;
	auto dcf = DCF::function(style);

	std::size_t const lastIdx = bondLattice.timeDimension() - 1;
	std::size_t const revertBranchesSize = timeIdx;

	Node dt{};
	Node coupon{};
	TimeAxis date{};
	std::tuple<Node, Node, Node> prob;


	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n >= timeIdx; --n) {
		nodesSize = bondLattice.nodesAtIdx(n).size();
		dt = DT::deltaTime(n, deltaTime);
		coupon = 0.0;
		date = bondLattice.timeAt(n);
		if (couponDates.find(date) != couponDates.end()) {
			coupon = nominal * couponRate * dt;
		}
		prob = generator.nodeRiskNeutralProb(revertBranchesSize, nodesSize, 0, dt);
		bondLattice(n, 0) = (std::get<0>(prob)*(bondLattice(n + 1, 0) + coupon) +
			std::get<1>(prob)*(bondLattice(n + 1, 1) + coupon) +
			std::get<2>(prob)*(bondLattice(n + 1, 2) + coupon))*
			dcf(calibratedRateLattice(n, 0), dt);

		for (auto l = 1; l < nodesSize - 1; ++l) {
			prob = generator.nodeRiskNeutralProb(revertBranchesSize, nodesSize, l, dt);
			bondLattice(n, l) = (std::get<0>(prob)*(bondLattice(n + 1, l - 1) + coupon) +
				std::get<1>(prob)*(bondLattice(n + 1, l) + coupon) +
				std::get<2>(prob)*(bondLattice(n + 1, l + 1) + coupon))*
				dcf(calibratedRateLattice(n, l), dt);

		}

		prob = generator.nodeRiskNeutralProb(revertBranchesSize, nodesSize, nodesSize - 1, dt);
		bondLattice(n, nodesSize - 1) = (std::get<0>(prob)*(bondLattice(n + 1, nodesSize - 3) + coupon) +
			std::get<1>(prob)*(bondLattice(n + 1, nodesSize - 2) + coupon) +
			std::get<2>(prob)*(bondLattice(n + 1, nodesSize - 1) + coupon))*
			dcf(calibratedRateLattice(n, nodesSize - 1), dt);
	}
}

template<typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, DeltaTime, Node>::
_buildReverting(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
	Node nominal, Node couponRate, std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
	DeltaTime const &deltaTime, lattice_types::DiscountingStyle style) {

	// If couponDates is an empty set -> ignore possible non-zero couponRate variable
	// and  compute pure discount bond tree
	if (couponDates.empty()) {
		_buildRevertingPureDiscountBondTree(timeIdx, bondLattice, calibratedRateLattice,
			generator, nominal, deltaTime, style);
	}
	else {
		_buildRevertingCouponBondTree(timeIdx, bondLattice, calibratedRateLattice,
			generator, nominal, couponRate, couponDates, deltaTime, style);
	}
}

template<typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, DeltaTime, Node>::
_buildBondTree(LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
	Node nominal, Node couponRate,Node lastCoupon, std::set<typename LatticeObject::TimeAxis_type> const &couponDates,
	DeltaTime const &deltaTime, lattice_types::DiscountingStyle style) {

	const std::size_t firstRevertIdx = bondLattice.firstRevertingIdx();
	const std::size_t treeSize = bondLattice.timeDimension();

	const std::size_t lastIdx = treeSize - 1;
	const std::size_t lastNodesSize = bondLattice.nodesAtIdx(lastIdx).size();

	for (auto i = 0; i < lastNodesSize; ++i) {
		bondLattice(lastIdx, i) = nominal + lastCoupon;
	}

	if (firstRevertIdx == 0) {
		// This trinomial tree does not have reverting property:
		_buildNormal(lastIdx, bondLattice, calibratedRateLattice, generator, nominal, couponRate,couponDates, deltaTime, style);
	}
	else {
		// This trinomial tree does have reverting property:
		_buildReverting(firstRevertIdx - 1, bondLattice, calibratedRateLattice, generator, nominal, couponRate, couponDates, deltaTime, style);
		_buildNormal(firstRevertIdx - 1, bondLattice, calibratedRateLattice, generator, nominal, couponRate, couponDates, deltaTime, style);

	}
}







#endif ///_LATTICE_BOND_BUILDERS