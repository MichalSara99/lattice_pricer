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
		template<typename LatticeObject,typename Generator>
		static void _buildBondTree(LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
			Node nominal, Node couponRate, DeltaTime const &deltaTime, DiscountingStyle style);

	public:
		template<typename LatticeObject,typename Generator>
		void operator()(LatticeObject &bondLattice, LatticeObject const &calibartedRateLattice, Generator const &generator,
			Node nominal, Node couponRate,DeltaTime const &deltaTime, DiscountingStyle style = DiscountingStyle::Continuous) {
			LASSERT(bondLattice.timeDimension() == calibartedRateLattice.timeDimension(),
				"Lattices must have the same dimension.");
			LASSERT(bondLattice.type() == generator.latticeType(),
				"Mismatch between lattice types");
			LASSERT(Generator::assetClass() == AssetClass::InterestRate,
				"The passed model must be interest-rate model");
			_buildBondTree(bondLattice, calibartedRateLattice, generator, nominal, couponRate, deltaTime, style);
		}

	};


	template<typename DeltaTime,
			typename Node>
	class BondBuilder<LatticeType::Trinomial, DeltaTime, Node> {
	private:
		template<typename LatticeObject,typename Generator>
		static void _buildNormal(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
			Node nominal, Node couponRate, DeltaTime const &deltaTime, DiscountingStyle style);

		template<typename LatticeObject,typename Generator>
		static void _buildReverting(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
			Node nominal, Node couponRate, DeltaTime const &deltaTime, DiscountingStyle style);

		template<typename LatticeObject, typename Generator>
		static void _buildBondTree(LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
			Node nominal, Node couponRate, DeltaTime const &deltaTime, DiscountingStyle style);

	public:
		template<typename LatticeObject, typename Generator>
		void operator()(LatticeObject &bondLattice, LatticeObject const &calibartedRateLattice, Generator const &generator,
			Node nominal, Node couponRate, DeltaTime const &deltaTime, DiscountingStyle style = DiscountingStyle::Continuous) {
			LASSERT(bondLattice.timeDimension() == calibartedRateLattice.timeDimension(),
				"Lattices must have the same dimension.");
			LASSERT(bondLattice.type() == generator.latticeType(),
				"Mismatch between lattice types");
			LASSERT(Generator::assetClass() == AssetClass::InterestRate,
				"The passed model must be interest-rate model");
			_buildBondTree(bondLattice, calibartedRateLattice, generator, nominal, couponRate, deltaTime, style);
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
void lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Binomial,DeltaTime,Node>::
_buildBondTree(LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
	Node nominal, Node couponRate, DeltaTime const &deltaTime, lattice_types::DiscountingStyle style) {

	typedef DiscountingFactor<Node> DCF;
	typedef DeltaTimeHolder<DeltaTime> DT;
	auto dcf = DCF::function(style);

	const std::size_t lastIdx = bondLattice.timeDimension() - 1;
	const std::size_t lastNodesSize = bondLattice.nodesAtIdx(lastIdx).size();
	Node dt{};
	Node coupon{};
	Node prob = generator.nodeRiskNeutralProb();

	for (auto i = 0; i < lastNodesSize; ++i) {
		bondLattice(lastIdx, i) = nominal;
	}

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n > 0; --n) {
		nodesSize = bondLattice.nodesAtIdx(n).size();
		dt = DT::deltaTime(n, deltaTime);
		coupon = nominal * couponRate * dt;
		for (auto i = 0; i < nodesSize; ++i) {
			bondLattice(n, i) =
				prob * (bondLattice(n + 1, i) + bondLattice(n + 1, i + 1) + 2.0 * coupon)*
				dcf(calibratedRateLattice(n, i), dt);
		}
	}

	dt = DT::deltaTime(0, deltaTime);
	bondLattice(0, 0) = prob * (bondLattice(1, 0) + bondLattice(1, 1) + 2.0*coupon)*
		dcf(calibratedRateLattice(0, 0), dt);
}


template<typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, DeltaTime, Node>::
_buildNormal(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
	Node nominal, Node couponRate, DeltaTime const &deltaTime, lattice_types::DiscountingStyle style) {

	typedef DiscountingFactor<Node> DCF;
	typedef DeltaTimeHolder<DeltaTime> DT;
	auto dcf = DCF::function(style);

	std::size_t const lastIdx = bondLattice.timeDimension() - 1;
	std::size_t const revertBranchesSize = timeIdx - 1;

	std::tuple<Node, Node, Node> prob;
	Node dt{};
	Node coupon{};

	std::size_t nodesSize{ 0 };
	for (auto n = timeIdx - 1; n > 0; --n) {
		nodesSize = bondLattice.nodesAtIdx(n).size();
		dt = DT::deltaTime(n, deltaTime);
		coupon = nominal * couponRate * dt;
		for (auto l = 0; l < nodesSize; ++l) {
			prob = generator.nodeRiskNeutralProb(revertBranchesSize, nodesSize, l, dt);
			bondLattice(n, l) = (std::get<0>(prob)*(bondLattice(n + 1, l) + coupon) +
				std::get<1>(prob)*(bondLattice(n + 1, l + 1) + coupon) +
				std::get<2>(prob)*(bondLattice(n + 1, l + 2) + coupon))*
				dcf(calibratedRateLattice(n, l), dt);

		}
	}

	dt = DT::deltaTime(0, deltaTime);
	bondLattice(0, 0) = (std::get<0>(prob)*(bondLattice(1, 0) + coupon) +
		std::get<1>(prob)*(bondLattice(1, 1) + coupon) +
		std::get<2>(prob)*(bondLattice(1, 2) + coupon))*
		dcf(calibratedRateLattice(0, 0), dt);
}


template<typename DeltaTime,
	typename Node>
	template<typename LatticeObject, typename Generator>
void lattice_bond_builders::BondBuilder<lattice_types::LatticeType::Trinomial, DeltaTime, Node>::
_buildReverting(std::size_t timeIdx, LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
	Node nominal, Node couponRate, DeltaTime const &deltaTime, lattice_types::DiscountingStyle style) {

	typedef DiscountingFactor<Node> DCF;
	typedef DeltaTimeHolder<DeltaTime> DT;
	auto dcf = DCF::function(style);

	std::size_t const lastIdx = bondLattice.timeDimension() - 1;
	std::size_t const revertBranchesSize = timeIdx - 1;

	std::tuple<Node, Node, Node> prob;
	Node dt{};
	Node coupon{};

	std::size_t nodesSize{ 0 };
	for (auto n = lastIdx - 1; n >= timeIdx; --n) {
		nodesSize = bondLattice.nodesAtIdx(n).size();
		dt = DT::deltaTime(n, deltaTime);
		coupon = nominal * couponRate * dt;
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
_buildBondTree(LatticeObject &bondLattice, LatticeObject const &calibratedRateLattice, Generator const &generator,
	Node nominal, Node couponRate, DeltaTime const &deltaTime, lattice_types::DiscountingStyle style) {

	const std::size_t firstRevertIdx = calibratedRateLattice.firstRevertingIdx();
	const std::size_t treeSize = calibratedRateLattice.timeDimension();

	const std::size_t lastIdx = treeSize - 1;
	const std::size_t lastNodesSize = calibratedRateLattice.nodesAtIdx(lastIdx).size();

	for (auto i = 0; i < lastNodesSize; ++i) {
		bondLattice(lastIdx, i) = nominal;
	}

	if (firstRevertIdx == 0) {
		// This trinomial tree does not have reverting property:
		_buildNormal(lastIdx, bondLattice, calibratedRateLattice, generator, nominal, couponRate, deltaTime, style);
	}
	else {
		// This trinomial tree does have reverting property:
		_buildReverting(firstRevertIdx - 1, bondLattice, calibratedRateLattice, generator, nominal, couponRate, deltaTime, style);
		_buildNormal(firstRevertIdx - 1, bondLattice, calibratedRateLattice, generator, nominal, couponRate, deltaTime, style);

	}
}







#endif ///_LATTICE_BOND_BUILDERS