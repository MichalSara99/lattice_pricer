#pragma once
#if! defined(_LATTICE_CALIBRATOR_IR)
#define _LATTICE_CALIBRATOR_IR

#include<memory>
#include<cmath>
#include"lattice_types.h"
#include"lattice_calibrator_results.h"


#include"optimization_methods/unconstrained_optimization.h"

namespace lattice_calibrator_ir {

	using lattice_types::LatticeType;
	using lattice_calibrator_results::CalibratorResults;
	using lattice_types::AssetClass;
	using optimization_odm::GoldenSectionMethod;
	using optimization_odm::Range;

	// ==============================================================================
	// ========================== CalibratorIR ======================================
	// ==============================================================================

	template<LatticeType Type,
		typename TimeAxis,
		typename DeltaTime,
		typename DiscountCurve,
		typename Node>
	struct CalibratorIR {};


	template<typename TimeAxis,
			typename DeltaTime,
			typename DiscountCurve,
			typename Node>
	struct CalibratorIR<LatticeType::Binomial, TimeAxis, DeltaTime, DiscountCurve, Node> {
	private:
		static Node const _discount_impl(std::size_t idx, DiscountCurve const &discountCurve, std::true_type) {
			return discountCurve.at(idx);
		}

		static Node const _discount_impl(std::size_t idx, DiscountCurve const &discountCurve, std::false_type) {
			return discountCurve;
		}

		static Node const _discount(std::size_t idx, DiscountCurve const &discountCurve) {
			return _discount_impl(idx, discountCurve, std::is_compound<DiscountCurve>());
		}

		template<typename LatticeObject, typename Generator>
		static std::shared_ptr<CalibratorResults<AssetClass::InterestRate, LatticeObject, Node>> const
			_calibrate_impl(LatticeObject &rateLattice, Generator &&generator,
				DeltaTime const &deltaTime, DiscountCurve const &discountCurve, std::true_type);

		template<typename LatticeObject, typename Generator>
		static std::shared_ptr<CalibratorResults<AssetClass::InterestRate, LatticeObject, Node>> const
			_calibrate_impl(LatticeObject &rateLattice, Generator &&generator,
				DeltaTime const &deltaTime, DiscountCurve const &discountCurve, std::false_type);

	public:
		template<typename LatticeObject, typename Generator>
		static std::shared_ptr<CalibratorResults<AssetClass::InterestRate, LatticeObject, Node>> const
			calibrate(LatticeObject &rateLattice, Generator &&generator,
						DeltaTime const &deltaTime, DiscountCurve const &discountCurve) {
			return _calibrate_impl(rateLattice, std::forward<Generator>(generator),
				deltaTime, discountCurve, std::is_compound<DeltaTime>());
		}





	};




}

template<typename TimeAxis,
	typename DeltaTime,
	typename DiscountCurve,
	typename Node>
template<typename LatticeObject, typename Generator>
std::shared_ptr<lattice_calibrator_results::CalibratorResults<lattice_types::AssetClass::InterestRate, LatticeObject, Node>> const
lattice_calibrator_ir::CalibratorIR<lattice_types::LatticeType::Binomial,TimeAxis, DeltaTime, DiscountCurve, Node>::
	_calibrate_impl(LatticeObject &rateLattice, Generator &&generator,DeltaTime const &deltaTime,
		DiscountCurve const &discountCurve, std::true_type) {

	// get the objective function for calibration:
	auto calibFun = generator.calibrationObjective();
	// get calibration forward generator:
	auto forwardGen = generator.calibrationForwardGenerator();
	// instantiate the optimization:
	auto range = Range<Node>{ 0.005,1.0 };
	auto tol = 10e-7;
	GoldenSectionMethod<Node> gsm{ range,tol };
	// prepare container for optimizers:
	std::vector<std::tuple<Node, Node, Node, std::size_t>> optimizers;


	const std::size_t latticeSize = rateLattice.timeDimension();
	Node theta{ 0.0 };
	std::tuple<Node, Node> tuple;
	std::size_t nodesSize{ 0 };

	LatticeObject arrowDebreuLattice(rateLattice);
	const Node rateApex = (-1.0*std::log(_discount(1, discountCurve)) / deltaTime[0]);
	rateLattice(0, 0) = rateApex;
	arrowDebreuLattice(0, 0) = 1.0;

	for (std::size_t t = 1; t < latticeSize; ++t) {
		arrowDebreuLattice(t, 0) = (0.5 / (1.0 + rateLattice(t - 1, 0)*deltaTime[t - 1]))*arrowDebreuLattice(t - 1, 0);
		if (t >= 2) {
			nodesSize = rateLattice.nodesAtIdx(t).size();
			for (std::size_t l = 1; l < nodesSize - 1; ++l) {

				arrowDebreuLattice(t, l) = (0.5 / (1.0 + rateLattice(t - 1, l)*deltaTime[t - 1]))*arrowDebreuLattice(t - 1, l) +
					(0.5 / (1.0 + rateLattice(t - 1, l - 1)*deltaTime[t - 1]))*arrowDebreuLattice(t - 1, l - 1);

			}
		}
		arrowDebreuLattice(t, t) = (0.5 / (1.0 + rateLattice(t - 1, t - 1)*deltaTime[t - 1]))*arrowDebreuLattice(t - 1, t - 1);

		// one-dimensional optimization:
		auto objective = std::bind(calibFun, std::placeholders::_1,
			deltaTime[t - 1], _discount(t, discountCurve), arrowDebreuLattice.nodesAtIdx(t));
		auto optimizer = gsm(objective);
		
		// forward induction: 
		for (std::size_t l = 0; l < rateLattice.nodesAtIdx(t - 1).size(); ++l) {
			tuple = forwardGen(std::get<0>(optimizer), deltaTime[t - 1], l);
			rateLattice(t, l) = std::get<0>(tuple);
			rateLattice(t, l + 1) = std::get<1>(tuple);
		}

		optimizers.emplace_back(optimizer);
	}
	return std::shared_ptr<lattice_calibrator_results::
		CalibratorResults<lattice_types::AssetClass::InterestRate,
		LatticeObject, Node>>{new lattice_calibrator_results::
		CalibratorResults<lattice_types::AssetClass::InterestRate,
		LatticeObject, Node>(arrowDebreuLattice, optimizers) };

}

template<typename TimeAxis,
	typename DeltaTime,
	typename DiscountCurve,
	typename Node>
	template<typename LatticeObject, typename Generator>
std::shared_ptr<lattice_calibrator_results::CalibratorResults<lattice_types::AssetClass::InterestRate, LatticeObject, Node>> const
lattice_calibrator_ir::CalibratorIR<lattice_types::LatticeType::Binomial, TimeAxis, DeltaTime, DiscountCurve, Node>::
_calibrate_impl(LatticeObject &rateLattice, Generator &&generator, DeltaTime const &deltaTime,
	DiscountCurve const &discountCurve, std::false_type) {

	// get the objective function for calibration:
	auto calibFun = generator.calibrationObjective();
	// get calibration forward generator:
	auto forwardGen = generator.calibrationForwardGenerator();
	// instantiate the optimization:
	auto range = Range<Node>{ 0.005,1.0 };
	auto tol = 10e-7;
	GoldenSectionMethod<Node> gsm{ range,tol };
	// prepare container for optimizers:
	std::vector<std::tuple<Node, Node, Node, std::size_t>> optimizers;


	const std::size_t latticeSize = rateLattice.timeDimension();
	Node theta{ 0.0 };
	std::tuple<Node, Node> tuple;
	std::size_t nodesSize{ 0 };

	LatticeObject arrowDebreuLattice(rateLattice);
	const Node rateApex = (-1.0*std::log(_discount(1, discountCurve)) / deltaTime);
	rateLattice(0, 0) = rateApex;
	arrowDebreuLattice(0, 0) = 1.0;

	for (std::size_t t = 1; t < latticeSize; ++t) {
		arrowDebreuLattice(t, 0) = (0.5 / (1.0 + rateLattice(t - 1, 0)*deltaTime))*arrowDebreuLattice(t - 1, 0);
		if (t >= 2) {
			nodesSize = rateLattice.nodesAtIdx(t).size();
			for (std::size_t l = 1; l < nodesSize - 1; ++l) {

				arrowDebreuLattice(t, l) = (0.5 / (1.0 + rateLattice(t - 1, l)*deltaTime))*arrowDebreuLattice(t - 1, l) +
					(0.5 / (1.0 + rateLattice(t - 1, l - 1)*deltaTime))*arrowDebreuLattice(t - 1, l - 1);

			}
		}
		arrowDebreuLattice(t, t) = (0.5 / (1.0 + rateLattice(t - 1, t - 1)*deltaTime))*arrowDebreuLattice(t - 1, t - 1);

		// one-dimensional optimization:
		auto objective = std::bind(calibFun, std::placeholders::_1,
			deltaTime, _discount(t, discountCurve), arrowDebreuLattice.nodesAtIdx(t));
		auto optimizer = gsm(objective);

		// forward induction: 
		for (std::size_t l = 0; l < rateLattice.nodesAtIdx(t - 1).size(); ++l) {
			tuple = forwardGen(std::get<0>(optimizer), deltaTime, l);
			rateLattice(t, l) = std::get<0>(tuple); //down l
			rateLattice(t, l + 1) = std::get<1>(tuple); //up l + 1  
		}

		optimizers.emplace_back(optimizer);
	}
	return std::shared_ptr<lattice_calibrator_results::
		CalibratorResults<lattice_types::AssetClass::InterestRate,
		LatticeObject, Node>>{new lattice_calibrator_results::
		CalibratorResults<lattice_types::AssetClass::InterestRate,
		LatticeObject, Node>(arrowDebreuLattice, optimizers) };

}





#endif ///_LATTICE_CALIBRATOR_IR