#pragma once
#if! defined(_LATTICE_CALIBRATOR_IR)
#define _LATTICE_CALIBRATOR_IR

#include<memory>
#include<cmath>
#include"lattice_types.h"
#include"lattice_utility.h"
#include"lattice_calibrator_results.h"
#include"lattice_miscellaneous.h"


#include"optimization_methods/unconstrained_optimization.h"

namespace lattice_calibrator_ir {

	using lattice_types::LatticeType;
	using lattice_calibrator_results::CalibratorResults;
	using lattice_types::AssetClass;
	using lattice_types::DiscountingStyle;
	using lattice_types::MinimizerMethod;
	using optimization_odm::GoldenSectionMethod;
	using optimization_odm::PowellMethod;
	using optimization_odm::Range;
	using lattice_miscellaneous::lerp;
	using lattice_utility::DeltaTimeHolder;


	// ==============================================================================
	// ========================== DiscountingFactor =================================
	// ==============================================================================

	template<typename T>
	struct DiscountingFactor {

		static std::function<T(T,T)> function(DiscountingStyle style) {
			if (style == DiscountingStyle::Continuous) {
				return [=](T rate, T delta)->T {
					return std::exp(-1.0*rate*delta);
				};
			}
			return [=](T rate, T delta)->T {
				return (1.0 / (1.0 + rate * delta));
			};
		}


	};



	// ==============================================================================
	// ==================== DiscountCurveHolder =====================================
	// ==============================================================================

	template<typename DiscountCurve>
	struct DiscountCurveHolder {
	private:
		static typename DiscountCurve::value_type const _discount_impl(std::size_t idx, DiscountCurve const &discountCurve, std::true_type) {
			return discountCurve.at(idx);
		}

		static DiscountCurve const _discount_impl(std::size_t idx, DiscountCurve const &discountCurve, std::false_type) {
			return discountCurve;
		}
	public:
		static auto const discount(std::size_t idx, DiscountCurve const &discountCurve) {
			return _discount_impl(idx, discountCurve, std::is_compound<DiscountCurve>());
		}

	};





	// ==============================================================================
	// ========================== CalibratorIR ======================================
	// ==============================================================================

	template<LatticeType Type,
		typename TimeAxis,
		typename DeltaTime,
		typename DiscountCurve>
	struct CalibratorIR {};


	template<typename TimeAxis,
			typename DeltaTime,
			typename DiscountCurve>
	struct CalibratorIR<LatticeType::Binomial, TimeAxis, DeltaTime, DiscountCurve> {
	private:
		template<typename LatticeObject, typename Generator>
		static std::shared_ptr<CalibratorResults<AssetClass::InterestRate, LatticeObject>> const
			_calibrate_impl(LatticeObject &rateLattice, Generator &&generator,
				DeltaTime const &deltaTime, DiscountCurve const &discountCurve);

	public:
		template<typename LatticeObject, typename Generator>
		static std::shared_ptr<CalibratorResults<AssetClass::InterestRate, LatticeObject>> const
			calibrate(LatticeObject &rateLattice, Generator &&generator,
						DeltaTime const &deltaTime, DiscountCurve const &discountCurve) {
			return _calibrate_impl(rateLattice, std::forward<Generator>(generator),
				deltaTime, discountCurve);
		}

	};

	template<typename TimeAxis,
		typename DeltaTime,
		typename DiscountCurve>
	struct CalibratorIR<LatticeType::Trinomial, TimeAxis, DeltaTime, DiscountCurve> {
	private:
		template<typename LatticeObject, typename Generator>
		static std::shared_ptr<CalibratorResults<AssetClass::InterestRate, LatticeObject>> const
			_calibrate_impl(LatticeObject &rateLattice, Generator &&generator,
				DeltaTime const &deltaTime, DiscountCurve const &discountCurve);


	public:
		template<typename LatticeObject, typename Generator>
		static std::shared_ptr<CalibratorResults<AssetClass::InterestRate, LatticeObject>> const
			calibrate(LatticeObject &rateLattice, Generator &&generator,
				DeltaTime const &deltaTime, DiscountCurve const &discountCurve) {
			return _calibrate_impl(rateLattice, std::forward<Generator>(generator),
				deltaTime, discountCurve);
		}
	
	
	};


}





template<typename TimeAxis,
	typename DeltaTime,
	typename DiscountCurve>
template<typename LatticeObject, typename Generator>
std::shared_ptr<lattice_calibrator_results::CalibratorResults<lattice_types::AssetClass::InterestRate, LatticeObject>> const
lattice_calibrator_ir::CalibratorIR<lattice_types::LatticeType::Binomial,TimeAxis, DeltaTime, DiscountCurve>::
	_calibrate_impl(LatticeObject &rateLattice, Generator &&generator,DeltaTime const &deltaTime,
		DiscountCurve const &discountCurve) {

	// typedef node type:
	typedef typename LatticeObject::Node_type Node;
	// typedef discount curve holder:
	typedef DiscountCurveHolder<DiscountCurve> DC;
	// typedef discounting factor:
	typedef DiscountingFactor<Node> DCF;
	// typedef deltaTime holder:
	typedef DeltaTimeHolder<DeltaTime> DT;


	// get the objective function for calibration:
	auto calibFun = generator.calibrationObjective();
	// get calibration forward generator:
	auto forwardGen = generator.calibrationForwardGenerator();
	// get correct discounting factor:
	auto dcf = DCF::function(DiscountingStyle::Continuous);

	// instantiate the optimization:
	auto range = Range<Node>{ -1.0,1.0 };
	auto tol = 10e-7;
	GoldenSectionMethod<Node> gsm{ range,tol };
	// prepare container for optimizers:
	std::vector<std::tuple<Node, Node, Node, std::size_t>> optimizers;


	const std::size_t latticeSize = rateLattice.timeDimension();
	Node theta{ 0.0 };
	std::tuple<Node, Node> tuple;
	std::size_t nodesSize{ 0 };

	LatticeObject arrowDebreuLattice(rateLattice);
	// Take 10 % of the deltaTime:
	const Node dt_0 = 0.01 * DT::deltaTime(0, deltaTime);
	// Linear interpolation between today discount factor (which is 1) and next discount factor 
	// with maturity dt_0:
	const Node d_0 = lerp(DC::discount(0, discountCurve), DC::discount(1, discountCurve), dt_0);
	// get the near-today short rate:
	const Node rateApex = (-1.0*std::log(d_0) / dt_0);
	rateLattice(0, 0) = rateApex;
	arrowDebreuLattice(0, 0) = 1.0;
	Node dt{};

	for (std::size_t t = 1; t < latticeSize; ++t) {
		dt = DT::deltaTime(t - 1, deltaTime);
		arrowDebreuLattice(t, 0) = (0.5 *dcf(rateLattice(t - 1, 0), dt))*arrowDebreuLattice(t - 1, 0);
		if (t >= 2) {
			nodesSize = rateLattice.nodesAtIdx(t).size();
			for (std::size_t l = 1; l < nodesSize - 1; ++l) {

				arrowDebreuLattice(t, l) = (0.5 * dcf(rateLattice(t - 1, l), dt))*arrowDebreuLattice(t - 1, l) +
					(0.5 * dcf(rateLattice(t - 1, l - 1), dt))*arrowDebreuLattice(t - 1, l - 1);

			}
		}
		arrowDebreuLattice(t, t) = (0.5 * dcf(rateLattice(t - 1, t - 1), dt))*arrowDebreuLattice(t - 1, t - 1);

		// one-dimensional optimization:
		auto objective = std::bind(calibFun, std::placeholders::_1,
			dt, DC::discount(t, discountCurve), rateLattice.nodesAtIdx(t - 1),
			arrowDebreuLattice.nodesAtIdx(t), dcf);
		auto optimizer = gsm(objective);
		
		// forward induction: 
		for (std::size_t l = 0; l < rateLattice.nodesAtIdx(t - 1).size(); ++l) {
			tuple = forwardGen(std::get<0>(optimizer), rateLattice(t - 1, l), dt, l);
			rateLattice(t, l) = std::get<0>(tuple);
			rateLattice(t, l + 1) = std::get<1>(tuple);
		}

		optimizers.emplace_back(optimizer);
	}
	return std::shared_ptr<lattice_calibrator_results::
		CalibratorResults<lattice_types::AssetClass::InterestRate,
		LatticeObject>>{new lattice_calibrator_results::
		CalibratorResults<lattice_types::AssetClass::InterestRate,
		LatticeObject>(arrowDebreuLattice, optimizers) };
}




template<typename TimeAxis,
	typename DeltaTime,
	typename DiscountCurve>
	template<typename LatticeObject, typename Generator>
std::shared_ptr<lattice_calibrator_results::CalibratorResults<lattice_types::AssetClass::InterestRate, LatticeObject>> const
lattice_calibrator_ir::CalibratorIR<lattice_types::LatticeType::Trinomial, TimeAxis, DeltaTime, DiscountCurve>::
_calibrate_impl(LatticeObject &rateLattice, Generator &&generator, DeltaTime const &deltaTime,
	DiscountCurve const &discountCurve) {

	// typedef node type:
	typedef typename LatticeObject::Node_type Node;
	// typedef discount curve holder:
	typedef DiscountCurveHolder<DiscountCurve> DC;
	// typedef discounting factor:
	typedef DiscountingFactor<Node> DCF;
	// typedef deltaTime holder:
	typedef DeltaTimeHolder<DeltaTime> DT;


	// get the objective function for calibration:
	auto calibFun = generator.calibrationObjective();
	// get calibration forward generator:
	auto forwardGen = generator.calibrationForwardGenerator();
	// get correct discounting factor:
	auto dcf = DCF::function(DiscountingStyle::Continuous);

	// instantiate the optimization:
	auto range = Range<Node>{ -1.0,1.0 };
	auto tol = 10e-7;
	GoldenSectionMethod<Node> gsm{ range,tol };
	// prepare container for optimizers:
	std::vector<std::tuple<Node, Node, Node, std::size_t>> optimizers;


	const std::size_t latticeSize = rateLattice.timeDimension();
	Node theta{ 0.0 };
	std::tuple<Node, Node> tuple;
	std::size_t nodesSize{ 0 };

	LatticeObject arrowDebreuLattice(rateLattice);
	// Take 10 % of the deltaTime:
	const Node dt_0 = 0.01 * DT::deltaTime(0, deltaTime);
	// Linear interpolation between today discount factor (which is 1) and next discount factor 
	// with maturity dt_0:
	const Node d_0 = lerp(DC::discount(0, discountCurve), DC::discount(1, discountCurve), dt_0);
	// get the near-today short rate:
	const Node rateApex = (-1.0*std::log(d_0) / dt_0);
	rateLattice(0, 0) = rateApex;
	arrowDebreuLattice(0, 0) = 1.0;
	Node dt{};

	for (std::size_t t = 1; t < latticeSize; ++t) {
		dt = DT::deltaTime(t - 1, deltaTime);
		arrowDebreuLattice(t, 0) = (0.5 *dcf(rateLattice(t - 1, 0), dt))*arrowDebreuLattice(t - 1, 0);
		if (t >= 2) {
			nodesSize = rateLattice.nodesAtIdx(t).size();
			for (std::size_t l = 1; l < nodesSize - 1; ++l) {

				arrowDebreuLattice(t, l) = (0.5 * dcf(rateLattice(t - 1, l), dt))*arrowDebreuLattice(t - 1, l) +
					(0.5 * dcf(rateLattice(t - 1, l - 1), dt))*arrowDebreuLattice(t - 1, l - 1);

			}
		}
		arrowDebreuLattice(t, t) = (0.5 * dcf(rateLattice(t - 1, t - 1), dt))*arrowDebreuLattice(t - 1, t - 1);

		// one-dimensional optimization:
		auto objective = std::bind(calibFun, std::placeholders::_1,
			dt, DC::discount(t, discountCurve), rateLattice.nodesAtIdx(t - 1),
			arrowDebreuLattice.nodesAtIdx(t), dcf);
		auto optimizer = gsm(objective);

		// forward induction: 
		for (std::size_t l = 0; l < rateLattice.nodesAtIdx(t - 1).size(); ++l) {
			tuple = forwardGen(std::get<0>(optimizer), rateLattice(t - 1, l), dt, l);
			rateLattice(t, l) = std::get<0>(tuple);
			rateLattice(t, l + 1) = std::get<1>(tuple);
		}

		optimizers.emplace_back(optimizer);
	}
	return std::shared_ptr<lattice_calibrator_results::
		CalibratorResults<lattice_types::AssetClass::InterestRate,
		LatticeObject>>{new lattice_calibrator_results::
		CalibratorResults<lattice_types::AssetClass::InterestRate,
		LatticeObject>(arrowDebreuLattice, optimizers) };

}






#endif ///_LATTICE_CALIBRATOR_IR