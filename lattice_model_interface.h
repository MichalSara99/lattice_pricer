#pragma once
#if !defined(_LATTICE_MODEL_INTERFACE)
#define _LATTICE_MODEL_INTERFACE


#include<cmath>
#include"lattice_types.h"
#include"lattice_utility.h"
#include"lattice_miscellaneous.h"

namespace lattice_model {

	using lattice_types::LeafForwardGenerator;
	using lattice_types::LeafBackwardGenerator;
	using lattice_types::LatticeType;
	using lattice_types::AssetClass;
	using lattice_types::MinimizerMethod;
	using lattice_types::BranchingStyle;
	using lattice_utility::sign;

	template<std::size_t FactorCount, typename T>
	class BinomialModel {};

	template<std::size_t FactorCount, typename T>
	class TrinomialModel {};

	template<typename T>
	class BinomialModel<1, T> {
	public:
		// Forward generator:
		virtual std::tuple<T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting)const = 0;

		// Backward generator:
		virtual T operator()(T currValue, T upValue, T downValue, T dt) = 0;

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
			LeafForwardGenerator<T, T, T>> forwardGenerator() const = 0;

		// Factor count:
		enum { FactorCount = 2 };

		// LatticeType:
		LatticeType latticeType()const { return LatticeType::Binomial; }

	};

	template<typename T>
	class TrinomialModel<1, T> {
	public:
		// Forward generator:
		virtual std::tuple<T, T, T> operator()(T value, T dt, std::size_t leafIdx, std::size_t timeIdx, bool isMeanReverting)const = 0;

		// Backward generator:
		virtual T operator()(T currValue, T upValue, T midValue, T downValue, T dt,
			std::size_t revertBranchesSize, std::size_t nodesSize, std::size_t leafIdx) = 0;

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
			LeafForwardGenerator<T, T, T, T>> forwardGenerator()const = 0;

		// Forward generator 2:
		virtual std::pair<LeafBackwardGenerator<T, T, T, T, T, T>,
			LeafBackwardGenerator<T, T, T, T, T, T>> backwardGenerator()const = 0;

		// Factor count:
		enum { FactorCount = 2 };

		// LatticeType:
		LatticeType latticeType()const { return LatticeType::Trinomial; }

	};




}





#endif ///_LATTICE_MODEL_INTERFACE