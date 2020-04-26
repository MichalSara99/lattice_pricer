#pragma once
#if !defined(_LATTICE_UTILITY)
#define _LATTICE_UTILITY

#include"lattice_structure.h"
#include"lattice_traits.h"
#include<tuple>

namespace lattice_utility {

	using lattice_types::LatticeType;
	using lattice_traits::PrintTraits;

	// ==============================================================================
	// ================================== print =====================================
	// ==============================================================================

	template<typename LatticeObject,
		typename cIter,
		typename Traits = PrintTraits<typename LatticeObject::Node_type>>
	void _print_impl(LatticeObject const &lattice, cIter cbegin, cIter cend,std::true_type) {

		for (auto itr = cbegin; itr != cend; ++itr) {
			std::cout << "["
				<< Traits::printLine(*itr);
			std::cout << "]\n";
		}
	}

	template<typename LatticeObject,
		typename cIter,
		typename Traits = PrintTraits<typename LatticeObject::Node_type>>
	void _print_impl(LatticeObject const &lattice, cIter cbegin,cIter cend, std::false_type) {
		for (auto itr = cbegin; itr != cend; ++itr) {
			std::cout << "(" << (*itr).first << "): ["
				<< Traits::printLine((*itr).second);
			std::cout << "]\n";
		}
	}


	template<typename LatticeObject,
		typename cIter,
		typename Traits = PrintTraits<typename LatticeObject::Node_type>>
		void print(LatticeObject const &lattice, cIter cbegin,cIter cend) {
		_print_impl(lattice, cbegin, cend, std::is_integral<typename LatticeObject::TimeAxis_type>());
	}



	// ==============================================================================
	// ======================== DeltaTimeHolder =====================================
	// ==============================================================================

	template<typename DeltaTime>
	struct DeltaTimeHolder {
	private:
		static auto const _deltaTime_impl(std::size_t idx, DeltaTime const &deltaTime, std::true_type) {
			return deltaTime.at(idx);
		}

		static DeltaTime const _deltaTime_impl(std::size_t idx, DeltaTime const &deltaTime, std::false_type) {
			return deltaTime;
		}
	public:
		static auto const deltaTime(std::size_t idx, DeltaTime const &deltaTime) {
			return _deltaTime_impl(idx, deltaTime, std::is_compound<DeltaTime>());
		}

	};


	// ==============================================================================
	// =================================== sign =====================================
	// ==============================================================================

	template<typename T>
	T sign(T x) {
		if (x < 0)
			return -1.0;
		else if (x > 0)
			return 1.0;
		else
			return 0.0;
	}

};



#endif ///_LATTICE_UTILITY