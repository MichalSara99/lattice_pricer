#pragma once
#if !defined(_LATTICE_UTILITY)
#define _LATTICE_UTILITY

#include"lattice_structure.h"
#include"lattice_traits.h"
#include<tuple>

namespace lattice_utility {

	using lattice_types::LatticeType;
	using lattice_traits::PrintTraits;


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


};



#endif ///_LATTICE_UTILITY