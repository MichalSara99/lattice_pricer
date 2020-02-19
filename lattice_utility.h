#pragma once
#if !defined(_LATTICE_UTILITY)
#define _LATTICE_UTILITY

#include"lattice_structure.h"
#include<tuple>

namespace lattice_utility {


	using lattice_structure::GeneralLattice;
	using lattice_types::LatticeType;


	namespace {
		template<LatticeType Type,
			typename Node,
			typename TimeAxis,
			typename NodeContainerType>
			void _print_lattice_impl(GeneralLattice<Type, Node, TimeAxis, NodeContainerType> const &lattice,
				typename GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::Iterator_type begin,
				typename GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::Iterator_type end,
				std::true_type) {
			for (auto itr = begin; itr != end; ++itr) {
				std::cout << "[";
				for (auto const &e : (*itr)) {
					std::cout << e << ",";
				}
				std::cout << "]\n";
			}
			std::cout << "\n";
		}

		template<LatticeType Type,
			typename Node,
			typename TimeAxis,
			typename NodeContainerType>
			void _print_lattice_impl(GeneralLattice<Type, Node, TimeAxis, NodeContainerType> const &lattice,
				typename GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::Iterator_type begin,
				typename GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::Iterator_type end,
				std::false_type) {
			for (auto itr = begin; itr != end; ++itr) {
				std::cout << "(" << (*itr).first << "): [";
				for (auto const &e : (*itr).second) {
					std::cout << e << ",";
				}
				std::cout << "]\n";
			}
			std::cout << "\n";
		}


		template<LatticeType Type,
			typename Node,
			typename TimeAxis,
			typename NodeContainerType>
			void _print_lattice_c_impl(GeneralLattice<Type, Node, TimeAxis, NodeContainerType> const &lattice,
				typename GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::Const_iterator_type cbegin,
				typename GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::Const_iterator_type cend,
				std::true_type) {
			for (auto itr = cbegin; itr != cend; ++itr) {
				std::cout << "[";
				for (auto const &e : (*itr)) {
					std::cout << e << ",";
				}
				std::cout << "]\n";
			}
			std::cout << "\n";
		}

		template<LatticeType Type,
			typename Node,
			typename TimeAxis,
			typename NodeContainerType>
			void _print_lattice_c_impl(GeneralLattice<Type, Node, TimeAxis, NodeContainerType> const &lattice,
				typename GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::Const_iterator_type cbegin,
				typename GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::Const_iterator_type cend,
				std::false_type) {
			for (auto itr = cbegin; itr != cend; ++itr) {
				std::cout << "(" << (*itr).first << "): [";
				for (auto const &e : (*itr).second) {
					std::cout << e << ",";
				}
				std::cout << "]\n";
			}
			std::cout << "\n";
		}

	};

	template<LatticeType Type,
			typename Node,
			typename TimeAxis,
			typename NodeContainerType>
	void print(GeneralLattice<Type, Node, TimeAxis, NodeContainerType> const &lattice,
		typename GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::Iterator_type begin,
		typename GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::Iterator_type end) {
		_print_lattice_impl(lattice, begin, end, std::is_integral<TimeAxis>());
	}


	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		void print(GeneralLattice<Type, Node, TimeAxis, NodeContainerType> const &lattice,
			typename GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::Const_iterator_type cbegin,
			typename GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::Const_iterator_type cend) {
		_print_lattice_c_impl(lattice, cbegin, cend, std::is_integral<TimeAxis>());
	}


};



#endif ///_LATTICE_UTILITY