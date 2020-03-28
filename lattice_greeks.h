#pragma once
#if !defined(_LATTICE_GREEKS)
#define _LATTICE_GREEKS


#include"lattice_structure.h"
#include<numeric>
#include<cassert>

namespace lattice_greeks {

	using lattice_types::LatticeType;
	using lattice_types::LeafForwardGenerator;
	using lattice_types::LeafBackwardGenerator;
	using lattice_types::Payoff;
	using lattice_types::PayoffAdjuster;
	using lattice_structure::IndexedLattice;
	using lattice_structure::Lattice;


	// implementation details of delta:

	namespace {

		template<typename Node>
		Node _delta_indexed_binomial_impl(IndexedLattice<LatticeType::Binomial, Node> const &underlyingLattice,
											IndexedLattice<LatticeType::Binomial, Node> const &priceLattice) {
			assert(priceLattice.timeDimension() >= 2);
			assert(underlyingLattice.timeDimension() >= 2);
			auto p_first = priceLattice.cbegin();
			auto u_first = underlyingLattice.cbegin();
			auto p_itr = std::next(p_first);
			auto u_itr = std::next(u_first);
			auto u_change = (u_itr->at(1) - u_itr->at(0));
			auto p_change = (p_itr->at(1) - p_itr->at(0));
			return (p_change / u_change);
		}

		template<typename Node>
		Node _delta_indexed_trinomial_impl(IndexedLattice<LatticeType::Trinomial, Node> const &underlyingLattice,
											IndexedLattice<LatticeType::Trinomial, Node> const &priceLattice) {
			throw std::exception("Not yet implemented.");
		}


		template<typename Node, typename TimeAxis>
		Node _delta_lattice_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis> const& underlyingLattice,
											Lattice<LatticeType::Binomial, Node, TimeAxis> const &priceLattice) {
			assert(priceLattice.timeDimension() >= 2);
			assert(underlyingLattice.timeDimension() >= 2);
			auto p_first = priceLattice.cbegin();
			auto u_first = underlyingLattice.cbegin();
			auto p_itr = std::next(p_first);
			auto u_itr = std::next(u_first);
			auto u_change = (u_itr->at(1) - u_itr->at(0));
			auto p_change = (p_itr->at(1) - p_itr->at(0));
			return (p_change / u_change);
		}

		template<typename Node, typename TimeAxis>
		Node _delta_lattice_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>const& underlyingLattice,
											Lattice<LatticeType::Trinomial, Node, TimeAxis>const& priceLattice) {
			throw std::exception("Not yet implemented.");
		}


	}




	template<typename Node>
	Node delta(IndexedLattice<LatticeType::Binomial, Node> const &underlyingLattice,
		IndexedLattice<LatticeType::Binomial, Node> const &priceLattice) {
		return _delta_indexed_binomial_impl(underlyingLattice, priceLattice);
	}

	template<typename Node>
	Node delta(IndexedLattice<LatticeType::Trinomial, Node> const &underlyingLattice,
		IndexedLattice<LatticeType::Trinomial, Node> const &priceLattice) {
		return _delta_indexed_trinomial_impl(underlyingLattice, priceLattice);
	}


	template<typename Node, typename TimeAxis>
	Node delta(Lattice<LatticeType::Binomial, Node, TimeAxis> const& underlyingLattice,
		Lattice<LatticeType::Binomial,Node,TimeAxis> const &priceLattice) {
		return _delta_lattice_binomial_impl(underlyingLattice, priceLattice);
	}

	template<typename Node, typename TimeAxis>
	Node delta(Lattice<LatticeType::Trinomial, Node, TimeAxis>const& underlyingLattice,
		Lattice<LatticeType::Trinomial, Node, TimeAxis>const& priceLattice) {
		return _delta_lattice_trinomial_impl(underlyingLattice, priceLattice);
	}

	// implementation details of gamma:

	namespace {

		template<typename Node>
		Node _gamma_indexed_binomial_impl(IndexedLattice<LatticeType::Binomial, Node> const &underlyingLattice,
			IndexedLattice<LatticeType::Binomial, Node> const &priceLattice) {
			throw std::exception("Not yet implemented.");
		}

		template<typename Node>
		Node _gamma_indexed_trinomial_impl(IndexedLattice<LatticeType::Trinomial, Node> const &underlyingLattice,
			IndexedLattice<LatticeType::Trinomial, Node> const &priceLattice) {
			throw std::exception("Not yet implemented.");
		}


		template<typename Node, typename TimeAxis>
		Node _gamma_lattice_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis> const& underlyingLattice,
			Lattice<LatticeType::Binomial, Node, TimeAxis> const &priceLattice) {
			throw std::exception("Not yet implemented.");
		}

		template<typename Node, typename TimeAxis>
		Node _gamma_lattice_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>const& underlyingLattice,
			Lattice<LatticeType::Trinomial, Node, TimeAxis>const& priceLattice) {
			throw std::exception("Not yet implemented.");
		}

	}


	template<typename Node>
	Node gamma(IndexedLattice<LatticeType::Binomial, Node> const &underlyingLattice,
		IndexedLattice<LatticeType::Binomial, Node> const &priceLattice) {
		return _gamma_indexed_binomial_impl(underlyingLattice, priceLattice);
	}

	template<typename Node>
	Node gamma(IndexedLattice<LatticeType::Trinomial, Node> const &underlyingLattice,
		IndexedLattice<LatticeType::Trinomial, Node> const &priceLattice) {
		return _gamma_indexed_trinomial_impl(underlyingLattice, priceLattice);
	}


	template<typename Node, typename TimeAxis>
	Node gamma(Lattice<LatticeType::Binomial, Node, TimeAxis> const& underlyingLattice,
		Lattice<LatticeType::Binomial, Node, TimeAxis> const &priceLattice) {
		return _gamma_lattice_binomial_impl(underlyingLattice, priceLattice);
	}

	template<typename Node, typename TimeAxis>
	Node gamma(Lattice<LatticeType::Trinomial, Node, TimeAxis>const& underlyingLattice,
		Lattice<LatticeType::Trinomial, Node, TimeAxis>const& priceLattice) {
		return _gamma_lattice_trinomial_impl(underlyingLattice, priceLattice);
	}


	// implementation details of theta:

	namespace {

		template<typename Node>
		Node _theta_indexed_binomial_impl(IndexedLattice<LatticeType::Binomial, Node> const &underlyingLattice,
			IndexedLattice<LatticeType::Binomial, Node> const &priceLattice) {
			throw std::exception("Not yet implemented.");
		}

		template<typename Node>
		Node _theta_indexed_trinomial_impl(IndexedLattice<LatticeType::Trinomial, Node> const &underlyingLattice,
			IndexedLattice<LatticeType::Trinomial, Node> const &priceLattice) {
			throw std::exception("Not yet implemented.");
		}


		template<typename Node, typename TimeAxis>
		Node _theta_lattice_binomial_impl(Lattice<LatticeType::Binomial, Node, TimeAxis> const& underlyingLattice,
			Lattice<LatticeType::Binomial, Node, TimeAxis> const &priceLattice) {
			throw std::exception("Not yet implemented.");
		}

		template<typename Node, typename TimeAxis>
		Node _theta_lattice_trinomial_impl(Lattice<LatticeType::Trinomial, Node, TimeAxis>const& underlyingLattice,
			Lattice<LatticeType::Trinomial, Node, TimeAxis>const& priceLattice) {
			throw std::exception("Not yet implemented.");
		}

	}


	template<typename Node>
	Node theta(IndexedLattice<LatticeType::Binomial, Node> const &underlyingLattice,
		IndexedLattice<LatticeType::Binomial, Node> const &priceLattice) {
		return _theta_indexed_binomial_impl(underlyingLattice, priceLattice);
	}

	template<typename Node>
	Node theta(IndexedLattice<LatticeType::Trinomial, Node> const &underlyingLattice,
		IndexedLattice<LatticeType::Trinomial, Node> const &priceLattice) {
		return _theta_indexed_trinomial_impl(underlyingLattice, priceLattice);
	}


	template<typename Node, typename TimeAxis>
	Node theta(Lattice<LatticeType::Binomial, Node, TimeAxis> const& underlyingLattice,
		Lattice<LatticeType::Binomial, Node, TimeAxis> const &priceLattice) {
		return _theta_lattice_binomial_impl(underlyingLattice, priceLattice);
	}

	template<typename Node, typename TimeAxis>
	Node theta(Lattice<LatticeType::Trinomial, Node, TimeAxis>const& underlyingLattice,
		Lattice<LatticeType::Trinomial, Node, TimeAxis>const& priceLattice) {
		return _theta_lattice_trinomial_impl(underlyingLattice, priceLattice);
	}



}



#endif ///_LATTICE_GREEKS