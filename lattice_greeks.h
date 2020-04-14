#pragma once
#if !defined(_LATTICE_GREEKS)
#define _LATTICE_GREEKS


#include"lattice_structure.h"
#include"lattice_policy.h"
#include<numeric>
#include<cassert>

namespace lattice_greeks {

	using lattice_types::LatticeType;
	using lattice_policy::greeks::DeltaPolicy;
	using lattice_policy::greeks::GammaPolicy;
	using lattice_policy::greeks::ThetaPolicy;

	// ==============================================================================
	// ================================= Greeks =====================================
	// ==============================================================================


	template<typename TimeAxis>
	class Greeks {
	private:
		template<typename LatticeObject,
			typename Itr1,
			typename Itr2,
			typename Policy = DeltaPolicy<LatticeObject::type()>>
		static auto _delta_impl(LatticeObject const &underlyingLattice, LatticeObject const &priceLattice,
			Itr1 uItr,Itr2 pItr,std::true_type) {
			return Policy::sensitivity<typename LatticeObject::Node_type>(uItr, pItr);
		}
		template<typename LatticeObject,
			typename Itr1,
			typename Itr2,
			typename Policy = DeltaPolicy<LatticeObject::type()>>
		static auto _delta_impl(LatticeObject const &underlyingLattice, LatticeObject const &priceLattice,
			Itr1 uItr, Itr2 pItr, std::false_type) {
			auto u_itr = &uItr->second;
			auto p_itr = &pItr->second;
			return Policy::sensitivity<typename LatticeObject::Node_type>(u_itr,p_itr);
		}

		template<typename LatticeObject,
			typename Itr1,
			typename Itr2,
			typename Policy = GammaPolicy<LatticeObject::type()>>
			static auto _gamma_impl(LatticeObject const &underlyingLattice, LatticeObject const &priceLattice,
				Itr1 uItr, Itr2 pItr, std::true_type) {
			return Policy::sensitivity<typename LatticeObject::Node_type>(uItr, pItr);
		}
		template<typename LatticeObject,
			typename Itr1,
			typename Itr2,
			typename Policy = GammaPolicy<LatticeObject::type()>>
			static auto _gamma_impl(LatticeObject const &underlyingLattice, LatticeObject const &priceLattice,
				Itr1 uItr, Itr2 pItr, std::false_type) {
			auto u_itr = &uItr->second;
			auto p_itr = &pItr->second;
			return Policy::sensitivity<typename LatticeObject::Node_type>(u_itr, p_itr);
		}

		template<typename LatticeObject,
			typename Itr,
			typename Node,
			typename Policy = ThetaPolicy<LatticeObject::type()>>
			static auto _theta_impl_(LatticeObject const &underlyingLattice, LatticeObject const &priceLattice,
				Itr pItr, Node delta1, Node delta2,Node price, std::true_type) {
			auto mid = pItr->at(1);
			return Policy::sensitivity<typename LatticeObject::Node_type>(delta1, delta2, mid, price);
		}

		template<typename LatticeObject,
			typename Itr,
			typename Node,
			typename Policy = ThetaPolicy<LatticeObject::type()>>
			static auto _theta_impl_(LatticeObject const &underlyingLattice, LatticeObject const &priceLattice,
				Itr pItr, Node delta1, Node delta2,Node price ,std::false_type) {
			auto mid = (pItr->second).at(1);
			return Policy::sensitivity<typename LatticeObject::Node_type>(delta1, delta2, mid, price);
		}

		template<typename LatticeObject,
			typename DeltaTime,
			typename Policy = ThetaPolicy<LatticeObject::type()>>
			static auto _theta_impl(LatticeObject const &underlyingLattice, LatticeObject const &priceLattice, DeltaTime const &deltaTime, std::true_type) {
			LASSERT(underlyingLattice.type() == priceLattice.type(), "Mismatch of lattice types");
			LASSERT(underlyingLattice.timeDimension() >= 3, "Lattice must have at least three periods");
			LASSERT(priceLattice.timeDimension() >= 3, "Lattice must have at least three periods");
			auto p_first = priceLattice.cbegin();
			auto price = priceLattice.apex();
			std::size_t offSet = 1;
			if (LatticeObject::type() == LatticeType::Binomial) {
				offSet = 2;
			}
			auto p_itr = std::next(p_first, offSet);
			auto delta_1 = deltaTime.at(0);
			auto delta_2 = deltaTime.at(1);

			return _theta_impl_(underlyingLattice, priceLattice, p_itr, delta_1, delta_2, price,
				std::is_integral<typename LatticeObject::TimeAxis_type>());
		}

		template<typename LatticeObject,
			typename DeltaTime,
			typename Policy = ThetaPolicy<LatticeObject::type()>>
			static auto _theta_impl(LatticeObject const &underlyingLattice, LatticeObject const &priceLattice, DeltaTime const &deltaTime, std::false_type) {
			LASSERT(underlyingLattice.type() == priceLattice.type(), "Mismatch of lattice types");
			LASSERT(underlyingLattice.timeDimension() >= 3, "Lattice must have at least three periods");
			LASSERT(priceLattice.timeDimension() >= 3, "Lattice must have at least three periods");
			auto p_first = priceLattice.cbegin();
			auto price = priceLattice.apex();
			std::size_t offSet = 1;
			if (LatticeObject::type() == LatticeType::Binomial) {
				offSet = 2;
			}
			auto p_itr = std::next(p_first, offSet);
			return _theta_impl_(underlyingLattice, priceLattice, p_itr, deltaTime, deltaTime, price,
				std::is_integral<typename LatticeObject::TimeAxis_type>());
		}

	public:
		template<typename LatticeObject>
			static auto delta(LatticeObject const &underlyingLattice, LatticeObject const &priceLattice) {
			LASSERT(underlyingLattice.type() == priceLattice.type(), "Mismatch of lattice types");
			LASSERT(underlyingLattice.timeDimension() >= 2, "Lattice must have at least three periods");
			LASSERT(priceLattice.timeDimension() >= 2, "Lattice must have at least three periods");
			auto p_first = priceLattice.cbegin();
			auto u_first = underlyingLattice.cbegin();
			auto p_itr = std::next(p_first);
			auto u_itr = std::next(u_first);
			return _delta_impl(underlyingLattice, priceLattice, u_itr, p_itr,
				std::is_integral<typename LatticeObject::TimeAxis_type>());
		}

		template<typename LatticeObject>
			static auto gamma(LatticeObject const &underlyingLattice, LatticeObject const &priceLattice) {
			LASSERT(underlyingLattice.type() == priceLattice.type(), "Mismatch of lattice types");
			LASSERT(underlyingLattice.timeDimension() >= 3, "Lattice must have at least three periods");
			LASSERT(priceLattice.timeDimension() >= 3, "Lattice must have at least three periods");
			auto p_first = priceLattice.cbegin();
			auto u_first = underlyingLattice.cbegin();
			std::size_t offSet = 1;
			if (LatticeObject::type() == LatticeType::Binomial) {
				offSet = 2;
			}
			auto p_itr = std::next(p_first, offSet);
			auto u_itr = std::next(u_first, offSet);
			return _gamma_impl(underlyingLattice, priceLattice, u_itr, p_itr,
				std::is_integral<typename LatticeObject::TimeAxis_type>());
		}

		template<typename LatticeObject,
			typename DeltaTime>
			static auto theta(LatticeObject const &underlyingLattice, LatticeObject const &priceLattice, DeltaTime const &deltaTime) {
			return _theta_impl(underlyingLattice, priceLattice, deltaTime, std::is_compound<DeltaTime>());
		}


	};



}



#endif ///_LATTICE_GREEKS