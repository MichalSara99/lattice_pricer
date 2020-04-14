#pragma once
#if !defined(_LATTICE_POLICY)
#define _LATTICE_POLICY

#include"lattice_miscellaneous.h"
#include"lattice_types.h"
#include<memory>

namespace lattice_policy{

	using lattice_types::LatticeType;

	namespace greeks {

		template<LatticeType Type>
		struct DeltaPolicy{};

		template<>
		struct DeltaPolicy<LatticeType::Binomial> {

			template<typename Node,typename I1,typename I2>
			static const Node sensitivity(I1 uIter, I2 pIter) {
				Node u_change = (uIter->at(1) - uIter->at(0));
				auto p_change = (pIter->at(1) - pIter->at(0));
				return (p_change / u_change);
			}


		};

		template<>
		struct DeltaPolicy<LatticeType::Trinomial> {

			template<typename Node, typename I1, typename I2>
			static const Node sensitivity(I1 uIter, I2 pIter) {
				auto u_change = ((uIter->at(2) - uIter->at(1)) + (uIter->at(1) - uIter->at(0)));
				auto p_change = (pIter->at(2) - pIter->at(0));
				return (p_change / u_change);
			}
		};

		template<LatticeType Type>
		struct GammaPolicy {};

		template<>
		struct GammaPolicy<LatticeType::Binomial> {

			template<typename Node, typename I1, typename I2>
			static const Node sensitivity(I1 uIter, I2 pIter) {
				auto p_uu = pIter->at(2);
				auto p_ud = pIter->at(1);
				auto p_dd = pIter->at(0);
				auto u_uu = uIter->at(2);
				auto u_ud = uIter->at(1);
				auto u_dd = uIter->at(0);
				return (((p_uu - p_ud) / (u_uu - u_ud)) - ((p_ud - p_dd) / (u_ud - u_dd))) / (0.5*(u_uu - u_dd));;
			}


		};

		template<>
		struct GammaPolicy<LatticeType::Trinomial> {

			template<typename Node, typename I1, typename I2>
			static const Node sensitivity(I1 uIter, I2 pIter) {
				auto p_uu = pIter->at(2);
				auto p_ud = pIter->at(1);
				auto p_dd = pIter->at(0);
				auto u_uu = uIter->at(2);
				auto u_ud = uIter->at(1);
				auto u_dd = uIter->at(0);
				return ((p_uu - 2.0*p_ud + p_dd) / ((u_uu - u_ud)*(u_ud - u_dd)));
			}
		};

		template<LatticeType Type>
		struct ThetaPolicy {};

		template<>
		struct ThetaPolicy<LatticeType::Binomial> {

			template<typename Node, typename Time>
			static const Node sensitivity(Time delta1, Time delta2, Node mid, Node price) {
				return ((mid - price) / (delta1 + delta2));
			}

		};

		template<>
		struct ThetaPolicy<LatticeType::Trinomial> {

			template<typename Node, typename Time>
			static const Node sensitivity(Time delta1, Time delta2, Node mid, Node price) {
				return ((mid - price) / (delta1));
			}
		};



	}



}




#endif ///_LATTICE_POLICY