#pragma once
#if !defined(_LATTICE_MODEL_COMPONENTS)
#define _LATTICE_MODEL_COMPONENTS

#include<cassert>
#include"lattice_types.h"
#include"lattice_miscellaneous.h"

namespace lattice_model_components {


	namespace leisen_reimer_inversion {

		using lattice_miscellaneous::sign;

		template<typename T=double>
		struct PeizerPrattFirstInversion {
		private:
			std::size_t n_;

		public:
			explicit PeizerPrattFirstInversion(std::size_t numberTimePoints)
				:n_{ numberTimePoints } {
				// numberTimePoints must be always odd number:
				assert(n_ % 2 != 0);
			}

			T operator()(T x)const {
				auto firstBra = (x / (n_ + (1.0 / 3.0)));
				auto secondBra = (n_ + (1.0 / 6.0));
				auto sqrt = std::sqrt(0.25 - 0.25 * std::exp(-1.0 * firstBra * firstBra * secondBra));
				return (0.5 + sign(x) * sqrt);
			}

		};

		template<typename T = double>
		struct PeizerPrattSecondInversion {
		private:
			std::size_t n_;

		public:
			explicit PeizerPrattSecondInversion(std::size_t numberTimePoints)
				:n_{ numberTimePoints } {
				// numberTimePoints must be always odd number:
				assert(n_ % 2 != 0);
			}

			T operator()(T x)const {
				auto firstBra = (x / (n_ + (1.0 / 3.0) + (0.1 / (n_ + 1.0))));
				auto secondBra = (n_ + (1.0 / 6.0));
				auto sqrt = std::sqrt(0.25 - 0.25 * std::exp(-1.0 * firstBra * firstBra * secondBra));
				return (0.5 + sign(x) * sqrt);
			}
		};





	}





};


#endif //_LATTICE_MODEL_COMPONENTS