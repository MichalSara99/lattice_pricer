#pragma once
#if !defined(_LATTICE_UTILITY)
#define _LATTICE_UTILITY

#include"lattice_structure.h"
#include"lattice_traits.h"
#include<tuple>

namespace lattice_utility {

	using lattice_types::LatticeType;
	using lattice_traits::PrintTraits;
	using lattice_types::DiscountingStyle;

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
	// ========================== DiscountingFactor =================================
	// ==============================================================================

	template<typename T>
	struct DiscountingFactor {

		static std::function<T(T, T)> function(DiscountingStyle style) {
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

	// ==============================================================================
	// ============================= LinearInterpolator =============================
	// ==============================================================================


	template<typename Container>
	class LinearInterpolator{
	private:
		Container x_;
		Container y_;

	public:
		typedef Container Container_type;
		typedef typename Container::value_type Container_value_type;

		LinearInterpolator()
			:x_(),y_(){}
		LinearInterpolator(LinearInterpolator const &other)
			:x_(other.x_),y_(other.y_){}
		LinearInterpolator(LinearInterpolator &&other)noexcept
			:x_(std::move(other.x_)),y_(std::move(other.y_)){}
		~LinearInterpolator(){}

		LinearInterpolator& operator=(LinearInterpolator const &other){
			if (this != &other) {
				x_ = other.x_;
				y_ = other.y_;
			}
			return *this;
		}
		LinearInterpolator& operator=(LinearInterpolator &&other)noexcept {
			if (this != &other) {
				x_ = std::move(other.x_);
				y_ = std::move(other.y_);
			}
			return *this;
		}

		void setPoints(Container const &xpoints, Container const &ypoints,bool sortXPoints = false) {
			x_ = xpoints;
			y_ = ypoints;

			if(sortXPoints)
				std::sort(x_.begin(), x_.end());

			for (std::size_t i = 0; i < x_.size(); ++i) {
				for (std::size_t j = 0; j < x_.size(); ++j) {
					if (x_[i] == xpoints[j]) {
						y_[i] = ypoints[j];
						break;
					}
				}
			}
		}

		Container_value_type const getValue(Container_value_type const &x)const {
			Container_value_type x0{}, y0{}, x1{}, y1{};
			if ((x < x_[0]) || (x > x_[x_.size() - 1])) {
				return 0.0; // outside of domain
			}
			for (std::size_t i = 0; i < x_.size(); ++i) {
				if (x_[i] < x) {
					x0 = x_[i];
					y0 = y_[i];
				}
				else if (x_[i] >= x) {
					x1 = x_[i];
					y1 = y_[i];
					break;
				}
			}
			return ((y0*(x - x1) / (x0 - x1)) + (y1 * (x - x0) / (x1 - x0)));
		}

	};



};




#endif ///_LATTICE_UTILITY