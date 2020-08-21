#pragma once
#if !defined(_LATTICE_UTILITY)
#define _LATTICE_UTILITY

#include"lattice_structure.h"
#include"lattice_traits.h"
#include"lattice_types.h"
#include<tuple>
#include<iostream>


namespace lattice_utility {

	using lattice_types::LatticeType;
	using lattice_traits::PrintTraits;
	using lattice_types::DiscountingStyle;
	using lattice_types::BarrierType;

	// ==============================================================================
	// ===================================== Logger =================================
	// ==============================================================================

	class Logger {
	protected:
		explicit Logger() {};

		void inline what(std::ostream &out, std::string text)const {
			out << text;
		}

	public:
		virtual ~Logger() {}

		Logger(Logger const &) = delete;
		Logger(Logger &&) = delete;
		Logger& operator=(Logger const &) = delete;
		Logger& operator=(Logger &&) = delete;

		static inline Logger& get() {
			static Logger logger;
			return logger;
		}

		void inline warning(std::ostream &out, std::string text)const {
			this->what(out, std::move(std::string{ "WARNING: " + text }));
		}

		void inline critical(std::ostream &out, std::string text)const {
			this->what(out, std::move(std::string{ "CRITICAL: " + text }));
		}
	};

	// ==============================================================================
	// ================================ BarrierComparer =============================
	// ==============================================================================

	template<typename T>
	struct BarrierComparer {
		static std::function<bool(T, T)> const comparer(BarrierType BType) {
			switch (BType)
			{
			case lattice_types::BarrierType::DownAndOut:
				return [](T stock, T barrier)->bool {
					return (stock > barrier);
				};
				break;
			case lattice_types::BarrierType::UpAndIn:
				return [](T stock, T barrier)->bool {
					return (stock >= barrier);
				};
				break;
			case lattice_types::BarrierType::DownAndIn:
				return [](T stock, T barrier)->bool {
					return (stock <= barrier);
				};
				break;
			case lattice_types::BarrierType::UpAndOut:
				return [](T stock, T barrier)->bool{
					return (stock < barrier);
				};
				break;
			}
		}
	};

	// ==============================================================================
	// ================================ DermanKaniErgenerAdjuster ===================
	// ==============================================================================

	template<typename T>
	using CheckAdjusterPair = std::pair<std::function<bool(T, T, T, T)>, std::function<T(T, T, T, T, T, T)>>;

	template<typename T>
	struct DermanKaniErgenerAdjuster {
		static CheckAdjusterPair<T> const adjuster(lattice_types::BarrierType BType) {
			auto cmp = BarrierComparer<T>::comparer(BType);
			switch (BType)
			{
				case lattice_types::BarrierType::DownAndOut:
				case lattice_types::BarrierType::UpAndIn:
					{
						auto checker =  [=](T stock, T stockDown,T stockUp, T barrier)->bool {
							return (cmp(stock,barrier) && (stockDown <= barrier));
						};

						auto adjuster = [](T stock,T stockDown,T stockUp,T barrier,T rebate,T optionPrice)->T {
							return (std::max(rebate - optionPrice, optionPrice - rebate) / (stockDown - stock))*(barrier - stock);
						};
						return std::make_pair(checker, adjuster);
					}
					break;
				case lattice_types::BarrierType::DownAndIn:
				case lattice_types::BarrierType::UpAndOut:
					{
						auto checker = [=](T stock, T stockDown, T stockUp, T barrier)->bool {
							return (cmp(stock, barrier) && (stockUp >= barrier));
						};

						auto adjuster = [](T stock, T stockDown, T stockUp, T barrier, T rebate, T optionPrice)->T {
							return (std::max(rebate - optionPrice, optionPrice - rebate) / (stockUp - stock))*(barrier - stock);
						};
						return std::make_pair(checker, adjuster);
					}
					break;
			}
		}
	};

	// ==============================================================================
	// ================================== print =====================================
	// ==============================================================================

	template<typename LatticeObject,
		typename cIter,
		typename Traits = PrintTraits<typename LatticeObject::Node_type, LatticeObject::type()>>
	void _print_impl(LatticeObject const &lattice, cIter cbegin, cIter cend,std::ostream &out, std::true_type) {

		for (auto itr = cbegin; itr != cend; ++itr) {
			out << "["
				<< Traits::printLine(*itr);
			out << "]\n";
		}
	}

	template<typename LatticeObject,
		typename cIter,
		typename Traits = PrintTraits<typename LatticeObject::Node_type, LatticeObject::type()>>
	void _print_impl(LatticeObject const &lattice, cIter cbegin,cIter cend, std::ostream &out, std::false_type) {
		std::string bl{ "" };
		if (LatticeObject::type() == LatticeType::TwoVariableBinomial)
			bl = "\n";
		for (auto itr = cbegin; itr != cend; ++itr) {
			out << "(" << (*itr).first << "):" << bl << "["
				<< Traits::printLine((*itr).second);
			out << "]\n";
		}
	}


	template<typename LatticeObject,
		typename cIter,
		typename Traits = PrintTraits<typename LatticeObject::Node_type, LatticeObject::type()>>
		void print(LatticeObject const &lattice, cIter cbegin,cIter cend,std::ostream &out = std::cout) {
		_print_impl(lattice, cbegin, cend, out, std::is_integral<typename LatticeObject::TimeAxis_type>());
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
	// ===================== RiskFreeRateHolder =====================================
	// ==============================================================================

	template<typename RiskFreeRate>
	struct RiskFreeRateHolder {
	private:
		static auto const _rate_impl(std::size_t idx, RiskFreeRate const &riskFreeRate, std::true_type) {
			return riskFreeRate.at(idx);
		}

		static auto const _rate_impl(std::size_t idx, RiskFreeRate const &riskFreeRate, std::false_type) {
			return riskFreeRate;
		}
	public:
		static auto const rate(std::size_t idx, RiskFreeRate const &riskFreeRate) {
			return _rate_impl(idx, riskFreeRate, std::is_compound<RiskFreeRate>());
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
	// ================================ probFloorCapper =============================
	// ==============================================================================

	/* If this function is used it should be logged so one can see when the probability 
	is not within <0.0,1.0> */
	
	template<typename T>
	T probFloorCapper(T x) {
		if (x < 0.0 || x > 1.0) {
			std::stringstream ss{};
			ss << "Probability value " << x << " has been modified to fit <0,1>\n";
			Logger::get().critical(std::cout, ss.str());
		}
		return std::min(1.0, std::max(0.0, x));
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
			if (x < x_[0]) {
				std::stringstream ss{};
				ss << "Lower constant extrapolation for x ( " << x << " ) occured.\n";
				Logger::get().warning(std::cout, ss.str());
				return y_[0];
			}
			if ((x > x_[x_.size() - 1])) {
				std::stringstream ss{};
				ss << "Upper constant extrapolation for x ( " << x << " ) occured.\n";
				Logger::get().warning(std::cout, ss.str());
				return y_[x_.size() - 1];
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