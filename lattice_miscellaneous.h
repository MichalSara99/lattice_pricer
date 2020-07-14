#pragma once
#if !defined(_LATTICE_MISCELLANEOUS)
#define _LATTICE_MISCELLANEOUS

#include<type_traits>
#include<vector>
#include<thread>

namespace lattice_miscellaneous {

	static std::string const clipComma(std::string str) {
		const std::string lstr{ str };
		const std::size_t sz = lstr.size();
		return str.substr(0, sz - 2);
	}


	template<typename T>
	T sign(T x) {
		if (x < 0)
			return -1;
		else if (x > 0)
			return 1;
		return 0;
	}

	template<typename T>
	T lerp(T y0, T y1, T t) {
		return std::fma(t, (y1 - y0), y0);
	}

	/* TO BE DELETED SOON */
	template<typename T,
			typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	struct OptionData {
		T Underlying;
		T Strike;
		T Volatility;
		T RiskFreeRate;
		T DividentRate;
		T ReversionSpeed;
		//T Correlation;
	};
	/* TO BE DELETED SOON */

	template<typename T>
		struct MeanRevertingParams {
		T ReversionSpeed;
	};



	class scoped_thread {
	private:
		std::thread t_;
	public:
		explicit scoped_thread(std::thread thread)
			:t_{std::move(thread)}{
			if (!t_.joinable()) {
				throw std::logic_error("No thread\n");
			}
		}

		~scoped_thread() {
			t_.join();
		}

		scoped_thread(scoped_thread const &other) = delete;

		scoped_thread &operator=(scoped_thread const &other) = delete;

	};



	namespace map_traits {

		template<typename ...>
		struct voider { using type = void; };

		template<typename ...T>
		using void_t = typename voider<T...>::type;

		template<typename T, typename U = void>
		struct is_map_impl :std::false_type {};

		template<typename T>
		struct is_map_impl<T, void_t<typename T::key_type,
			typename T::mapped_type,
			decltype(std::declval<T&>()[std::declval<const typename T::key_type&>()])>>
			:std::true_type{};
	}

	template<typename T>
	struct is_map :map_traits::is_map_impl<T>::type {};



}




#endif //_LATTICE_MISCELLANEOUS