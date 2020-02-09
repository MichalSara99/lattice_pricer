#pragma once
#if !defined(_LATTICE_MISCELLANEOUS)
#define _LATTICE_MISCELLANEOUS

#include<type_traits>
#include<vector>

namespace lattice_miscellaneous {

	template<typename T = std::enable_if<std::is_arithmetic<T>::value>::type>
	struct OptionData {
		T Underlying;
		T Strike;
		T Volatility;
		T RiskFreeRate;
		T DividentRate;
	};


}




#endif //_LATTICE_MISCELLANEOUS