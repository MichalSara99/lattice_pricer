#pragma once
#if !defined(_LATTICE_TYPES)
#define _LATTICE_TYPES

#include<tuple>
#include<functional>

namespace lattice_types {


	enum class LatticeType { Binomial, Trinomial };

	template<typename Node,typename ...Nodes>
	using LeafForwardGenerator = std::function<std::tuple<Nodes...>(Node,Node)>;

	template<typename Node,typename ...Nodes>
	using LeafBackwardGenerator = std::function<Node(Nodes...)>;

	template<typename ReturnType,typename ...Args>
	using Payoff = std::function<ReturnType(Args...)>;

	template<typename ...Args>
	using PayoffAdjuster = std::function<void(Args...)>;

}




#endif ///_LATTICE_TYPES




