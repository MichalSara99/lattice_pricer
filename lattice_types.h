#pragma once
#if !defined(_LATTICE_TYPES)
#define _LATTICE_TYPES

#include<tuple>
#include<functional>

namespace lattice_types {


	enum class LaunchPolicy { Sequential, Parallel };

	enum class LatticeType { Binomial, Trinomial, TwoVariableBinomial };

	enum class AssetClass { InterestRate, Equity };

	enum class DiscountingStyle { Continuous, Discrete };

	enum class MinimizerMethod { Analytic, Numeric };

	enum class LatticeClass { Normal, MeanReverting };

	enum class BranchingStyle { Normal, Reverting };


	template<typename Node,typename ...Nodes>
	using LeafForwardGenerator = std::function<std::tuple<Nodes...>(Node,Node,std::size_t, std::size_t,bool)>;

	template<typename Node,typename ...Nodes>
	using LeafBackwardGenerator = std::function<Node(Nodes...)>;

}




#endif ///_LATTICE_TYPES




