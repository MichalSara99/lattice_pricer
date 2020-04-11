#pragma once
#if !defined(_LATTICE_TRAITS)
#define _LATTICE_TRAITS

#include<tuple>

namespace lattice_traits {

	template<std::size_t N,typename Node>
	struct MergeTraits{};

	template<typename Node>
	struct MergeTraits<0, Node> {
		using NodeHolder = Node;
	};
	template<typename Node>
	struct MergeTraits<1, Node> {
		using NodeHolder = std::tuple<Node, Node>;
	};
	template<typename Node>
	struct MergeTraits<2, Node> {
		using NodeHolder = std::tuple<Node, Node, Node>;
	};
	template<typename Node>
	struct MergeTraits<3, Node> {
		using NodeHolder = std::tuple<Node, Node, Node, Node>;
	};
	template<typename Node>
	struct MergeTraits<4, Node> {
		using NodeHolder = std::tuple<Node, Node, Node, Node, Node>;
	};
	template<typename Node>
	struct MergeTraits<5, Node> {
		using NodeHolder = std::tuple<Node, Node, Node, Node, Node, Node>;
	};







}




#endif ///_LATTICE_TRAITS