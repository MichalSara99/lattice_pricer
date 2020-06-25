#pragma once
#if !defined(_LATTICE_MULTIDIMENSIONAL)
#define _LATTICE_MULTIDIMENSIONAL

#include"lattice_types.h"
#include<array>

namespace lattice_multidimensional {

	using lattice_types::LatticeType;

	// ==============================================================================
	// ============================= MultidimGeneralLattice =========================
	// ==============================================================================

	template<LatticeType Type,
			typename Node,
			typename TimeAxis,
			typename NodeContainerType,
			template<LatticeType,typename...> typename ...LatticeObjects>
	class MultidimGeneralLattice {
	
	
	};


	// ==============================================================================
	// ============================= MultidimIndexedLattice =========================
	// ==============================================================================


	template<LatticeType Type,
			typename Node,
			template<LatticeType,typename> typename ...IndexLattices>
	class MultidimIndexedLattice :public MultidimGeneralLattice<Type, Node, std::size_t, std::vector<Node>,IndexLattices ...> {



	};


	// ==============================================================================
	// =================================== MultidimLattice ==========================
	// ==============================================================================
	 
	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		template<LatticeType, typename, typename> typename ...Lattices>
	class MultidimLattice :public MultidimGeneralLattice<Type, Node, TimeAxis, std::vector<Node>,Lattices ...> {



	};




}
















#endif ///_LATTICE_MULTIDIMENSIONAL