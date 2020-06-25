#pragma once
#if !defined(_LATTICE_MULTIDIMENSIONAL_T)
#define _LATTICE_MULTIDIMENSIONAL_T


#include"lattice_types.h"
#include"lattice_structure.h"
#include"lattice_multidimensional.h"


void testCreate2DIndexedLattice(){


	using lattice_multidimensional::MultidimIndexedLattice;
	using lattice_structure::IndexedLattice;
	using lattice_types::LatticeType;


	MultidimIndexedLattice<LatticeType::Binomial, double, IndexedLattice, IndexedLattice> index2Dtree;







}











#endif ///_LATTICE_MULTIDIMENSIONAL_T