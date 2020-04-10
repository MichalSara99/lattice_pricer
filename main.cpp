#include<iostream>
#include<string>
#include<set>

#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>

//#include"lattice_structure.h"
#include"lattice_algorithms_t.h"
//#include"lattice_utility_t.h"
#include"lattice_model_t.h"
#include"lattice_examples.h"
//#include"lattice_greeks_t.h"

using namespace boost::gregorian;



int main(int argc, char const *argv[]) {

	auto today = date(day_clock::local_day());
	auto dd = date_duration{ 2 };
	auto today2 = today + dd;
	std::cout << "today: " << today << "\n";
	std::cout << "today + 2: " << today2 << "\n";
	std::cout << "duration: " << (today2 - today)<< "\n";
	


	std::cout << "===============================================\n";

	//crrBinomialLatticeParallelPricing();
	//mcrrBinomialLatticeParallelPricingScoped();
	//jrBinomialLatticeParallelPricingScoped();
	//trimBinomialLatticeParallelPricingScoped();
	//tmBinomialLatticeParallelPricingScoped();
	//lrBinomialLatticeParallelPricingScoped();
	//pascalTriangleIndexedTest();
	//pascalTriangleTest();

	//mergeIndexedBinomial();
	//mergeBinomial();


	//crrIndexedLattice();
	//mcrrIndexedLattice();
	//jrIndexedLattice();
	//trimIndexedLattice();
	//tmIndexedLattice();
	//lrIndexedLattice();
	//bmIndexedLattice();

	//crrLattice();
	//mcrrLattice();
	//jrLattice();
	//trimLattice();
	//tmLattice();
	//lrLattice();
	//bmLattice();

	//indexedLatticeBinomialDelta();
	//crrIndexedLatticeGreeks();
	//mcrrIndexedLatticeGreeks();
	//jrIndexedLatticeGreeks();
	//trimIndexedLatticeGreeks();
	//tmIndexedLatticeGreeks();
	//lrIndexedLatticeGreeks();
	//bmIndexedLatticeGreeks();

	//indexedLatticeTrinomialForwardInduction();
	//indexedLatticeTrinomialForwardInductionNew();
	//latticeBinomialForwardInduction();
	//latticeBinomialForwardInductionNew();

	//latticeTrinomialForwardInduction();
	//latticeTrinomialForwardInductionNew();
	//indexedLatticeBinomialForwardInductionDividend();
	//indexedLatticeBinomialForwardInductionDividendNew();
	
	//indexedLatticeTrinomialForwardInductionDividends();
	//indexedLatticeTrinomialForwardInductionDividendsNew();

	//latticeBinomialForwardInductionDividends();
	//latticeBinomialForwardInductionDividendsNew();


	//latticeTrinomialForwardInductionDividend();
	//latticeTrinomialForwardInductionDividendNew();
	

	//indexedLatticeBinomialBackwardInduction();
	//indexedLatticeBinomialBackwardInductionNew();

	//indexedLatticeTrinomialBackwardInduction();
	//indexedLatticeTrinomialBackwardInductionNew();

	//latticeBinomialBackwardInduction();
	//latticeBinomialBackwardInductionNew();


	//latticeTrinomialBackwardInduction();
	//latticeTrinomialBackwardInductionNew();

	//crrIndexedLattice();
	//crrIndexedLatticeNew();

	//mcrrIndexedLattice();
	//mcrrIndexedLatticeNew();

	//jrIndexedLattice();
	//jrIndexedLatticeNew();

	//trimIndexedLattice();
	//trimIndexedLatticeNew();

	//tmIndexedLattice();
	//tmIndexedLatticeNew();

	//lrIndexedLattice();
	//lrIndexedLatticeNew();

	//bmIndexedLattice();
	//bmIndexedLatticeNew();

	//crrIndexedLatticeAmerican();
	//crrIndexedLatticeAmericanNew();

	//mcrrIndexedLatticeAmerican();
	//mcrrIndexedLatticeAmericanNew();


	//jrIndexedLatticeAmerican();
	//jrIndexedLatticeAmericanNew();

	//trimIndexedLatticeAmerican();
	//trimIndexedLatticeAmericanNew();

	//tmIndexedLatticeAmerican();
	//tmIndexedLatticeAmericanNew();

	//lrIndexedLatticeAmerican();
	//lrIndexedLatticeAmericanNew();

	//crrLattice();
	//crrLatticeNew();

	//mcrrLattice();
	//mcrrLatticeNew();

	//jrLattice();
	//jrLatticeNew();

	//trimLattice();
	//trimLatticeNew();

	//tmLattice();
	//tmLatticeNew();

	//lrLattice();
	//lrLatticeNew();

	//bmLattice();
	//bmLatticeNew();

	//crrLatticeAmerican();
	//crrLatticeAmericanNew();


	//crrBinomialLatticeParallelPricing();
	//crrBinomialLatticeParallelPricingNew();

	//crrBinomialLatticeParallelPricingScoped();
	//crrBinomialLatticeParallelPricingScopedNew();

	//mcrrBinomialLatticeParallelPricingScoped();
	//mcrrBinomialLatticeParallelPricingScopedNew();
	

	//jrBinomialLatticeParallelPricingScoped();
	//jrBinomialLatticeParallelPricingScopedNew();

	//trimBinomialLatticeParallelPricingScoped();
	//trimBinomialLatticeParallelPricingScopedNew();


	//tmBinomialLatticeParallelPricingScoped();
	//tmBinomialLatticeParallelPricingScopedNew();


	lrBinomialLatticeParallelPricingScoped();
	lrBinomialLatticeParallelPricingScopedNew();


	std::cout << "\n";
	std::cin.get();
	std::cin.get();
	return 0;
}