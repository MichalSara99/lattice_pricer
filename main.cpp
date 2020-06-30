#include<iostream>
#include<string>
#include<set>
#include<tuple>
#include<array>
#include<chrono>

#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>

#include"lattice_structure_t.h"
#include"lattice_algorithms_t.h"
#include"lattice_model_t.h"
#include"lattice_examples.h"
#include"lattice_greeks_t.h"
#include"lattice_utility_t.h"
#include"lattice_calibrator_ir_t.h"
#include"lattice_examples_calibrate_price_ir.h"

#include"lattice_calibrator_equity_t.h"
#include"lattice_multidimensional_t.h"
#include"lattice_builder_t.h"

using namespace boost::gregorian;



int main(int argc, char const *argv[]) {

	auto today = date(day_clock::local_day());
	auto dd = date_duration{ 2 };
	auto today2 = today + dd;
	std::cout << "today: " << today << "\n";
	std::cout << "today + 2: " << today2 << "\n";
	std::cout << "duration: " << (today2 - today)<< "\n";
	
	

	std::cout << "===============================================\n";

	// ===============================================
	// ======== lattice_calibrator_equity_t.h ========
	// ===============================================

	// testIndexedImpliedStatePriceLattice();
	// testImpliedStatePriceLattice();

	// testIndexedImpliedProbabilities();

	// testIndexedSPLImpliedBinomial();

	// ===============================================

	// ===============================================
	// === lattice_examples_calibrate_price_ir.h =====
	// ===============================================

	// testIndexedBinomialPureDiscountBondLattice();
	// testIndexedBinomialPureDiscountBondLatticeVariable(1);
	// testBinomialPureDicsountBondLattice();
	// testBinomialPureDicsountBondLatticeVariable(1);
	// testIndexedTrinomialPureDiscountBondLattice();
	// testIndexedTrinomialPureDiscountBondLatticeVariable(20);
	// testTrinomialPureDiscountBondLattice();
	// testTrinomialPureDiscountBondLatticeVariable(6);
	// testIndexedBinomialCouponBondLattice();
	// testIndexedBinomialCouponBondLatticeVariable(10);
	// testBinomialCouponBondLattice();
	// testBinomialCouponBondLatticeVariable(11);
	// testIndexedTrinomialCouponBondLattice();
	// testIndexedTrinomialCouponBondLatticeVariable(10);
	// testTrinomialCouponBondLattice();
	// testTrinomialCouponBondLatticeVariable(11);

	// testIndexedBinomialEuropeanOptionOnCouponBondLattice();
	// testBinomialEuropeanOptionOnCouponBondLattice();
	// testIndexedTrinomialEuropeanOptionOnCouponBondLattice();
	// testTrinomialEuropeanOptionOnCouponBondLattice();

	// ===============================================

	// ===================================
	// ==== lattice_calibrator_ir_t.h ====
	// ===================================

	// testIndexedBDT();
	// testBDT();
	// testIndexedHL();
	// testHL();
	// testIndexedHW();
	// testHW();
	// testIndexedBK();
	// testBK();

	// ===================================

	// ===============================
	// ===== lattice_greeks_t.h ======
	// ===============================

	// testIndexedLatticesGreeks();
	// testLatticesGreeks();

	// indexedLatticeBinomialDelta();
	// crrIndexedLatticeGreeks();
	// mcrrIndexedLatticeGreeks();
	// jrIndexedLatticeGreeks();
	// trimIndexedLatticeGreeks();
	// tmIndexedLatticeGreeks();
	// lrIndexedLatticeGreeks();
	// bmIndexedLatticeGreeks();

	// crrLatticeGreeks();
	// mcrrLatticeGreeks();
	// jrLatticeGreeks();
	// trimLatticeGreeks();
	// tmLatticeGreeks();
	// lrLatticeGreeks();
	// bmLatticeGreeks();

	// ===============================


	// ===============================
	// ==== lattice_structure_t.h ====
	// ===============================

	// testIndexedLatticeCreation();
	// testLatticeCreation();

	// ===============================


	// ===============================
	// === lattice_algorithms_t.h ====
	// ===============================

	// testIndexedForwardInduction();
	// testForwardInduction();

	// testIndexedForwardInductionDividends();
	// testForwardInductionDividends();

	// testIndexedBackwardInductionDividends();
	// testBackwardInductionDividends();

	// testIndexedMerge();
	// testMerge();

	// ===============================


	// ===============================
	// ====== lattice_model_t.h ======
	// ===============================

	// testIndexedEuropeanLattices();
	// testIndexedAmericanLattices();
	// testEuropeanLattices();
	// testAmericanLattices();

	// ===============================

	// =================================
	// ====== lattice_builder_t.h ======
	// =================================

	// testCreateIndexBasedLattice();
	// testCreateLattice();

	// =================================

	// ==========================================
	// ====== lattice_multidimensional_t.h ======
	// ==========================================

    // testCreateMultidimensionalLattice();
	// testCreateMultidimensionalMRLattice();

	// ==========================================

	// ==========================================
	// =========== lattice_utility_t.h ==========
	// ==========================================

	// utilityTest3();

	 // ==========================================


	std::cout << "\n";
	std::cin.get();
	std::cin.get();
	return 0;
}