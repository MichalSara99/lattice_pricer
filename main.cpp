#include<iostream>
#include<string>
#include<set>

#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>

#include"lattice_structure.h"
#include"lattice_algorithms_t.h"
#include"lattice_utility_t.h"
#include"lattice_model_t.h"
#include"lattice_examples.h"

using namespace boost::gregorian;



int main(int argc, char const *argv[]) {

	auto today = date(day_clock::local_day());
	auto dd = date_duration{ 2 };
	auto today2 = today + dd;
	std::cout << "today: " << today << "\n";
	std::cout << "today + 2: " << today2 << "\n";
	std::cout << "duration: " << (today2 - today)<< "\n";
	


	std::cout << "===============================================\n";

	crrBinomialLatticeParallelPricing();
	mcrrBinomialLatticeParallelPricingScoped();
	jrBinomialLatticeParallelPricingScoped();
	trimBinomialLatticeParallelPricingScoped();
	tmBinomialLatticeParallelPricingScoped();
	lrBinomialLatticeParallelPricingScoped();


	//mergeIndexedBinomial();
	//mergeBinomial();

	//crrLattice();
	//mcrrLattice();
	//jrLattice();
	//trimLattice();
	//tmLattice();
	//lrLattice();
	
	std::cout << "\n";
	std::cin.get();
	std::cin.get();
	return 0;
}