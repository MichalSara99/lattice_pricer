#include<iostream>
#include<string>
#include<set>

#include<boost/date_time/gregorian/gregorian.hpp>
#include<boost/date_time/gregorian/gregorian_types.hpp>

#include"lattice_structure.h"
#include"lattice_algorithms_t.h"
#include"lattice_utility_t.h"

using namespace boost::gregorian;



int main(int argc, char const *argv[]) {


	std::cout << "===============================================\n";

	utilityTest1();
	utilityTest2();

	
	std::cout << "\n";
	std::cin.get();
	std::cin.get();
	return 0;
}