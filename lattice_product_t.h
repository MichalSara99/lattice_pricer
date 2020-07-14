#pragma once
#if !defined(_LATTICE_PRODUCT_H)
#define _LATTICE_PRODUCT_H


#include"lattice_product.h"
#include"lattice_product_builder.h"



void testOptionBuilder() {
	using lattice_product::Option;
	using lattice_product_builder::OptionBuilder;

	auto option = OptionBuilder<double>()
		.withName("plain option")
		.withDividend(12)
		.withMaturity(1.0)
		.withSpot(60.0)
		.build();

	auto params = option.modelParams();

	std::cout << "option builder:\n";
	std::cout << "option div: " << option.dividend() << "\n";
	std::cout << "option name: " << option.name() << "\n";
	std::cout << "option maturity: " << option.maturity() << "\n";
}

void testBarrierOptionBuilder() {
	using lattice_product::BarrierOption;
	using lattice_product_builder::BarrierOptionBuilder;

	auto option = BarrierOptionBuilder<double>()
		.withName("plain barrier option")
		.withDividend(12)
		.withMaturity(1.0)
		.build();

	auto params = option.modelParams();

	std::cout << "option builder:\n";
	std::cout << "option div: " << option.dividend() << "\n";
	std::cout << "option name: " << option.name() << "\n";
	std::cout << "option maturity: " << option.maturity() << "\n";

}


void testPureDiscountBondBuilder() {
	using lattice_product::PureDiscountBond;
	using lattice_product_builder::PureDiscountBondBuilder;

	auto option = PureDiscountBondBuilder<double>()
		.withName("pure discount option")
		.withNominal(100)
		.build();

	std::cout << "option builder:\n";
	std::cout << "option nominal: " << option.nominal() << "\n";
	std::cout << "option name: " << option.name() << "\n";
}



#endif ///_LATTICE_PRODUCT_H