#pragma once
#if !defined(_LATTICE_TRAITS)
#define _LATTICE_TRAITS

#include<tuple>
#include<array>
#include<sstream>
#include"lattice_miscellaneous.h"
#include"lattice_types.h"

namespace lattice_traits {

	// ===================================================================
	// =========================== MergeTraits ===========================
	// ===================================================================

	template<std::size_t N,typename Node>
	struct MergeTraits{
		using NodeHolder = std::array<Node, N>;

		template<typename T,typename ...Ts>
		static std::array<Node, N> holder(T t,Ts... ts) {
			return std::array<Node, N>{t,ts...};
		}
	};

	template<typename Node>
	struct MergeTraits<1, Node> {
		using NodeHolder = std::tuple<Node, Node>;

		template<typename T,typename ...Ts>
		static std::tuple<Node,Node> holder(T t,Ts... ts) {
			return std::make_tuple(t, ts...);
		}
	};
	template<typename Node>
	struct MergeTraits<2, Node> {
		using NodeHolder = std::tuple<Node, Node, Node>;

		template<typename T, typename ...Ts>
		static std::tuple<Node,Node, Node> holder(T t, Ts... ts) {
			return std::make_tuple(t, ts...);
		}
	};
	template<typename Node>
	struct MergeTraits<3, Node> {
		using NodeHolder = std::tuple<Node, Node, Node, Node>;

		template<typename T, typename ...Ts>
		static std::tuple<Node, Node,Node, Node> holder(T t, Ts... ts) {
			return std::make_tuple(t, ts...);
		}
	};



	// ===================================================================
	// =========================== PrintTraits ===========================
	// ===================================================================
	

	using lattice_miscellaneous::clipComma;
	using lattice_types::LatticeType;

	template<typename Node,LatticeType Type>
	struct PrintTraits{
		static const std::size_t size = 1;

		template<typename Container>
		static std::string printLine(Container &cont) {
			std::ostringstream ss;
			for (auto e : cont) {
				ss << e << ", ";
			}
			return clipComma(std::move(ss.str()));
		}
	};

	template<typename Node>
	struct PrintTraits<Node, LatticeType::TwoVariableBinomial> {
		static const std::size_t size = 1;

		template<typename Container>
		static std::string printLine(Container &cont) {
			std::ostringstream ss;
			std::size_t const side = static_cast<std::size_t>(std::sqrt(cont.size()));
			for (std::size_t i = 0; i < cont.size(); ++i) {
				if ((i > 0) && (i%side == 0))
					ss << "\n ";
				ss << cont[i] << ", ";
			}
			return clipComma(std::move(ss.str()));
		}
	};


	template<typename Node,LatticeType Type>
	struct PrintTraits<std::tuple<Node, Node>, Type> {
		static const std::size_t size = 2;

		template<typename Container>
		static std::string printLine(Container &cont) {
			std::ostringstream ss;
			for (auto e : cont) {
				ss << "(" << std::get<0>(e) << ","
					<< std::get<1>(e) << "), ";
			}
			return clipComma(std::move(ss.str()));
		}
	};

	template<typename Node>
	struct PrintTraits<std::tuple<Node, Node>, LatticeType::TwoVariableBinomial> {
		static const std::size_t size = 1;

		template<typename Container>
		static std::string printLine(Container &cont) {
			std::ostringstream ss;
			std::size_t const side = static_cast<std::size_t>(std::sqrt(cont.size()));
			std::tuple<Node, Node> val{};
			for (std::size_t i = 0; i < cont.size(); ++i) {
				if ((i > 0) && (i%side == 0))
					ss << "\n ";
				val = cont[i];
				ss << "(" << std::get<0>(val) << ","
					<< std::get<1>(val) << "), ";
			}
			return clipComma(std::move(ss.str()));
		}
	};

	template<typename Node,LatticeType Type>
	struct PrintTraits<std::tuple<Node,Node, Node>, Type> {
		static const std::size_t size = 3;

		template<typename Container>
		static std::string printLine(Container &cont) {
			std::ostringstream ss;
			for (auto e : cont) {
				ss << "(" << std::get<0>(e) << ","
					<< std::get<1>(e) << ","
					<< std::get<2>(e) << "), ";
			}
			return clipComma(std::move(ss.str()));
		}
	};

	template<typename Node>
	struct PrintTraits<std::tuple<Node, Node, Node>, LatticeType::TwoVariableBinomial> {
		static const std::size_t size = 1;

		template<typename Container>
		static std::string printLine(Container &cont) {
			std::ostringstream ss;
			std::size_t const side = static_cast<std::size_t>(std::sqrt(cont.size()));
			std::tuple<Node, Node, Node> val{};
			for (std::size_t i = 0; i < cont.size(); ++i) {
				if ((i > 0) && (i%side == 0))
					ss << "\n ";
				val = cont[i];
				ss << "(" << std::get<0>(val) << ","
					<< std::get<1>(val) << ","
					<< std::get<2>(val) << "), ";
			}
			return clipComma(std::move(ss.str()));
		}
	};

	template<typename Node,LatticeType Type>
	struct PrintTraits<std::tuple<Node, Node, Node, Node>, Type> {
		static const std::size_t size = 4;

		template<typename Container>
		static std::string printLine(Container &cont) {
			std::ostringstream ss;
			for (auto e : cont) {
				ss << "(" << std::get<0>(e) << ","
					<< std::get<1>(e) << ","
					<< std::get<2>(e) << ","
					<< std::get<3>(e) << "), ";
			}
			return clipComma(std::move(ss.str()));
		}
	};

	template<typename Node>
	struct PrintTraits<std::tuple<Node, Node, Node, Node>, LatticeType::TwoVariableBinomial> {
		static const std::size_t size = 1;

		template<typename Container>
		static std::string printLine(Container &cont) {
			std::ostringstream ss;
			std::size_t const side = static_cast<std::size_t>(std::sqrt(cont.size()));
			std::tuple<Node, Node, Node, Node> val{};
			for (std::size_t i = 0; i < cont.size(); ++i) {
				if ((i > 0) && (i%side == 0))
					ss << "\n ";
				val = cont[i];
				ss << "(" << std::get<0>(val) << ","
					<< std::get<1>(val) << ","
					<< std::get<2>(val) << ","
					<< std::get<3>(val) << "), ";
			}
			return clipComma(std::move(ss.str()));
		}
	};

	template<typename Node,LatticeType Type,std::size_t N>
	struct PrintTraits<std::array<Node,N>, Type> {
		static const std::size_t size = N;

		template<typename Container>
		static std::string printLine(Container &cont) {
			std::ostringstream ss;
			for (auto e : cont) {
				ss << "(";
				for (auto i = 0; i < N; ++i) {
					ss << e[i] << ",";
				}
				ss << "), ";
			}
			return clipComma(std::move(ss.str()));
		}
	};

	template<typename Node, std::size_t N>
	struct PrintTraits<std::array<Node, N>, LatticeType::TwoVariableBinomial> {
		static const std::size_t size = N;

		template<typename Container>
		static std::string printLine(Container &cont) {
			std::ostringstream ss;
			std::size_t const side = static_cast<std::size_t>(std::sqrt(cont.size()));
			std::array<Node, N> val{};
			for (std::size_t i = 0; i < cont.size(); ++i) {
				if ((i > 0) && (i%side == 0))
					ss << "\n ";
				val = cont[i];
				ss << "(";
				for (std::size_t j = 0; j < N; ++j) {
					ss << val[j] << ",";
				}
				ss << "), ";
			}
			return clipComma(std::move(ss.str()));
		}
	};


}




#endif ///_LATTICE_TRAITS