#pragma once
#if !defined(_LATTICE_STRUCTURE)
#define _LATTICE_STRUCTURE

#include<vector>
#include<map>
#include<set>
#include<type_traits>
#include<string>
#include<initializer_list>
#include"lattice_types.h"

namespace lattice_structure {

	using lattice_types::LatticeType;

	template<LatticeType Type,
			typename Node,
			typename TimeAxis,
			typename NodeContainerType>
	class GeneralLattice {
	private:
		template<typename Arg>
		void buildTree_impl(Arg const &arg, std::true_type);

		template<typename Arg>
		void buildTree_impl(Arg const &arg, std::false_type);

		Node apex_impl(std::true_type)const;

		Node apex_impl(std::false_type)const;

	public:
		typedef typename std::conditional<std::is_integral<TimeAxis>::value,
								std::vector<NodeContainerType>,
								std::map<TimeAxis,NodeContainerType>>::type TreeType;

		typedef typename TreeType::iterator Iterator_type;

	protected:
		TreeType tree_;

		template<typename Arg>
		void buildTree(Arg const &arg);

	public:
		TreeType const &tree()const { return this->tree_; }

		Node const &operator()(TimeAxis timeIdx, std::size_t leafIdx)const {
			return this->tree_[timeIdx][leafIdx];
		}

		Node &operator()(TimeAxis timeIdx, std::size_t leafIdx) {
			return this->tree_[timeIdx][leafIdx];
		}
			
		NodeContainerType const &nodesAt(TimeAxis timeIdx)const {
			return this->tree_[timeIdx];
		}

		NodeContainerType &nodesAt(TimeAxis timeIdx) {
			return this->tree_[timeIdx];
		}

		Node apex()const;

		Iterator_type begin(){
			return this->tree_.begin();
		}

		Iterator_type end() {
			return this->tree_.end();
		}
	};


	template<LatticeType Type,
			typename Node,
			typename TimeAxis,
			typename NodeContainerType>
	template<typename Arg>
	void GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::buildTree_impl(Arg const &arg, std::true_type) {
		auto numberNodes = [](std::size_t nodeIdx) {
			auto intRep = static_cast<std::underlying_type<LatticeType>::type>(Type);
			return ((intRep + 1)*(nodeIdx + 1) - (intRep));
		};
		for (std::size_t t = 0; t < arg.size(); ++t) {
			tree_[arg[t]] = std::move(NodeContainerType(numberNodes(t)));
		}
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
	template<typename Arg>
	void GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::buildTree_impl(Arg const &arg, std::false_type) {
		auto numberNodes = [](std::size_t nodeIdx) {
			auto intRep = static_cast<std::underlying_type<LatticeType>::type>(Type);
			return ((intRep + 1)*(nodeIdx + 1) - (intRep));
		};
		tree_.reserve(arg);
		for (std::size_t t = 0; t <= arg; ++t) {
			tree_.emplace_back(NodeContainerType(numberNodes(t)));
		}
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
	template<typename Arg>
	void GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::buildTree(Arg const &arg) {
		buildTree_impl(arg,std::is_compound<Arg>());
	}


	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		Node GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::apex_impl(std::false_type) const{
		return *((*this->tree_.begin()).second.begin());
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		Node GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::apex_impl(std::true_type) const{
		return *((*this->tree_.begin()).begin());
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		Node GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::apex()const {
		return apex_impl(std::is_integral<TimeAxis>());
	}


	// IndexedLattice:
	template<LatticeType Type,
		typename Node>
	class IndexedLattice :public GeneralLattice<Type, Node, std::size_t, std::vector<Node>> {
	private:
		std::size_t minIndex_{ 0 };
		std::size_t maxIndex_;
	public:
		IndexedLattice(std::size_t numberPeriods)
			:maxIndex_{ numberPeriods } {
			this->buildTree(numberPeriods);
		}
		virtual ~IndexedLattice(){}


		IndexedLattice(IndexedLattice<Type, Node> const &other)
			:minIndex_{ other.minIndex_ },
			maxIndex_{other.maxIndex_} {
			this->tree_ = other.tree_;
		}

		IndexedLattice(IndexedLattice<Type, Node> &&other)noexcept
			:minIndex_{ std::move(other.minIndex_) },
			maxIndex_{ std::move(other.maxIndex_) } {
			this->tree_ = std::move(other.tree_);
		}

		IndexedLattice &operator=(IndexedLattice<Type, Node> const &other) {
			if (this != &other) {
				minIndex_ = other.minIndex_;
				maxIndex_ = other.maxIndex_;
				this->tree_ = other.tree_;
			}
			return *this;
		}

		IndexedLattice &operator=(IndexedLattice<Type, Node> &&other) noexcept{
			if (this != &other) {
				minIndex_ = std::move(other.minIndex_);
				maxIndex_ = std::move(other.maxIndex_);
				this->tree_ = std::move(other.tree_);
			}
			return *this;
		}

		std::size_t minIndex()const { return minIndex_; }
		std::size_t maxIndex()const { return maxIndex_; }
	};


	template<LatticeType Type,
			typename Node,
			typename TimeAxis>
	class Lattice :public GeneralLattice<Type, Node, TimeAxis, std::vector<Node>> {
	private:
		std::set<TimeAxis> fixingDatesSet_;
		std::vector<TimeAxis> fixingDates_;
	
	public:
		Lattice(std::set<TimeAxis> const &fixingDatesSet)
			:fixingDatesSet_{ fixingDatesSet } {
			fixingDates_.reserve(fixingDatesSet_.size());
			for (auto e : fixingDatesSet_) {
				fixingDates_.emplace_back(std::move(e));
			}
			this->buildTree(fixingDates_);
		}

		Lattice(std::initializer_list<TimeAxis> const &fixingDates)
			:fixingDatesSet_{fixingDates} {
			fixingDates_.reserve(fixingDatesSet_.size());
			for (auto e : fixingDatesSet_) {
				fixingDates_.emplace_back(std::move(e));
			}
			this->buildTree(fixingDates_);
		}

		virtual ~Lattice() {}

		Lattice(Lattice<Type, Node, TimeAxis> const &other)
			:fixingDates_{other.fixingDates_} {
			this->tree_ = other.tree_;
		}

		Lattice(Lattice<Type, Node, TimeAxis> &&other)noexcept 
			:fixingDates_{std::move(other.fixingDates_)} {
			this->tree_ = std::move(this->tree_);
		}

		Lattice &operator=(Lattice<Type, Node, TimeAxis> const &other) {
			if (this != &other) {
				fixingDates_ = other.fixingDates_;
				this->tree = other.tree_;
			}
			return *this;
		}

		Lattice &operator=(Lattice<Type, Node, TimeAxis> &&other)noexcept {
			if (this != &other) {
				fixingDates_ = std::move(other.fixingDates_);
				this->tree = std::move(other.tree_);
			}
			return *this;
		}

		std::vector<TimeAxis> fixingDates()const { return fixingDates_; }


	};


}




#endif //_LATTICE_STRUCTURE