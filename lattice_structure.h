#pragma once
#if !defined(_LATTICE_STRUCTURE)
#define _LATTICE_STRUCTURE

#include<vector>
#include<map>
#include<set>
#include<type_traits>
#include<string>
#include<numeric>
#include<initializer_list>
#include"lattice_types.h"
#include"lattice_miscellaneous.h"
#include"lattice_macros.h"

namespace lattice_structure {

	using lattice_types::LatticeType;
	using lattice_types::LatticeClass;
	using lattice_miscellaneous::MeanRevertingParams;
	using lattice_miscellaneous::is_map;


	// =======================================================================================
	// =============================== GeneralLattice ========================================
	// =======================================================================================

	template<LatticeType Type,
			typename Node,
			typename TimeAxis,
			typename NodeContainerType>
	class GeneralLattice {
	private:
		std::size_t firstRevertIdx_;

		template<typename Arg>
		void buildTree_impl(Arg const &arg, std::true_type);

		template<typename Arg>
		void buildTree_impl(Arg const &arg, std::false_type);

		template<typename Arg,typename DeltaTime>
		void buildRevertingTree_impl(Arg const &arg,MeanRevertingParams<Node> const &params,DeltaTime const &deltaTime, std::true_type);

		template<typename Arg,typename DeltaTime>
		void buildRevertingTree_impl(Arg const &arg, MeanRevertingParams<Node> const &params, DeltaTime const &deltaTime, std::false_type);

		std::size_t indexOf_impl(TimeAxis time, std::true_type)const;

		std::size_t indexOf_impl(TimeAxis time, std::false_type)const;

		Node apex_impl(std::true_type)const;

		Node apex_impl(std::false_type)const;

		Node const &at_impl(TimeAxis timeIdx, std::size_t leafIdx, std::true_type)const;

		Node const &at_impl(TimeAxis timeIdx, std::size_t leafIdx, std::false_type)const;

		Node const &operator_const_impl(std::size_t timeIdx, std::size_t leafIdx, std::true_type)const;

		Node const &operator_const_impl(std::size_t timeIdx, std::size_t leafIdx, std::false_type)const;

		Node &operator_impl(std::size_t timeIdx, std::size_t leafIdx, std::true_type);

		Node &operator_impl(std::size_t timeIdx, std::size_t leafIdx, std::false_type);

		NodeContainerType nodesAt_impl(TimeAxis timeIdx, std::true_type)const;

		NodeContainerType nodesAt_impl(TimeAxis timeIdx, std::false_type)const;

		NodeContainerType nodesAtIdx_impl(std::size_t idx, std::true_type)const;

		NodeContainerType nodesAtIdx_impl(std::size_t idx, std::false_type)const;

		TimeAxis const &timeAt_impl(std::size_t idx, std::true_type)const;

		TimeAxis const &timeAt_impl(std::size_t idx, std::false_type)const;

		template<typename DeltaTime>
		typename DeltaTime::value_type const getMinDeltaTime_impl(DeltaTime const &deltaTime, std::true_type)const {
			auto itr = std::min_element(deltaTime.cbegin(), deltaTime.cend());
			return *itr;
		}

		template<typename DeltaTime>
		DeltaTime const getMinDeltaTime_impl(DeltaTime const &deltaTime, std::false_type)const {
			return deltaTime;
		}

	public:
		typedef typename std::conditional<std::is_integral<TimeAxis>::value,
								std::vector<NodeContainerType>,
								std::map<TimeAxis,NodeContainerType>>::type TreeType;

		typedef typename TreeType::iterator Iterator_type;
		typedef typename TreeType::const_iterator Const_iterator_type;
		typedef Node Node_type;
		typedef TimeAxis TimeAxis_type;
		typedef NodeContainerType NodeContainerType_type;

	protected:
		TreeType tree_;

		void _firstRevertingIdx();

		template<typename DeltaTime>
		auto const getMinDeltaTime(DeltaTime const &deltaTime)const  {
			return getMinDeltaTime_impl(deltaTime, std::is_compound<DeltaTime>());
		}

		template<typename Arg>
		void buildTree(Arg const &arg);

		template<typename Arg,typename DeltaTime>
		void buildRevertingTree(Arg const &arg, MeanRevertingParams<Node> const &params, DeltaTime const &deltaTime);

		std::size_t numberNodes(std::size_t timeIdx)const;

	public:
		TreeType const &tree()const { return this->tree_; }

		static constexpr LatticeType type(){ return Type; }

		std::size_t constexpr factors()const { return 1; }

		constexpr std::size_t timeDimension()const { return std::distance(tree_.begin(), tree_.end()); }
		

		Node const &at(TimeAxis timeIdx, std::size_t leafIdx)const {
			return at_impl(timeIdx, leafIdx, is_map<TreeType>());
		}

		Node &at(TimeAxis timeIdx, std::size_t leafIdx) {
			return (this->tree_[timeIdx][leafIdx]);
		}

		Node const &operator()(std::size_t timeIdx, std::size_t leafIdx)const {
			return operator_const_impl(timeIdx, leafIdx, is_map<TreeType>());
		}

		Node &operator()(std::size_t timeIdx, std::size_t leafIdx) {
			return operator_impl(timeIdx, leafIdx, is_map<TreeType>());
		}
			
		NodeContainerType nodesAt(TimeAxis timeIdx)const {
			return nodesAt_impl(timeIdx, is_map<TreeType>());
		}

		NodeContainerType &nodesAt(TimeAxis timeIdx) {
			return this->tree_[timeIdx];
		}

		NodeContainerType nodesAtIdx(std::size_t idx)const {
			return nodesAtIdx_impl(idx, is_map<TreeType>());
		}

		std::size_t indexOf(TimeAxis time)const {
			return indexOf_impl(time, is_map<TreeType>());
		}

		TimeAxis const &timeAt(std::size_t idx)const {
			return timeAt_impl(idx, is_map<TreeType>());
		}

		Node apex()const;

		Const_iterator_type cbegin()const noexcept{
			return this->tree_.cbegin();
		}

		Const_iterator_type cend() const noexcept{
			return this->tree_.cend();
		}

		Iterator_type begin(){
			return this->tree_.begin();
		}

		Iterator_type end() {
			return this->tree_.end();
		}

		std::size_t const firstRevertingIdx()const { return this->firstRevertIdx_; }

		bool isFirstReverting(std::size_t timeIdx) { return (timeIdx <= 0) ?  false :  timeIdx == firstRevertIdx_; }
	};


	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
	void GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::_firstRevertingIdx(){
		std::vector<int> states;
		const std::size_t treeSize = timeDimension();
		for (std::size_t t = 0; t < treeSize; ++t) {
			states.push_back(nodesAtIdx(t).size());
		}
		std::adjacent_difference(states.begin(), states.end(), states.begin());
		auto zeroItr = std::find_if(states.begin(), states.end(), [](int val) {return (val == 0); });
		if (zeroItr == states.end()) {
			this->firstRevertIdx_ = 0;
		}
		else {
			this->firstRevertIdx_ = std::distance(states.begin(), zeroItr);
		}
	}


	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		NodeContainerType GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::nodesAtIdx_impl(std::size_t idx, std::true_type)const {
		LASSERT(idx >= 0, "Index must be nonegative");
		auto first = this->tree_.cbegin();
		if (idx == 0) {
			return first->second;
		}
		else {
			auto itr = std::next(first, idx);
			return itr->second;
		}
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		NodeContainerType GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::nodesAtIdx_impl(std::size_t idx, std::false_type)const {
		LASSERT(idx >= 0, "Index must be nonegative");
		return this->tree_[idx];
	}



	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		TimeAxis const &GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::timeAt_impl(std::size_t idx, std::true_type)const {
		LASSERT(idx >= 0, "Index must be nonegative");
		auto first = this->tree_.cbegin();
		if (idx == 0) {
			return first->first;
		}
		else {
			auto itr = std::next(first, idx);
			return itr->first;
		}
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		TimeAxis const &GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::timeAt_impl(std::size_t idx, std::false_type)const {
		return idx;
	}


	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		std::size_t GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::indexOf_impl(TimeAxis timeIdx, std::true_type)const {
		typename TreeType::const_iterator citer(this->tree_.find(timeIdx));
		return ((citer != this->tree_.cend()) ? (std::distance(this->tree_.cbegin(), citer) ) : throw std::out_of_range("Error: timeIdx out of range.\n"));
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		std::size_t GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::indexOf_impl(TimeAxis timeIdx, std::false_type)const {
		return timeIdx;
	}


	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		Node const &GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::at_impl(TimeAxis timeIdx, std::size_t leafIdx, std::true_type)const {
		typename TreeType::const_iterator citer(this->tree_.find(timeIdx));
		return ((citer != this->tree_.end()) ? (citer->second[leafIdx]) : throw std::out_of_range("Error: timeIdx out of range.\n"));
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		Node const &GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::at_impl(TimeAxis timeIdx, std::size_t leafIdx, std::false_type)const {
		return (this->tree_[timeIdx][leafIdx]);
	}

	template<LatticeType Type,
			typename Node,
			typename TimeAxis,
			typename NodeContainerType>
	Node &GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::operator_impl(std::size_t timeIdx, std::size_t leafIdx, std::true_type) {
		LASSERT(timeIdx >= 0, "Index must be nonegative");
		auto first = this->tree_.begin();
		if (timeIdx == 0) {
			return first->second[leafIdx];
		}
		else {
			auto itr = std::next(first, timeIdx);
			return itr->second[leafIdx];
		}
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
	Node &GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::operator_impl(std::size_t timeIdx, std::size_t leafIdx, std::false_type) {
		LASSERT(timeIdx >= 0, "Index must be nonegative");
		auto first = this->tree_.begin();
		if (timeIdx == 0) {
			return (*first)[leafIdx];
		}
		else {
			auto itr = std::next(first, timeIdx);
			return (*itr)[leafIdx];
		}
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		Node const &GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::operator_const_impl(std::size_t timeIdx, std::size_t leafIdx, std::true_type)const {
		LASSERT(timeIdx >= 0, "Index must be nonegative");
		auto first = this->tree_.cbegin();
		if (timeIdx == 0) {
			return first->second[leafIdx];
		}
		else {
			auto itr = std::next(first, timeIdx);
			return itr->second[leafIdx];
		}
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		Node const &GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::operator_const_impl(std::size_t timeIdx, std::size_t leafIdx, std::false_type)const {
		LASSERT(timeIdx >= 0, "Index must be nonegative");
		auto first = this->tree_.cbegin();
		if (timeIdx == 0) {
			return (*first)[leafIdx];
		}
		else {
			auto itr = std::next(first, timeIdx);
			return (*itr)[leafIdx];
		}
	}

	template<LatticeType Type,
		typename Node,
	typename TimeAxis,
	typename NodeContainerType>
	NodeContainerType GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::nodesAt_impl(TimeAxis timeIdx, std::true_type)const {
		typename TreeType::const_iterator citer(this->tree_.find(timeIdx));
		return ((citer != this->tree_.end()) ? (citer->second) : throw std::out_of_range("Error: timeIdx out of range.\n"));
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
	NodeContainerType GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::nodesAt_impl(TimeAxis timeIdx, std::false_type)const {
		return (this->tree_[timeIdx]);
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
	std::size_t GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::numberNodes(std::size_t timeIdx)const {
		if (Type == LatticeType::TwoVariableBinomial) {
			return ((timeIdx + 1)*(timeIdx + 1));
		}
		auto intRep = static_cast<std::underlying_type<LatticeType>::type>(Type);
		return ((intRep + 1)*(timeIdx + 1) - (intRep));
	};


	template<LatticeType Type,
			typename Node,
			typename TimeAxis,
			typename NodeContainerType>
	template<typename Arg>
	void GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::buildTree_impl(Arg const &arg, std::true_type) {
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
		template<typename Arg,typename DeltaTime>
	void GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::buildRevertingTree_impl(Arg const &arg,MeanRevertingParams<Node> const &params,DeltaTime const &deltaTime, std::true_type) {
		auto const dt = getMinDeltaTime(deltaTime);
		std::size_t const maxStatesUp = static_cast<std::size_t>(std::floor(1 / (2.0*params.ReversionSpeed*dt))) + 1;
		auto const intRep = static_cast<std::underlying_type<LatticeType>::type>(LatticeType::Trinomial);

		auto numberNodes = [&](std::size_t timeIdx) {
			std::size_t normalNodesSize = ((intRep + 1)*(timeIdx + 1) - (intRep));
			return std::min(normalNodesSize, 2 * maxStatesUp + 1);
		};
		for (std::size_t t = 0; t < arg.size(); ++t) {
			tree_[arg[t]] = std::move(NodeContainerType(numberNodes(t)));
		}
	}

	template<LatticeType Type,
		typename Node,
		typename TimeAxis,
		typename NodeContainerType>
		template<typename Arg, typename DeltaTime>
	void GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::buildRevertingTree_impl(Arg const &arg, MeanRevertingParams<Node> const &params, DeltaTime const &deltaTime, std::false_type) {
		
		auto const dt = getMinDeltaTime(deltaTime);
		std::size_t const maxStatesUp = static_cast<std::size_t>(std::floor(1 / (2.0*params.ReversionSpeed*dt))) + 1;
		auto const intRep = static_cast<std::underlying_type<LatticeType>::type>(LatticeType::Trinomial);
		
		auto numberNodes = [&](std::size_t timeIdx) {
			std::size_t normalNodesSize = ((intRep + 1)*(timeIdx + 1) - (intRep));
			return std::min(normalNodesSize, 2 * maxStatesUp + 1);
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
		template<typename Arg,typename DeltaTime>
	void GeneralLattice<Type, Node, TimeAxis, NodeContainerType>::buildRevertingTree(Arg const &arg,MeanRevertingParams<Node> const &params,DeltaTime const &deltaTime) {
		buildRevertingTree_impl(arg, params, deltaTime, std::is_compound<Arg>());
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


	// =======================================================================================
	// =============================== IndexedLattice ========================================
	// =======================================================================================


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

		static LatticeClass const latticeClass() { return LatticeClass::Normal; }

		std::size_t minIndex()const { return minIndex_; }
		std::size_t maxIndex()const { return maxIndex_; }
	};

	

	// =======================================================================================
	// =================================== Lattice ===========================================
	// =======================================================================================

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
			this->tree_ = std::move(other.tree_);
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

		static LatticeClass const latticeClass() { return LatticeClass::Normal; }

		std::vector<TimeAxis> fixingDates()const { return fixingDates_; }


	};


	// =======================================================================================
	// ========================= MeanRevertingIndexedLattice =================================
	// = Make sense only for trinomial trees so far											 =
	// =======================================================================================

	template<typename Node>
		class MeanRevertingIndexedLattice 
			:public GeneralLattice<LatticeType::Trinomial, Node, std::size_t, std::vector<Node>> {
		private:
			std::size_t minIndex_{ 0 };
			std::size_t maxIndex_;
		public:
			template<typename DeltaTime>
			MeanRevertingIndexedLattice(std::size_t numberPeriods, MeanRevertingParams<Node> const &params, DeltaTime const &deltaTime)
				:maxIndex_{ numberPeriods } {
				this->buildRevertingTree(numberPeriods, params, deltaTime);
				this->_firstRevertingIdx();
			}
			virtual ~MeanRevertingIndexedLattice() {}


			MeanRevertingIndexedLattice(MeanRevertingIndexedLattice<Node> const &other)
				:minIndex_{ other.minIndex_ },
				maxIndex_{ other.maxIndex_ } {
				this->tree_ = other.tree_;
			}

			MeanRevertingIndexedLattice(MeanRevertingIndexedLattice<Node> &&other)noexcept
				:minIndex_{ std::move(other.minIndex_) },
				maxIndex_{ std::move(other.maxIndex_) } {
				this->tree_ = std::move(other.tree_);
			}

			MeanRevertingIndexedLattice &operator=(MeanRevertingIndexedLattice<Node> const &other) {
				if (this != &other) {
					minIndex_ = other.minIndex_;
					maxIndex_ = other.maxIndex_;
					this->tree_ = other.tree_;
				}
				return *this;
			}

			MeanRevertingIndexedLattice &operator=(MeanRevertingIndexedLattice<Node> &&other) noexcept {
				if (this != &other) {
					minIndex_ = std::move(other.minIndex_);
					maxIndex_ = std::move(other.maxIndex_);
					this->tree_ = std::move(other.tree_);
				}
				return *this;
			}

			static LatticeClass const latticeClass() { return LatticeClass::MeanReverting; }

			std::size_t minIndex()const { return minIndex_; }
			std::size_t maxIndex()const { return maxIndex_; }
	};


	// =======================================================================================
	// ============================== MeanRevertingLattice ===================================
	// = Make sense only for trinomial trees so far											 =
	// =======================================================================================

	template<typename Node,
		typename TimeAxis>
		class MeanRevertingLattice 
			:public GeneralLattice<LatticeType::Trinomial, Node, TimeAxis, std::vector<Node>> {

		private:
			std::set<TimeAxis> fixingDatesSet_;
			std::vector<TimeAxis> fixingDates_;

		public:
			template<typename DeltaTime>
			MeanRevertingLattice(std::set<TimeAxis> const &fixingDatesSet, MeanRevertingParams<Node> const &params, DeltaTime const &deltaTime)
				:fixingDatesSet_{ fixingDatesSet } {
				fixingDates_.reserve(fixingDatesSet_.size());
				for (auto e : fixingDatesSet_) {
					fixingDates_.emplace_back(std::move(e));
				}
				this->buildRevertingTree(fixingDates_, params, deltaTime);
				this->_firstRevertingIdx();
			}

			virtual ~MeanRevertingLattice() {}

			MeanRevertingLattice(MeanRevertingLattice<Node, TimeAxis> const &other)
				:fixingDates_{ other.fixingDates_ } {
				this->tree_ = other.tree_;
			}

			MeanRevertingLattice(MeanRevertingLattice<Node, TimeAxis> &&other)noexcept
				:fixingDates_{ std::move(other.fixingDates_) } {
				this->tree_ = std::move(other.tree_);
			}

			MeanRevertingLattice &operator=(MeanRevertingLattice<Node, TimeAxis> const &other) {
				if (this != &other) {
					fixingDates_ = other.fixingDates_;
					this->tree = other.tree_;
				}
				return *this;
			}

			MeanRevertingLattice &operator=(MeanRevertingLattice<Node, TimeAxis> &&other)noexcept {
				if (this != &other) {
					fixingDates_ = std::move(other.fixingDates_);
					this->tree = std::move(other.tree_);
				}
				return *this;
			}

			static LatticeClass const latticeClass() { return LatticeClass::MeanReverting; }

			std::vector<TimeAxis> fixingDates()const { return fixingDates_; }
	};


}




#endif //_LATTICE_STRUCTURE