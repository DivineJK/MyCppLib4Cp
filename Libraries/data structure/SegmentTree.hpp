//
//  SegmentTree.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef SegmentTree_hpp
#define SegmentTree_hpp

template <typename T>
class SegmentTree {
private:
	function<T(T, T)> func;
	vector<T> tree;
private:
	static constexpr uint32_t getBinMin(uint32_t n) {
		uint32_t ret = 1;
		while (ret < n) { ret <<= 1; }
		return ret;
	}
public:
	SegmentTree(uint32_t n, const T& idv, const function<T(T, T)>& aFunc) {
		tree = vector<T>(getBinMin(n) << 1, idv);
		func = aFunc;
	}
	SegmentTree(const vector<T>& vec, const T& idv, const function<T(T, T)>& aFunc) {
		uint32_t n = (int)vec.size();
		uint32_t b = getBinMin(n);
		tree = vector<T>(b << 1, idv);
		func = aFunc;
	}
};

#endif /* SegmentTree_hpp */
