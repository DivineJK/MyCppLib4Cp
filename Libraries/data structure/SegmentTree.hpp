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
	function<T(const T&, const T&)> func;
	T unit;
	vector<T> tree;
private:
	static constexpr uint32_t getBinMin(uint32_t n) {
		uint32_t ret = 1;
		while (ret < n) { ret <<= 1; }
		return ret;
	}
	void setArray(const vector<T>& vec) {
		uint32_t n = (int)vec.size();
		uint32_t b = getBinMin(n);
		tree = vector<T>(b << 1, unit);
		for (uint32_t i = 0; i < n; i++) { tree[i + b] = vec[i]; }
		for (uint32_t i = b; i > 0; i--) {
			tree[i - 1] = func(tree[(i - 1) << 1], tree[((i - 1) << 1) | 1]);
		}
	}
	void extendBottom(uint32_t s) {
		uint32_t b = size() >> 1;
		if (s <= b) { return; }
		uint32_t nb = getBinMin(s);
		tree.resize(nb << 1);
		for (uint32_t i = 0; i < b; i++) { tree[i + nb] = tree[i + b]; }
		for (uint32_t i = b; i < nb; i++) { tree[i + nb] = unit; }
		for (uint32_t i = nb; i > 0; i--) {
			tree[i - 1] = func(tree[(i - 1) << 1], tree[((i - 1) << 1) | 1]);
		}
	}
public:
	SegmentTree(uint32_t n, const T& aUnit, const function<T(const T&, const T&)>& aFunc) {
		unit = aUnit;
		func = aFunc;
		tree = vector<T>(getBinMin(n) << 1, unit);
	}
	SegmentTree(const vector<T>& vec, const T& aUnit, const function<T(const T&, const T&)>& aFunc) {
		unit = aUnit;
		func = aFunc;
		setArray(vec);
	}
	SegmentTree(const SegmentTree& segtree) { operator=(segtree); }
	SegmentTree& operator=(const SegmentTree& segtree) {
		unit = segtree.unit;
		func = segtree.func;
		tree = segtree.tree;
	}
	uint32_t size() const { return (uint32_t)tree.size(); }
	void setValue(uint32_t i, const T& v) {
		uint32_t b = size() >> 1;
		if (i >= b) { extendBottom(i + 1); }
		b = size() >> 1;
		tree[i + b] = v;
		uint32_t idx = i + b;
		while (idx) {
			idx >>= 1;
			tree[idx] = func(tree[idx << 1], tree[(idx << 1) | 1]);
		}
	}
	T getValue(uint32_t i) const {
		uint32_t b = size() >> 1;
		if (i >= b) { return unit; }
		return tree[i + b];
	}
	T getSegment(uint32_t l, uint32_t r) const {
		if (l >= r) { return unit; }
		uint32_t b = size() >> 1;
		if (l >= b) { return unit; }
		uint32_t rl = l + b;
		uint32_t rr = r + b;
		T retl = unit, retr = unit;
		while (rl < rr) {
			if (rl & 1) { retl = func(retl, tree[rl++]); }
			if (rr & 1) { retr = func(tree[--rr], retr); }
			rl >>= 1;
			rr >>= 1;
		}
		return func(retl, retr);
	}
	T getSum(uint32_t n) const { return getSegment(0, n); }
};

#endif /* SegmentTree_hpp */
