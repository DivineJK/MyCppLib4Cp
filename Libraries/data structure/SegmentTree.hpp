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
	T identity;
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
		tree = vector<T>(b << 1, identity);
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
		for (uint32_t i = b; i < nb; i++) { tree[i + nb] = identity; }
		for (uint32_t i = nb; i > 0; i--) {
			tree[i - 1] = func(tree[(i - 1) << 1], tree[((i - 1) << 1) | 1]);
		}
	}
public:
	SegmentTree(uint32_t n, const T& idv, const function<T(T, T)>& aFunc) {
		tree = vector<T>(getBinMin(n) << 1, idv);
		identity = idv;
		func = aFunc;
	}
	SegmentTree(const vector<T>& vec, const T& idv, const function<T(T, T)>& aFunc) {
		identity = idv;
		func = aFunc;
		setArray(vec);
	}
	SegmentTree(const SegmentTree& segtree) { operator=(segtree); }
	SegmentTree& operator=(const SegmentTree& segtree) {
		tree = segtree.tree;
		func = segtree.func;
		identity = segtree.identity;
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
		if (i >= b) { return identity; }
		return tree[i + b];
	}
	T getSegment(uint32_t l, uint32_t r) const {
		if (l >= r) { return identity; }
		uint32_t b = size() >> 1;
		if (l >= b) { return identity; }
		uint32_t rl = l + b;
		uint32_t rr = r + b;
		vector<uint32_t> ls, rs;
		while (rl < rr) {
			if (rl & 1) {
				ls.push_back(rl);
				rl++;
			}
			if (rr & 1) {
				rs.push_back(rr - 1);
				rr--;
			}
			rl >>= 1;
			rr >>= 1;
		}
		T ret = identity;
		vector<uint32_t>::const_iterator it;
		for (it = ls.cbegin(); it != ls.cend(); ++it) {
			uint32_t idx = *it;
			ret = func(ret, tree[idx]);
		}
		vector<uint32_t>::const_reverse_iterator rit;
		for (rit = rs.crbegin(); rit != rs.crend(); ++rit) {
			ret = func(ret, tree[*rit]);
		}
		return ret;
	}
	T getSum(uint32_t n) const {
		uint32_t b = size() >> 1;
		uint32_t idx = (n >= b) ? b << 1 : n + b;
		vector<uint32_t> ids;
		while (idx) {
			if (idx & 1) { ids.push_back(idx - 1); }
			idx >>= 1;
		}
		T ret = identity;
		vector<uint32_t>::const_reverse_iterator rit;
		for (rit = ids.crbegin(); rit != ids.crend(); ++rit) {
			ret = func(ret, tree[*rit]);
		}
		return ret;
	}
};

#endif /* SegmentTree_hpp */
