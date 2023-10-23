//
//  UnionFind.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef UnionFind_hpp
#define UnionFind_hpp

class UnionFind {
private:
	vector<uint32_t> roots;
	vector<uint32_t> sizes;
	uint32_t groupCount = 0;
public:
	UnionFind(uint32_t n) {
		roots.resize(n);
		sizes.resize(n);
		groupCount = n;
		for (uint32_t i = 0; i < n; i++) {
			roots[i] = i;
			sizes[i] = 1;
		}
	}
	uint32_t getVerticeCount() const { return (uint32_t)roots.size(); }
	uint32_t getRoot(uint32_t x) {
		assert(0 <= x && x < getVerticeCount());
		while (x != roots[x]) { x = roots[x]; }
		return x;
	}
	bool isSame(uint32_t x, uint32_t y) { return getRoot(x) == getRoot(y); }
	uint32_t getGroupCount() const { return groupCount; }
	void unite(uint32_t x, uint32_t y) {
		uint32_t n = getVerticeCount();
		assert(0 <= x && x < n && 0 <= y && y < n);
		x = getRoot(x);
		y = getRoot(y);
		if (x == y) { return; }
		if (sizes[x] < sizes[y]) {
			uint32_t t = x;
			x = y;
			y = t;
		}
		sizes[x] += sizes[y];
		roots[y] = x;
		groupCount--;
	}
};

#endif /* UnionFind_hpp */
