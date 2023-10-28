//
//  BinaryIndexedTree.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef BinaryIndexedTree_hpp
#define BinaryIndexedTree_hpp

template <typename T>
class BinaryIndexedTree {
private:
	function<T(T, T)> func;
	function<T(T)> inv;
	T identity = 0;
	vector<T> arr;
private:
	void setArray(const vector<T>& aArr) {
		uint32_t n = (uint32_t)aArr.size();
		arr = vector<T>(n + 1, identity);
		for (uint32_t i = 0; i < n; i++) {
			uint32_t p = i + 1;
			while (p <= n) {
				arr[p] = func(arr[p], aArr[i]);
				p += p ^ (p & (p - 1));
			}
		}
	}
	void extend(uint32_t newLen) {
		uint32_t oldLen = (uint32_t)arr.size();
		if (oldLen >= newLen) { return; }
		arr.resize(newLen);
		for (uint32_t i = oldLen; i < newLen; i++) {
			arr[i] = identity;
			if (i & 1) { continue; }
			uint32_t p = i ^ (i & (i - 1));
			if (i - p + 1 >= oldLen) { continue; }
			uint32_t t = i - p;
			p >>= 1;
			t += p;
			while (p) {
				arr[i] = func(arr[i], arr[t]);
				p >>= 1;
				t += p;
			}
		}
	}
public:
	BinaryIndexedTree() {
		identity = 0;
		func = [](T a, T b) { return a + b; };
		inv = [](T a) { return -a; };
		arr = vector<T>(1, identity);
	}
	BinaryIndexedTree(uint32_t n) {
		identity = 0;
		func = [](T a, T b) { return a + b; };
		inv = [](T a) { return -a; };
		arr = vector<T>(n + 1, identity);
	}
	BinaryIndexedTree(const vector<T>& aArr) { operator=(aArr); }
	BinaryIndexedTree(const T& aId, const function<T(T, T)>& aFunc, const function<T(T)>& aInv) {
		identity = aId;
		func = aFunc;
		inv = aInv;
		arr = vector<T>(1, aId);
	}
	BinaryIndexedTree(uint32_t n, const T& aId, const function<T(T, T)>& aFunc, const function<T(T)>& aInv) {
		identity = aId;
		func = aFunc;
		inv = aInv;
		arr = vector<T>(n + 1, aId);
	}
	BinaryIndexedTree(const vector<T>& aArr, const T& aId, const function<T(T, T)>& aFunc, const function<T(T)>& aInv) {
		identity = aId;
		func = aFunc;
		inv = aInv;
		setArray(aArr);
	}
	BinaryIndexedTree(const BinaryIndexedTree& other) { operator=(other); }
	BinaryIndexedTree& operator=(const BinaryIndexedTree& other) {
		func = other.func;
		inv = other.inv;
		identity = other.identity;
		arr = other.arr;
		return *this;
	}
	BinaryIndexedTree& operator=(const vector<T>& aArr) {
		identity = 0;
		func = [](T a, T b) { return a + b; };
		inv = [](T a) { return -a; };
		setArray(aArr);
		return *this;
	}
	void add(uint32_t i, const T& v) {
		uint32_t n = (uint32_t)arr.size() - 1;
		if (i + 1 > n) {
			extend(i + 2);
			n = i + 1;
		}
		uint32_t p = i + 1;
		while (p <= n) {
			arr[p] = func(arr[p], v);
			p += p ^ (p & (p - 1));
		}
	}
	T getSum(uint32_t i) const {
		uint32_t n = (uint32_t)arr.size() - 1;
		T ret = identity;
		uint32_t p = (i >= n) ? n : i;
		while (p) {
			ret = func(ret, arr[p]);
			p &= p - 1;
		}
		return ret;
	}
	T getSegment(uint32_t l, uint32_t r) const {
		if (l > r) { return identity; }
		return func(getSum(r), inv(getSum(l)));
	}
	T getValue(uint32_t i) const {
		return func(getSum(i + 1), inv(getSum(i)));
	}
};

#endif /* BinaryIndexedTree_hpp */
