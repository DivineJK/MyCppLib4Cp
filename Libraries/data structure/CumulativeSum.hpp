//
//  CumulativeSum.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/11/05.
//

#ifndef CumulativeSum_h
#define CumulativeSum_h

template <typename T>
class CumulativeSum {
private:
	T unit = 0;
	vector<T> arr;
public:
	CumulativeSum() { arr = vector<T>(1, unit); }
	explicit CumulativeSum(const vector<T>& aArr, const T& aUnit = 0) {
		unit = aUnit;
		size_t n = aArr.size();
		arr = vector<T>(n + 1, unit);
		for (size_t i = 0; i < n; i++) { arr[i + 1] = arr[i] + aArr[i]; }
	}
	CumulativeSum(const CumulativeSum& other) { operator=(other); }
	CumulativeSum& operator=(const CumulativeSum& other) {
		unit = other.unit;
		arr = other.arr;
	}
	size_t size() const { return arr.size() - 1; }
	void setValue(size_t i, const T& v) { arr[i] = v; }
	void flush() {
		size_t n = size();
		for (size_t i = 0; i < n; i++) { arr[i + 1] += arr[i]; }
	}
	T getSum(size_t n) const { return arr[min(n, size())]; }
	T getSum(size_t l, size_t r) const {
		size_t n = size();
		size_t regl = (l >= n) ? n : l;
		size_t regr = (r >= n) ? n : r;
		if (regl >= regr) { return unit; }
		return arr[regr] - arr[regl];
	}
};

#endif /* CumulativeSum_h */
