//
//  LIS.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2024/02/10.
//

#ifndef LIS_h
#define LIS_h

template <typename T>
vector<size_t> LIS(const vector<T>& v) {
	size_t n = v.size();
	if (n == 0) { return {}; }
	vector<size_t> min_back(n + 1, n);
	vector<size_t> prev(n, n);
	min_back[0] = n + 1;
	size_t lm = 0;
	for (size_t i = 0; i < n; i++) {
		size_t l = 0, r = i + 1;
		size_t d = i >> 1;
		while (r - l > 1) {
			if (min_back[d] == n + 1 || (min_back[d] < n && v[min_back[d]] < v[i])) {
				l = d;
			} else {
				r = d;
			}
			d = l + ((r - l) >> 1);
		}
		if (min_back[d + 1] == n || v[i] < v[min_back[d + 1]]) {
			min_back[d + 1] = i;
			prev[i] = min_back[d];
			lm = max(lm, d + 1);
		}
	}
	vector<size_t> ret;
	ret.reserve(lm);
	size_t idx = min_back[lm];
	for (size_t i = lm; i > 0; i--) {
		ret.push_back(idx);
		idx = prev[idx];
	}
	reverse(ret.begin(), ret.end());
	return ret;
}

#endif /* LIS_h */
