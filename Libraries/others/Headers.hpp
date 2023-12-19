//
//  Headers.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/11/10.
//

#ifndef Headers_h
#define Headers_h

#include <algorithm>
#include <atomic>
#include <bit>
#include <bitset>
#include <cassert>
#include <cstdarg>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <err.h>
#include <errno.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <locale>
#include <map>
#include <memory>
#include <new>
#include <numbers>
#include <numeric>
#include <optional>
#include <queue>
#include <random>
#include <ranges>
#include <regex>
#include <set>
#include <span>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <wchar.h>

//#include <bits/stdc++.h>
using namespace std;
template <typename T1, typename T2>
ostream& operator<<(ostream& os, const pair<T1, T2>& p) {
	os << p.first << ' ' << p.second;
	return os;
}
template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v) {
	size_t n = v.size();
	if (n == 0) { return os; }
	os << v[0];
	if (n == 1) { return os; }
	for (size_t i = 1; i < n; i++) { os << ' ' << v[i]; }
	return os;
}
template <typename T>
ostream& operator<<(ostream& os, const vector<vector<T>>& mat) {
	size_t n = mat.size();
	if (n == 0) { return os; }
	size_t m0 = mat[0].size();
	if (m0 > 0) {
		os << mat[0][0];
		for (size_t i = 1; i < m0; i++) { os << ' ' << mat[0][i]; }
	}
	for (size_t i = 1; i < n; i++) {
		os << '\n';
		size_t m = mat[i].size();
		if (m == 0) { continue; }
		os << mat[i][0];
		for (size_t j = 1; j < m; j++) { os << ' ' << mat[i][j]; }
	}
	return os;
}
template <typename T>
void readVector(vector<T>* vec) {
	for (size_t i = 0; i < vec->size(); i++) { cin >> (*vec)[i]; }
}
template <typename T>
void printVector(const vector<T>& vec) {
	if (vec.size() == 0) {
		cout << '\n';
		return;
	}
	for (size_t i = 0; i < vec.size(); i++) {
		if (i > 0) { cout << " "; }
		cout << vec[i];
	}
	cout << '\n';
}
template <typename T>
void printArray(T* arr, size_t s) {
	if (s == 0) {
		cout << '\n';
		return;
	}
	for (size_t i = 0; i < s; i++) {
		if (i > 0) { cout << " "; }
		cout << arr[i];
	}
	cout << '\n';
}

/* Headers end */

#endif /* Headers_h */
