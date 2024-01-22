//
//  CountSubsequences.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2024/01/21.
//

#ifndef CountSubsequences_h
#define CountSubsequences_h

template <typename S, typename T>
S countSubsequences(const vector<T>& v) {
	S ret = 1, tmp, x;
	unordered_map<T, S> cm;
	cm.reserve(v.size());
	for (size_t i = 0; i < v.size(); i++) {
		S tmp = ret;
		S x = (cm.count(v[i])) ? cm[v[i]] : 0;
		ret += ret - x;
		cm[v[i]] = tmp;
	}
	return ret;
}

template <typename S>
S countSubsequences(const string& s) {
	S ret = 1;
	unordered_map<char, S> cm;
	cm.reserve(s.size());
	for (size_t i = 0; i < s.size(); i++) {
		S tmp = ret;
		S x = (cm.count(s[i])) ? cm[s[i]] : 0;
		ret += ret - x;
		cm[s[i]] = tmp;
	}
	return ret;
}

#endif /* CountSubsequences_h */
