//
//  DebugUtil.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/23.
//

#ifndef DebugUtil_h
#define DebugUtil_h

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

#endif /* DebugUtil_h */
