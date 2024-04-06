//
//  BinarySearch.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2024/02/11.
//

#ifndef BinarySearch_h
#define BinarySearch_h

template <typename T>
size_t bisect_increase_leftmost_greq(const vector<T>& vec, T val) {
	size_t l = 0, r = vec.size();
	size_t d = r >> 1;
	while (r - l > 1) {
		if (vec[d] <= val) {
			l = d;
		} else {
			r = d;
		}
		d = l + ((r - l) >> 1);
	}
	return d;
}

#endif /* BinarySearch_h */
