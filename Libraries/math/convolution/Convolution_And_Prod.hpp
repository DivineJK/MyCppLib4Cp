//
//  Convolution_And_Prod.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef Convolution_And_Prod_hpp
#define Convolution_And_Prod_hpp

template <typename T>
void transform_upper_subset(vector<T>* a) {
	size_t n = a->size();
	if (n == 0) { return; }
	assert(!(n & (n - 1)));
	for (size_t i = 1; i < n; i <<= 1) {
		for (size_t j = 0; j < n; j += (i << 1)) {
			for (size_t k = j; k < j + i; k++) { (*a)[k] += (*a)[k + i]; }
		}
	}
}

template <typename T>
void inv_transform_upper_subset(vector<T>* a) {
	size_t n = a->size();
	if (n == 0) { return; }
	assert(!(n & (n - 1)));
	for (size_t i = 1; i < n; i <<= 1) {
		for (size_t j = 0; j < n; j += (i << 1)) {
			for (size_t k = j; k < j + i; k++) { (*a)[k] -= (*a)[k + i]; }
		}
	}
}

template <typename T>
void convolve_and_prod(vector<T>* a, const vector<T>& b) {
	size_t n = max(a->size(), b.size());
	if (n == 0) { return; }
	if (n & (n - 1)) {
		size_t m = 1;
		while (m < n) { m <<= 1; }
		n = m;
		a->resize(n);
	}
	vector<T> y(n);
	for (size_t i = 0; i < b.size(); i++) { y[i] = b[i]; }
	transform_upper_subset(a);
	transform_upper_subset(&y);
	for (size_t i = 0; i < n; i++) { (*a)[i] *= y[i]; }
	inv_transform_upper_subset(a);
}

#endif /* Convolution_And_Prod_hpp */
