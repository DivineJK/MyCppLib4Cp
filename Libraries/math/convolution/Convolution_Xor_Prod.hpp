//
//  Convolution_Xor_Prod.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef Convolution_Xor_Prod_hpp
#define Convolution_Xor_Prod_hpp

template <typename T>
void FHT(vector<T>* a) {
	size_t n = a->size();
	if (n <= 1) { return; }
	size_t m = 1;
	if (n & (n - 1)) {
		while (m < n) { m <<= 1; }
		a->resize(m);
	} else {
		m = n;
	}
	for (size_t i = 1; i <= (m >> 1); i <<= 1) {
		for (size_t j = 0; j < m; j += (i << 1)) {
			for (size_t k = j; k < j + i; k++) {
				T x = (*a)[k], y = (*a)[k + i];
				(*a)[k] = x + y;
				(*a)[k + i] = x - y;
			}
		}
	}
}

template <typename T>
void convolve_xor_prod(vector<T>* a, const vector<T>& b) {
	if (a->size() == 0 || b.size() == 0) {
		a->clear();
		return;
	}
	size_t n = max(a->size(), b.size());
	if (n & (n - 1)) {
		size_t m = 1;
		while (m < n) { m <<= 1; }
		n = m;
	}
	a->resize(n);
	vector<T> y = b;
	y.resize(n);
	FHT(a);
	FHT(&y);
	for (size_t i = 0; i < n; i++) { (*a)[i] *= y[i]; }
	FHT(a);
	for (size_t i = 0; i < n; i++) { (*a)[i] /= (T)n; }
}

#endif /* Convolution_Xor_Prod_hpp */
