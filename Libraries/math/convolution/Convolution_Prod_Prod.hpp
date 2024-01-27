//
//  Convolution_Prod_Prod.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef Convolution_Prod_Prod_hpp
#define Convolution_Prod_Prod_hpp

template <typename T>
void convolve_prod_prod(vector<T>* a, const vector<T>& b, size_t lim = UINT64_MAX) {
	size_t m = min(a->size() * b.size(), lim);
	a->resize(m);
	vector<T> bck = *a;
	*a = vector<T>(m);
	for (size_t i = 1; i <= m; i++) {
		T val = bck[i - 1];
		for (size_t j = 1; i * j <= m && j <= b.size(); j++) {
			(*a)[i * j - 1] += val * b[j - 1];
		}
	}
}

#endif /* Convolution_Prod_Prod_hpp */
