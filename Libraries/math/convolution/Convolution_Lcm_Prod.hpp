//
//  Convolution_Lcm_Prod.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2024/01/27.
//

#ifndef Convolution_Lcm_Prod_h
#define Convolution_Lcm_Prod_h

namespace conv_lcm_internal {
static inline vector<size_t>& persistent_sieve(size_t n) {
	static vector<bool> is_p = {false, false, true};
	static vector<size_t> primes = {2};
	if (n >= is_p.size()) {
		size_t old = is_p.size();
		is_p.resize(n + 1);
		for (size_t i = old; i <= n; i++) { is_p[i] = true; }
		for (size_t i = ((old + 1) >> 1) << 1; i <= n; i += 2) { is_p[i] = false; }
		size_t m = 3;
		while (m * m <= n) {
			if (!is_p[m]) {
				m += 2;
				continue;
			}
			size_t iv = max(m * m, (old + m - 1) / m * m);
			for (size_t i = iv; i <= n; i += m) { is_p[i] = false; }
			m += 2;
		}
		for (size_t i = old; i <= n; i++) {
			if (is_p[i]) { primes.push_back(i); }
		}
	}
	return primes;
}
template <typename T>
void FZT(vector<T>* a) {
	size_t n = a->size();
	vector<size_t> primes = persistent_sieve(n);
	for (size_t p : primes) {
		for (size_t i = 1; i * p <= n; i++) { (*a)[i * p - 1] += (*a)[i - 1]; }
	}
}
template <typename T>
void FMT(vector<T>* a) {
	size_t n = a->size();
	vector<size_t> primes = persistent_sieve(n);
	for (size_t p : primes) {
		for (size_t i = n / p; i > 0; i--) { (*a)[i * p - 1] -= (*a)[i - 1]; }
	}
}
} // namespace conv_lcm_internal

template <typename T>
void convolve_lcm_prod(vector<T>* a, const vector<T>& b, size_t lim = UINT64_MAX) {
	if (a->size() == 0 || b.size() == 0) {
		a->clear();
		return;
	}
	size_t n = min(a->size() * b.size(), lim);
	a->resize(n);
	vector<T> y(n);
	for (size_t i = 0; i < b.size(); i++) { y[i] = b[i]; }
	conv_lcm_internal::FZT(a);
	conv_lcm_internal::FZT(&y);
	for (size_t i = 0; i < n; i++) { (*a)[i] *= y[i]; }
	conv_lcm_internal::FMT(a);
}

#endif /* Convolution_Lcm_Prod_h */
