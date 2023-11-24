//
//  Integer_All.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef Integer_All_hpp
#define Integer_All_hpp

class Integer_All {
public:
	template <typename T, enable_if_t<is_integral_v<T>>* = nullptr>
	static void safety_divide(T a, T b, T* d, T* m) {
		if (b == 0) {
			*d = (a == 0) ? NAN : ((a < 0) ? -INFINITY : INFINITY);
			*m = (a == 0) ? NAN : INFINITY;
			return;
		}
		*d = a / b;
		*m = a % b;
		if (*m < 0) {
			*m += (b < 0) ? -b : b;
			*d += (b < 0) ? 1 : -1;
		}
	}
	template <typename T, enable_if_t<is_integral_v<T>>* = nullptr>
	static bool solveLinearEquation(T a, T b, T c, pair<T, T>* sol) {
		if (a == 0 && b == 0) {
			*sol = make_pair(0, 0);
			return c == 0;
		}
		if (a == 0) {
			*sol = make_pair(0, c / b);
			return c % b == 0;
		}
		if (b == 0) {
			*sol = make_pair(c / a, 0);
			return c % a == 0;
		}
	}
	template <typename T, enable_if_t<is_integral_v<T>>* = nullptr,
		typename S, enable_if_t<is_integral_v<S>>* = nullptr>
	static T modinv(T a, S m) {
		
	}
	template <typename T, enable_if_t<is_integral_v<T>>* = nullptr,
		typename S, enable_if_t<is_integral_v<S>>* = nullptr>
	static T modpow(T a, S n, T m) {
		
	}
	static uint64_t floorSum(uint64_t n, uint64_t m, uint64_t a, uint64_t b) {
		uint64_t ret = 0;
		bool sgn = false;
		while (m != 0 && n != 0) {
			ret = (sgn) ? ret - n * (n - 1) * (a / m) / 2 - (b / m) * n : ret + n * (n - 1) * (a / m) / 2 + (b / m) * n;
			a %= m;
			b %= m;
			uint64_t t = (a * (n - 1) + b) / m;
			ret = (sgn) ? ret - t * (n - 1) : ret + t * (n - 1);
			sgn ^= true;
			n = t;
			uint64_t tmp = m;
			m = a;
			a = tmp;
			b = tmp - b - 1;
		}
		return ret;
	}
};

#endif /* Integer_All_hpp */
