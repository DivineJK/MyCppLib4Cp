//
//  Integer_All.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef Integer_All_hpp
#define Integer_All_hpp

uint64_t floorSum(uint64_t n, uint64_t m, uint64_t a, uint64_t b) {
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

#endif /* Integer_All_hpp */
