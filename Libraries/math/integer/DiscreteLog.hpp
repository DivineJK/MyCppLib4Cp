//
//  DiscreteLog.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/09.
//

#ifndef DiscreteLog_hpp
#define DiscreteLog_hpp

class DiscreteLog {
private:
	static int modpow(int n, int m, int mod) {
		int ret = 1, b = n;
		while (m) {
			if (m & 1) { ret = (int64_t)ret * b % mod; }
			b = (int64_t)b * b % mod;
			m >>= 1;
		}
		return ret;
	}
	static unordered_map<int, int> factorize(int n) {
		unordered_map<int, int> ret;
		int a = n, c = 0;
		while (!(a & 1)) {
			a >>= 1;
			c++;
		}
		if (c) { ret[2] = c; }
		int p = 3;
		while (a != 1) {
			c = 0;
			while (a % p == 0) {
				a /= p;
				c++;
			}
			if (c) { ret[p] = c; }
			p += 2;
			if (a != 1 && p * p > n) {
				ret[a] = 1;
				break;
			}
		}
		return ret;
	}
	static int getTotient(int n) {
		int ret = (n & 1) ? n : n >> 1;
		int a = ret, p = 3;
		while (!(a & 1)) { a >>= 1; }
		while (a > 1) {
			if (a % p == 0) {
				ret /= p;
				ret *= p - 1;
				while (a % p == 0) { a /= p; }
			}
			p += 2;
			if (a > 1 && p * p > n) {
				ret /= a;
				ret *= a - 1;
				break;
			}
		}
		return ret;
	}
	static vector<int> getDivisors(int n) {
		int p = 1;
		vector<int> ret;
		while (p * p <= n) {
			if (n % p == 0) { ret.push_back(p); }
			p++;
		}
		size_t s = ret.size();
		if (s == 0) { return ret; }
		if (ret[s - 1] * ret[s - 1] == n) { s--; }
		for (size_t i = s; i > 0; i--) { ret.push_back(n / ret[i - 1]); }
		return ret;
	}
	static int gcd(int a, int b) {
		int x = a, y = b, tmp;
		while (y) {
			tmp = x;
			x = y;
			y = tmp % y;
		}
		return x;
	}
public:
	static int solve(int x, int y, int m) {
		if (m <= 0 || y >= m || y < 0) { return -1; }
		x = (x < 0) ? x % m + m : x % m;
		if (x == 0) {
			if (m == 1) {
				return 0;
			} else {
				return (y > 1) ? -1 : 1 - y;
			}
		} else if (y == 1) { return 0; }
		unordered_map<int, int> xf = factorize(x);
		int mp = m;
		int t = 0;
		for (const pair<const int, int>& p : xf) {
			int c = 0;
			while (mp % p.first == 0) {
				mp /= p.first;
				c++;
			}
			c = (c + p.second - 1) / p.second;
			if (t < c) { t = c; }
		}
		int v = 1;
		for (int i = 0; i < t; i++) {
			if (v == y) { return i; }
			v = (int64_t)v * x % m;
		}
		if (y % gcd(v, m)) { return -1; }
		mp = getTotient(mp);
		vector<int> divs = getDivisors(mp);
		int l = mp;
		for (size_t i = 0; i < divs.size(); i++) {
			if ((int64_t)v * modpow(x, divs[i], m) % m == v) {
				l = divs[i];
				break;
			}
		}
		int ix = modpow(x, l - 1, m) % m, lsq = l >> 1;
		int lv = 0, rv = l;
		while (rv - lv > 1) {
			if ((int64_t)lsq * lsq <= l) {
				lv = lsq;
			} else {
				rv = lsq;
			}
			lsq = lv + ((rv - lv) >> 1);
		}
		if (lsq * lsq < l) { lsq++; }
		unordered_map<int, int> baby_dict;
		int f = v;
		for (int i = 0; i < lsq; i++) {
			baby_dict[f] = i;
			f = (int64_t)f * x % m;
		}
		ix = modpow(ix, lsq, m);
		int g = y;
		for (int i = 0; i < lsq; i++) {
			if (baby_dict.count(g)) { return i * lsq + baby_dict[g] + t; }
			g = (int64_t)g * ix % m;
		}
		return -1;
	}
};

#endif /* DiscreteLog_hpp */
