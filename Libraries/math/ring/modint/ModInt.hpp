//
//  ModInt.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef ModInt_hpp
#define ModInt_hpp

namespace modint_internal {
struct internal_primitive_bases {
private:
	uint32_t mod = 1;
	int base_size = 0;
	uint32_t roots[2]; // primitive root, inverse root
	vector<uint32_t> bases[3]; // primitive bases, inverse bases, cumulative bases
private:
	internal_primitive_bases() = delete;
	bool is_primitive(uint32_t n) const {
		uint32_t t = n;
		for (int i = 0; i < base_size; i++) {
			if (t == 1) { return false; }
			t = (uint32_t)((uint64_t)t * t % mod);
		}
		return t == 1;
	}
	int modinv(uint32_t n) const {
		uint32_t p = n, t = 1;
		while (p > 1) {
			t = (uint64_t)t * (mod - (mod / p)) % mod;
			p = mod % p;
		}
		return t;
	}
	static constexpr bool is_prime(uint32_t n) {
		if (n <= 1) { return false; }
		if (n != 2 && !(n & 1)) { return false; }
		uint32_t p = 3;
		while (p * p <= n) {
			if (n % p == 0) { return false; }
			p += 2;
		}
		return true;
	}
public:
	internal_primitive_bases(uint32_t aMod) {
		assert(is_prime(aMod));
		mod = aMod;
		make_base();
	}
	void make_base() {
		int r = 1, p = mod - 1;
		base_size = 0;
		while (!(p & 1)) {
			p >>= 1;
			base_size++;
		}
		for (int i = 0; i < 3; i++) { bases[i] = vector<uint32_t>(base_size + 1); }
		while (true) {
			if (is_primitive(r)) {
				roots[0] = r;
				roots[1] = modinv(r);
				bases[0][base_size] = r;
				bases[1][base_size] = roots[1];
				for (int i = base_size; i > 0; i--) {
					bases[0][i - 1] = (uint64_t)bases[0][i] * bases[0][i] % mod;
					bases[1][i - 1] = (uint64_t)bases[1][i] * bases[1][i] % mod;
				}
				int ie = 1;
				for (int i = 0; i < base_size - 1; i++) {
					bases[2][i] = (uint64_t)ie * bases[0][i + 2] % mod;
					ie = (uint64_t)ie * bases[1][i + 2] % mod;
				}
				break;
			}
			r++;
		}
	}
	uint32_t get_base(int t, int idx) const {
		assert(0 <= t || t < 3);
		return (idx < 0 || idx > base_size) ? 1 : bases[t][idx];
	}
	uint32_t get_root(int t) const {
		assert(t == 0 || t == 1);
		return roots[t];
	}
	int get_base_size() const { return base_size; }
};
static unordered_map<uint32_t, unique_ptr<internal_primitive_bases>> primitive_bases;
}

template <int MOD>
class ModInt {
	using mint_prim_base = modint_internal::internal_primitive_bases;
private:
	uint32_t x = 0;
	static constexpr int PRE_MOD_21 = 924844033;
	static constexpr int PRE_MOD_22 = 985661441;
	static constexpr int PRE_MOD_23 = 998244353;
	static constexpr int PRE_MOD_24 = 754974721;
	static constexpr int PRE_MOD_25 = 167772161;
	static constexpr int PRE_MOD_26 = 469762049;
	using mint21 = ModInt<PRE_MOD_21>;
	using mint22 = ModInt<PRE_MOD_22>;
	using mint23 = ModInt<PRE_MOD_23>;
	using mint24 = ModInt<PRE_MOD_24>;
	using mint25 = ModInt<PRE_MOD_25>;
	using mint26 = ModInt<PRE_MOD_26>;
private:
	static constexpr int mod() { return MOD; }
	static constexpr uint32_t umod() { return (uint32_t)MOD; }
	static constexpr bool is_prime(uint32_t n) {
		if (n <= 1) { return false; }
		if (n != 2 && !(n & 1)) { return false; }
		uint32_t p = 3;
		while (p * p <= n) {
			if (n % p == 0) { return false; }
			p += 2;
		}
		return true;
	}
	static void make_base(uint32_t m) { modint_internal::primitive_bases[m] = make_unique<mint_prim_base>(m); }
	static mint_prim_base* get_bases(uint32_t m) {
		if (modint_internal::primitive_bases.find(m) == modint_internal::primitive_bases.end()) {
			make_base(m);
		}
		return modint_internal::primitive_bases[m].get();
	}
	static int get_base(int m, int t, int idx) { return get_bases(m)->get_base(t, idx); }
	static void ntt(vector<ModInt>* arr) {
		size_t n = arr->size();
		if (n <= 1) { return; }
		assert(!(n & (n - 1)));
		int m = 0;
		while ((1 << m) < n) { m++; }
		int p = 1, t = 1 << (m - 1);
		for (int i = 0; i < m; i++) {
			ModInt g = 1;
			for (int j = 0; j < p; j++) {
				ModInt u, v;
				int offset = j << (m - i);
				for (int k = 0; k < t; k++) {
					u = (*arr)[k + offset], v = (*arr)[k + offset + t] * g;
					(*arr)[k + offset] = u + v;
					(*arr)[k + offset + t] = u - v;
				}
				int w = (j + 1) ^ (j & (j + 1)), z = -1;
				while (w) {
					z++;
					w >>= 1;
				}
				g *= get_base(umod(), 2, z);
			}
			p <<= 1;
			t >>= 1;
		}
	}
	static void intt(vector<ModInt>* arr) {
		size_t n = arr->size();
		if (n <= 1) { return; }
		ModInt invn = inv((int)n);
		assert(!(n & (n - 1)));
		int m = 0;
		while ((1 << m) < n) { m++; }
		for (int i = 0; i < m; i++) {
			ModInt g = 1, s = get_base(umod(), 1, i + 1);
			int offset = 1 << i;
			for (int k = 0; k < offset; k++) {
				for (int j = k; j < n; j += (1 << (i + 1))) {
					ModInt x = (*arr)[j], y = (*arr)[j + offset] * g;
					(*arr)[j] = x + y;
					(*arr)[j + offset] = x - y;
				}
				g *= s;
			}
		}
		for (int i = 0; i < n; i++) {
			(*arr)[i] *= invn;
		}
	}
	static constexpr bool is_prepared_base() {
		return MOD == PRE_MOD_21 || MOD == PRE_MOD_22 || MOD == PRE_MOD_23
			|| MOD == PRE_MOD_24 || MOD == PRE_MOD_25 || MOD == PRE_MOD_26;
	}
	static bool is_prepared() {
		if (is_prepared_base()) { return true; }
		if (modint_internal::primitive_bases.find(MOD) != modint_internal::primitive_bases.end()) {
			return true;
		}
		if (!is_prime(MOD)) { return false; }
		int p = MOD - 1;
		for (int i = 0; i < 20; i++) {
			if (p & 1) { return false; }
			p >>= 1;
		}
		return true;
	}
	static void fht(vector<ModInt>* a) {
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
					ModInt x = (*a)[k], y = (*a)[k + i];
					(*a)[k] = x + y;
					(*a)[k + i] = x - y;
				}
			}
		}
	}
public:
	ModInt() : x(0) {}
	ModInt(int a) { x = (a % mod() < 0) ? a % mod() + umod() : a % mod(); }
	ModInt(int64_t a) { x = (a % mod() < 0) ? a % mod() + umod() : a % mod(); }
	ModInt(uint32_t a) { x = a % umod(); }
	ModInt(uint64_t a) { x = a % umod(); }
	ModInt(const ModInt& other) { operator=(other); }
	uint32_t val() const { return x; }
	ModInt inv() const {
		if (x == 0) { return ModInt(); }
		uint32_t a = x, t = 1;
		while (a > 1) {
			t = (uint64_t)t * (umod() - umod() / a) % umod();
			a = umod() % a;
		}
		return ModInt(t);
	}
	ModInt extended_inv() const {
		if (x == 0) { return 0; }
		int a = 1, b = 0, c = 0, d = 1, k = x, l = mod();
		while (l > 0) {
			uint32_t t1 = a, t2 = b, t3 = k;
			a = c;
			b = d;
			c = t1 - c * (k / l);
			d = t2 - d * (k / l);
			k = l;
			l = t3 % l;
		}
		return a;
	}
	static ModInt inv(const ModInt& v) {
		if (v.x == 0) { return v; }
		uint32_t a = v.x, t = 1;
		while (a > 1) {
			t = (uint64_t)t * (umod() - umod() / a) % umod();
			a = umod() % a;
		}
		return ModInt(t);
	}
	static ModInt extended_inv(const ModInt& v) {
		return v.extended_inv();
	}
	ModInt& operator=(const ModInt& other) {
		x = other.x;
		return *this;
	}
	friend ostream& operator<<(ostream& os, const ModInt& other) {
		os << other.x;
		return os;
	}
	friend istream& operator>>(istream& ist, ModInt& other) {
		int64_t a;
		ist >> a;
		int v = a % mod();
		other.x = (v < 0) ? v + umod() : v;
		return ist;
	}
	explicit operator bool() const { return x != 0; }
	bool operator!() const { return x == 0; }
	bool operator==(const ModInt& other) const { return x == other.x; }
	bool operator!=(const ModInt& other) const { return x != other.x; }
	ModInt operator+() const { return ModInt(*this); }
	ModInt operator-() const {
		ModInt ret = ModInt(*this);
		ret.x = (ret.x == 0) ? 0 : umod() - ret.x;
		return ret;
	}
	ModInt& operator++() {
		x = (x + 1 == umod()) ? 0 : x + 1;
		return *this;
	}
	ModInt operator++(int) {
		ModInt ret = ModInt(*this);
		++*this;
		return ret;
	}
	ModInt& operator--() {
		x = (x == 0) ? umod() - 1 : x - 1;
		return *this;
	}
	ModInt operator--(int) {
		ModInt ret = ModInt(*this);
		--*this;
		return ret;
	}
	ModInt& operator+=(const ModInt& other) {
		x = (x + other.x >= umod()) ? x + other.x - umod() : x + other.x;
		return *this;
	}
	ModInt& operator-=(const ModInt& other) {
		x = (x < other.x) ? x + umod() - other.x : x - other.x;
		return *this;
	}
	ModInt& operator*=(const ModInt& other) {
		x = (uint32_t)((uint64_t)x * other.x % umod());
		return *this;
	}
	ModInt& operator/=(const ModInt& other) { return *this *= inv(other); }
	ModInt operator+(const ModInt& other) const { return ModInt(*this) += other; }
	ModInt operator-(const ModInt& other) const { return ModInt(*this) -= other; }
	ModInt operator*(const ModInt& other) const { return ModInt(*this) *= other; }
	ModInt operator/(const ModInt& other) const { return ModInt(*this) /= other; }
	template <typename T>
	friend ModInt operator+(T a, const ModInt& other) { return ModInt(a) += other; }
	template <typename T>
	friend ModInt operator-(T a, const ModInt& other) { return ModInt(a) -= other; }
	template <typename T>
	friend ModInt operator*(T a, const ModInt& other) { return ModInt(a) *= other; }
	template <typename T>
	friend ModInt operator/(T a, const ModInt& other) { return ModInt(a) /= other; }
	template <typename T>
	static ModInt pow(const ModInt& n, T m) {
		ModInt ret = 1, v;
		if (m < 0) {
			m = -m;
			v = inv(n);
		} else { v = n; }
		while (m) {
			if (m & 1) { ret *= v; }
			v *= v;
			m >>= 1;
		}
		return ret;
	}
	template <typename T>
	ModInt operator^(T m) const { return pow(*this, m); }
	static int sqrt(const ModInt& n) {
		if (n.val() == 0 || n.val() == 1) { return n.val(); }
		if (pow(n, umod() >> 1) + 1 == 0) { return -1; }
		int m = mod() - 1, q = 0;
		while (!(m & 1)) {
			q++;
			m >>= 1;
		}
		ModInt z = 1, r = pow(n, (m + 1) >> 1), t = pow(n, m);
		while (pow(z, umod() >> 1) == 1) { z++; }
		z = pow(z, m);
		while (t != 1) {
			ModInt c = t;
			int i = 0;
			while (c != 1) {
				i++;
				c *= c;
			}
			ModInt b = pow(z, 1 << (q - i - 1));
			q = i;
			z = b * b;
			t *= z;
			r *= b;
		}
		return (r.x > mod() - r.x) ? mod() - r.x : r.x;
	}
	static void convolve_one_base(vector<ModInt>* a, const vector<ModInt>& b) {
		if (a->size() == 0 || b.size() == 0) {
			*a = {};
			return;
		}
		size_t n = a->size(), m = b.size(), bas = 1;
		while (bas < n + m) { bas <<= 1; }
		a->resize(bas);
		vector<ModInt> v(bas);
		for (size_t i = 0; i < m; i++) { v[i] = b[i]; }
		ntt(a);
		ntt(&v);
		for (size_t i = 0; i < bas; i++) { (*a)[i] *= v[i]; }
		intt(a);
	}
	static void convolve_one(vector<ModInt>* a, const vector<ModInt>& b) {
		if (is_prepared()) {
			convolve_one_base(a, b);
			return;
		}
		size_t n = a->size(), m = b.size();
		vector<mint23> x23(n);
		vector<mint24> x24(n);
		vector<mint26> x26(n);
		vector<mint23> y23(m);
		vector<mint24> y24(m);
		vector<mint26> y26(m);
		for (int i = 0; i < n; i++) {
			x23[i] = (*a)[i].val();
			x24[i] = (*a)[i].val();
			x26[i] = (*a)[i].val();
		}
		for (int i = 0; i < m; i++) {
			y23[i] = b[i].val();
			y24[i] = b[i].val();
			y26[i] = b[i].val();
		}
		mint23::convolve_one_base(&x23, y23);
		mint24::convolve_one_base(&x24, y24);
		mint26::convolve_one_base(&x26, y26);
		size_t r = x23.size();
		a->resize(r);
		function<int(int, int)> rev = [](int m, int a) {
			int p = a, ret = 1;
			while (p > 1) {
				ret = (int64_t)ret * (m - m / p) % m;
				p = m % p;
			}
			return ret;
		};
		function<int(int, int, int, int, int, int, int)> rem;
		rem = [&rev](int m, int v0, int v1, int v2, int m0, int m1, int m2) {
			int m01 = rev(m1, m0), m02 = rev(m2, m0), m12 = rev(m2, m1);
			int a0 = v0;
			int a1 = (int64_t)m01 * (v1 - a0) % m1;
			if (a1 < 0) { a1 += m1; }
			int a2 = (int64_t)m02 * (v2 - a0) % m2;
			if (a2 < 0) { a2 += m2; }
			a2 = (int64_t)m12 * (a2 - a1) % m2;
			if (a2 < 0) { a2 += m2; }
			a1 = (int64_t)a1 * m0 % m;
			a2 = (int64_t)a2 * m0 % m;
			a2 = (int64_t)a2 * m1 % m;
			return ((a0 + a1) % m + a2) % m;
		};
		for (int i = 0; i < r; i++) {
			(*a)[i] = rem(MOD, x26[i].val(), x24[i].val(), x23[i].val(),
						  PRE_MOD_26, PRE_MOD_24, PRE_MOD_23);
		}
	}
	static void convolve(vector<ModInt>* a, const vector<ModInt>& b) {
		if (a->size() == 0 || b.size() == 0) {
			*a = {};
			return;
		}
		int n = (int)a->size(), m = (int)b.size();
		int bas = 1, c = 0;
		while (bas < n + m) {
			bas <<= 1;
			c++;
		}
		int bs = 21;
		if (is_prepared_base()) {
			if (mod() == PRE_MOD_21) { bs = 19; }
			else if (mod() == PRE_MOD_22) { bs = 20; }
		}
		if (bas <= (1 << bs)) {
			convolve_one(a, b);
			return;
		}
		int lim = 1 << 18;
		int s = (n + lim - 1) / lim, t = (m + lim - 1) / lim;
		a->resize(s * lim);
		vector<ModInt> y(t * lim);
		for (int i = 0; i < m; i++) { y[i] = b[i]; }
		vector<ModInt> ret((s + t) * lim);
		for (int i = 0; i < s; i++) {
			vector<ModInt> u(lim);
			for (int j = 0; j < lim; j++) { u[j] = (*a)[lim * i + j]; }
			for (int j = 0; j < t; j++) {
				vector<ModInt> res(lim);
				for (int k = 0; k < lim; k++) { res[k] = y[lim * j + k]; }
				convolve_one(&u, res);
				for (int k = 0; k < u.size(); k++) { ret[(i + j) * lim + k] += u[k]; }
			}
		}
		*a = ret;
	}
	static void get_fact_table(int n, vector<ModInt>* fact,
							   vector<ModInt>* invf) {
		assert(n >= 0);
		ModInt t = 1;
		if (fact) {
			fact->resize(n + 1);
			(*fact)[0] = 1;
			for (int i = 0; i < n; i++) { (*fact)[i + 1] = (i + 1) * (*fact)[i]; }
		} else { for (int i = 0; i < n; i++) { t *= i + 1; } }
		if (invf) {
			invf->resize(n + 1);
			(*invf)[n] = (fact) ? (*fact)[n].inv() : t.inv();
			for (int i = n; i > 0; i--) { (*invf)[i - 1] = i * (*invf)[i]; }
		}
	}
	static void get_inv_table(int n, vector<ModInt>* invn) {
		assert(invn && n >= 0);
		invn->resize(n + 1);
		for (int i = 0; i < 2; i++) {
			(*invn)[i] = 1;
			if (n == i) { return; }
		}
		for (int i = 2; i <= n; i++) {
			int t = umod() % i;
			(*invn)[i] = (umod() - umod() / i) * (*invn)[t];
		}
	}
	static vector<ModInt> getMontmort(size_t n) {
		if (n == 0) { return {}; }
		vector<ModInt> ret(n);
		ret[0] = 0;
		for (size_t i = 1; i < n; i++) {
			ret[i] = ret[i - 1] * (i + 1) + ((i & 1) ? 1 : -1);
		}
		return ret;
	}
	static vector<ModInt> getPowerTable(size_t n, uint32_t k) {
		vector<ModInt> ret(n + 1, 1);
		if (k == 0) {
			for (int i = 0; i <= n; i++) { ret[i] = 1; }
			return ret;
		}
		if (k == 1) {
			for (int i = 0; i <= n; i++) { ret[i] = i; }
			return ret;
		}
		ret[0] = 0;
		if (n == 0) { return ret; }
		vector<bool> flg(n + 1, false);
		for (int i = 2; i <= n; i++) {
			if (i % mod() == 0 || flg[i]) { continue; }
			ret[i] = pow(i, k);
			flg[i] = true;
			for (int j = 2; j * i <= n && j <= i; j++) {
				ret[i * j] = ret[i] * ret[j];
				flg[i * j] = true;
			}
		}
		return ret;
	}
	static ModInt getStirlingFirst(uint64_t n, uint64_t k) {
		static vector<vector<ModInt>> seeds(mod(), vector<ModInt>(mod(), 0));
		static bool needCalc = true;
		if (needCalc) {
			needCalc = false;
			seeds[0][0] = 1;
			for (int i = 1; i < mod(); i++) {
				seeds[i][0] = 0;
				for (int j = 1; j < mod(); j++) {
					seeds[i][j] = seeds[i - 1][j - 1] - (i - 1) * seeds[i - 1][j];
				}
			}
		}
		if (k < n / mod()) { return 0; }
		uint64_t q = n / mod();
		if ((k - q) % (mod() - 1) > n % mod()) { return 0; }
		ModInt ret = 0;
		function<ModInt(uint64_t, uint64_t)> vf = [n, q](uint64_t l, uint64_t m) {
			ModInt s = ((q & 1) == (l & 1)) ? 1 : -1;
			return getCombinationOne(q, l) * s * seeds[n % mod()][m];
		};
		if ((n + 1) % mod() == 0 && (k - q) % (mod() - 1) == 0 && (k - q) + 1 >= mod()) {
			ret = vf((k - q + 1 - mod()) / (mod() - 1), mod() - 1);
		}
		return ret + vf((k - q) / (mod() - 1), (k - q) % (mod() - 1));
	}
	static ModInt getStirlingSecond(uint64_t n, uint64_t k) {
		static vector<vector<ModInt>> seeds(mod(), vector<ModInt>(mod(), 0));
		static bool needCalc = true;
		if (needCalc) {
			needCalc = false;
			seeds[0][0] = 1;
			for (int i = 1; i < mod(); i++) {
				seeds[i][0] = 0;
				for (int j = 1; j < mod(); j++) {
					seeds[i][j] = seeds[i - 1][j - 1] + j * seeds[i - 1][j];
				}
			}
		}
		if (n < k) { return 0; }
		if (n < mod()) { return seeds[n][k]; }
		uint64_t q0 = k / mod(), r0 = k % mod();
		uint64_t q1 = (n - 1 - q0) / (mod() - 1), r1 = (n - 1 - q0) % (mod() - 1);
		ModInt ret = getCombinationOne(q1, q0) * seeds[r1 + 1][r0];
		if (r1 == mod() - 2 && q0 > 0) {
			ret += getCombinationOne(q1, q0 - 1) * seeds[0][r0];
		}
		return ret;
	}
	static ModInt getCombinationOne(uint64_t n, uint64_t k) {
		assert(is_prime(mod()));
		static vector<ModInt> fact, invf;
		static bool needCalc = true;
		if (needCalc) {
			needCalc = false;
			get_fact_table(mod() - 1, &fact, &invf);
		}
		ModInt ret = 1;
		while (n) {
			int n0 = n % mod(), k0 = k % mod();
			if (n0 < k0 || n0 < 0 || k0 < 0) { return 0; }
			ret *= fact[n0] * invf[k0] * invf[n0 - k0];
			n /= mod();
			k /= mod();
		}
		return ret;
	}
	static void convolve_xor_prod(vector<ModInt>* a,
									 const vector<ModInt>& b) {
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
		vector<ModInt> y = b;
		y.resize(n);
		fht(a);
		fht(&y);
		for (size_t i = 0; i < n; i++) { (*a)[i] *= y[i]; }
		fht(a);
		ModInt iv = (uint64_t)n;
		iv = iv.inv();
		for (size_t i = 0; i < n; i++) { (*a)[i] *= iv; }
	}
};

#endif /* ModInt_hpp */
