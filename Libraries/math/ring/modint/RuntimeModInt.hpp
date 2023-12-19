//
//  DynamicModInt.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef RuntimeModInt_hpp
#define RuntimeModInt_hpp

namespace dmint_internal {
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

using dmint_prim_base = dmint_internal::internal_primitive_bases;

class RuntimeModInt {
private:
	uint32_t x = 0;
	static constexpr int PRE_MOD_21 = 924844033;
	static constexpr int PRE_MOD_22 = 985661441;
	static constexpr int PRE_MOD_23 = 998244353;
	static constexpr int PRE_MOD_24 = 754974721;
	static constexpr int PRE_MOD_25 = 167772161;
	static constexpr int PRE_MOD_26 = 469762049;
private:
	static int mod_val() { return mod(); }
	static uint32_t umod_val() { return (uint32_t)mod(); }
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
	static void make_base(uint32_t m) {
		dmint_internal::primitive_bases[m] = make_unique<dmint_prim_base>(m);
	}
	static dmint_prim_base* get_bases(uint32_t m) {
		if (dmint_internal::primitive_bases.find(m) == dmint_internal::primitive_bases.end()) {
			make_base(m);
		}
		return dmint_internal::primitive_bases[m].get();
	}
	static int get_base(int m, int t, int idx) { return get_bases(m)->get_base(t, idx); }
	static uint32_t inv_raw(int md, uint32_t v) {
		if (v == 0) { return v; }
		uint32_t a = v, t = 1;
		while (a > 1) {
			t = (uint64_t)t * ((uint32_t)md - (uint32_t)md / a) % (uint32_t)md;
			a = (uint32_t)md % a;
		}
		return t;
	}
	static void ntt(int md, vector<uint32_t>* arr) {
		size_t n = arr->size();
		if (n <= 1) { return; }
		assert(!(n & (n - 1)));
		int m = 0;
		while ((1 << m) < n) { m++; }
		int p = 1, t = 1 << (m - 1);
		for (int i = 0; i < m; i++) {
			uint32_t g = 1;
			for (int j = 0; j < p; j++) {
				uint32_t u, v;
				int offset = j << (m - i);
				for (int k = 0; k < t; k++) {
					u = (*arr)[k + offset] % md;
					v = (uint64_t)(*arr)[k + offset + t] * g % (uint32_t)md;
					(*arr)[k + offset] = (u + v >= md) ? u + v - md : u + v;
					(*arr)[k + offset + t] = (u < v) ? md + u - v : u - v;
				}
				int w = (j + 1) ^ (j & (j + 1)), z = -1;
				while (w) {
					z++;
					w >>= 1;
				}
				g = (uint64_t)g * get_base(md, 2, z) % (uint32_t)md;
			}
			p <<= 1;
			t >>= 1;
		}
	}
	static void ntt(vector<RuntimeModInt>* arr) {
		size_t n = arr->size();
		if (n <= 1) { return; }
		assert(!(n & (n - 1)));
		int m = 0;
		while ((1 << m) < n) { m++; }
		int p = 1, t = 1 << (m - 1);
		for (int i = 0; i < m; i++) {
			RuntimeModInt g = 1;
			for (int j = 0; j < p; j++) {
				RuntimeModInt u, v;
				int offset = j << (m - i);
				for (int k = 0; k < t; k++) {
					u = (*arr)[k + offset];
					v = (*arr)[k + offset + t] * g;
					(*arr)[k + offset] = u + v;
					(*arr)[k + offset + t] = u - v;
				}
				int w = (j + 1) ^ (j & (j + 1)), z = -1;
				while (w) {
					z++;
					w >>= 1;
				}
				g *= get_base(umod_val(), 2, z);
			}
			p <<= 1;
			t >>= 1;
		}
	}
	static void intt(int md, vector<uint32_t>* arr) {
		size_t n = arr->size();
		if (n <= 1) { return; }
		uint32_t invn = inv_raw(md, (uint32_t)n);
		assert(!(n & (n - 1)));
		int m = 0;
		while ((1 << m) < n) { m++; }
		for (int i = 0; i < m; i++) {
			uint32_t g = 1, s = get_base(md, 1, i + 1);
			int offset = 1 << i;
			for (int k = 0; k < offset; k++) {
				for (int j = k; j < n; j += (1 << (i + 1))) {
					uint32_t x = (*arr)[j] % md;
					uint32_t y = (uint64_t)(*arr)[j + offset] * g % (uint32_t)md;
					(*arr)[j] = (x + y >= md) ? x + y - md : x + y;
					(*arr)[j + offset] = (x < y) ? md + x - y : x - y;
				}
				g = (uint64_t)g * s % (uint32_t)md;
			}
		}
		for (int i = 0; i < n; i++) {
			(*arr)[i] = (uint64_t)invn * (*arr)[i] % (uint32_t)md;
		}
	}
	static void intt(vector<RuntimeModInt>* arr) {
		size_t n = arr->size();
		if (n <= 1) { return; }
		RuntimeModInt invn = inv((int)n);
		assert(!(n & (n - 1)));
		int m = 0;
		while ((1 << m) < n) { m++; }
		for (int i = 0; i < m; i++) {
			RuntimeModInt g = 1, s = get_base(umod_val(), 1, i + 1);
			int offset = 1 << i;
			for (int k = 0; k < offset; k++) {
				for (int j = k; j < n; j += (1 << (i + 1))) {
					RuntimeModInt x = (*arr)[j], y = (*arr)[j + offset] * g;
					(*arr)[j] = x + y;
					(*arr)[j + offset] = x - y;
				}
				g *= s;
			}
		}
		for (int i = 0; i < n; i++) { (*arr)[i] *= invn; }
	}
	static void convolve_one_base_raw(int md, vector<uint32_t>* a,
									  const vector<uint32_t>& b) {
		if (a->size() == 0 || b.size() == 0) {
			*a = {};
			return;
		}
		size_t n = a->size(), m = b.size(), bas = 1;
		while (bas < n + m) { bas <<= 1; }
		a->resize(bas);
		vector<uint32_t> v(bas);
		for (size_t i = 0; i < m; i++) { v[i] = b[i]; }
		ntt(md, a);
		ntt(md, &v);
		for (size_t i = 0; i < bas; i++) { (*a)[i] = (uint64_t)(*a)[i] * v[i] % (uint32_t)md; }
		intt(md, a);
	}
	static bool is_prepared_base() {
		return mod_val() == PRE_MOD_21 || mod_val() == PRE_MOD_22 || mod_val() == PRE_MOD_23
		|| mod_val() == PRE_MOD_24 || mod_val() == PRE_MOD_25 || mod_val() == PRE_MOD_26;
	}
	static bool is_prepared() {
		if (is_prepared_base()) { return true; }
		if (dmint_internal::primitive_bases.find(mod_val()) != dmint_internal::primitive_bases.end()) {
			return true;
		}
		if (!is_prime(mod_val())) { return false; }
		int p = mod_val() - 1;
		for (int i = 0; i < 20; i++) {
			if (p & 1) { return false; }
			p >>= 1;
		}
		return true;
	}
public:
	static int& mod() {
		static int m = 0;
		return m;
	}
	RuntimeModInt() : x(0) {}
	template <typename T, enable_if_t<is_integral_v<T>>* = nullptr>
	RuntimeModInt(T v) { operator=(v); }
	RuntimeModInt(const RuntimeModInt& other) { operator=(other); }
	uint32_t val() const { return x; }
	RuntimeModInt inv() const {
		if (x == 0) { return RuntimeModInt(mod_val()); }
		uint32_t a = x, t = 1;
		while (a > 1) {
			t = (uint64_t)t * (umod_val() - umod_val() / a) % umod_val();
			a = umod_val() % a;
		}
		return RuntimeModInt(t);
	}
	RuntimeModInt extended_inv() const {
		if (x == 0) { return RuntimeModInt(0); }
		int a = 1, b = 0, c = 0, d = 1, k = x, l = mod_val();
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
	static RuntimeModInt inv(const RuntimeModInt& v) {
		if (v.x == 0) { return v; }
		uint32_t a = v.x, t = 1;
		while (a > 1) {
			t = (uint64_t)t * (umod_val() - umod_val() / a) % umod_val();
			a = umod_val() % a;
		}
		return RuntimeModInt(t);
	}
	static RuntimeModInt extended_inv(const RuntimeModInt& v) {
		return v.extended_inv();
	}
	template <typename T, enable_if_t<is_integral_v<T> && is_signed_v<T>>* = nullptr>
	RuntimeModInt& operator=(T aX) {
		x = (aX % mod_val() < 0) ? aX % mod_val() + umod_val() : aX % mod_val();
		return *this;
	}
	template <typename T, enable_if_t<is_integral_v<T> && is_unsigned_v<T>>* = nullptr>
	RuntimeModInt& operator=(T aX) {
		x = aX % umod_val();
		return *this;
	}
	RuntimeModInt& operator=(const RuntimeModInt& other) {
		x = other.x;
		return *this;
	}
	friend ostream& operator<<(ostream& os, const RuntimeModInt& other) {
		os << other.x;
		return os;
	}
	friend istream& operator>>(istream& ist, RuntimeModInt& other) {
		uint64_t a;
		ist >> a;
		other = a;
		return ist;
	}
	explicit operator bool() const { return x != 0; }
	bool operator!() const { return x == 0; }
	bool operator==(const RuntimeModInt& other) const { return x == other.x; }
	bool operator!=(const RuntimeModInt& other) const { return x != other.x; }
	RuntimeModInt operator+() const { return RuntimeModInt(*this); }
	RuntimeModInt operator-() const {
		RuntimeModInt ret = RuntimeModInt(*this);
		ret.x = (ret.x == 0) ? 0 : umod_val() - ret.x;
		return ret;
	}
	RuntimeModInt& operator++() {
		x = (x + 1 == umod_val()) ? 0 : x + 1;
		return *this;
	}
	RuntimeModInt operator++(int) {
		RuntimeModInt ret = RuntimeModInt(*this);
		++*this;
		return ret;
	}
	RuntimeModInt& operator--() {
		x = (x == 0) ? umod_val() - 1 : x - 1;
		return *this;
	}
	RuntimeModInt operator--(int) {
		RuntimeModInt ret = RuntimeModInt(*this);
		--*this;
		return ret;
	}
	RuntimeModInt& operator+=(const RuntimeModInt& other) {
		x = (x + other.x >= umod_val()) ? x + other.x - umod_val() : x + other.x;
		return *this;
	}
	RuntimeModInt& operator-=(const RuntimeModInt& other) {
		x = (x < other.x) ? x + umod_val() - other.x : x - other.x;
		return *this;
	}
	RuntimeModInt& operator*=(const RuntimeModInt& other) {
		x = (uint32_t)((uint64_t)x * other.x % umod_val());
		return *this;
	}
	RuntimeModInt& operator/=(const RuntimeModInt& other) { return *this *= inv(other); }
	RuntimeModInt operator+(const RuntimeModInt& other) const { return RuntimeModInt(*this) += other; }
	RuntimeModInt operator-(const RuntimeModInt& other) const { return RuntimeModInt(*this) -= other; }
	RuntimeModInt operator*(const RuntimeModInt& other) const { return RuntimeModInt(*this) *= other; }
	RuntimeModInt operator/(const RuntimeModInt& other) const { return RuntimeModInt(*this) /= other; }
	template <typename T>
	friend RuntimeModInt operator+(T a, const RuntimeModInt& other) { return RuntimeModInt(a) += other; }
	template <typename T>
	friend RuntimeModInt operator-(T a, const RuntimeModInt& other) { return RuntimeModInt(a) -= other; }
	template <typename T>
	friend RuntimeModInt operator*(T a, const RuntimeModInt& other) { return RuntimeModInt(a) *= other; }
	template <typename T>
	friend RuntimeModInt operator/(T a, const RuntimeModInt& other) { return RuntimeModInt(a) /= other; }
	template <typename T>
	static RuntimeModInt pow(const RuntimeModInt& n, T m) {
		RuntimeModInt ret = 1, v;
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
	RuntimeModInt operator^(T m) const { return pow(*this, m); }
	static int sqrt(const RuntimeModInt& n) {
		if (n.val() == 0 || n.val() == 1) { return n.val(); }
		if (pow(n, umod_val() >> 1) + 1 == 0) { return -1; }
		int m = mod_val() - 1, q = 0;
		while (!(m & 1)) {
			q++;
			m >>= 1;
		}
		RuntimeModInt z = 1, r = pow(n, (m + 1) >> 1), t = pow(n, m);
		while (pow(z, umod_val() >> 1) == 1) { z++; }
		z = pow(z, m);
		while (t != 1) {
			RuntimeModInt c = t;
			int i = 0;
			while (c != 1) {
				i++;
				c *= c;
			}
			RuntimeModInt b = pow(z, 1 << (q - i - 1));
			q = i;
			z = b * b;
			t *= z;
			r *= b;
		}
		return (r.x > mod_val() - r.x) ? mod_val() - r.x : r.x;
	}
	static void convolve_one_base(vector<RuntimeModInt>* a, const vector<RuntimeModInt>& b) {
		if (a->size() == 0 || b.size() == 0) {
			*a = {};
			return;
		}
		size_t n = a->size(), m = b.size(), bas = 1;
		while (bas < n + m) { bas <<= 1; }
		a->resize(bas);
		vector<RuntimeModInt> v(bas);
		for (size_t i = 0; i < m; i++) { v[i] = b[i]; }
		ntt(a);
		ntt(&v);
		for (size_t i = 0; i < bas; i++) { (*a)[i] *= v[i]; }
		intt(a);
	}
	static void convolve_one(vector<RuntimeModInt>* a, const vector<RuntimeModInt>& b) {
		if (is_prepared()) {
			convolve_one_base(a, b);
			return;
		}
		size_t n = a->size(), m = b.size();
		vector<uint32_t> x23(n);
		vector<uint32_t> x24(n);
		vector<uint32_t> x26(n);
		vector<uint32_t> y23(m);
		vector<uint32_t> y24(m);
		vector<uint32_t> y26(m);
		for (int i = 0; i < n; i++) {
			x23[i] = (*a)[i].val() % PRE_MOD_23;
			x24[i] = (*a)[i].val() % PRE_MOD_24;
			x26[i] = (*a)[i].val() % PRE_MOD_26;
		}
		for (int i = 0; i < m; i++) {
			y23[i] = b[i].val() % PRE_MOD_23;
			y24[i] = b[i].val() % PRE_MOD_24;
			y26[i] = b[i].val() % PRE_MOD_26;
		}
		convolve_one_base_raw(PRE_MOD_23, &x23, y23);
		convolve_one_base_raw(PRE_MOD_24, &x24, y24);
		convolve_one_base_raw(PRE_MOD_26, &x26, y26);
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
			(*a)[i] = rem(mod_val(), x26[i], x24[i], x23[i],
						  PRE_MOD_26, PRE_MOD_24, PRE_MOD_23);
		}
	}
	static void convolve(vector<RuntimeModInt>* a, const vector<RuntimeModInt>& b) {
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
			if (mod_val() == PRE_MOD_21) { bs = 19; }
			else if (mod_val() == PRE_MOD_22) { bs = 20; }
		}
		if (bas <= (1 << bs)) {
			convolve_one(a, b);
			return;
		}
		int lim = 1 << 18;
		int s = (n + lim - 1) / lim, t = (m + lim - 1) / lim;
		a->resize(s * lim);
		vector<RuntimeModInt> y(t * lim);
		for (int i = 0; i < m; i++) { y[i] = b[i]; }
		vector<RuntimeModInt> ret((s + t) * lim);
		for (int i = 0; i < s; i++) {
			vector<RuntimeModInt> u(lim);
			for (int j = 0; j < lim; j++) { u[j] = (*a)[lim * i + j]; }
			for (int j = 0; j < t; j++) {
				vector<RuntimeModInt> res(lim);
				for (int k = 0; k < lim; k++) { res[k] = y[lim * j + k]; }
				convolve_one(&u, res);
				for (int k = 0; k < u.size(); k++) { ret[(i + j) * lim + k] += u[k]; }
			}
		}
		*a = ret;
	}
	static void get_fact_table(int n, vector<RuntimeModInt>* fact,
							   vector<RuntimeModInt>* invf) {
		assert(n >= 0);
		RuntimeModInt t = 1;
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
	static void get_inv_table(int n, vector<RuntimeModInt>* invn) {
		assert(invn && n >= 0);
		invn->resize(n + 1);
		for (int i = 0; i < 2; i++) {
			(*invn)[i] = 1;
			if (n == i) { return; }
		}
		for (int i = 2; i <= n; i++) {
			int t = umod_val() % i;
			(*invn)[i] = (umod_val() - umod_val() / i) * (*invn)[t];
		}
	}
	static vector<RuntimeModInt> getMontmort(size_t n) {
		if (n == 0) { return {}; }
		vector<RuntimeModInt> ret(n);
		ret[0] = 0;
		for (size_t i = 1; i < n; i++) {
			ret[i] = ret[i - 1] * (i + 1) + ((i & 1) ? 1 : -1);
		}
		return ret;
	}
	static vector<RuntimeModInt> getPowerTable(size_t n, uint32_t k) {
		vector<RuntimeModInt> ret(n + 1, 1);
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
	static RuntimeModInt getStirlingFirst(uint64_t n, uint64_t k) {
		static vector<vector<RuntimeModInt>> seeds(mod(), vector<RuntimeModInt>(mod(), 0));
		static int preCalcMod = mod();
		if (preCalcMod != mod()) {
			preCalcMod = mod();
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
		RuntimeModInt ret = 0;
		function<RuntimeModInt(uint64_t, uint64_t)> vf = [n, k, q](uint64_t l, uint64_t m) {
			RuntimeModInt s = ((q - l) & 1) ? -1 : 1;
			return getCombinationOne(q, l) * s * seeds[n % mod()][m];
		};
		if (n % mod() == -1 && (k - q) % (mod() - 1) == 0 && (k - q) + 1 >= mod()) {
			ret = vf((k - q + 1 - mod()) / (mod() - 1), mod() - 1);
		}
		return ret + vf((k - q) / (mod() - 1), (k - q) % (mod() - 1));
	}
	static RuntimeModInt getStirlingSecond(uint64_t n, uint64_t k) {
		static vector<vector<RuntimeModInt>> seeds(mod(), vector<RuntimeModInt>(mod(), 0));
		static int preCalcMod = mod();
		if (preCalcMod != mod()) {
			preCalcMod = mod();
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
		if (r1 < mod() - 2) {
			return getCombinationOne(q1, q0) * seeds[1 + r1][r0];
		} else if (r1 == mod() - 2) {
			if (r1 != 0) {
				return getCombinationOne(q1, q0 - 1);
			} else {
				return getCombinationOne(q1, q0) * seeds[mod() - 1][r1];
			}
		}
		assert(false);
		return 0;
	}
	static RuntimeModInt getCombinationOne(uint64_t n, uint64_t k) {
		assert(is_prime(mod()));
		static vector<RuntimeModInt> fact, invf;
		static int preCalcMod = mod();
		if (preCalcMod != mod()) {
			preCalcMod = mod();
			get_fact_table(mod(), &fact, &invf);
		}
		RuntimeModInt ret = 1;
		while (n) {
			int n0 = n % mod(), k0 = k % mod();
			ret *= fact[n0] * invf[k0] * invf[n0 - k0];
			n /= mod();
			k /= mod();
		}
		return ret;
	}
};

#endif /* RuntimeModInt_hpp */
