//
//  BigInteger.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef BigInteger_hpp
#define BigInteger_hpp

namespace bint_internal {
static constexpr int MOD_0 = 880803841;
static constexpr int MOD_1 = 897581057;
static constexpr int MOD_2 = 998244353; // 998244353 = 119 * 2 ^ 23 + 1
static constexpr int BASE_SIZE = 23;
static constexpr int get_inv_primitive(int a) {
	for (int i = 1; i < a; i++) {
		int b = i;
		for (int j = 0; j < BASE_SIZE; j++) {
			if (b == 1) {
				b = 0;
				break;
			}
			b = (int64_t)b * b % a;
		}
		if (b == 1) { return i; }
	}
	return a;
}
static constexpr int get_primitive(int a) {
	int b = get_inv_primitive(a);
	int t = a - 2;
	int ret = 1;
	int bas = b;
	while (t) {
		if (t & 1) { ret = (int64_t)ret * bas % a; }
		bas = (int64_t)bas * bas % a;
		t >>= 1;
	}
	return ret;
}
static constexpr int PRI_0 = get_primitive(MOD_0);
static constexpr int PRI_1 = get_primitive(MOD_1);
static constexpr int PRI_2 = get_primitive(MOD_2);
static constexpr int INV_0 = get_inv_primitive(MOD_0);
static constexpr int INV_1 = get_inv_primitive(MOD_1);
static constexpr int INV_2 = get_inv_primitive(MOD_2);
static int PRS_0[BASE_SIZE + 1];
static int PRS_1[BASE_SIZE + 1];
static int PRS_2[BASE_SIZE + 1];
static int IRS_0[BASE_SIZE + 1];
static int IRS_1[BASE_SIZE + 1];
static int IRS_2[BASE_SIZE + 1];
static int CRS_0[BASE_SIZE - 1];
static int CRS_1[BASE_SIZE - 1];
static int CRS_2[BASE_SIZE - 1];
template <int M>
struct internal_mint {
private:
	uint32_t x;
private:
	static void make_base() {
		PRS_0[BASE_SIZE] = PRI_0;
		PRS_1[BASE_SIZE] = PRI_1;
		PRS_2[BASE_SIZE] = PRI_2;
		IRS_0[BASE_SIZE] = INV_0;
		IRS_1[BASE_SIZE] = INV_1;
		IRS_2[BASE_SIZE] = INV_2;
		for (int i = BASE_SIZE; i > 0; i--) {
			PRS_0[i - 1] = (int64_t)PRS_0[i] * PRS_0[i] % MOD_0;
			PRS_1[i - 1] = (int64_t)PRS_1[i] * PRS_1[i] % MOD_1;
			PRS_2[i - 1] = (int64_t)PRS_2[i] * PRS_2[i] % MOD_2;
			IRS_0[i - 1] = (int64_t)IRS_0[i] * IRS_0[i] % MOD_0;
			IRS_1[i - 1] = (int64_t)IRS_1[i] * IRS_1[i] % MOD_1;
			IRS_2[i - 1] = (int64_t)IRS_2[i] * IRS_2[i] % MOD_2;
		}
		int ie0 = 1, ie1 = 1, ie2 = 1;
		for (int i = 0; i < BASE_SIZE - 1; i++) {
			CRS_0[i] = (int64_t)ie0 * PRS_0[i + 2] % MOD_0;
			CRS_1[i] = (int64_t)ie1 * PRS_1[i + 2] % MOD_1;
			CRS_2[i] = (int64_t)ie2 * PRS_2[i + 2] % MOD_2;
			ie0 = (int64_t)ie0 * IRS_0[i + 2] % MOD_0;
			ie1 = (int64_t)ie1 * IRS_1[i + 2] % MOD_1;
			ie2 = (int64_t)ie2 * IRS_2[i + 2] % MOD_2;
		}
	}
	static void select_roots(const int** pr, const int** iv, int** prs, int** irs, int** crs) {
		if (M == MOD_0) {
			*pr = &PRI_0;
			*iv = &INV_0;
			*prs = &PRS_0[0];
			*irs = &IRS_0[0];
			*crs = &CRS_0[0];
		} else if (M == MOD_1) {
			*pr = &PRI_1;
			*iv = &INV_1;
			*prs = &PRS_1[0];
			*irs = &IRS_1[0];
			*crs = &CRS_1[0];
		} else if (M == MOD_2) {
			*pr = &PRI_2;
			*iv = &INV_2;
			*prs = &PRS_2[0];
			*irs = &IRS_2[0];
			*crs = &CRS_2[0];
		} else {
			*pr = nullptr;
			*iv = nullptr;
			*prs = nullptr;
			*irs = nullptr;
			*crs = nullptr;
		}
	}
	static void ntt(vector<internal_mint>* a) {
		int n = (int)a->size();
		if (a->size() <= 1) { return; }
		static_assert(M == MOD_0 || M == MOD_1 || M == MOD_2,
					  "M is not NTT-friendly modulo");
		static bool need_base = true;
		if (need_base) {
			need_base = false;
			make_base();
		}
		const int* pr;
		const int* iv;
		int* prs;
		int* irs;
		int* crs;
		select_roots(&pr, &iv, &prs, &irs, &crs);
		assert(!(n & (n - 1)));
		int m = 0;
		while ((1 << m) < n) { m++; }
		int p = 1, t = 1 << (m - 1);
		for (int i = 0; i < m; i++) {
			internal_mint g = 1;
			for (int j = 0; j < p; j++) {
				internal_mint u, v;
				int offset = j << (m - i);
				for (int k = 0; k < t; k++) {
					u = (*a)[k + offset];
					v = (*a)[k + offset + t] * g;
					(*a)[k + offset] = u + v;
					(*a)[k + offset + t] = u - v;
				}
				int w = (j + 1) ^ (j & (j + 1)), z = -1;
				while (w) {
					z++;
					w >>= 1;
				}
				g *= crs[z];
			}
			p <<= 1;
			t >>= 1;
		}
	}
	static void intt(vector<internal_mint>* a) {
		int n = (int)a->size();
		if (n <= 1) { return; }
		static_assert(is_prepared(), "M is not NTT-friendly modulo");
		static bool need_base = true;
		if (need_base) {
			need_base = false;
			make_base();
		}
		const int* pr;
		const int* iv;
		int* prs = nullptr;
		int* irs = nullptr;
		int* crs = nullptr;
		select_roots(&pr, &iv, &prs, &irs, &crs);
		internal_mint invn = inv(n);
		assert(!(n & (n - 1)));
		int m = 0;
		while ((1 << m) < n) { m++; }
		for (int i = 0; i < m; i++) {
			internal_mint g = 1, s = irs[i + 1];
			int offset = 1 << i;
			for (int k = 0; k < offset; k++) {
				for (int j = k; j < n; j += (1 << (i + 1))) {
					internal_mint x = (*a)[j], y = (*a)[j + offset] * g;
					(*a)[j] = x + y;
					(*a)[j + offset] = x - y;
				}
				g *= s;
			}
		}
		for (int i = 0; i < n; i++) { (*a)[i] *= invn; }
	}
	static constexpr int is_prepared() {
		return M == MOD_0 || M == MOD_1 || M == MOD_2;
	}
public:
	static constexpr int mod() { return M; }
	static constexpr uint32_t umod() { return M; }
	internal_mint() : x(0) {}
	template <typename T>
	internal_mint(T v) { operator=(v); }
	template <typename T, enable_if_t<is_integral_v<T> && is_signed_v<T>>* = nullptr>
	internal_mint& operator=(T v) {
		x = (v % mod() < 0) ? mod() - v % mod() : v % mod();
		return *this;
	}
	template <typename T, enable_if_t<is_integral_v<T> && is_unsigned_v<T>>* = nullptr>
	internal_mint& operator=(T v) {
		x = v % umod();
		return *this;
	}
	internal_mint& operator=(const internal_mint& other) {
		x = other.x;
		return *this;
	}
	explicit operator bool() const { return x != 0; }
	bool operator!() const { return x == 0; }
	bool operator==(const internal_mint& other) const { return x == other.x; }
	bool operator!=(const internal_mint& other) const { return x != other.x; }
	internal_mint operator+() const { return internal_mint(*this); }
	internal_mint operator-() const {
		internal_mint ret = *this;
		ret.x = (ret.x > 0) ? umod() - ret.x : 0;
		return ret;
	}
	internal_mint& operator++() {
		if (++x >= umod()) { x = 0; }
		return *this;
	}
	internal_mint operator++(int) {
		internal_mint ret = *this;
		++*this;
		return ret;
	}
	internal_mint& operator--() {
		x = (x == 0) ? umod() - 1 : x - 1;
		return *this;
	}
	internal_mint operator--(int) {
		internal_mint ret = *this;
		--*this;
		return ret;
	}
	internal_mint& operator+=(const internal_mint& other) {
		x = (x + other.x >= umod()) ? x + other.x - umod() : x + other.x;
		return *this;
	}
	internal_mint& operator-=(const internal_mint& other) {
		x = (x < other.x) ? x + umod() - other.x : x - other.x;
		return *this;
	}
	internal_mint& operator*=(const internal_mint& other) {
		x = (uint64_t)x * other.x % umod();
		return *this;
	}
	internal_mint& operator/=(const internal_mint& other) {
		return *this *= other.inv;
	}
	internal_mint operator+(const internal_mint& other) const {
		return internal_mint(*this) += other;
	}
	internal_mint operator-(const internal_mint& other) const {
		return internal_mint(*this) -= other;
	}
	internal_mint operator*(const internal_mint& other) const {
		return internal_mint(*this) *= other;
	}
	internal_mint operator/(const internal_mint& other) const {
		return internal_mint(*this) /= other;
	}
	friend ostream& operator<<(ostream& os, const internal_mint& mint) {
		os << mint.x;
		return os;
	}
	friend istream& operator<<(istream& ist, internal_mint& mint) {
		ist >> mint.x;
		return ist;
	}
	uint32_t val() const { return x; }
	internal_mint inv() const { return inv(*this); }
	static internal_mint inv(const internal_mint& mint) {
		if (mint.x == 0) { return mint; }
		uint32_t a = mint.x, t = 1;
		while (a > 1) {
			t = (uint64_t)t * (mod() - mod() / a) % mod();
			a = mod() % a;
		}
		return internal_mint(t);
	}
	static void convolve_one(vector<internal_mint>* a, const vector<internal_mint>& b) {
		size_t n = a->size(), m = b.size();
		size_t bs = 1;
		while (bs < n + m) { bs <<= 1; }
		a->resize(bs);
		vector<internal_mint> y(bs);
		for (size_t i = 0; i < m; i++) { y[i] = b[i]; }
		ntt(a);
		ntt(&y);
		for (size_t i = 0; i < bs; i++) { (*a)[i] *= y[i]; }
		intt(a);
	}
	static void convolve_prepared(vector<internal_mint>* a,
								  const vector<internal_mint>& b) {
		int n = (int)a->size(), m = (int)b.size();
		if (n == 0 || m == 0) {
			*a = {};
			return;
		}
		if (n + m <= 256) {
			vector<internal_mint> temp = *a;
			a->clear();
			a->resize(n + m - 1);
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) { (*a)[i + j] += temp[i] * b[j]; }
			}
			return;
		}
		int bs = 1;
		while (bs < n + m) { bs <<= 1; }
		if (bs < (1 << (BASE_SIZE - 2))) {
			convolve_one(a, b);
			return;
		}
		int lim = 1 << (BASE_SIZE - 5);
		int s = (n + lim - 1) / lim, t = (m + lim - 1) / lim;
		vector<internal_mint> ret((s + t + 1) * lim);
		for (int i = 0; i < s; i++) {
			vector<internal_mint> u(lim);
			int p = i * lim;
			for (int j = 0; j < lim; j++) { u[j] = (p >= n) ? 0 : (*a)[p++]; }
			for (int j = 0; j < t; j++) {
				vector<internal_mint> v(lim);
				int q = j * lim;
				for (int k = 0; k < lim; k++) { v[k] = (q >= m) ? 0 : b[q++]; }
				convolve_one(&v, u);
				int r = (i + j) * lim;
				for (int i = 0; i < v.size(); i++) { ret[r + i] += v[i]; }
			}
		}
		*a = ret;
	}
};
} // namespace bint_internal

class BigInteger {
	template <int M> using mint = bint_internal::internal_mint<M>;
	using mint0 = mint<bint_internal::MOD_0>;
	using mint1 = mint<bint_internal::MOD_1>;
	using mint2 = mint<bint_internal::MOD_2>;
private:
	enum state_type : uint8_t {
		state_valid = 0,
		state_infinity,
		state_invalid,
		state_count,
	};
	static constexpr int DIGITS_BASE = 100000000;
	static constexpr int DIGITS_NUM = 8;
	static constexpr int BASE = 10;
	state_type state = state_valid;
	bool is_negative = false;
	vector<uint32_t> digits;
private:
	void set_invalid() {
		state = state_invalid;
		digits.clear();
		is_negative = false;
	}
	void shrink_digits() {
		if (state != state_valid) {
			digits.clear();
			return;
		}
		shrink_vector(&digits);
		if (digits.size() == 0) { is_negative = false; }
	}
	static void shrink_vector(vector<uint32_t>* v) {
		size_t s = v->size();
		while (s > 0) {
			if ((*v)[s - 1]) { break; }
			s--;
		}
		v->resize(s);
	}
	template <typename T>
	void assign_value(T v) {
		state = state_valid;
		digits.clear();
		while (v) {
			digits.push_back(v % DIGITS_BASE);
			v /= DIGITS_BASE;
		}
	}
	static void increment_digits(vector<uint32_t>* v) {
		size_t p = 0;
		while (true) {
			if (p >= v->size()) { v->resize(p + 1); }
			(*v)[p]++;
			if ((*v)[p] < DIGITS_BASE) { break; }
			(*v)[p] = 0;
			p++;
		}
	}
	static void decrement_digits(vector<uint32_t>* v) {
		size_t p = 0;
		assert(v->size() > 0);
		while (true) {
			if ((*v)[p]) {
				(*v)[p]--;
				break;
			}
			(*v)[p] = DIGITS_BASE - 1;
			p++;
		}
	}
	static void add_digits(vector<uint32_t>* u, const vector<uint32_t>& v, size_t st) {
		size_t n = u->size(), m = v.size();
		size_t ml = std::max(n, m + st);
		if (u->size() < ml) { u->resize(ml); }
		for (size_t i = 0; i < m; i++) {
			if (i + st >= u->size()) { u->resize(i + st + 1); }
			(*u)[i + st] += (i >= m) ? 0 : v[i];
			if ((*u)[i + st] >= DIGITS_BASE) {
				if (i + st + 1 >= u->size()) { u->resize(i + st + 2); }
				(*u)[i + st + 1]++;
				(*u)[i + st] -= DIGITS_BASE;
			}
		}
		size_t p = m + st;
		while (p < u->size() && (*u)[p] >= DIGITS_BASE) {
			if (p + 1 >= u->size()) { u->resize(p + 2); }
			(*u)[p + 1]++;
			(*u)[p] -= DIGITS_BASE;
			p++;
		}
	}
	static bool subtract_digits(vector<uint32_t>* u, const vector<uint32_t>& v, size_t st) {
		size_t n = u->size(), m = v.size();
		size_t ml = std::max(n, m + st);
		u->resize(ml);
		for (size_t i = 0; i < ml; i++) {
			if (i + st >= u->size()) { u->resize(i + st + 1); }
			(*u)[i + st] += (i >= m) ? 0 : DIGITS_BASE - 1 - v[i];
			if ((*u)[i + st] >= DIGITS_BASE) {
				if (i + st == ml - 1) { u->push_back(0); }
				(*u)[i + st + 1]++;
				(*u)[i + st] -= DIGITS_BASE;
			}
		}
		increment_digits(u);
		shrink_vector(u);
		if (u->size() > m + st) {
			size_t p = m + st;
			while (true) {
				if ((*u)[p]) {
					(*u)[p]--;
					break;
				}
				(*u)[p] = DIGITS_BASE - 1;
				++p;
			}
			shrink_vector(u);
			return false;
		} else {
			u->resize(m + st);
			for (size_t i = 0; i < m + st; i++) { (*u)[i] = DIGITS_BASE - 1 - (*u)[i]; }
			increment_digits(u);
			shrink_vector(u);
			return true;
		}
	}
	static vector<uint32_t> restore_digits(const mint0& v0, const mint1& v1, const mint2& v2) {
		static bool is_calc = false;
		static mint0 di0;
		static mint1 di1, m01;
		static mint2 di2, m02, m12;
		if (!is_calc) {
			is_calc = true;
			di0 = mint0::inv(DIGITS_BASE);
			di1 = mint1::inv(DIGITS_BASE);
			di2 = mint2::inv(DIGITS_BASE);
			m01 = mint1::inv(bint_internal::MOD_0);
			m02 = mint2::inv(bint_internal::MOD_0);
			m12 = mint2::inv(bint_internal::MOD_1);
		}
		vector<uint32_t> temp;
		mint0 a0 = v0;
		mint1 a1 = v1;
		mint2 a2 = v2;
		while (a0 != 0 || a1 != 0 || a2 != 0) {
			mint0 u0 = a0;
			mint1 u1 = (a1 - u0.val()) * m01;
			mint2 u2 = (a2 - u0.val() - mint2(u1.val()) * bint_internal::MOD_0) * m02 * m12;
			uint32_t z = (uint32_t)(u2.val() % DIGITS_BASE);
			z = (uint32_t)(((uint64_t)z * bint_internal::MOD_1 + u1.val()) % DIGITS_BASE);
			z = (uint32_t)(((uint64_t)z * bint_internal::MOD_0 + u0.val()) % DIGITS_BASE);
			temp.push_back(z);
			a0 = di0 * (a0 - z);
			a1 = di1 * (a1 - z);
			a2 = di2 * (a2 - z);
		}
		return temp;
	}
	static void carry_convolve_naive(vector<uint32_t>* a, const vector<uint32_t>& b) {
		size_t n = a->size(), m = b.size();
		vector<mint<DIGITS_BASE>> temp(n);
		for (size_t i = 0; i < n; i++) { temp[i] = (*a)[i]; }
		a->clear();
		a->resize(n + m - 1);
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < m; j++) {
				uint64_t t = (uint64_t)temp[i].val() * b[j];
				size_t k = i + j;
				while (t) {
					if (k >= a->size()) { a->resize(k + 1); }
					(*a)[k] += t % DIGITS_BASE;
					if ((*a)[k] >= DIGITS_BASE) {
						if (k + 1 >= a->size()) { a->resize(k + 2); }
						(*a)[k] -= DIGITS_BASE;
						(*a)[k + 1]++;
					}
					t /= DIGITS_BASE;
					k++;
				}
				if (k < a->size() && (*a)[k] >= DIGITS_BASE) {
					if (k + 1 >= a->size()) { a->resize(k + 2); }
					(*a)[k] -= DIGITS_BASE;
					(*a)[k + 1]++;
				}
			}
		}
	}
	static void carry_convolve(vector<uint32_t>* a, const vector<uint32_t>& b) {
		size_t n = a->size(), m = b.size();
		if (n == 0 || m == 0) {
			a->clear();
			return;
		}
		if (n + m <= 256) {
			carry_convolve_naive(a, b);
			return;
		}
		vector<mint0> x0(n), y0(m);
		vector<mint1> x1(n), y1(m);
		vector<mint2> x2(n), y2(m);
		for (size_t i = 0; i < n; i++) {
			x0[i] = (*a)[i];
			x1[i] = (*a)[i];
			x2[i] = (*a)[i];
		}
		for (size_t i = 0; i < m; i++) {
			y0[i] = b[i];
			y1[i] = b[i];
			y2[i] = b[i];
		}
		mint0::convolve_prepared(&x0, y0);
		mint1::convolve_prepared(&x1, y1);
		mint2::convolve_prepared(&x2, y2);
		a->clear();
		a->reserve(n + m + 19);
		for (size_t i = 0; i < x0.size(); i++) {
			vector<uint32_t> temp = restore_digits(x0[i], x1[i], x2[i]);
			if (temp.size() == 0) { continue; }
			add_digits(a, temp, i);
		}
	}
	static void divide_digits(vector<uint32_t>* a, const vector<uint32_t>& b) {
		assert(b.size() != 0);
		if (a->size() == 0) { return; }
		size_t n = a->size(), m = b.size();
		vector<uint32_t> y(n);
		for (uint32_t i = 0; i < n; i++) { y[i] = (*a)[i]; }
		size_t da = DIGITS_NUM * (n - 1);
		size_t db = DIGITS_NUM * (m - 1);
		uint32_t la = (*a)[n - 1], lb = b[m - 1];
		while (la) {
			da++;
			la /= BASE;
		}
		while (lb) {
			db++;
			lb /= BASE;
		}
		if (da < db) {
			a->clear();
			return;
		}
		uint32_t iv = divide_naive(b, DIGITS_NUM);
		if (iv == DIGITS_BASE) {
			divide_vector_ten(a, db - 1);
			return;
		}
		vector<uint32_t> x;
		x.push_back(iv);
		size_t e = DIGITS_NUM, f = DIGITS_NUM;
		vector<uint32_t> bc(m);
		for (uint32_t i = 0; i < m; i++) { bc[i] = b[i]; }
		while (((f + 1) >> 1) + db < da) {
			f = (f << 1) - 1;
			size_t g = (f + 1 > db) ? f + 1 - db : 0;
			vector<uint32_t> bx(x.size());
			for (size_t i = 0; i < x.size(); i++) { bx[i] = x[i]; }
			if (f + 1 < db) { divide_vector_ten(&bc, db - f - 1); }
			carry_convolve(&bx, bc);
			if (f + 1 < db) {
				bc.resize(m);
				for (uint32_t i = 0; i < m; i++) { bc[i] = b[i]; }
			}
			size_t qf = (e + f - g) / DIGITS_NUM;
			size_t mf = (e + f - g) % DIGITS_NUM;
			uint32_t t = 2;
			for (uint32_t i = 0; i < mf; i++) { t *= BASE; }
			bx.resize(qf + 1);
			for (size_t i = 0; i < qf; i++) { bx[i] = DIGITS_BASE - 1 - bx[i]; }
			bx[qf] = t - 1 - bx[qf];
			increment_digits(&bx);
			carry_convolve(&x, bx);
			e = (e << 1) + f - g;
			divide_vector_ten(&x, e - f);
			e = f;
		}
		carry_convolve(a, x);
		divide_vector_ten(a, e + db - 1);
		shrink_vector(a);
		vector<uint32_t> bx(a->size());
		for (size_t i = 0; i < a->size(); i++) { bx[i] = (*a)[i]; }
		carry_convolve(&bx, b);
		shrink_vector(&bx);
		if (bx.size() > y.size()) {
			decrement_digits(a);
			shrink_vector(a);
			return;
		} else if (bx.size() == y.size()) {
			bool f = false;
			for (size_t i = bx.size(); i > 0; i--) {
				if (y[i - 1] != bx[i - 1]) {
					f = y[i - 1] < bx[i - 1];
					break;
				}
			}
			if (f) {
				decrement_digits(a);
				shrink_vector(a);
				return;
			}
		}
		add_digits(&bx, b, 0);
		if (y.size() > bx.size()) { increment_digits(a); }
		else if (y.size() == bx.size()) {
			bool f = true;
			for (size_t i = bx.size(); i > 0; i--) {
				if (y[i - 1] != bx[i - 1]) {
					f = y[i - 1] > bx[i - 1];
					break;
				}
			}
			if (f) { increment_digits(a); }
		}
	}
	static uint32_t divide_naive(const vector<uint32_t>& b, uint32_t prec) {
		assert(1 <= prec && prec <= DIGITS_NUM);
		size_t m = b.size();
		assert(m > 0);
		uint32_t lim = 1;
		uint32_t v = 1;
		while (v * BASE <= b[m - 1]) { v *= BASE; }
		for (uint32_t i = 0; i < prec; i++) { lim *= BASE; }
		uint32_t l = 0, r = lim + 1;
		uint32_t d = lim / 2;
		uint64_t t = (uint64_t)lim * v;
		while (r - l > 1) {
			uint32_t c = 0;
			bool flg = true;
			for (size_t i = 0; i < m - 1; i++) {
				uint64_t x = (uint64_t)b[i] * d + c;
				c = (uint32_t)(x / DIGITS_BASE);
				flg = (x % DIGITS_BASE > 0) ? false : flg;
			}
			uint64_t val = (uint64_t)b[m - 1] * d + c;
			if (val < t || (val == t && flg)) { l = d; }
			else { r = d; }
			d = l + (r - l) / 2;
		}
		return d;
	}
	static void multiple_vector_ten(vector<uint32_t>* v, size_t n) {
		size_t m = v->size();
		if (m == 0 || n == 0) { return; }
		size_t bs = n % DIGITS_NUM, gs = n / DIGITS_NUM;
		v->resize(m + gs + 1);
		uint32_t b = 1;
		for (size_t i = 0; i < bs; i++) { b *= BASE; }
		uint32_t cb = DIGITS_BASE / b;
		for (size_t i = m; i > 0; i--) {
			(*v)[i + gs] += (*v)[i - 1] / cb;
			(*v)[i + gs - 1] = (*v)[i - 1] % cb * b;
		}
		for (size_t i = gs; i > 0; i--) { (*v)[i - 1] = 0; }
		shrink_vector(v);
	}
	static void divide_vector_ten(vector<uint32_t>* v, size_t n) {
		size_t m = v->size();
		if (m == 0) { return; }
		size_t bs = n % DIGITS_NUM, gs = n / DIGITS_NUM;
		uint32_t bs_base = 1;
		for (size_t i = 0; i < bs; i++) { bs_base *= BASE; }
		uint32_t bs_c_base = DIGITS_BASE / bs_base;
		if (gs > 0) {
			for (size_t i = 0; i + gs < m; i++) {
				(*v)[i] = (*v)[i + gs];
			}
			if (m >= gs) { v->resize(m - gs); }
			else { v->clear(); }
		}
		if (bs > 0 && v->size() > 0) {
			(*v)[0] /= bs_base;
			for (size_t i = 1; i + gs < m; i++) {
				(*v)[i - 1] += (*v)[i] % bs_base * bs_c_base;
				(*v)[i] /= bs_base;
			}
		}
		shrink_vector(v);
	}
public:
	BigInteger() {}
	template <typename T>
	BigInteger(T other) { operator=(other); }
	BigInteger& operator=(bool b) {
		state = state_valid;
		is_negative = false;
		digits.clear();
		if (b) { digits.push_back(1); }
		return *this;
	}
	BigInteger& operator=(int8_t v) {
		uint8_t r = 0;
		if (v < 0) {
			is_negative = true;
			if (v == INT8_MIN) { r = (uint8_t)(-(v + 1)) + 1; }
			else { r = (uint8_t)-v; }
		} else { r = (uint8_t)v; }
		return operator=(r);
	}
	BigInteger& operator=(int16_t v) {
		uint16_t r = 0;
		if (v < 0) {
			is_negative = true;
			if (v == INT16_MIN) { r = (uint16_t)(-(v + 1)) + 1; }
			else { r = (uint16_t)-v; }
		} else { r = (uint16_t)v; }
		return operator=(r);
	}
	BigInteger& operator=(int32_t v) {
		uint32_t r = 0;
		if (v < 0) {
			is_negative = true;
			if (v == INT32_MIN) { r = (uint32_t)(-(v + 1)) + 1; }
			else { r = (uint32_t)-v; }
		} else { r = (uint32_t)v; }
		return operator=(r);
	}
	BigInteger& operator=(int64_t v) {
		uint64_t r = 0;
		if (v < 0) {
			is_negative = true;
			if (v == INT64_MIN) { r = (uint64_t)(-(v + 1)) + 1; }
			else { r = (uint64_t)-v; }
		} else { r = (uint64_t)v; }
		return operator=(r);
	}
	BigInteger& operator=(uint8_t v) {
		assign_value(v);
		return *this;
	}
	BigInteger& operator=(uint16_t v) {
		assign_value(v);
		return *this;
	}
	BigInteger& operator=(uint32_t v) {
		assign_value(v);
		return *this;
	}
	BigInteger& operator=(uint64_t v) {
		assign_value(v);
		return *this;
	}
	BigInteger& operator=(const string& s) {
		if (s[0] == '-') { is_negative = true; }
		digits.clear();
		state = state_valid;
		uint32_t b = 1;
		digits.resize(s.length() / 2 + 1);
		size_t p = 0;
		for (size_t i = s.length(); i > (is_negative) ? 1 : 0; i--) {
			if (s[i - 1] < '0' || '9' < s[i - 1]) {
				state = state_invalid;
				return *this;
			}
			digits[p] += b * (s[i - 1] ^ '0');
			b *= 10;
			if (b >= DIGITS_BASE) {
				p++;
				b = 1;
			}
		}
		shrink_digits();
		return *this;
	}
	BigInteger& operator=(const BigInteger& other) {
		state = other.state;
		is_negative = other.is_negative;
		digits = other.digits;
		return *this;
	}
	explicit operator bool() const { return is_inf() || (!is_nan() && !digits.empty()); }
	bool operator!() const { return !is_inf() && (is_nan() || digits.empty()); }
	bool operator==(const BigInteger& other) const {
		if (is_nan() && other.is_nan()) { return true; }
		if (is_inf() && other.is_inf()) { return is_negative == other.is_negative; }
		return state == other.state && is_negative == other.is_negative && digits == other.digits;
	}
	bool operator!=(const BigInteger& other) const {
		if (is_nan() && other.is_nan()) { return false; }
		if (is_inf() && other.is_inf()) { return is_negative != other.is_negative; }
		return state != other.state || is_negative != other.is_negative || digits != other.digits;
	}
	bool operator<=(const BigInteger& other) const {
		if (is_nan() || other.is_nan()) { return false; }
		if (is_inf() && other.is_inf()) { return is_negative || !other.is_negative; }
		if (is_inf() && !other.is_inf()) { return is_negative; }
		if (!is_inf() && other.is_inf()) { return !other.is_negative; }
		if (is_negative != other.is_negative) { return is_negative; }
		size_t ts = digits.size(), os = other.digits.size();
		if (ts != os) { return (is_negative && ts > os) || (!is_negative && ts < os); }
		for (size_t i = ts; i > 0; i--) {
			uint32_t ln = digits[i - 1], rn = other.digits[i - 1];
			if (ln != rn) { return (is_negative && ln > rn) || (!is_negative && ln < rn); }
		}
		return true;
	}
	bool operator<(const BigInteger& other) const {
		if (is_nan() || other.is_nan()) { return false; }
		if (is_inf() && other.is_inf()) { return is_negative && !other.is_negative; }
		if (is_inf() && !other.is_inf()) { return is_negative; }
		if (!is_inf() && other.is_inf()) { return !other.is_negative; }
		if (is_negative != other.is_negative) { return is_negative; }
		size_t ts = digits.size(), os = other.digits.size();
		if (ts != os) { return (is_negative && ts > os) || (!is_negative && ts < os); }
		for (size_t i = ts; i > 0; i--) {
			uint32_t ln = digits[i - 1], rn = other.digits[i - 1];
			if (ln != rn) { return (is_negative && ln > rn) || (!is_negative && ln < rn); }
		}
		return false;
	}
	bool operator>=(const BigInteger& other) const { return other <= *this; }
	bool operator>(const BigInteger& other) const { return other < *this; }
	BigInteger operator+() const { return BigInteger(*this); }
	BigInteger operator-() const {
		BigInteger ret = *this;
		if (ret.is_inf() || !ret.digits.empty()) { ret.is_negative = !ret.is_negative; }
		return ret;
	}
	BigInteger& operator++() {
		if (is_inf() || is_nan()) { return *this; }
		if (is_negative) { decrement_digits(&digits); }
		else { increment_digits(&digits); }
		shrink_digits();
		return *this;
	}
	BigInteger operator++(int) {
		BigInteger ret = *this;
		++*this;
		return ret;
	}
	BigInteger& operator--() {
		if (is_inf() || is_nan()) { return *this; }
		if (digits.size() == 0) { is_negative = true; }
		if (is_negative) { increment_digits(&digits); }
		else { decrement_digits(&digits); }
		shrink_digits();
		return *this;
	}
	BigInteger operator--(int) {
		BigInteger ret = *this;
		--*this;
		return ret;
	}
	BigInteger& operator+=(const BigInteger& other) {
		if (is_nan() || is_inf()) { return *this; }
		if (other.is_nan()) {
			set_invalid();
			return *this;
		}
		if (other.is_inf()) {
			state = state_infinity;
			is_negative = other.is_negative;
			digits.clear();
			return *this;
		}
		if (!other) { return *this; }
		if (!*this) { return *this = other; }
		if (is_negative == other.is_negative) { add_digits(&digits, other.digits, 0); }
		else {
			bool ch = subtract_digits(&digits, other.digits, 0);
			if (ch) { is_negative = !is_negative; }
		}
		shrink_digits();
		return *this;
	}
	BigInteger& operator-=(const BigInteger& other) {
		if (is_nan() || is_inf()) { return *this; }
		if (other.is_nan()) {
			set_invalid();
			return *this;
		}
		if (other.is_inf()) {
			state = state_infinity;
			is_negative = other.is_negative;
			digits.clear();
			return *this;
		}
		if (!other) { return *this; }
		if (!*this) { return *this = -other; }
		if (is_negative != other.is_negative) { add_digits(&digits, other.digits, 0); }
		else {
			bool ch = subtract_digits(&digits, other.digits, 0);
			if (ch) { is_negative = !is_negative; }
		}
		shrink_digits();
		return *this;
	}
	BigInteger& operator*=(const BigInteger& other) {
		if (is_nan()) { return *this; }
		if (other.is_nan()) {
			set_invalid();
			return *this;
		}
		is_negative ^= other.is_negative;
		if (is_inf() || other.is_inf()) {
			state = state_infinity;
			return *this;
		}
		carry_convolve(&digits, other.digits);
		shrink_digits();
		return *this;
	}
	BigInteger& operator/=(const BigInteger& other) {
		if (is_nan()) { return *this; }
		if (other.is_nan()) {
			set_invalid();
			return *this;
		}
		if (is_inf() || other.is_inf()) {
			state = state_infinity;
			is_negative ^= other.is_negative;
			return *this;
		}
		pair<BigInteger, BigInteger> bp = div_mod(*this, other);
		return *this = bp.first;
	}
	BigInteger& operator%=(const BigInteger& other) {
		BigInteger q = *this;
		q /= other;
		return *this -= q * other;
	}
	BigInteger operator+(const BigInteger& other) const { return BigInteger(*this) += other; }
	BigInteger operator-(const BigInteger& other) const { return BigInteger(*this) -= other; }
	BigInteger operator*(const BigInteger& other) const { return BigInteger(*this) *= other; }
	BigInteger operator/(const BigInteger& other) const { return BigInteger(*this) /= other; }
	BigInteger operator%(const BigInteger& other) const { return BigInteger(*this) %= other; }
	template <typename T>
	BigInteger operator^(const T& m) const { return pow(*this, m); }
	template <typename T>
	friend BigInteger operator+(const T& lhs, const BigInteger& rhs) {
		return BigInteger(lhs) + rhs;
	}
	template <typename T>
	friend BigInteger operator-(const T& lhs, const BigInteger& rhs) {
		return BigInteger(lhs) - rhs;
	}
	template <typename T>
	friend BigInteger operator*(const T& lhs, const BigInteger& rhs) {
		return BigInteger(lhs) * rhs;
	}
	template <typename T>
	friend BigInteger operator/(const T& lhs, const BigInteger& rhs) {
		return BigInteger(lhs) / rhs;
	}
	template <typename T>
	friend BigInteger operator%(const T& lhs, const BigInteger& rhs) {
		return BigInteger(lhs) % rhs;
	}
	BigInteger& operator>>=(int a) { return *this /= pow(BigInteger(2), a); }
	BigInteger operator>>(int a) const { return BigInteger(*this) >>= a; }
	BigInteger& operator<<=(int a) { return *this *= pow(BigInteger(2), a); }
	BigInteger operator<<(int a) const { return BigInteger(*this) <<= a; }
	friend ostream& operator<<(ostream& os, const BigInteger& other) {
		os << other.to_string();
		return os;
	}
	friend istream& operator>>(istream& ist, BigInteger& other) {
		string s;
		ist >> s;
		other = s;
		return ist;
	}
	static void read(BigInteger* bint) {
		string s;
		char c;
		c = getchar();
		while (c <= ' ') { c = getchar(); }
		while (c > ' ') {
			s += c;
			c = getchar();
		}
		*bint = s;
	}
	static void write(const BigInteger& bint) { printf("%s", bint.to_string().c_str()); }
	static void writeln(const BigInteger& bint) { printf("%s\n", bint.to_string().c_str()); }
	bool is_nan() const { return state != state_valid && state != state_infinity; }
	bool is_inf() const { return state == state_infinity; }
	bool is_positive_inf() const { return state == state_infinity && !is_negative; }
	bool is_negative_inf() const { return state == state_infinity && is_negative; }
	string to_string() const {
		if (is_nan()) { return "nan"; }
		if (is_inf()) { return (is_negative) ? "-inf" : "inf"; }
		if (digits.size() == 0) { return "0"; }
		string s;
		if (is_negative) { s += '-'; }
		size_t sz = digits.size();
		uint32_t v = DIGITS_BASE;
		while (v > 1) {
			v /= BASE;
			if (digits[sz - 1] / v == 0) { continue; }
			s += '0' ^ digits[sz - 1] / v % BASE;
		}
		for (size_t i = sz - 1; i > 0; i--) {
			v = DIGITS_BASE;
			while (v > 1) {
				v /= BASE;
				s += '0' ^ digits[i - 1] / v % BASE;
			}
		}
		return s;
	}
	uint64_t digit_sum() const {
		if (is_nan()) { return numeric_limits<uint64_t>::quiet_NaN(); }
		if (is_inf()) { return numeric_limits<uint64_t>::infinity(); }
		uint64_t ret = 0;
		for (uint32_t i : digits) {
			while (i) {
				ret += i % BASE;
				i /= BASE;
			}
		}
		return ret;
	}
	size_t digit_count() const {
		if (is_nan()) { return numeric_limits<uint64_t>::quiet_NaN(); }
		if (is_inf()) { return numeric_limits<uint64_t>::infinity(); }
		size_t n = digits.size();
		if (n == 0) { return 0; }
		size_t ret = DIGITS_NUM * (n - 1);
		uint32_t v = 1;
		while (v <= digits[n - 1]) {
			ret++;
			v *= BASE;
		}
		return ret;
	}
	template <typename T>
	static BigInteger pow(const BigInteger& bint, T m) {
		assert(m >= 0);
		BigInteger ret = 1, b = bint;
		while (m) {
			if (m % 2 == 1) { ret *= b; }
			b *= BigInteger(b);
			m >>= 1;
		}
		return ret;
	}
	template <typename T, typename S>
	static BigInteger pow(const BigInteger& bint, const T& m, const S& mod) {
		assert(m >= 0);
		BigInteger ret = 1, b = bint;
		while (m) {
			if (m % 2 == 1) {
				ret *= b;
				ret %= mod;
			}
			b *= BigInteger(b);
			b %= mod;
			m >>= 1;
		}
		return ret;
	}
	static BigInteger fact(uint32_t n) {
		BigInteger ret = 1;
		for (uint32_t i = 0; i < n; i++) { ret *= i + 1; }
		return ret;
	}
	static BigInteger gcd(const BigInteger& x, const BigInteger& y) {
		BigInteger k = x, l = y;
		while (l) {
			BigInteger tmp = k;
			k = l;
			l = tmp % l;
		}
		return k;
	}
	static BigInteger lcm(const BigInteger& x, const BigInteger& y) {
		return x * y / gcd(x, y);
	}
	static BigInteger abs(const BigInteger& x) {
		return (x.is_negative) ? -x : x;
	}
	static BigInteger min(const BigInteger& x, const BigInteger& y) {
		return (x <= y) ? x : y;
	}
	static BigInteger max(const BigInteger& x, const BigInteger& y) {
		return (x >= y) ? x : y;
	}
	static pair<BigInteger, BigInteger> div_mod(const BigInteger& x, const BigInteger& y) {
		BigInteger q, m;
		if (y == 0) {
			q.state = state_infinity;
			m.state = state_infinity;
			return make_pair(q, m);
		}
		if (x == 0) {
			return make_pair(0, 0);
		}
		q = x;
		divide_digits(&q.digits, y.digits);
		q.is_negative ^= y.is_negative;
		m = x - q * y;
		if (m < 0) {
			if (x < 0) {
				--q;
				m += y;
			} else {
				++q;
				m -= y;
			}
		}
		return make_pair(q, m);
	}
};

#endif /* BigInteger_hpp */
