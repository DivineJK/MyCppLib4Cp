//
//  FormalPowerSeries.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/09.
//

#ifndef FormalPowerSeries_hpp
#define FormalPowerSeries_hpp

namespace fps_internal {
constexpr int INTERNAL_MOD = 998244353; // 998244353 = 119 * 2 ^ 23 + 1
constexpr int get_base_size(int m) {
	int v = m - 1;
	int r = 0;
	while (!(v & 1)) {
		r++;
		v >>= 1;
	}
	return r;
}
static constexpr int modinv(uint32_t n) {
	uint32_t p = n, t = 1;
	while (p > 1) {
		t = (uint64_t)t * (INTERNAL_MOD - (INTERNAL_MOD / p)) % INTERNAL_MOD;
		p = INTERNAL_MOD % p;
	}
	return t;
}
constexpr int INTERNAL_BASE_SIZE = get_base_size(INTERNAL_MOD);
static_assert(INTERNAL_BASE_SIZE >= 15, "INTERNAL_MOD is not suitable for NTT");
static constexpr bool is_primitive(int i) {
	int v = i;
	for (int i = 0; i < INTERNAL_BASE_SIZE; i++) {
		if (v == 1) { return false; }
		v = (int64_t)v * v % INTERNAL_MOD;
	}
	return v == 1;
}
static constexpr int get_primitive_base() {
	for (int i = 1; i < INTERNAL_MOD; i++) {
		if (is_primitive(i)) { return i; }
	}
	return 0;
}
constexpr int INTERNAL_PRI_ROOT = get_primitive_base();
constexpr int INTERNAL_INV_ROOT = modinv(INTERNAL_PRI_ROOT);
static bool is_calc_bases = false;
static int internal_pri_bases[INTERNAL_BASE_SIZE + 1];
static int internal_inv_bases[INTERNAL_BASE_SIZE + 1];
static int internal_cml_bases[INTERNAL_BASE_SIZE + 1];
class internal_modint {
	using modint = internal_modint;
private:
	uint32_t x = 0;
private:
	static void initialize() {
		if (!is_calc_bases) {
			make_base();
			is_calc_bases = true;
		}
	}
	static void make_base() {
		internal_pri_bases[INTERNAL_BASE_SIZE] = INTERNAL_PRI_ROOT;
		internal_inv_bases[INTERNAL_BASE_SIZE] = INTERNAL_INV_ROOT;
		for (int i = INTERNAL_BASE_SIZE; i > 0; i--) {
			int pv = internal_pri_bases[i];
			int iv = internal_inv_bases[i];
			internal_pri_bases[i - 1] = (int64_t)pv * pv % INTERNAL_MOD;
			internal_inv_bases[i - 1] = (int64_t)iv * iv % INTERNAL_MOD;
		}
		int ie = 1;
		for (int i = 0; i < INTERNAL_BASE_SIZE - 1; i++) {
			internal_cml_bases[i] = (int64_t)ie * internal_pri_bases[i + 2] % INTERNAL_MOD;
			ie = (int64_t)ie * internal_inv_bases[i + 2] % INTERNAL_MOD;
		}
	}
	static void ntt(vector<modint>* arr) {
		size_t n = arr->size();
		if (n <= 1) { return; }
		assert(!(n & (n - 1)));
		int m = 0;
		while ((1 << m) < n) { m++; }
		int p = 1, t = 1 << (m - 1);
		for (int i = 0; i < m; i++) {
			modint g = 1;
			for (int j = 0; j < p; j++) {
				modint u, v;
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
				g *= internal_cml_bases[z];
			}
			p <<= 1;
			t >>= 1;
		}
	}
	static void intt(vector<modint>* arr) {
		size_t n = arr->size();
		if (n <= 1) { return; }
		modint invn = inv((int)n);
		assert(!(n & (n - 1)));
		int m = 0;
		while ((1 << m) < n) { m++; }
		for (int i = 0; i < m; i++) {
			modint g = 1, s = internal_inv_bases[i + 1];
			int offset = 1 << i;
			for (int k = 0; k < offset; k++) {
				for (int j = k; j < n; j += (1 << (i + 1))) {
					modint x = (*arr)[j], y = (*arr)[j + offset] * g;
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
public:
	internal_modint() : x(0) { initialize(); }
	internal_modint(bool b) {
		x = (b) ? 1 : 0;
		initialize();
	}
	internal_modint(int a) {
		int v = a % INTERNAL_MOD;
		x = (v < 0) ? INTERNAL_MOD + v : v;
		initialize();
	}
	internal_modint(int64_t a) {
		int v = a % INTERNAL_MOD;
		x = (v < 0) ? INTERNAL_MOD + v : v;
		initialize();
	}
	internal_modint(uint32_t a) {
		x = a % INTERNAL_MOD;
		initialize();
	}
	internal_modint(uint64_t a) {
		x = a % INTERNAL_MOD;
		initialize();
	}
	internal_modint(const modint& mint) { operator=(mint); }
	uint32_t val() const { return x; }
	modint& operator=(const modint& mint) {
		x = mint.x;
		initialize();
		return *this;
	}
	explicit operator bool() const { return x != 0; }
	bool operator!() const { return x == 0; }
	friend ostream& operator<<(ostream& os, const modint& mint) {
		os << mint.x;
		return os;
	}
	friend istream& operator>>(istream& ist, modint& mint) {
		int a;
		ist >> a;
		mint = a;
		return ist;
	}
	bool operator==(const modint& mint) const { return x == mint.x; }
	bool operator!=(const modint& mint) const { return x != mint.x; }
	modint inv() const {
		if (x == 0) { return modint(); }
		uint32_t a = x, t = 1;
		while (a > 1) {
			t = (uint64_t)t * (INTERNAL_MOD - INTERNAL_MOD / a) % INTERNAL_MOD;
			a = INTERNAL_MOD % a;
		}
		return modint(t);
	}
	static modint inv(const modint& mint) {
		if (mint.x == 0) { return mint; }
		uint32_t a = mint.x, t = 1;
		while (a > 1) {
			t = (uint64_t)t * (INTERNAL_MOD - INTERNAL_MOD / a) % INTERNAL_MOD;
			a = INTERNAL_MOD % a;
		}
		return modint(t);
	}
	modint operator+() const { return modint(*this); }
	modint operator-() const {
		modint ret = modint(*this);
		ret.x = (ret.x == 0) ? ret.x : INTERNAL_MOD - ret.x;
		return ret;
	}
	modint& operator++() {
		x++;
		if (x == INTERNAL_MOD) { x = 0; }
		return *this;
	}
	modint operator++(int) {
		modint ret = modint(*this);
		++*this;
		return ret;
	}
	modint& operator--() {
		x = (x == 0) ? INTERNAL_MOD - 1 : x - 1;
		return *this;
	}
	modint operator--(int) {
		modint ret = modint(*this);
		--*this;
		return ret;
	}
	modint& operator+=(const modint& mint) {
		x = (x + mint.x >= INTERNAL_MOD) ? x + mint.x - INTERNAL_MOD : x + mint.x;
		return *this;
	}
	modint& operator-=(const modint& mint) {
		x = (x < mint.x) ? INTERNAL_MOD + x - mint.x : x - mint.x;
		return *this;
	}
	modint& operator*=(const modint& mint) {
		x = (uint64_t)x * mint.x % INTERNAL_MOD;
		return *this;
	}
	modint& operator/=(const modint& mint) {
		*this *= inv(mint);
		return *this;
	}
	modint operator+(const modint& mint) const { return modint(*this) += mint; }
	modint operator-(const modint& mint) const { return modint(*this) -= mint; }
	modint operator*(const modint& mint) const { return modint(*this) *= mint; }
	modint operator/(const modint& mint) const { return modint(*this) /= mint; }
	template <typename T>
	friend modint operator+(const T& lhs, const modint& mint) { return modint(lhs) += mint; }
	template <typename T>
	friend modint operator-(const T& lhs, const modint& mint) { return modint(lhs) -= mint; }
	template <typename T>
	friend modint operator*(const T& lhs, const modint& mint) { return modint(lhs) *= mint; }
	template <typename T>
	friend modint operator/(const T& lhs, const modint& mint) { return modint(lhs) /= mint; }
	template <typename T>
	static modint pow(const modint& a, T m) {
		modint ret = 1, v;
		if (m < 0) {
			m = -m;
			v = inv(a);
		} else { v = a; }
		while (m) {
			if (m & 1) { ret *= v; }
			v *= v;
			m >>= 1;
		}
		return ret;
	}
	template <typename T>
	modint operator^(T m) const { return pow(*this, m); }
	static int sqrt(const modint& n) {
		if (n.val() == 0 || n.val() == 1) { return n.val(); }
		if (pow(n, INTERNAL_MOD >> 1) + 1 == 0) { return -1; }
		int m = INTERNAL_MOD - 1, q = 0;
		while (!(m & 1)) {
			q++;
			m >>= 1;
		}
		modint z = 1, r = pow(n, (m + 1) >> 1), t = pow(n, m);
		while (pow(z, INTERNAL_MOD >> 1) == 1) { z++; }
		z = pow(z, m);
		while (t != 1) {
			modint c = t;
			int i = 0;
			while (c != 1) {
				i++;
				c *= c;
			}
			modint b = pow(z, 1 << (q - i - 1));
			q = i;
			z = b * b;
			t *= z;
			r *= b;
		}
		return (r.x > INTERNAL_MOD - r.x) ? INTERNAL_MOD - r.x : r.x;
	}
	static void convolve_one(vector<modint>* a, const vector<modint>& b) {
		if (a->size() == 0 || b.size() == 0) {
			*a = {};
			return;
		}
		size_t n = a->size(), m = b.size(), bas = 1;
		while (bas < n + m) { bas <<= 1; }
		a->resize(bas);
		vector<modint> v(bas);
		for (size_t i = 0; i < m; i++) { v[i] = b[i]; }
		ntt(a);
		ntt(&v);
		for (size_t i = 0; i < bas; i++) { (*a)[i] *= v[i]; }
		intt(a);
	}
	static void convolve(vector<modint>* a, const vector<modint>& b) {
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
		if (bas <= (1 << bs)) {
			convolve_one(a, b);
			return;
		}
		int lim = 1 << 18;
		int s = (n + lim - 1) / lim, t = (m + lim - 1) / lim;
		a->resize(s * lim);
		vector<modint> y(t * lim);
		for (int i = 0; i < m; i++) { y[i] = b[i]; }
		vector<modint> ret((s + t) * lim);
		for (int i = 0; i < s; i++) {
			vector<modint> u(lim);
			for (int j = 0; j < lim; j++) { u[j] = (*a)[lim * i + j]; }
			for (int j = 0; j < t; j++) {
				vector<modint> res(lim);
				for (int k = 0; k < lim; k++) { res[k] = y[lim * j + k]; }
				convolve_one(&u, res);
				for (int k = 0; k < u.size(); k++) { ret[(i + j) * lim + k] += u[k]; }
			}
		}
		*a = ret;
	}
	static void get_fact_table(int n, vector<modint>* fact,
							   vector<modint>* invf) {
		assert(n >= 0);
		modint t = 1;
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
	static void get_inv_table(int n, vector<modint>* invn) {
		assert(invn && n >= 0);
		invn->resize(n + 1);
		for (int i = 0; i < 2; i++) {
			(*invn)[i] = 1;
			if (n == i) { return; }
		}
		for (int i = 2; i <= n; i++) {
			int t = INTERNAL_MOD % i;
			(*invn)[i] = (INTERNAL_MOD - INTERNAL_MOD / i) * (*invn)[t];
		}
	}
};

}	// namespace fps_internal

static constexpr int MOD_FOR_FPS = fps_internal::INTERNAL_MOD;
using modint_for_fps = fps_internal::internal_modint;

class FormalPowerSeries {
private:
	// state ==  0: valid
	// state == -1: infinity by divide by zero
	// state == -2: invalid
	int state = 0;
	vector<modint_for_fps> poly;
private:
	static constexpr uint32_t getMinBinary(uint32_t v) {
		uint32_t ret = 1;
		while (ret < v) { ret <<= 1; }
		return ret;
	}
	static constexpr int getIntLog2(int n) {
		if (n <= 0) { return -1; }
		for (int i = 0; i < 32; i++) {
			if ((n >> i) == 1) { return i; }
		}
		return 32;
	}
	static FormalPowerSeries babystep_compose(const FormalPowerSeries& f,
											  const FormalPowerSeries& g) {
		int n = f.getDegree();
		uint32_t lg = getMinBinary(n + 1);
		uint32_t b = lg >> 1;
		vector<FormalPowerSeries> fpss(lg);
		for (int i = 0; i <= n; i++) { fpss[i] = f.getCoeff(i); }
		FormalPowerSeries h = g;
		while (b) {
			for (int i = 0; i < b; i++) {
				fpss[i] = fpss[i << 1] + h * fpss[(i << 1) | 1];
				fpss[i].crop(n + 1);
			}
			h *= h;
			h.crop(n + 1);
			b >>= 1;
		}
		return fpss[0];
	}
	static FormalPowerSeries compose_internal_lower_zero(const FormalPowerSeries& f,
											  const FormalPowerSeries& q, int k) {
		int n = f.getDegree();
		FormalPowerSeries ret = 0, rq = 1;
		for (int i = 0; k * i <= n; i++) {
			ret += f.getCoeff(i) * (rq << (i * k));
			ret.crop(n + 1);
			rq *= q;
			rq.setDegree(n - k * i);
		}
		return ret;
	}
	static FormalPowerSeries invert_internal(const FormalPowerSeries& f) {
		FormalPowerSeries ret = f.getCoeff(0).inv(), g = f.getCoeff(0);
		uint32_t correct = 1;
		int n = f.getDegree();
		uint32_t lg = getMinBinary(n + 1);
		while (correct < lg) {
			g.setDegree((correct << 1) - 1);
			for (int i = correct; i < (correct << 1); i++) { g[i] = f.getCoeff(i); }
			correct <<= 1;
			FormalPowerSeries h = 2 - ret * g;
			h.crop(correct);
			ret *= h;
			ret.crop(correct);
		}
		ret.setDegree(n);
		return ret;
	}
	static FormalPowerSeries pow_internal(const FormalPowerSeries& f, uint64_t m,
										  int modDeg = INT32_MAX) {
		FormalPowerSeries lps = f;
		lps.setDegree((int)modDeg - 1);
		lps = m * log(lps);
		lps = exp(lps);
		lps.crop(modDeg);
		return lps;
	}
	static void getDivMod_internal(const FormalPowerSeries& f, const FormalPowerSeries& g,
						  FormalPowerSeries* outDiv, FormalPowerSeries* outRem) {
		int n = f.getDegree(), m = g.getDegree();
		if (n < m) {
			*outDiv = 0;
			*outRem = f;
			return;
		}
		FormalPowerSeries rf = f.getReversed(), rg = g.getReversed();
		rg.setDegree(n - m);
		rg.invert();
		rf *= rg;
		rf.setDegree(n - m);
		*outDiv = rf.getReversed();
		if (outRem) { *outRem = f - g * (*outDiv); }
	}
	static FormalPowerSeries sqrt_internal(const FormalPowerSeries& f) {
		int z = modint_for_fps::sqrt(f.getCoeff(0));
		if (z == -1) {
			FormalPowerSeries ret;
			ret.state = -2;
			return ret;
		}
		modint_for_fps mz = z;
		int n = f.getDegree(), b = 1;
		uint32_t lg = getMinBinary(n + 1);
		FormalPowerSeries ret = mz, h = f.getCoeff(0);
		while (b <= lg) {
			ret += h * ret.getInverted();
			ret /= 2;
			h.setDegree((b << 1) - 1);
			for (int i = b; i < (b << 1); i++) { h[i] = f.getCoeff(i); }
			b <<= 1;
			ret.setDegree(b - 1);
		}
		ret.setDegree(n);
		return ret;
	}
	static FormalPowerSeries getPositiveBinomialExpansion(uint32_t n, int lim) {
		if (lim == 0) { return FormalPowerSeries(); }
		uint32_t deg = (lim < 0) ? n : lim - 1;
		deg = min(deg, n);
		FormalPowerSeries ret;
		ret.setDegree(deg);
		ret[0] = 1;
		if (n == 0) { return ret; }
		vector<modint_for_fps> invn;
		modint_for_fps::get_inv_table(deg, &invn);
		for (int i = 1; i <= deg; i++) {
			ret[i] = ret[i - 1] * (n - i + 1) * invn[i];
		}
		return ret;
	}
	static FormalPowerSeries getNegativeBinomialExpansion(uint32_t n, int lim) {
		assert(lim >= 0);
		if (lim == 0) { return FormalPowerSeries(); }
		FormalPowerSeries ret;
		ret.setDegree(lim - 1);
		ret[0] = 1;
		if (n == 0) { return ret; }
		vector<modint_for_fps> invn;
		modint_for_fps::get_inv_table(lim, &invn);
		for (int i = 1; i < lim; i++) { ret[i] = -ret[i - 1] * (i + n - 1) * invn[i]; }
		return ret;
	}
	static void hgcd(const FormalPowerSeries& a, const FormalPowerSeries& b,
					 FormalPowerSeries* p, FormalPowerSeries* q,
					 FormalPowerSeries* r, FormalPowerSeries* s) {
		int n = a.getDegree();
		int m = (n + 1) / 2;
		if (b.getDegree() < m) {
			*p = 1;
			*q = 0;
			*r = 0;
			*s = 1;
			return;
		}
		FormalPowerSeries a0 = a >> m, b0 = b >> m;
		FormalPowerSeries r0, r1, r2, r3;
		hgcd(a0, b0, &r0, &r1, &r2, &r3);
		a0 = r3 * a - r1 * b;
		b0 = r0 * b - r2 * a;
		a0.regularize();
		b0.regularize();
		if (b0.getDegree() < m) {
			*p = r0;
			*q = r1;
			*r = r2;
			*s = r3;
			return;
		}
		int l = b0.getDegree();
		FormalPowerSeries ad, d;
		getDivMod(a0, b0, &ad, &d);
		int k = 2 * m - l;
		FormalPowerSeries c0 = b0 >> k, d0 = d >> k;
		FormalPowerSeries s0, s1, s2, s3;
		hgcd(c0, d0, &s0, &s1, &s2, &s3);
		r1 += r0 * ad;
		r3 += r2 * ad;
		*p = r1 * s0 + r0 * s2;
		*q = r1 * s1 + r0 * s3;
		*r = r3 * s0 + r2 * s2;
		*s = r3 * s1 + r2 * s3;
		p->regularize();
		q->regularize();
		r->regularize();
		s->regularize();
	}
public:
	FormalPowerSeries() {}
	template <typename T>
	FormalPowerSeries(const T& fps) { operator=(fps); }
	FormalPowerSeries& operator=(int val) { return operator=((modint_for_fps)val); }
	FormalPowerSeries& operator=(uint32_t val) { return operator=((modint_for_fps)val); }
	FormalPowerSeries& operator=(int64_t val) { return operator=((modint_for_fps)val); }
	FormalPowerSeries& operator=(uint64_t val) { return operator=((modint_for_fps)val); }
	FormalPowerSeries& operator=(const modint_for_fps& mint) {
		state = 0;
		poly = {mint};
		regularize();
		return *this;
	}
	template <typename S>
	FormalPowerSeries& operator=(const vector<S>& arr) {
		state = 0;
		poly.reserve(arr.size());
		for (int i = 0; i < arr.size(); i++) { poly.push_back((modint_for_fps)arr[i]); }
		regularize();
		return *this;
	}
	FormalPowerSeries& operator=(const FormalPowerSeries& fps) {
		state = fps.state;
		poly = fps.poly;
		regularize();
		return *this;
	}
	friend ostream& operator<<(ostream& os, const FormalPowerSeries& rhs) {
		if (rhs.state < 0) {
			os << -1;
			return os;
		}
		if (rhs.poly.size() == 0) {
			os << 0;
			return os;
		}
		for (int i = 0; i < rhs.poly.size(); i++) {
			if (i > 0) { os << " "; }
			os << rhs.poly[i];
		}
		return os;
	}
	friend istream& operator>>(istream& ist, FormalPowerSeries& rhs) {
		int n = rhs.getDegree();
		for (int i = 0; i <= n; i++) {
			modint_for_fps v;
			ist >> rhs[i];
		}
		return ist;
	}
	explicit operator bool() const {
		if (state < 0) { return state == -1; }
		for (const modint_for_fps& v : poly) {
			if (v) { return true; }
		}
		return false;
	}
	bool operator!() const {
		if (state < 0) { return state == -2; }
		for (const modint_for_fps& v : poly) {
			if (v) { return false; }
		}
		return true;
	}
	bool operator==(const FormalPowerSeries& other) const {
		if (state < 0) {
			if (state == -1) { return other.state == -1; }
			else if (state == -2) { return false; }
		}
		int n = max(getDegree(), other.getDegree());
		for (int i = 0; i <= n; i++) {
			if (getCoeff(i) != other.getCoeff(i)) { return false; }
		}
		return true;
	}
	bool operator!=(const FormalPowerSeries& other) const {
		if (state < 0) {
			if (state == -1) { return other.state != -1; }
			else if (state == -2) { return false; }
		}
		int n = max(getDegree(), other.getDegree());
		for (int i = 0; i <= n; i++) {
			if (getCoeff(i) != other.getCoeff(i)) { return true; }
		}
		return false;
	}
	modint_for_fps& operator[](size_t i) {
		if (i >= poly.size()) { poly.resize(i + 1); }
		return poly[i];
	}
	const modint_for_fps& operator[](size_t i) const { return poly[i]; }
	FormalPowerSeries operator+() const { return FormalPowerSeries(*this); }
	FormalPowerSeries operator-() const {
		FormalPowerSeries ret = FormalPowerSeries(*this);
		for (int i = 0; i < poly.size(); i++) {
			ret.poly[i] = -poly[i];
		}
		return ret;
	}
	FormalPowerSeries& operator++() {
		if (poly.size() == 0) { poly = {0}; }
		++poly[0];
		regularize();
		return *this;
	}
	FormalPowerSeries operator++(int) {
		FormalPowerSeries ret = FormalPowerSeries(*this);
		++*this;
		return ret;
	}
	FormalPowerSeries& operator--() {
		if (poly.size() == 0) { poly = {0}; }
		--poly[0];
		regularize();
		return *this;
	}
	FormalPowerSeries operator--(int) {
		FormalPowerSeries ret = FormalPowerSeries(*this);
		--*this;
		return ret;
	}
	FormalPowerSeries& operator+=(const FormalPowerSeries& other) {
		int n = other.getDegree();
		for (int i = n; i >= 0; i--) { (*this)[i] += other.getCoeff(i); }
		regularize();
		return *this;
	}
	FormalPowerSeries& operator-=(const FormalPowerSeries& other) {
		int n = other.getDegree();
		for (int i = n; i >= 0; i--) { (*this)[i] -= other.getCoeff(i); }
		regularize();
		return *this;
	}
	FormalPowerSeries& operator*=(const FormalPowerSeries& other) {
		int n = getDegree(), m = other.getDegree();
		if (n == -1 || m == -1) {
			poly = {};
			return *this;
		}
		modint_for_fps::convolve(&poly, other.poly);
		setDegree(n + m);
		regularize();
		return *this;
	}
	FormalPowerSeries& operator*=(int v) { return operator*=((modint_for_fps)v); }
	FormalPowerSeries& operator*=(uint32_t v) { return operator*=((modint_for_fps)v); }
	FormalPowerSeries& operator*=(int64_t v) { return operator*=((modint_for_fps)v); }
	FormalPowerSeries& operator*=(uint64_t v) { return operator*=((modint_for_fps)v); }
	FormalPowerSeries& operator*=(const modint_for_fps& v) {
		for (int i = 0; i <= getDegree(); i++) { (*this)[i] *= v; }
		regularize();
		return *this;
	}
	FormalPowerSeries& operator/=(const FormalPowerSeries& other) {
		FormalPowerSeries div;
		getDivMod(*this, other, &div, nullptr);
		*this = div;
		return *this;
	}
	FormalPowerSeries& operator/=(int v) { return operator/=((modint_for_fps)v); }
	FormalPowerSeries& operator/=(uint32_t v) { return operator/=((modint_for_fps)v); }
	FormalPowerSeries& operator/=(int64_t v) { return operator/=((modint_for_fps)v); }
	FormalPowerSeries& operator/=(uint64_t v) { return operator/=((modint_for_fps)v); }
	FormalPowerSeries& operator/=(const modint_for_fps& v) {
		if (v == 0) {
			state = -1;
			return *this;
		}
		modint_for_fps iv = v.inv();
		for (int i = 0; i <= getDegree(); i++) { (*this)[i] *= iv; }
		regularize();
		return *this;
	}
	FormalPowerSeries& operator%=(const FormalPowerSeries& other) {
		FormalPowerSeries div, mod;
		getDivMod(*this, other, &div, &mod);
		*this = mod;
		return *this;
	}
	FormalPowerSeries& operator<<=(int n) {
		int m = getDegree();
		extend(n);
		for (int i = m; i >= 0; i--) { (*this)[i + n] = (*this)[i]; }
		for (int i = 0; i < n; i++) { (*this)[i] = 0; }
		return *this;
	}
	FormalPowerSeries& operator>>=(int n) {
		int m = getDegree();
		if (m < n) {
			poly = {};
			return *this;
		}
		for (int i = 0; i <= m - n; i++) { (*this)[i] = (*this)[i + n]; }
		setDegree(m - n);
		return *this;
	}
	FormalPowerSeries operator+(const FormalPowerSeries& other) const {
		return FormalPowerSeries(*this) += other;
	}
	FormalPowerSeries operator-(const FormalPowerSeries& other) const {
		return FormalPowerSeries(*this) -= other;
	}
	FormalPowerSeries operator*(const FormalPowerSeries& other) const {
		return FormalPowerSeries(*this) *= other;
	}
	FormalPowerSeries operator*(int v) {
		return FormalPowerSeries(*this) *= v;
	}
	FormalPowerSeries operator*(uint32_t v) {
		return FormalPowerSeries(*this) *= v;
	}
	FormalPowerSeries operator*(int64_t v) {
		return FormalPowerSeries(*this) *= v;
	}
	FormalPowerSeries operator*(uint64_t v) {
		return FormalPowerSeries(*this) *= v;
	}
	FormalPowerSeries operator*(const modint_for_fps& other) const {
		return FormalPowerSeries(*this) *= other;
	}
	FormalPowerSeries operator/(const FormalPowerSeries& other) const {
		FormalPowerSeries div;
		getDivMod(*this, other, &div, nullptr);
		div.state = state;
		return div;
	}
	FormalPowerSeries operator/(int v) const {
		return FormalPowerSeries(*this) /= v;
	}
	FormalPowerSeries operator/(uint32_t v) const {
		return FormalPowerSeries(*this) /= v;
	}
	FormalPowerSeries operator/(int64_t v) const {
		return FormalPowerSeries(*this) /= v;
	}
	FormalPowerSeries operator/(uint64_t v) const {
		return FormalPowerSeries(*this) /= v;
	}
	FormalPowerSeries operator/(const modint_for_fps& v) const {
		return FormalPowerSeries(*this) /= v;
	}
	FormalPowerSeries operator%(const FormalPowerSeries& other) const {
		FormalPowerSeries div, mod;
		getDivMod(*this, other, &div, &mod);
		mod.state = state;
		return mod;
	}
	template <typename T>
	friend FormalPowerSeries operator+(const T& lhs, const FormalPowerSeries& rhs) {
		return FormalPowerSeries(lhs) += rhs;
	}
	template <typename T>
	friend FormalPowerSeries operator-(const T& lhs, const FormalPowerSeries& rhs) {
		return FormalPowerSeries(lhs) -= rhs;
	}
	template <typename T>
	friend FormalPowerSeries operator*(const T& lhs, const FormalPowerSeries& rhs) {
		FormalPowerSeries ret = FormalPowerSeries(rhs);
		int n = ret.getDegree();
		for (int i = 0; i <= n; i++) { ret[i] *= (modint_for_fps)lhs; }
		return ret;
	}
	FormalPowerSeries operator<<(int n) const { return FormalPowerSeries(*this) <<= n; }
	FormalPowerSeries operator>>(int n) const { return FormalPowerSeries(*this) >>= n; }
	template <typename T>
	FormalPowerSeries operator^(const T& aExp) const { return pow(*this, aExp); }
	modint_for_fps getCoeff(size_t i) const {
		if (i >= poly.size()) { return 0; }
		return poly[i];
	}
	void setCoeff(size_t i, const modint_for_fps& val) { (*this)[i] = val; }
	int getDegree() const { return (int)poly.size() - 1; }
	void setDegree(int aDeg) { poly.resize(aDeg + 1); }
	void extend(int aDeg) {
		if (aDeg > 0) { setDegree(getDegree() + aDeg); }
	}
	void extendTo(int aDeg) {
		if (aDeg > getDegree()) { setDegree(aDeg); }
	}
	void regularize() {
		int s = (int)poly.size();
		for (int i = (int)poly.size(); i > 0; i--) {
			if (poly[i - 1]) { break; }
			s--;
		}
		poly.resize(s);
	}
	int getState() const { return state; }
	static FormalPowerSeries getInfPolynomial() {
		FormalPowerSeries ret;
		ret.state = -1;
		return ret;
	}
	static FormalPowerSeries getInvalidPolynomial() {
		FormalPowerSeries ret;
		ret.state = -2;
		return ret;
	}
	modint_for_fps evaluate(const modint_for_fps& v) const {
		modint_for_fps ret = 0;
		int n = getDegree();
		for (int i = n; i >= 0; i--) {
			ret *= v;
			ret += getCoeff(i);
		}
		return ret;
	}
	FormalPowerSeries& reverse() {
		int n = getDegree();
		for (int i = 0; 2 * i < n; i++) {
			modint_for_fps tmp = (*this)[i];
			(*this)[i] = (*this)[n - i];
			(*this)[n - i] = tmp;
		}
		return *this;
	}
	FormalPowerSeries& differentiate() {
		int n = getDegree();
		if (n == -1) { return *this; }
		for (int i = 1; i <= n; i++) { (*this)[i - 1] = i * (*this)[i]; }
		(*this)[n] = 0;
		regularize();
		return *this;
	}
	FormalPowerSeries& integrate() {
		regularize();
		int n = getDegree();
		if (n == -1) { return *this; }
		vector<modint_for_fps> invn;
		modint_for_fps::get_inv_table(n + 1, &invn);
		for (int i = n; i >= 0; i--) { (*this)[i + 1] = (*this)[i] * invn[i + 1]; }
		(*this)[0] = 0;
		return *this;
	}
	FormalPowerSeries& invert() {
		if (getCoeff(0) == 0) {
			state = -1;
			poly = {};
			return *this;
		}
		return *this = invert_internal(*this);
	}
	FormalPowerSeries& crop(int m) {
		if (getDegree() >= m) { setDegree(m - 1); }
		regularize();
		return *this;
	}
	FormalPowerSeries& compose(const FormalPowerSeries& fps) {
		int n = getDegree();
		if (n <= 0) { return *this; }
		if (fps.getDegree() < 0) { return *this = getCoeff(0); }
		int lg = getIntLog2(n + 1);
		int k = 0;
		while ((int64_t)lg * k * k <= n) { k++; }
		FormalPowerSeries p, q;
		p.setDegree(k - 1);
		q.setDegree(n - k);
		for (int i = 0; i < k; i++) { p[i] = fps.getCoeff(i); }
		for (int i = k; i <= n; i++) { q[i - k] = fps.getCoeff(i); }
		FormalPowerSeries bf = babystep_compose(*this, p);
		if (p.differentiate() == 0) { return *this = compose_internal_lower_zero(*this, q, k); }
		int b = 0;
		while (p.getCoeff(b) == 0) { b++; }
		p >>= b;
		vector<modint_for_fps> invf;
		modint_for_fps::get_fact_table(n, nullptr, &invf);
		FormalPowerSeries rq = 1;
		p.setDegree(n);
		p.invert();
		*this = 0;
		for (int i = 0; i * k <= n; i++) {
			*this += bf * rq * invf[i];
			crop(n + 1);
			bf.differentiate();
			bf *= p;
			bf >>= b;
			bf.crop(n + 1);
			rq *= q;
			rq <<= k;
			rq.crop(n + 1);
		}
		return *this;
	}
	FormalPowerSeries& shift(const modint_for_fps& v) {
		int n = getDegree();
		if (n < 0) { return *this; }
		vector<modint_for_fps> a(n + 1), c(n + 1), invn(n + 1);
		modint_for_fps iv = 1;
		for (int i = 0; i <= n; i++) {
			a[n - i] = getCoeff(i) * iv;
			iv *= i + 1;
		}
		invn[n] = iv.inv() * (n + 1);
		for (int i = n; i > 0; i--) { invn[i - 1] = invn[i] * i; }
		modint_for_fps cv = 1;
		for (int i = 0; i <= n; i++) {
			c[i] = cv * invn[i];
			cv *= v;
		}
		modint_for_fps::convolve(&a, c);
		extendTo(n);
		for (int i = 0; i <= n; i++) { (*this)[i] = a[n - i] * invn[i]; }
		return *this;
	}
	FormalPowerSeries getReversed() const { return FormalPowerSeries(*this).reverse(); }
	FormalPowerSeries getDifferentiated() const { return FormalPowerSeries(*this).differentiate(); }
	FormalPowerSeries getIntegrated() const { return FormalPowerSeries(*this).integrate(); }
	FormalPowerSeries getInverted() const {
		if (getCoeff(0) == 0) {
			FormalPowerSeries ret;
			ret.state = -1;
			ret.poly = {};
			return ret;
		}
		FormalPowerSeries ret = invert_internal(*this);
		ret.state = state;
		return ret;
	}
	FormalPowerSeries getPolynomialInverse(const FormalPowerSeries& g) const {
		return getPolynomialInverse(*this, g);
	}
	FormalPowerSeries getCropped(int m) const { return FormalPowerSeries(*this).crop(m); }
	FormalPowerSeries getComposition(const FormalPowerSeries& fps) const {
		return FormalPowerSeries(*this).compose(fps);
	}
	FormalPowerSeries getShifted(const modint_for_fps& v) const {
		return FormalPowerSeries(*this).shift(v);
	}
	vector<modint_for_fps> evaluateMultipoint(const vector<modint_for_fps>& vec) const {
		int n = (int)vec.size();
		uint32_t lg = getMinBinary(n);
		vector<FormalPowerSeries> fpss(lg << 1, 1);
		for (int i = 0; i < n; i++) {
			fpss[lg + i] = FormalPowerSeries(vector<modint_for_fps>({-vec[i], 1}));
		}
		for (int i = lg - 1; i > 0; i--) {
			fpss[i] = fpss[i << 1] * fpss[(i << 1) | 1];
		}
		fpss[1] = (*this) % fpss[1];
		for (int i = 2; i < (lg << 1); i++) {
			fpss[i] = fpss[i >> 1] % fpss[i];
		}
		vector<modint_for_fps> ret(n);
		for (int i = 0; i < n; i++) { ret[i] = fpss[lg + i].getCoeff(0); }
		return ret;
	}
	static void getDivMod(const FormalPowerSeries& f, const FormalPowerSeries& g,
						  FormalPowerSeries* outDiv, FormalPowerSeries* outRem) {
		if (f == 0) {
			*outDiv = 0;
			*outRem = 0;
			return;
		}
		int n = f.getDegree(), m = g.getDegree();
		if (m < 0) {
			outDiv->state = -2;
			outRem->state = -2;
			return;
		}
		modint_for_fps ft = f.getCoeff(n), gt = g.getCoeff(m);
		FormalPowerSeries rf, rg;
		if (ft == 0) { rf = f; }
		if (gt == 0) { rg = g; }
		getDivMod_internal((ft == 0) ? rf : f, (gt == 0) ? rg : g, outDiv, outRem);
	}
	static FormalPowerSeries log(const FormalPowerSeries& f) {
		FormalPowerSeries ret = f.getDifferentiated() * f.getInverted();
		ret.setDegree(f.getDegree() - 1);
		ret.integrate();
		return ret;
	}
	static FormalPowerSeries exp(const FormalPowerSeries& f) {
		int n = f.getDegree();
		uint32_t m = 1, lg = getMinBinary(n + 1);
		if (n == -1) { return FormalPowerSeries(1); }
		FormalPowerSeries ret = 1, g = 1, h = f.getCoeff(0), dh = 0, dr = 0;
		while (2 * m <= lg) {
			FormalPowerSeries r = getShiftedMonomialMod(ret * dh, m, 1);
			FormalPowerSeries s = getShiftedMonomialMod((dr - r) << 1, m, 1);
			s *= g;
			s.crop(m);
			s <<= m - 1;
			s.integrate();
			h.extendTo(2 * m - 1);
			dh.extendTo(2 * m - 2);
			for (int i = m; i < (m << 1); i++) {
				h[i] = f.getCoeff(i);
				dh[i - 1] = h[i] * i;
			}
			FormalPowerSeries u = (h - s) >> m;
			u *= ret;
			u.crop(m);
			ret += u << m;
			dr.extendTo(2 * m - 2);
			for (int i = m; i < (m << 1); i++) { dr[i - 1] = ret.getCoeff(i) * i; }
			g *= 2 - ret * g;
			m <<= 1;
			g.crop(m);
		}
		ret.setDegree(n);
		ret.regularize();
		return ret;
	}
	static FormalPowerSeries pow(const FormalPowerSeries& f, int m, int modDeg = INT32_MAX) {
		return pow(f, (int64_t)m, modDeg);
	}
	static FormalPowerSeries pow(const FormalPowerSeries& f, int64_t m, int modDeg = INT32_MAX) {
		if (f.getCoeff(0) == 0 && m < 0) {
			FormalPowerSeries ret;
			ret.state = -1;
			return ret;
		}
		if (m >= 0) { return pow(f, (uint64_t)m, modDeg); }
		FormalPowerSeries g = f.getInverted();
		int n = f.getDegree();
		g.setDegree(n);
		int md = (modDeg < 0) ? n + 1 : modDeg;
		modint_for_fps tv = f[0];
		modint_for_fps ptv = modint_for_fps::pow(tv, m);
		return ptv * pow_internal(g / tv, (m == INT64_MIN) ? (uint64_t)INT64_MAX + 1U : -m, md);
	}
	static FormalPowerSeries pow(const FormalPowerSeries& f, uint32_t m, int modDeg = INT32_MAX) {
		return pow(f, (uint64_t)m, modDeg);
	}
	static FormalPowerSeries pow(const FormalPowerSeries& f, uint64_t m, int modDeg = INT32_MAX) {
		int n = f.getDegree();
		int t = n + 1;
		uint64_t s = 0;
		modint_for_fps tv;
		for (int i = 0; i <= n; i++) {
			if (f[i] != 0) {
				t = i;
				tv = f[i];
				break;
			}
			s += m;
			if (s >= modDeg) { return FormalPowerSeries(); }
		}
		if (t == n + 1) {
			if (m == 0) { return FormalPowerSeries(1).crop(modDeg); }
			else { return FormalPowerSeries(); }
		}
		int md = (int)min((uint64_t)n * m + 1, (uint64_t)modDeg);
		int ls = (int)min(s, (uint64_t)INT32_MAX);
		modint_for_fps ptv = modint_for_fps::pow(tv, m);
		FormalPowerSeries fps = (f >> t) / tv;
		fps.setDegree(n - t);
		return (ptv * pow_internal(fps, m, md - ls)) << ls;
	}
	static FormalPowerSeries sqrt(const FormalPowerSeries& f) {
		int n = f.getDegree();
		int b = 0;
		while (b <= n) {
			if (f.getCoeff(b) != 0) { break; }
			b++;
		}
		if (b == n + 1) { return f; }
		if (b & 1) {
			FormalPowerSeries ret;
			ret.state = -2;
			return ret;
		}
		FormalPowerSeries rf = f >> b;
		rf.setDegree(n);
		rf = sqrt_internal(rf);
		if (rf.state < 0) { return rf; }
		return rf << (b >> 1);
	}
	static FormalPowerSeries interpolate(const vector<modint_for_fps>& x,
										 const vector<modint_for_fps>& y) {
		assert(x.size() == y.size());
		int n = (int)x.size();
		uint32_t lg = getMinBinary(n);
		vector<FormalPowerSeries> fpss1(lg << 1, 1);
		vector<FormalPowerSeries> fpss2(lg << 1);
		for (int i = 0; i < n; i++) {
			fpss1[i + lg] = FormalPowerSeries(vector<modint_for_fps>({-x[i], 1}));
		}
		for (int i = lg - 1; i >= 1; i--) {
			fpss1[i] = fpss1[i << 1] * fpss1[(i << 1) | 1];
		}
		fpss2[1] = fpss1[1].getDifferentiated();
		for (int i = 2; i < (lg << 1); i++) {
			fpss2[i] = fpss2[i >> 1] % fpss1[i];
		}
		for (int i = 0; i < n; i++) {
			fpss2[i + lg] = y[i] / fpss2[i + lg].getCoeff(0);
		}
		for (int i = lg - 1; i >= 1; i--) {
			fpss2[i] = fpss2[i << 1] * fpss1[(i << 1) | 1] + fpss2[(i << 1) | 1] * fpss1[i << 1];
		}
		fpss2[1].setDegree(n - 1);
		return fpss2[1];
	}
	static FormalPowerSeries getProduct(const vector<FormalPowerSeries>& polys) {
		int n = (int)polys.size();
		if (n == 0) {
			return 1;
		}
		if (n == 1) {
			return polys[0];
		}
		vector<FormalPowerSeries> temp(n);
		priority_queue<pair<int, int>> pq;
		int pos = 0;
		for (int i = 0; i < n; i++) { pq.push(make_pair(-polys[i].getDegree(), i)); }
		while (pq.size() > 1) {
			pair<int, int> r1 = pq.top();
			pq.pop();
			pair<int, int> r2 = pq.top();
			pq.pop();
			int idx1 = r1.second, idx2 = r2.second;
			int idx;
			if (idx1 < n && idx2 < n) {
				idx = pos;
				pos++;
			} else if (idx1 >= n) {
				idx = idx1 - n;
			} else {
				idx = idx2 - n;
			}
			const FormalPowerSeries* f1 = (idx1 < n) ? &polys[idx1] : &temp[idx1 - n];
			const FormalPowerSeries* f2 = (idx2 < n) ? &polys[idx2] : &temp[idx2 - n];
			temp[idx] = (*f1) * (*f2);
			pq.push(make_pair(r1.first + r2.first, idx + n));
		}
		pair<int, int> r = pq.top();
		return temp[r.second - n];
	}
	static FormalPowerSeries getShiftedMonomialDiv(const FormalPowerSeries& f,
												   int m, const modint_for_fps& a) {
		if (m <= 0) { return f / (-a); }
		int n = f.getDegree();
		if (n < m) { return FormalPowerSeries(); }
		FormalPowerSeries ret;
		ret.setDegree(n - m);
		for (int i = m; i <= n; i++) {
			int t = (i - m) / m, r = i % m;
			modint_for_fps c = 1;
			for (int j = t; j >= 0; j--) {
				ret[r + j * m] += c * f[i];
				c *= a;
			}
		}
		ret.regularize();
		return ret;
	}
	static FormalPowerSeries getShiftedMonomialMod(const FormalPowerSeries& f,
													   int m, const modint_for_fps& a) {
		int n = f.getDegree();
		if (m <= 0 || n < 0) { return FormalPowerSeries(); }
		vector<modint_for_fps> t(n / m + 1);
		t[0] = 1;
		for (int i = 1; i <= n / m; i++) { t[i] = t[i - 1] * a; }
		FormalPowerSeries ret;
		ret.setDegree(m - 1);
		for (int i = 0; i <= n; i++) { ret[i % m] += t[i / m] * f[i]; }
		ret.regularize();
		return ret;
	}
	static vector<modint_for_fps> getShiftedSamplingPoints(const vector<modint_for_fps>& val,
														   const modint_for_fps& c,
														   int len) {
		int n = (int)val.size();
		vector<modint_for_fps> fact, invf;
		modint_for_fps::get_fact_table(max(n, len), &fact, &invf);
		vector<modint_for_fps> prod(n, 1);
		modint_for_fps t = 1;
		for (int i = 1; i < n; i++) {
			prod[i] = t * (c - i + 1);
			t = prod[i];
			prod[i] *= invf[i];
		}
		vector<modint_for_fps> a(n), b(n);
		for (int i = 0; i < n; i++) {
			a[i] = val[i] * invf[i];
			b[i] = ((i & 1) ? -1 : 1) * invf[i];
		}
		modint_for_fps::convolve(&a, b);
		vector<modint_for_fps> p(n), q(len);
		for (int i = 0; i < n; i++) {
			p[i] = a[n - 1 - i] * fact[n - 1 - i];
		}
		modint_for_fps::convolve(&p, prod);
		for (int i = 0; i < min(n, len); i++) {
			q[i] = p[n - i - 1] * invf[i];
		}
		modint_for_fps::convolve(&q, invf);
		vector<modint_for_fps> ret(len);
		modint_for_fps v = 1;
		for (int i = 0; i < len; i++) {
			ret[i] = q[i] * v;
			v *= i + 1;
		}
		return ret;
	}
	static FormalPowerSeries getRegularExp(int deg) {
		vector<modint_for_fps> invf;
		modint_for_fps::get_fact_table(deg, nullptr, &invf);
		return FormalPowerSeries(invf);
	}
	static FormalPowerSeries getBinomialExpansion(int n, int lim = -1) {
		FormalPowerSeries ret;
		if (n >= 0) {
			ret = getPositiveBinomialExpansion(n, lim);
		} else {
			ret = getNegativeBinomialExpansion(-n, lim);
		}
		return ret;
	}
	static FormalPowerSeries getMonomial(int deg) {
		if (deg < 0) { return FormalPowerSeries(); }
		FormalPowerSeries ret;
		ret[deg] = 1;
		return ret;
	}
	static FormalPowerSeries getBernoulli(int deg) {
		FormalPowerSeries fps = getRegularExp(deg + 1) >> 1;
		fps.invert();
		fps.setDegree(deg);
		modint_for_fps v = 1;
		for (int i = 0; i <= deg; i++) {
			fps[i] *= v;
			v *= i + 1;
		}
		return fps;
	}
	static FormalPowerSeries getPartitions(int deg) {
		FormalPowerSeries fps;
		fps.setDegree(deg);
		vector<modint_for_fps> invn;
		modint_for_fps::get_inv_table(deg, &invn);
		for (int i = 1; i <= deg; i++) {
			for (int j = 0; i * j <= deg; j++) {
				fps[i * j] += invn[j];
			}
		}
		return exp(fps);
	}
	static FormalPowerSeries getStirlingFirst(int n) {
		uint32_t lg = getMinBinary(n);
		uint32_t b = lg >> 1;
		vector<FormalPowerSeries> fpss(lg, 1);
		for (int i = 0; i < n; i++) {
			fpss[i] = FormalPowerSeries(vector<int>({-i, 1}));
		}
		while (b > 0) {
			for (int i = 0; i < b; i++) {
				fpss[i] = fpss[i << 1] * fpss[(i << 1) | 1];
			}
			b >>= 1;
		}
		return fpss[0];
	}
	static FormalPowerSeries getStirlingSecond(int n) {
		vector<modint_for_fps> invf;
		modint_for_fps::get_fact_table(n, nullptr, &invf);
		FormalPowerSeries f(invf);
		f.setDegree(n);
		for (int i = 0; i <= n; i++) { f[i] *= modint_for_fps::pow(i, n); }
		FormalPowerSeries g(invf);
		g.setDegree(n);
		for (int i = 1; i <= n; i += 2) { g[i] = -g[i]; }
		f *= g;
		f.setDegree(n);
		return f;
	}
	static FormalPowerSeries countSubsetSum(const vector<uint32_t>& s, uint32_t lim) {
		int n = (int)s.size();
		vector<uint32_t> cnts(1);
		vector<modint_for_fps> invn;
		modint_for_fps::get_inv_table(lim, &invn);
		uint32_t mx = 0;
		for (int i = 0; i < n; i++) {
			if (mx < s[i]) {
				cnts.resize(s[i] + 1);
				mx = s[i];
			}
			cnts[s[i]]++;
		}
		FormalPowerSeries ret;
		ret.setDegree(lim);
		modint_for_fps c =  modint_for_fps::pow(2, cnts[0]);
		for (int i = 1; i <= mx; i++) {
			for (int j = 1; i * j <= lim; j++) {
				ret[i * j] += ((j & 1) ? 1 : -1) * invn[j] * cnts[i];
			}
		}
		return c * exp(ret);
	}
	static FormalPowerSeries gcd(const FormalPowerSeries& a, const FormalPowerSeries& b) {
		FormalPowerSeries x = a, y = b;
		if (x.getDegree() <= y.getDegree()) {
			FormalPowerSeries tmp = x;
			x = y;
			y = tmp % y;
		}
		while (y) {
			FormalPowerSeries p, q, r, s;
			hgcd(x, y, &p, &q, &r, &s);
			FormalPowerSeries tmp = x;
			x = s * x - q * y;
			y = p * y - r * tmp;
			x.regularize();
			y.regularize();
			if (y == 0) { break; }
			tmp = x;
			x = y;
			y = tmp % y;
		}
		return x;
	}
	static FormalPowerSeries getPolynomialInverse(const FormalPowerSeries& a,
												  const FormalPowerSeries& b) {
		FormalPowerSeries x = 1, y = 0, u = 0, v = 1, k = a, l = b;
		if (k.getDegree() <= l.getDegree()) {
			FormalPowerSeries xTmp = x, yTmp = y;
			FormalPowerSeries q, m;
			getDivMod(k, l, &q, &m);
			x = u;
			y = v;
			u = xTmp - u * q;
			v = yTmp - v * q;
			k = l;
			l = m;
		}
		while (l) {
			FormalPowerSeries p, q, r, s;
			hgcd(k, l, &p, &q, &r, &s);
			modint_for_fps det = (p * s - q * r).getCoeff(0);
			assert(det - 1 == 0 || det + 1 == 0);
			FormalPowerSeries xTmp = x, yTmp = y;
			x = s * x - q * u;
			y = s * y - q * v;
			u = p * u - r * xTmp;
			v = p * yTmp - r * v;
			FormalPowerSeries kTmp = k;
			k = s * k - q * l;
			l = p * l - r * kTmp;
			k.regularize();
			l.regularize();
			if (l == 0) { break; }
			FormalPowerSeries qu, md;
			getDivMod(k, l, &qu, &md);
			xTmp = x;
			yTmp = y;
			x = u;
			y = v;
			u = xTmp - u * qu;
			v = yTmp - v * md;
			k = l;
			l = md;
		}
		if (k.getDegree() > 0) {
			FormalPowerSeries ret;
			ret.state = -2;
			return ret;
		}
		modint_for_fps iv = k.getCoeff(0).inv();
		return x * iv;
	}
};

#endif /* FormalPowerSeries_hpp */
