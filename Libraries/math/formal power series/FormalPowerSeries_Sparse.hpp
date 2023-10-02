//
//  FormalPowerSeries_Sparse.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/30.
//

#ifndef FormalPowerSeries_Sparse_hpp
#define FormalPowerSeries_Sparse_hpp

namespace fpss_internal {
constexpr int INTERNAL_MOD = 998244353; // 998244353 = 119 * 2 ^ 23 + 1
class internal_modint {
	using modint = internal_modint;
private:
	uint32_t x = 0;
private:
public:
	internal_modint() : x(0) {}
	internal_modint(bool b) { x = (b) ? 1 : 0; }
	internal_modint(int a) {
		int v = a % INTERNAL_MOD;
		x = (v < 0) ? INTERNAL_MOD + v : v;
	}
	internal_modint(int64_t a) {
		int v = a % INTERNAL_MOD;
		x = (v < 0) ? INTERNAL_MOD + v : v;
	}
	internal_modint(uint32_t a) { x = a % INTERNAL_MOD; }
	internal_modint(uint64_t a) { x = a % INTERNAL_MOD; }
	internal_modint(const modint& mint) { operator=(mint); }
	uint32_t val() const { return x; }
	modint& operator=(const modint& mint) {
		x = mint.x;
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

}	// namespace formalpowerseries_sparse_internal

static constexpr int MOD_FOR_FPSS = fpss_internal::INTERNAL_MOD;
using modint_for_fpss = fpss_internal::internal_modint;

class FormalPowerSeries_Sparse {
private:
	// state ==  0: valid
	// state == -1: infinity by divide by zero
	// state == -2: invalid
	int state = 0;
	// degree_state == 0: dynamic
	// degree_state == 1: fixed
	// degree_state == 2: recalculate needed
	mutable int degree_state = 0;
	mutable int degree = -1;
	unordered_map<int, modint_for_fpss> coeffs;
private:
	void recalculateDegree() const {
		for (const pair<int, modint_for_fpss> c : coeffs) { degree = max(degree, c.first); }
		degree_state = 0;
	}
	void regularize() {
		if (degree_state == 1) {
			unordered_map<int, modint_for_fpss>::iterator it = coeffs.begin();
			while (it != coeffs.end()) {
				int idx = it->first;
				if (idx > degree) {
					it = coeffs.erase(it);
					continue;
				}
				++it;
			}
		} else if (degree_state == 0 || degree_state == 2) { recalculateDegree(); }
	}
	void initialize() {
		state = 0;
		degree_state = 0;
		degree = -1;
		coeffs.clear();
	}
	static FormalPowerSeries_Sparse pow_internal(const FormalPowerSeries_Sparse& fps,
												 uint64_t m, int modDeg = INT32_MAX) {
		modint_for_fpss b = fps.getCoeff(0);
		modint_for_fpss r0 = modint_for_fpss::pow(b, m);
		FormalPowerSeries_Sparse ret = r0;
		vector<modint_for_fpss> invn;
		modint_for_fpss::get_inv_table(modDeg, &invn);
		unordered_map<int, modint_for_fpss> temp;
		modint_for_fpss z = b.inv();
		for (int i = 1; i < modDeg; i++) {
			unordered_map<int, modint_for_fpss>::const_iterator it;
			modint_for_fpss v = 0;
			for (it = temp.begin(); it != temp.end(); ++it) {
				int idx = it->first;
				modint_for_fpss c = m;
				c *= idx;
				c += idx - i;
				v += it->second * ret.getCoeff(i - idx) * c * invn[i];
			}
			modint_for_fpss a = fps.getCoeff(i);
			if (a != 0) {
				v += a * r0 * m;
				temp[i] = a;
			}
			ret[i] = v * z;
		}
		return ret;
	}
	static FormalPowerSeries_Sparse sqrt_internal(const FormalPowerSeries_Sparse& fps) {
		modint_for_fpss b = fps.getCoeff(0);
		int z = modint_for_fpss::sqrt(b);
		FormalPowerSeries_Sparse ret;
		if (z == -1) {
			ret.state = -2;
			return ret;
		}
		int n = fps.getDegree();
		ret[0] = z;
		vector<modint_for_fpss> invn;
		modint_for_fpss::get_inv_table(n, &invn);
		unordered_map<int, modint_for_fpss> temp;
		modint_for_fpss i2 = modint_for_fpss::inv(2);
		modint_for_fpss ib = b.inv();
		for (int i = 1; i <= n; i++) {
			unordered_map<int, modint_for_fpss>::const_iterator it;
			modint_for_fpss v = 0;
			for (it = temp.begin(); it != temp.end(); ++it) {
				int idx = it->first;
				modint_for_fpss c = idx - (modint_for_fpss)2 * (i - idx);
				v += it->second * ret.getCoeff(i - idx) * c * invn[i];
			}
			modint_for_fpss a = fps.getCoeff(i);
			if (a != 0) {
				v += a * z;
				temp[i] = a;
			}
			ret[i] = v * i2 * ib;
		}
		return ret;
	}
public:
	FormalPowerSeries_Sparse() {}
	template <typename T>
	FormalPowerSeries_Sparse(const T& other) { operator=(other); }
	FormalPowerSeries_Sparse& operator=(int v) { return operator=((modint_for_fpss)v); }
	FormalPowerSeries_Sparse& operator=(uint32_t v) { return operator=((modint_for_fpss)v); }
	FormalPowerSeries_Sparse& operator=(int64_t v) { return operator=((modint_for_fpss)v); }
	FormalPowerSeries_Sparse& operator=(uint64_t v) { return operator=((modint_for_fpss)v); }
	FormalPowerSeries_Sparse& operator=(const unordered_map<int, modint_for_fpss>& aCoeffs) {
		initialize();
		coeffs.reserve(aCoeffs.size());
		for (const pair<int, modint_for_fpss> c : aCoeffs) {
			if (c.second != 0) { coeffs[c.first] = c.second; }
			degree = max(degree, c.first);
		}
		return *this;
	}
	FormalPowerSeries_Sparse& operator=(const modint_for_fpss& v) {
		initialize();
		if (v != 0) { coeffs[0] = v; }
		degree = 0;
		return *this;
	}
	template <typename T>
	FormalPowerSeries_Sparse& operator=(const vector<T>& v) {
		initialize();
		coeffs.reserve(v.size());
		for (int i = 0; i <= v.size(); i++) {
			if (v[i] == 0) { continue; }
			coeffs[i] = v[i];
			degree = i;
		}
		return *this;
	}
	FormalPowerSeries_Sparse& operator=(const FormalPowerSeries_Sparse& v) {
		state = v.state;
		degree_state = v.degree_state;
		degree = v.degree;
		coeffs = v.coeffs;
		return *this;
	}
	friend ostream& operator<<(ostream& os, const FormalPowerSeries_Sparse& fps) {
		if (fps.state < 0) {
			os << -1;
			return os;
		}
		int n = fps.getDegree();
		for (int i = 0; i <= n; i++) {
			if (i > 0) { os << " "; }
			os << fps.getCoeff(i);
		}
		return os;
	}
	friend istream& operator>>(istream& ist, FormalPowerSeries_Sparse& fps) {
		int i;
		modint_for_fpss v;
		ist >> i >> v;
		fps.setCoeff(i, v);
		return ist;
	}
	explicit operator bool() const {
		if (state < 0) { return state == -1; }
		return coeffs.size() > 0;
	}
	bool operator!() const {
		if (state < 0) { return state == -2; }
		return coeffs.size() == 0;
	}
	bool operator==(const FormalPowerSeries_Sparse& other) const {
		if (state < 0) {
			if (state == -1) { return other.state == -1; }
			else if (state == -2) { return false; }
		}
		unordered_set<int> s;
		s.reserve(coeffs.size());
		for (const pair<int, modint_for_fpss> p : coeffs) { s.insert(p.first); }
		for (const pair<int, modint_for_fpss> p : other.coeffs) {
			if (!s.count(p.first)) { return false; }
		}
		for (int i : s) {
			if (getCoeff(i) != other.getCoeff(i)) { return false; }
		}
		return true;
	}
	bool operator!=(const FormalPowerSeries_Sparse& other) const {
		if (state < 0) {
			if (state == -1) { return other.state != -1; }
			else if (state == -2) { return false; }
		}
		unordered_set<int> s;
		s.reserve(coeffs.size());
		for (const pair<int, modint_for_fpss> p : coeffs) { s.insert(p.first); }
		for (const pair<int, modint_for_fpss> p : other.coeffs) {
			if (!s.count(p.first)) { return true; }
		}
		for (int i : s) {
			if (getCoeff(i) != other.getCoeff(i)) { return true; }
		}
		return false;
	}
	FormalPowerSeries_Sparse operator+() const { return FormalPowerSeries_Sparse(*this); }
	FormalPowerSeries_Sparse operator-() const {
		FormalPowerSeries_Sparse ret = FormalPowerSeries_Sparse(*this);
		unordered_map<int, modint_for_fpss>::iterator it;
		for (it = ret.coeffs.begin(); it != ret.coeffs.end(); ++it) {
			it->second = -it->second;
		}
		return ret;
	}
	FormalPowerSeries_Sparse& operator++() {
		if (getCoeff(0) == 0) { coeffs[0] = 0; }
		++coeffs[0];
		regularize();
		return *this;
	}
	FormalPowerSeries_Sparse& operator--() {
		if (getCoeff(0) == 0) { coeffs[0] = 0; }
		--coeffs[0];
		regularize();
		return *this;
	}
	FormalPowerSeries_Sparse operator++(int) {
		FormalPowerSeries_Sparse ret = FormalPowerSeries_Sparse(*this);
		++*this;
		return ret;
	}
	FormalPowerSeries_Sparse operator--(int) {
		FormalPowerSeries_Sparse ret = FormalPowerSeries_Sparse(*this);
		--*this;
		return ret;
	}
	FormalPowerSeries_Sparse& operator+=(const FormalPowerSeries_Sparse& other) {
		unordered_map<int, modint_for_fpss>::const_iterator it;
		for (it = other.coeffs.begin(); it != other.coeffs.end(); ++it) {
			int idx = it->first;
			modint_for_fpss val = it->second;
			if (coeffs.count(idx)) { coeffs[idx] += val; }
			else { coeffs[idx] = val; }
			if (coeffs[idx] == 0) { coeffs.erase(idx); }
		}
		regularize();
		return *this;
	}
	FormalPowerSeries_Sparse& operator-=(const FormalPowerSeries_Sparse& other) {
		unordered_map<int, modint_for_fpss>::const_iterator it;
		for (it = other.coeffs.begin(); it != other.coeffs.end(); ++it) {
			int idx = it->first;
			modint_for_fpss val = it->second;
			if (coeffs.count(idx)) { coeffs[idx] -= val; }
			else { coeffs[idx] = -val; }
			if (coeffs[idx] == 0) { coeffs.erase(idx); }
		}
		regularize();
		return *this;
	}
	FormalPowerSeries_Sparse& operator*=(const FormalPowerSeries_Sparse& other) {
		unordered_map<int, modint_for_fpss> tmp = coeffs;
		coeffs.clear();
		for (pair<int, modint_for_fpss> p1 : tmp) {
			int idx1 = p1.first;
			modint_for_fpss val1 = p1.second;
			for (pair<int, modint_for_fpss> p2 : other.coeffs) {
				int idx2 = p2.first;
				modint_for_fpss val2 = p2.second;
				if (coeffs.count(idx1 + idx2)) { coeffs[idx1 + idx2] += val1 * val2; }
				else { coeffs[idx1 + idx2] = val1 * val2; }
				if (coeffs[idx1 + idx2] == 0) { coeffs.erase(idx1 + idx2); }
			}
		}
		regularize();
		return *this;
	}
	FormalPowerSeries_Sparse& operator/=(const FormalPowerSeries_Sparse& other);
	FormalPowerSeries_Sparse& operator%=(const FormalPowerSeries_Sparse& other);
	FormalPowerSeries_Sparse operator+(const FormalPowerSeries_Sparse& other) const {
		FormalPowerSeries_Sparse ret = FormalPowerSeries_Sparse(*this);
		ret += other;
		return ret;
	}
	FormalPowerSeries_Sparse operator-(const FormalPowerSeries_Sparse& other) const {
		FormalPowerSeries_Sparse ret = FormalPowerSeries_Sparse(*this);
		ret -= other;
		return ret;
	}
	FormalPowerSeries_Sparse operator*(const FormalPowerSeries_Sparse& other) const {
		FormalPowerSeries_Sparse ret = FormalPowerSeries_Sparse(*this);
		ret *= other;
		return ret;
	}
	FormalPowerSeries_Sparse operator/(const FormalPowerSeries_Sparse& other) const;
	FormalPowerSeries_Sparse operator%(const FormalPowerSeries_Sparse& other) const;
	FormalPowerSeries_Sparse& operator<<=(int k) {
		unordered_map<int, modint_for_fpss> tmp = coeffs;
		coeffs.clear();
		unordered_map<int, modint_for_fpss>::iterator it;
		for (it = tmp.begin(); it != tmp.end(); ++it) {
			int idx = it->first;
			coeffs[idx + k] = it->second;
		}
		if (degree_state == 0) { degree += k; }
		return *this;
	}
	FormalPowerSeries_Sparse& operator>>=(int k) {
		unordered_map<int, modint_for_fpss> tmp = coeffs;
		coeffs.clear();
		unordered_map<int, modint_for_fpss>::iterator it;
		for (it = tmp.begin(); it != tmp.end(); ++it) {
			int idx = it->first;
			if (idx < k) { continue; }
			coeffs[idx - k] = it->second;
		}
		if (degree_state == 0) { degree = max(-1, degree + k); }
		return *this;
	}
	FormalPowerSeries_Sparse operator<<(int k) const {
		FormalPowerSeries_Sparse ret;
		unordered_map<int, modint_for_fpss>::const_iterator it;
		for (it = coeffs.begin(); it != coeffs.end(); ++it) {
			int idx = it->first;
			ret[idx + k] = it->second;
		}
		ret.regularize();
		ret.state = state;
		return ret;
	}
	FormalPowerSeries_Sparse operator>>(int k) const {
		FormalPowerSeries_Sparse ret;
		unordered_map<int, modint_for_fpss>::const_iterator it;
		for (it = coeffs.begin(); it != coeffs.end(); ++it) {
			int idx = it->first;
			if (idx < k) { continue; }
			ret[idx - k] = it->second;
		}
		ret.regularize();
		ret.state = state;
		return ret;
	}
	template <typename T>
	FormalPowerSeries_Sparse operator^(T m) const { return pow(*this, m); }
	modint_for_fpss& operator[](int i) { return coeffs[i]; }
	const modint_for_fpss& operator[](int i) const { return coeffs.at(i); }
	void setDegree(int deg) {
		degree_state = 0;
		degree = deg;
	}
	void fixDegree(int deg) {
		degree_state = 1;
		degree = deg;
	}
	int getDegree() const {
		if (degree_state == 2) { recalculateDegree(); }
		return degree;
	}
	void setCoeff(int idx, const modint_for_fpss& v) {
		coeffs[idx] = v;
		if (v == 0) { coeffs.erase(idx); }
		if (degree_state == 0) { degree = max(idx, degree); }
	}
	modint_for_fpss getCoeff(int idx) const {
		if (!coeffs.count(idx)) { return 0; }
		return coeffs.at(idx);
	}
	modint_for_fpss evaluate(const modint_for_fpss& v) const {
		modint_for_fpss ret = 0;
		for (const pair<int, modint_for_fpss> p : coeffs) {
			ret += p.second * modint_for_fpss::pow(v, p.first);
		}
		return ret;
	}
	FormalPowerSeries_Sparse& differentiate() {
		unordered_map<int, modint_for_fpss> tmp = coeffs;
		coeffs.clear();
		for (const pair<int, modint_for_fpss> p : tmp) {
			int idx = p.first;
			modint_for_fpss val = p.second;
			if (idx * val == 0) { continue; }
			coeffs[idx - 1] = idx * val;
		}
		regularize();
		return *this;
	}
	FormalPowerSeries_Sparse& integrate() {
		unordered_map<int, modint_for_fpss> tmp = coeffs;
		coeffs.clear();
		for (const pair<int, modint_for_fpss> p : tmp) {
			int idx = p.first;
			modint_for_fpss val = p.second;
			if ((modint_for_fpss)idx + 1 == 0) { continue; }
			modint_for_fpss iv = modint_for_fpss::inv(idx + 1);
			if (iv * val == 0) { continue; }
			coeffs[idx + 1] = iv * val;
		}
		regularize();
		return *this;
	}
	FormalPowerSeries_Sparse& invert() {
		int n = getDegree();
		if (getCoeff(0) == 0) {
			state = -1;
			return *this;
		}
		modint_for_fpss z = getCoeff(0).inv();
		unordered_map<int, modint_for_fpss> tmp = coeffs;
		coeffs.clear();
		unordered_map<int, modint_for_fpss> refs;
		refs.reserve(tmp.size());
		coeffs[0] = z;
		for (int i = 1; i <= n; i++) {
			if (tmp.count(i)) { refs[i] = tmp[i]; }
			unordered_map<int, modint_for_fpss>::const_iterator it;
			modint_for_fpss val = 0;
			for (it = refs.begin(); it != refs.end(); ++it) {
				int idx = it->first;
				val += it->second * getCoeff(i - idx);
			}
			if (val == 0) { continue; }
			val *= -z;
			coeffs[i] = val;
		}
		return *this;
	}
	FormalPowerSeries_Sparse& crop(int n);
	FormalPowerSeries_Sparse getDifferentiated() const {
		return FormalPowerSeries_Sparse(*this).differentiate();
	}
	FormalPowerSeries_Sparse getIntegrated() const {
		return FormalPowerSeries_Sparse(*this).integrate();
	}
	FormalPowerSeries_Sparse getInverted() const {
		int n = getDegree();
		if (getCoeff(0) == 0) {
			FormalPowerSeries_Sparse ret;
			ret.state = -1;
			return ret;
		}
		modint_for_fpss z = getCoeff(0).inv();
		FormalPowerSeries_Sparse ret;
		ret.state = state;
		ret.fixDegree(n);
		unordered_map<int, modint_for_fpss> refs;
		refs.reserve(coeffs.size());
		ret[0] = z;
		for (int i = 1; i <= n; i++) {
			if (coeffs.count(i)) { refs[i] = getCoeff(i); }
			unordered_map<int, modint_for_fpss>::const_iterator it;
			modint_for_fpss val = 0;
			for (it = refs.begin(); it != refs.end(); ++it) {
				int idx = it->first;
				val += it->second * ret.getCoeff(i - idx);
			}
			if (val == 0) { continue; }
			val *= -z;
			ret[i] = val;
		}
		return ret;
	}
	FormalPowerSeries_Sparse getCropped(int n) const;
	static FormalPowerSeries_Sparse log(const FormalPowerSeries_Sparse& fps) {
		FormalPowerSeries_Sparse ret = fps.getDifferentiated() * fps.getInverted();
		ret.integrate();
		return ret;
	}
	static FormalPowerSeries_Sparse exp(const FormalPowerSeries_Sparse& fps) {
		FormalPowerSeries_Sparse ret = 1;
		int n = fps.getDegree();
		ret.fixDegree(n);
		vector<modint_for_fpss> invn;
		modint_for_fpss::get_inv_table(n, &invn);
		for (int i = 1; i <= n; i++) {
			modint_for_fpss v = 0;
			unordered_map<int, modint_for_fpss>::const_iterator it;
			for (it = fps.coeffs.begin(); it != fps.coeffs.end(); ++it) {
				int idx = it->first;
				if (idx > i) { continue; }
				v += idx * it->second * ret.getCoeff(i - idx);
			}
			ret[i] = v * invn[i];
		}
		return ret;
	}
	template <typename T>
	static FormalPowerSeries_Sparse pow(const FormalPowerSeries_Sparse& fps,
										T m, int modDeg = INT32_MAX) {
		int k = -1, mt = -1;
		unordered_map<int, modint_for_fpss>::const_iterator it;
		for (it = fps.coeffs.begin(); it != fps.coeffs.end(); ++it) {
			if (k == -1 || k > it->first) { k = it->first; }
			if (mt < it->first) { mt = it->first; }
		}
		if (k > 0 && m >= modDeg) { return 0; }
		assert(m >= 0 || k == 0);
		if (k >= modDeg) { return (m == 0) ? 1 : 0; }
		FormalPowerSeries_Sparse rf = fps;
		if (k > 0) { rf >>= k; }
		if (m < 0) { rf.invert(); }
		int md = modDeg;
		if (0 < m && m <= INT32_MAX) {
			int t = modDeg / m;
			if (mt <= t) { md = mt * (int)m + 1; }
		} else if (INT32_MIN < m && m < 0) {
			int t = modDeg / (-m);
			if (mt <= t) { md = mt * (int)-m + 1; }
		}
		if (md < k) { return (m == 0) ? 1 : 0; }
		rf = pow_internal(rf, (uint64_t)((m >= 0) ? m : (-m)), md - k);
		if (k > 0) { rf <<= k * (int)m; }
		rf.fixDegree(modDeg - 1);
		return rf;
	}
	static FormalPowerSeries_Sparse sqrt(const FormalPowerSeries_Sparse& fps) {
		int k = -1;
		unordered_map<int, modint_for_fpss>::const_iterator it;
		for (it = fps.coeffs.begin(); it != fps.coeffs.end(); ++it) {
			if (k == -1 || k > it->first) { k = it->first; }
		}
		int n = fps.getDegree();
		if (k == -1) { return 0; }
		if (k & 1) {
			FormalPowerSeries_Sparse ret;
			ret.state = -2;
			return ret;
		}
		FormalPowerSeries_Sparse rf = fps >> k;
		rf.fixDegree(n);
		return sqrt_internal(rf) << (k >> 1);
	}
};

#endif /* FormalPowerSeries_Sparse_hpp */
