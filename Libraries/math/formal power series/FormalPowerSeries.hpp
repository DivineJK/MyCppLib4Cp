//
//  FormalPowerSeries.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/09.
//

#ifndef FormalPowerSeries_hpp
#define FormalPowerSeries_hpp

namespace formalpowerseries_internal {
constexpr int INTERNAL_MOD = 998244353; // 998244353 = 119 * 2 ^ 23 + 1
class internal_modint {
	using modint = internal_modint;
private:
	uint32_t x = 0;
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
	modint& inv() {
		if (x == 0) { return *this; }
		uint32_t a = x;
		x = 1;
		while (a > 1) {
			x = (uint64_t)x * (INTERNAL_MOD - INTERNAL_MOD / a) % INTERNAL_MOD;
			a = INTERNAL_MOD % a;
		}
		return *this;
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
};

}	// namespace formalpowerseries_internal

static constexpr int MOD_FOR_FPS = formalpowerseries_internal::INTERNAL_MOD;
using modint_for_fps = formalpowerseries_internal::internal_modint;

class FormalPowerSeries {
private:
	vector<modint_for_fps> poly;
private:
	
public:
	FormalPowerSeries() {}
	FormalPowerSeries(const modint_for_fps& mint) {
		poly = (mint) ? vector<modint_for_fps>({mint}) : vector<modint_for_fps>();
	}
	FormalPowerSeries(const vector<modint_for_fps>& arr) { poly = arr; }
	FormalPowerSeries(const FormalPowerSeries& fps) { operator=(fps); }
	FormalPowerSeries& operator=(const FormalPowerSeries& fps) {
		poly = fps.poly;
		return *this;
	}
};

class FormalPowerSeries_Sparse {
	
};

#endif /* FormalPowerSeries_hpp */
