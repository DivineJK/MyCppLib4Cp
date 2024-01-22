//
//  GaussInteger.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/12/20.
//

#ifndef GaussInteger_h
#define GaussInteger_h

template <typename T, enable_if_t<is_integral_v<T>>* = nullptr>
class GaussInteger {
private:
	T x;
	T y;
private:
	static constexpr T half_div(T a, T b) {
		assert(b != 0);
		T q = a / b, r = a % b;
		T ab = (b < 0) ? -b : b;
		if (r < 0) {
			q += (b < 0) ? 1 : -1;
			r += ab;
		}
		if (ab - r < r) { q += (b < 0) ? -1 : 1; }
		return q;
	}
public:
	GaussInteger() : x(0), y(0) {}
	GaussInteger(T aX) : x(aX), y(0) {}
	GaussInteger(T aX, T aY) : x(aX), y(aY) {}
	template <typename S, enable_if_t<is_integral_v<S>>* = nullptr>
	GaussInteger& operator=(S other) {
		x = (T)other;
		y = 0;
		return *this;
	}
	GaussInteger& operator=(T other) {
		x = other;
		y = 0;
		return *this;
	}
	GaussInteger& operator=(const GaussInteger& other) {
		x = other.x;
		y = other.y;
		return *this;
	}
	friend ostream& operator<<(ostream& os, const GaussInteger& gint) {
		os << gint.real() << " " << gint.imag();
		return os;
	}
	friend istream& operator>>(istream& ist, GaussInteger& gint) {
		ist >> gint.x >> gint.y;
		return ist;
	}
	explicit operator bool() const { return x != 0 || y != 0; }
	bool operator!() const { return x == 0 && y == 0; }
	GaussInteger operator+() const { return GaussInteger(*this); }
	GaussInteger operator-() const { return GaussInteger(-x, -y); }
	GaussInteger& operator++() {
		++x;
		return *this;
	}
	GaussInteger operator++(int) {
		GaussInteger ret = *this;
		++*this;
		return ret;
	}
	GaussInteger& operator--() {
		--x;
		return *this;
	}
	GaussInteger operator--(int) {
		GaussInteger ret = *this;
		--*this;
		return ret;
	}
	GaussInteger& operator+=(const GaussInteger& other) {
		x += other.x;
		y += other.y;
		return *this;
	}
	GaussInteger& operator-=(const GaussInteger& other) {
		x -= other.x;
		y -= other.y;
		return *this;
	}
	GaussInteger& operator*=(const GaussInteger& other) {
		T s = x;
		x = s * other.x - y * other.y;
		y = s * other.y + y * other.x;
		return *this;
	}
	GaussInteger& operator/=(const T& v) {
		x = half_div(x, v);
		y = half_div(y, v);
		return *this;
	}
	GaussInteger& operator/=(const GaussInteger& other) {
		T n = other.norm();
		*this *= other.conj();
		*this /= n;
		return *this;
	}
	GaussInteger& operator%=(const GaussInteger& other) {
		*this -= other * (*this / other);
		return *this;
	}
	GaussInteger operator+(const GaussInteger& other) const { return GaussInteger(*this) += other;}
	GaussInteger operator-(const GaussInteger& other) const { return GaussInteger(*this) -= other;}
	GaussInteger operator*(const GaussInteger& other) const { return GaussInteger(*this) *= other;}
	GaussInteger operator/(const GaussInteger& other) const { return GaussInteger(*this) /= other;}
	GaussInteger operator%(const GaussInteger& other) const { return GaussInteger(*this) %= other;}
	T norm() const { return x * x + y * y; }
	T real() const { return x; }
	T imag() const { return y; }
	GaussInteger conj() const { return GaussInteger(x, -y); }
	static GaussInteger gcd(const GaussInteger& a, const GaussInteger& b) {
		GaussInteger x = a, y = b, tmp;
		while (y) {
			tmp = x;
			x = y;
			y = tmp % y;
		}
		return x;
	}
};

#endif /* GaussInteger_h */
