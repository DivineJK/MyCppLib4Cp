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
	explicit operator bool() const { return x != 0 || y != 0; }
	bool operator!() const { return x == 0 && y == 0; }
	GaussInteger operator+() const { return GaussInteger(*this); }
	GaussInteger operator~() const {
		GaussInteger ret = *this;
		ret.y = -y;
		return ret;
	}
};

#endif /* GaussInteger_h */
