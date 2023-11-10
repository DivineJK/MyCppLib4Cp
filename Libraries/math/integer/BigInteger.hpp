//
//  BigInteger.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef BigInt_hpp
#define BigInt_hpp

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
} // bint_internal

class BigInteger {
private:
	template <int M>
	struct internal_mint {
	private:
		uint32_t x;
	public:
		
	};
private:
	static constexpr int DIGITS_BASE = 1000000;
	vector<uint32_t> digits;
private:
	void shrink_digits();
public:
	BigInteger() {}
	BigInteger(bool b);
	BigInteger(int8_t v);
	BigInteger(int16_t v);
	BigInteger(int v);
	BigInteger(int64_t v);
	BigInteger(uint8_t v);
	BigInteger(uint16_t v);
	BigInteger(uint32_t v);
	BigInteger(uint64_t v);
	BigInteger(const string& s);
	BigInteger(const BigInteger& other);
	BigInteger& operator=(bool b);
	BigInteger& operator=(int8_t v);
	BigInteger& operator=(int16_t v);
	BigInteger& operator=(int v);
	BigInteger& operator=(int64_t v);
	BigInteger& operator=(uint8_t v);
	BigInteger& operator=(uint16_t v);
	BigInteger& operator=(uint32_t v);
	BigInteger& operator=(uint64_t v);
	BigInteger& operator=(const string& s);
	BigInteger& operator=(const BigInteger& other);
	explicit operator bool() const;
	bool operator!() const;
	bool operator==(const BigInteger& other) const;
	bool operator!=(const BigInteger& other) const;
	bool operator<=(const BigInteger& other) const;
	bool operator<(const BigInteger& other) const;
	bool operator>=(const BigInteger& other) const;
	bool operator>(const BigInteger& other) const;
	BigInteger operator+() const;
	BigInteger operator-() const;
	BigInteger& operator++();
	BigInteger operator++(int);
	BigInteger& operator--();
	BigInteger operator--(int);
	BigInteger& operator+=(const BigInteger& other);
	BigInteger& operator-=(const BigInteger& other);
	BigInteger& operator*=(const BigInteger& other);
	BigInteger& operator/=(const BigInteger& other);
	BigInteger& operator%=(const BigInteger& other);
	template <typename T>
	BigInteger& operator^=(const T& m);
	BigInteger operator+(const BigInteger& other) const;
	BigInteger operator-(const BigInteger& other) const;
	BigInteger operator*(const BigInteger& other) const;
	BigInteger operator/(const BigInteger& other) const;
	BigInteger operator%(const BigInteger& other) const;
	template <typename T>
	BigInteger operator^(const T& m) const;
	template <typename T>
	friend BigInteger operator+(const T& lhs, const BigInteger& rhs);
	template <typename T>
	friend BigInteger operator-(const T& lhs, const BigInteger& rhs);
	template <typename T>
	friend BigInteger operator*(const T& lhs, const BigInteger& rhs);
	template <typename T>
	friend BigInteger operator/(const T& lhs, const BigInteger& rhs);
	template <typename T>
	friend BigInteger operator%(const T& lhs, const BigInteger& rhs);
	BigInteger& operator>>=(int a);
	BigInteger operator>>(int a) const;
	BigInteger& operator<<=(int a);
	BigInteger operator<<(int a) const;
	friend ostream& operator<<(ostream& os, const BigInteger& other);
	friend istream& operator<<(istream& ist, BigInteger& other);
	string to_string() const;
	template <typename T>
	static BigInteger pow(const BigInteger& bint, const T& m);
	template <typename T, typename S>
	static BigInteger pow(const BigInteger& bint, const T& m, const S& mod);
	static BigInteger fact(int n);
	static BigInteger gcd(const BigInteger& x, const BigInteger& y);
	static BigInteger lcm(const BigInteger& x, const BigInteger& y);
	static BigInteger min(const BigInteger& x, const BigInteger& y);
	static BigInteger max(const BigInteger& x, const BigInteger& y);
	static pair<BigInteger, BigInteger> div_mod(const BigInteger& x, const BigInteger& y);
};

#endif /* BigInt_hpp */
