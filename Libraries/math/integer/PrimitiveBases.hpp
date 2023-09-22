//
//  PrimitiveBases.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/18.
//

#ifndef PrimitiveBases_hpp
#define PrimitiveBases_hpp

struct PrimitiveBases {
private:
	uint64_t mod = 1;
	int base_size = 0;
	uint64_t roots[2]; // primitive root, inverse root
	vector<uint64_t> bases[3]; // primitive bases, inverse bases, cumulative bases
private:
	PrimitiveBases() = delete;
	template <typename T>
	bool is_primitive(T n) const {
		T t = n;
		for (int i = 0; i < base_size; i++) {
			if (t == 1) { return false; }
			t = (T)((uint64_t)t * t % mod);
		}
		return t == 1;
	}
	template <typename T>
	T modinv(T n) const {
		T p = n, t = 1;
		while (p > 1) {
			t = (uint64_t)t * (mod - (mod / p)) % mod;
			p = mod % p;
		}
		return t;
	}
	template <typename T>
	static constexpr bool is_prime(T n) {
		if (n <= 1) { return false; }
		if (n != 2 && !(n & 1)) { return false; }
		T p = 3;
		while (p * p <= n) {
			if (n % p == 0) { return false; }
			p += 2;
		}
		return true;
	}
public:
	PrimitiveBases(uint32_t aMod, bool make_bas = true, bool chk_pr = true) {
		assert(!chk_pr || is_prime(aMod));
		mod = aMod;
		make_base(make_bas);
	}
	PrimitiveBases(uint64_t aMod, bool make_bas = true, bool chk_pr = true) {
		assert(!chk_pr || is_prime(aMod));
		mod = aMod;
		make_base(make_bas);
	}
	void make_base(bool make_bas) {
		uint64_t r = 1, p = mod - 1;
		base_size = 0;
		while (!(p & 1)) {
			p >>= 1;
			base_size++;
		}
		if (make_bas) {
			for (int i = 0; i < 3; i++) { bases[i] = vector<uint64_t>(base_size + 1); }
		}
		while (true) {
			if (is_primitive(r)) {
				roots[0] = r;
				roots[1] = modinv(r);
				if (!make_bas) { break; }
				bases[0][base_size] = r;
				bases[1][base_size] = roots[1];
				for (int i = base_size; i > 0; i--) {
					bases[0][i - 1] = (uint64_t)bases[0][i] * bases[0][i] % mod;
					bases[1][i - 1] = (uint64_t)bases[1][i] * bases[1][i] % mod;
				}
				uint64_t ie = 1;
				for (int i = 0; i < base_size - 1; i++) {
					bases[2][i] = (uint64_t)ie * bases[0][i + 2] % mod;
					ie = (uint64_t)ie * bases[1][i + 2] % mod;
				}
				break;
			}
			r++;
		}
	}
	uint64_t get_base(int t, int idx) const {
		assert(0 <= t || t < 3);
		return (idx < 0 || idx > base_size) ? 1 : bases[t][idx];
	}
	uint64_t get_root(int t) const {
		assert(t == 0 || t == 1);
		return roots[t];
	}
	int get_base_size() const { return base_size; }
};

#endif /* PrimitiveBases_hpp */
