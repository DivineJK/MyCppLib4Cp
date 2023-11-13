//
//  ModintMatrix.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/10/15.
//

#ifndef ModintMatrix_hpp
#define ModintMatrix_hpp

namespace mintmat_internal {
template <int MOD>
class internal_modint {
	using modint = internal_modint;
private:
	uint32_t x = 0;
private:
	static constexpr int mod() { return MOD; }
	static constexpr uint32_t umod() { return (uint32_t)MOD; }
public:
	internal_modint() : x(0) {}
	internal_modint(bool b) { x = (b) ? 1 : 0; }
	internal_modint(int a) {
		int v = a % mod();
		x = (v < 0) ? mod() + v : v;
	}
	internal_modint(int64_t a) {
		int v = a % mod();
		x = (v < 0) ? mod() + v : v;
	}
	internal_modint(uint32_t a) { x = a % umod(); }
	internal_modint(uint64_t a) { x = a % umod(); }
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
			t = (uint64_t)t * (umod() - umod() / a) % umod();
			a = umod() % a;
		}
		return modint(t);
	}
	modint extended_inv() const {
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
	static modint inv(const modint& mint) {
		if (mint.x == 0) { return mint; }
		uint32_t a = mint.x, t = 1;
		while (a > 1) {
			t = (uint64_t)t * (umod() - umod() / a) % umod();
			a = umod() % a;
		}
		return modint(t);
	}
	static modint extended_inv(const modint& v) {
		return v.extended_inv();
	}
	modint operator+() const { return modint(*this); }
	modint operator-() const {
		modint ret = modint(*this);
		ret.x = (ret.x == 0) ? ret.x : umod() - ret.x;
		return ret;
	}
	modint& operator++() {
		x++;
		if (x == umod()) { x = 0; }
		return *this;
	}
	modint operator++(int) {
		modint ret = modint(*this);
		++*this;
		return ret;
	}
	modint& operator--() {
		x = (x == 0) ? umod() - 1 : x - 1;
		return *this;
	}
	modint operator--(int) {
		modint ret = modint(*this);
		--*this;
		return ret;
	}
	modint& operator+=(const modint& mint) {
		x = (x + mint.x >= umod()) ? x + mint.x - umod() : x + mint.x;
		return *this;
	}
	modint& operator-=(const modint& mint) {
		x = (x < mint.x) ? umod() + x - mint.x : x - mint.x;
		return *this;
	}
	modint& operator*=(const modint& mint) {
		x = (uint64_t)x * mint.x % umod();
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
		if (pow(n, umod() >> 1) + 1 == 0) { return -1; }
		int m = umod() - 1, q = 0;
		while (!(m & 1)) {
			q++;
			m >>= 1;
		}
		modint z = 1, r = pow(n, (m + 1) >> 1), t = pow(n, m);
		while (pow(z, umod() >> 1) == 1) { z++; }
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
		return (r.x > umod() - r.x) ? umod() - r.x : r.x;
	}
};
} // namespace mintmat_internal

template <int MOD>
using modint_for_matrix = mintmat_internal::internal_modint<MOD>;

template <int MOD>
class ModintMatrix {
	using mint = modint_for_matrix<MOD>;
private:
	// state ==  0: valid
	// state == -1: infinity by divide by zero determinant
	// state == -2: invalid
	int state = 0;
	vector<vector<mint>> elms;
public:
	ModintMatrix() : state(0) {}
	ModintMatrix(uint32_t aRow, uint32_t aColumn)
	: state(0) {
		elms = vector<vector<mint>>(aRow, vector<mint>(aColumn, 0));
	}
	ModintMatrix(uint32_t aRow, uint32_t aColumn, const mint& v)
	: state(0) {
		elms = vector<vector<mint>>(aRow, vector<mint>(aColumn, v));
	}
	ModintMatrix(uint32_t aRow, uint32_t aColumn, const vector<mint>& aElms)
	: state(0) {
		assert(aElms.size() == aRow * aColumn);
		elms = vector<vector<mint>>(aRow, vector<mint>(aColumn));
		for (uint32_t i = 0; i < aRow; i++) {
			for (uint32_t j = 0; j < aColumn; j++) {
				elms[i][j] = aElms[i * aColumn + j];
			}
		}
	}
	ModintMatrix(const vector<vector<mint>>& aMat) { operator=(aMat); }
	ModintMatrix(const ModintMatrix& other) { operator=(other); }
	ModintMatrix& operator=(const vector<vector<mint>>& aMat) {
		state = 0;
		elms = aMat;
		return *this;
	}
	ModintMatrix& operator=(const ModintMatrix& other) {
		state = other.state;
		elms = other.elms;
		return *this;
	}
	friend ostream& operator<<(ostream& os, const ModintMatrix& mat) {
		if (mat.state < 0) {
			os << -1;
			return os;
		}
		uint32_t r = mat.getRow();
		uint32_t c = mat.getColumn();
		if (r == 0 || c == 0) { return os; }
		os << mat.elms[0][0];
		for (uint32_t j = 1; j < c; j++) { os << " " << mat.elms[0][j]; }
		for (uint32_t i = 1; i < r; i++) {
			os << '\n' << mat.elms[i][0];
			for (uint32_t j = 1; j < c; j++) { os << " " << mat.elms[i][j]; }
		}
		return os;
	}
	friend istream& operator>>(istream& ist, ModintMatrix& mat) {
		uint32_t r = mat.getRow();
		uint32_t c = mat.getColumn();
		for (uint32_t i = 0; i < r; i++) {
			for (uint32_t j = 0; j < c; j++) { ist >> mat.elms[i][j]; }
		}
		return ist;
	}
	vector<mint>& operator[](uint32_t i) {
		assert(0 <= i && i < getRow());
		return elms[i];
	}
	const vector<mint>& operator[](uint32_t i) const {
		assert(0 <= i && i < getRow());
		return elms[i];
	}
	bool operator==(const ModintMatrix& other) const {
		if (state != other.state) { return false; }
		return state < 0 || elms == other.elms;
	}
	bool operator!=(const ModintMatrix& other) const { return !(*this == other); }
	ModintMatrix operator+() const { return ModintMatrix(*this); }
	ModintMatrix operator-() const {
		ModintMatrix mat = ModintMatrix(*this);
		for (uint32_t i = 0; i < mat.getRow(); i++) {
			for (uint32_t j = 0; j < mat.getColumn(); j++) {
				mat.elms[i][j] = -mat.elms[i][j];
			}
		}
		return mat;
	}
	ModintMatrix& operator+=(const ModintMatrix& other) {
		assert(getRow() == other.getRow() && getColumn() == other.getColumn());
		uint32_t r = getRow();
		uint32_t c = getColumn();
		for (uint32_t i = 0; i < r; i++) {
			for (uint32_t j = 0; j < c; j++) { elms[i][j] += other.elms[i][j]; }
		}
		return *this;
	}
	ModintMatrix& operator-=(const ModintMatrix& other) {
		assert(getRow() == other.getRow() && getColumn() == other.getColumn());
		uint32_t r = getRow();
		uint32_t c = getColumn();
		for (uint32_t i = 0; i < r; i++) {
			for (uint32_t j = 0; j < c; j++) { elms[i][j] -= other.elms[i][j]; }
		}
		return *this;
	}
	ModintMatrix& operator*=(const ModintMatrix& other) {
		assert(getColumn() == other.getRow());
		if (state == -1 || other.state == -1) {
			state = -1;
			return *this;
		}
		if (state < 0 || other.state < 0) {
			state = -2;
			return *this;
		}
		uint32_t this_row = getRow();
		uint32_t this_column = getColumn();
		uint32_t other_column = other.getColumn();
		vector<mint> temp(this_column);
		for (uint32_t i = 0; i < this_row; i++) {
			for (uint32_t j = 0; j < this_column; j++) { temp[j] = elms[i][j]; }
			elms[i].resize(other_column);
			for (uint32_t j = 0; j < other_column; j++) {
				elms[i][j] = 0;
				for (uint32_t k = 0; k < this_column; k++) {
					elms[i][j] += temp[k] * other.elms[k][j];
				}
			}
		}
		return *this;
	}
	ModintMatrix& operator*=(const mint& other) {
		uint32_t r = getRow();
		uint32_t c = getColumn();
		for (uint32_t i = 0; i < r; i++) {
			for (uint32_t j = 0; j < c; j++) { elms[i] *= other; }
		}
		return *this;
	}
	ModintMatrix& operator/=(const ModintMatrix& other) { return *this *= other.getInverted(); }
	ModintMatrix& operator/=(const mint& other) {
		uint32_t r = getRow();
		uint32_t c = getColumn();
		mint iv = other.getInverted();
		for (uint32_t i = 0; i < r; i++) {
			for (uint32_t j = 0; j < c; j++) { elms[i] *= iv; }
		}
		return *this;
	}
	ModintMatrix operator+(const ModintMatrix& other) const { return ModintMatrix(*this) += other; }
	ModintMatrix operator-(const ModintMatrix& other) const { return ModintMatrix(*this) -= other; }
	ModintMatrix operator*(const ModintMatrix& other) const { return ModintMatrix(*this) *= other; }
	ModintMatrix operator*(const mint& other) const { return ModintMatrix(*this) *= other; }
	ModintMatrix operator/(const ModintMatrix& other) const { return Matrix(*this) /= other; }
	ModintMatrix operator/(const mint& other) const { return ModintMatrix(*this) /= other; }
	friend ModintMatrix operator*(const mint& lhs, const ModintMatrix& rhs) { return rhs * lhs; }
	template <typename S>
	ModintMatrix operator^(S aExp) const { return pow(*this, aExp); }
	uint32_t getRow() const { return (uint32_t)elms.size(); }
	uint32_t getColumn() const {
		return ((uint32_t)elms.size() == 0) ? 0 : (uint32_t)elms[0].size();
	}
	void setRow(uint32_t aR) {
		uint32_t r = getRow();
		uint32_t c = getColumn();
		elms.resize(aR);
		if (r < aR) {
			for (uint32_t i = r; i < aR; i++) { elms[i] = vector<mint>(c, 0); }
		}
	}
	void setColumn(uint32_t aC) {
		uint32_t r = getRow();
		for (uint32_t i = 0; i < r; i++) { elms[i].resize(aC); }
	}
	void setSize(uint32_t aR, uint32_t aC) {
		setRow(aR);
		setColumn(aC);
	}
	void setElement(uint32_t i, uint32_t j, const mint& v) {
		if (i >= getRow()) { setRow(i + 1); }
		if (j >= getColumn()) { setColumn(j + 1); }
		elms[i][j] = v;
	}
	void setRowElements(uint32_t r, const vector<mint>& v) {
		assert(v.size() == getColumn());
		if (r >= getRow()) { setRow(r + 1); }
		elms[r] = v;
	}
	void setColumnElements(uint32_t c, const vector<mint>& v) {
		uint32_t r = getRow();
		assert(v.size() == r);
		if (c >= getColumn()) { setColumn(c + 1); }
		for (uint32_t i = 0; i < r; i++) { elms[i][c] = v[i]; }
	}
	mint getElement(uint32_t i, uint32_t j) const {
		assert(0 <= i && i < getRow() && 0 <= j && j < getColumn());
		return elms[i][j];
	}
	vector<mint> getRowElements(uint32_t r) const {
		assert(0 <= r && r < getRow());
		return elms[r];
	}
	vector<mint> getColumnElements(uint32_t c) const {
		assert(0 <= c && c < getColumn());
		uint32_t r = getRow();
		vector<mint> ret(r);
		for (uint32_t i = 0; i < r; i++) { ret[i] = elms[i][c]; }
		return ret;
	}
	void swapRow(uint32_t i, uint32_t j) {
		assert(0 <= i && i < getRow() && 0 <= j && j < getRow());
		if (i == j) { return; }
		uint32_t c = getColumn();
		for (uint32_t k = 0; k < c; k++) {
			mint t = elms[i][k];
			elms[i][k] = elms[j][k];
			elms[j][k] = t;
		}
	}
	void swapColumn(uint32_t i, uint32_t j) {
		assert(0 <= i && i < getColumn()() && 0 <= j && j < getColumn());
		if (i == j) { return; }
		uint32_t r = getRow();
		for (uint32_t k = 0; k < r; k++) {
			mint t = elms[k][i];
			elms[k][i] = elms[k][j];
			elms[k][j] = t;
		}
	}
	ModintMatrix& transpose() {
		uint32_t r = getRow();
		uint32_t c = getColumn();
		if (r < c) {
			setRow(c);
			for (uint32_t i = 0; i < r; i++) {
				for (uint32_t j = r; j < c; j++) { elms[j][i] = elms[i][j]; }
			}
		}
		else if (r > c) {
			setColumn(r);
			for (uint32_t j = 0; j < c; j++) {
				for (uint32_t i = c; i < r; i++) { elms[j][i] = elms[i][j]; }
			}
		}
		uint32_t m = min(r, c);
		for (uint32_t i = 0; i < m; i++) {
			for (uint32_t j = i + 1; j < m; j++) {
				mint tmp = elms[i][j];
				elms[i][j] = elms[j][i];
				elms[j][i] = tmp;
			}
		}
		if (r < c) { setColumn(r); }
		else if (r > c) { setRow(c); }
		return *this;
	}
	ModintMatrix getTransposed() const {
		uint32_t r = getRow();
		uint32_t c = getColumn();
		ModintMatrix ret = ModintMatrix(c, r);
		ret.state = state;
		for (uint32_t i = 0; i < r; i++) {
			for (uint32_t j = 0; j < c; j++) { ret.elms[j][i] = elms[i][j]; }
		}
		return ret;
	}
	ModintMatrix& invert() {
		if (state < 0) { return *this; }
		assert(getRow() == getColumn());
		uint32_t n = getRow();
		setColumn(n << 1);
		uint32_t en = n << 1;
		for (uint32_t i = 0; i < n; i++) { setElement(i, i + n, 1); }
		for (uint32_t i = 0; i < n; i++) {
			if ((*this)[i][i] == 0) {
				uint32_t j = i;
				while (j < n && (*this)[j][i] == 0) { j++; }
				if (j == n) {
					state = -1;
					setColumn(n);
					return *this;
				}
				swapRow(i, j);
			}
			mint v = (*this)[i][i];
			mint iv = v.inv();
			(*this)[i][i] = 1;
			for (uint32_t j = i + 1; j < en; j++) { (*this)[i][j] *= iv; }
			for (uint32_t j = i + 1; j < n; j++) {
				mint pv = (*this)[j][i];
				if (pv == 0) { continue; }
				for (uint32_t k = i + 1; k < en; k++) {
					(*this)[j][k] -= (*this)[i][k] * pv;
				}
			}
		}
		for (uint32_t i = n; i > 0; i--) {
			for (uint32_t k = i - 1; k > 0; k--) {
				mint pv = (*this)[k - 1][i - 1];
				if (pv == 0) { continue; }
				for (uint32_t j = n; j < en; j++) {
					(*this)[k - 1][j] -= (*this)[i - 1][j] * pv;
				}
			}
			for (uint32_t j = 0; j < n; j++) { (*this)[i - 1][j] = (*this)[i - 1][j + n]; }
		}
		setColumn(n);
		return *this;
	}
	ModintMatrix getInverted() const { return ModintMatrix(*this).invert(); }
	bool solveLinearEquationsSystem(const vector<mint>& b, vector<mint>* sol,
									vector<vector<mint>>* kernel) const {
		return solveLinearEquationsSystem(*this, b, sol, kernel);
	}
	static bool solveLinearEquationsSystem(const ModintMatrix& mat, const vector<mint>& b,
										   vector<mint>* sol,
										   vector<vector<mint>>* kernel) {
		assert(mat.getRow() == b.size());
		uint32_t r = mat.getRow();
		uint32_t c = mat.getColumn();
		ModintMatrix extMat = ModintMatrix(mat);
		extMat.setColumn(c + 1);
		for (uint32_t i = 0; i < r; i++) { extMat[i][c] = b[i]; }
		uint32_t cr = 0;
		for (uint32_t i = 0; i < c; i++) {
			if (extMat[cr][i] == 0) {
				uint32_t j = cr;
				while (j < r && extMat[j][i] == 0) { j++; }
				if (j == r) { continue; }
				extMat.swapRow(cr, j);
			}
			mint iv = extMat[cr][i].inv();
			extMat[cr][i] = 1;
			for (uint32_t j = i + 1; j <= c; j++) { extMat[cr][j] *= iv; }
			for (uint32_t j = cr + 1; j < r; j++) {
				mint pv = extMat[j][i];
				if (pv == 0) { continue; }
				extMat[j][i] = 0;
				for (uint32_t k = i + 1; k <= c; k++) {
					extMat[j][k] -= extMat[cr][k] * pv;
				}
			}
			cr++;
			if (cr == r) { break; }
		}
		vector<uint32_t> lm(r, c);
		uint32_t idx = 0;
		for (uint32_t i = 0; i < r; i++) {
			while (idx < c && extMat[i][idx] == 0) { idx++; }
			if (idx == c) {
				if (extMat[i][c] != 0) { return false; }
				continue;
			}
			lm[i] = idx;
		}
		*sol = vector<mint>(c, 0);
		vector<bool> lk(c, false);
		for (uint32_t i = r; i > 0; i--) {
			if (lm[i - 1] == c) { continue; }
			uint32_t pos = lm[i - 1];
			for (uint32_t j = 0; j < i - 1; j++) {
				mint v = extMat[j][pos];
				extMat[j][pos] = 0;
				for (uint32_t k = pos + 1; k <= c; k++) { extMat[j][k] -= extMat[i - 1][k] * v; }
			}
			lk[pos] = true;
			(*sol)[pos] = extMat[i - 1][c];
		}
		kernel->clear();
		kernel->reserve(c);
		for (uint32_t i = 0; i < c; i++) {
			if (lk[i]) { continue; }
			vector<mint> ker(c, 0);
			ker[i] = 1;
			for (uint32_t j = 0; j < r; j++) {
				uint32_t lmv = lm[j];
				if (i < lmv) { break; }
				ker[lmv] = -extMat[j][i];
			}
			kernel->push_back(ker);
		}
		return true;
	}
	static ModintMatrix getIdentity(uint32_t s) {
		ModintMatrix ret(s, s);
		for (uint32_t i = 0; i < s; i++) { ret.elms[i][i] = 1; }
		return ret;
	}
	static mint det(const ModintMatrix& mat) {
		assert(mat.getRow() == mat.getColumn());
		uint32_t n = (uint32_t)mat.getRow();
		ModintMatrix m = mat;
		mint ret = 1;
		vector<int> rows(mat.getRow());
		for (uint32_t i = 0; i < n; i++) {
			if (m[i][i] == 0) {
				ret = -ret;
				uint32_t j = i;
				while (j < n && m[j][i] == 0) { j++; }
				if (j == n) { return 0; }
				m.swapRow(i, j);
				
			}
			mint iv = m[i][i].inv();
			ret *= m[i][i];
			for (uint32_t j = i + 1; j < n; j++) {
				mint pv = m[j][i];
				if (pv == 0) { continue; }
				for (uint32_t k = i + 1; k < n; k++) {
					m[j][k] -= m[i][k] * pv * iv;
				}
			}
		}
		return ret;
	}
	template <typename S>
	static ModintMatrix pow(const ModintMatrix& mat, S aExp) {
		assert(mat.getRow() == mat.getColumn());
		ModintMatrix ret = getIdentity(mat.getRow());
		ModintMatrix b = (aExp >= 0) ? mat : mat.getInverted();
		while (aExp) {
			if (aExp & 1) { ret *= b; }
			b *= ModintMatrix(b);
			aExp >>= 1;
		}
		return ret;
	}
	static ModintMatrix getKroneckerProduct(const ModintMatrix& mat1, const ModintMatrix& mat2) {
		uint32_t r1 = mat1.getRow();
		uint32_t r2 = mat2.getRow();
		uint32_t c1 = mat1.getColumn();
		uint32_t c2 = mat2.getColumn();
		ModintMatrix ret(r1 * r2, c1 * c2);
		for (uint32_t i = 0; i < r1; i++) {
			for (uint32_t j = 0; j < r2; j++) {
				for (uint32_t k = 0; k < c1; k++) {
					for (uint32_t l = 0; l < c2; l++) {
						ret[i * r2 + j][k * c2 + l] = mat1[i][k] * mat2[j][l];
					}
				}
			}
		}
		return ret;
	}
};

#endif /* ModintMatrix_hpp */
