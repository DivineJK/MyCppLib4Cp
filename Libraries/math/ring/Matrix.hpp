//
//  Matrix.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/09.
//

#ifndef Matrix_hpp
#define Matrix_hpp

namespace matrix_internal {
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
} // namespace matrix_internal

template <int MOD>
using modint_for_matrix = matrix_internal::internal_modint<MOD>;

template <typename T>
class LinearEquationSolution {
public:
	bool is_valid = true;
	vector<T> primalSolution;
	vector<vector<T>> kernelBases;
};

template <typename T>
class Matrix {
private:
	// state ==  0: valid
	// state == -1: infinity by divide by zero determinant
	// state == -2: invalid
	int state = 0;
	vector<vector<T>> elms_2d;
public:
	Matrix() : state(0) {}
	Matrix(uint32_t aRow, uint32_t aColumn)
	: state(0) {
		elms_2d = vector<vector<T>>(aRow, vector<T>(aColumn, (T)0));
	}
	Matrix(uint32_t aRow, uint32_t aColumn, const T& v)
	: state(0) {
		elms_2d = vector<vector<T>>(aRow, vector<T>(aColumn, v));
	}
	Matrix(uint32_t aRow, uint32_t aColumn, const vector<T>& aElms)
	: state(0) {
		assert(aElms.size() == aRow * aColumn);
		elms_2d = vector<vector<T>>(aRow, vector<T>(aColumn));
		for (uint32_t i = 0; i < aRow; i++) {
			for (uint32_t j = 0; j < aColumn; j++) {
				elms_2d[i][j] = aElms[i * aColumn + j];
			}
		}
	}
	Matrix(const vector<vector<T>>& aMat) { operator=(aMat); }
	Matrix(const Matrix& other) { operator=(other); }
	Matrix& operator=(const vector<vector<T>>& aMat) {
		state = 0;
		elms_2d = aMat;
		return *this;
	}
	Matrix& operator=(const Matrix& other) {
		state = other.state;
		elms_2d = other.elms_2d;
		return *this;
	}
	friend ostream& operator<<(ostream& os, const Matrix& mat) {
		if (mat.state < 0) {
			os << -1;
			return os;
		}
		uint32_t r = mat.getRow();
		uint32_t c = mat.getColumn();
		if (r == 0 || c == 0) { return os; }
		os << mat.elms_2d[0][0];
		for (uint32_t j = 1; j < c; j++) { os << " " << mat.elms_2d[0][j]; }
		for (uint32_t i = 1; i < r; i++) {
			os << '\n' << mat.elms_2d[i][0];
			for (uint32_t j = 1; j < c; j++) { os << " " << mat.elms_2d[i][j]; }
		}
		return os;
	}
	friend istream& operator>>(istream& ist, Matrix& mat) {
		uint32_t r = mat.getRow();
		uint32_t c = mat.getColumn();
		for (uint32_t i = 0; i < r; i++) {
			for (uint32_t j = 0; j < c; j++) { ist >> mat.elms_2d[i][j]; }
		}
		return ist;
	}
	vector<T>& operator[](uint32_t i) {
		assert(0 <= i && i < getRow());
		return elms_2d[i];
	}
	const vector<T>& operator[](uint32_t i) const {
		assert(0 <= i && i < getRow());
		return elms_2d[i];
	}
	bool operator==(const Matrix& other) const {
		if (state != other.state) { return false; }
		return state < 0 || elms_2d == other.elms_2d;
	}
	bool operator!=(const Matrix& other) const { return !(*this == other); }
	Matrix operator+() const { return Matrix(*this); }
	Matrix operator-() const {
		Matrix mat = Matrix(*this);
		for (uint32_t i = 0; i < mat.getRow(); i++) {
			for (uint32_t j = 0; j < mat.getColumn(); j++) {
				mat.elms_2d[i][j] = -mat.elms_2d[i][j];
			}
		}
		return mat;
	}
	Matrix& operator+=(const Matrix& other) {
		assert(getRow() == other.getRow() && getColumn() == other.getColumn());
		uint32_t r = getRow();
		uint32_t c = getColumn();
		for (uint32_t i = 0; i < r; i++) {
			for (uint32_t j = 0; j < c; j++) { elms_2d[i][j] += other.elms_2d[i][j]; }
		}
		return *this;
	}
	Matrix& operator-=(const Matrix& other) {
		assert(getRow() == other.getRow() && getColumn() == other.getColumn());
		uint32_t r = getRow();
		uint32_t c = getColumn();
		for (uint32_t i = 0; i < r; i++) {
			for (uint32_t j = 0; j < c; j++) { elms_2d[i][j] -= other.elms_2d[i][j]; }
		}
		return *this;
	}
	Matrix& operator*=(const Matrix& other) {
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
		vector<T> temp(this_column);
		for (uint32_t i = 0; i < this_row; i++) {
			for (uint32_t j = 0; j < this_column; j++) { temp[j] = elms_2d[i][j]; }
			elms_2d[i].resize(other_column);
			for (uint32_t j = 0; j < other_column; j++) {
				elms_2d[i][j] = 0;
				for (uint32_t k = 0; k < this_column; k++) {
					elms_2d[i][j] += temp[k] * other.elms_2d[k][j];
				}
			}
		}
		return *this;
	}
	Matrix& operator*=(const T& other) {
		uint32_t r = getRow();
		uint32_t c = getColumn();
		for (uint32_t i = 0; i < r; i++) {
			for (uint32_t j = 0; j < c; j++) { elms_2d[i] *= other; }
		}
		return *this;
	}
	Matrix& operator/=(const Matrix& other) { return *this *= other.getInverted(); }
	Matrix& operator/=(const T& other) {
		uint32_t r = getRow();
		uint32_t c = getColumn();
		for (uint32_t i = 0; i < r; i++) {
			for (uint32_t j = 0; j < c; j++) { elms_2d[i] /= other; }
		}
		return *this;
	}
	Matrix operator+(const Matrix& other) const { return Matrix(*this) += other; }
	Matrix operator-(const Matrix& other) const { return Matrix(*this) -= other; }
	Matrix operator*(const Matrix& other) const { return Matrix(*this) *= other; }
	Matrix operator*(const T& other) const { return Matrix(*this) *= other; }
	Matrix operator/(const Matrix& other) const { return Matrix(*this) /= other; }
	Matrix operator/(const T& other) const { return Matrix(*this) /= other; }
	friend Matrix operator*(const T& lhs, const Matrix& rhs) { return rhs * lhs; }
	template <typename S>
	Matrix operator^(S aExp) const { return pow(*this, aExp); }
	uint32_t getRow() const { return (uint32_t)elms_2d.size(); }
	uint32_t getColumn() const {
		return ((uint32_t)elms_2d.size() == 0) ? 0 : (uint32_t)elms_2d[0].size();
	}
	void setRow(uint32_t aR) {
		uint32_t r = getRow();
		uint32_t c = getColumn();
		elms_2d.resize(aR);
		if (r < aR) {
			for (uint32_t i = r; i < aR; i++) { elms_2d[i] = vector<T>(c, 0); }
		}
	}
	void setColumn(uint32_t aC) {
		uint32_t r = getRow();
		for (uint32_t i = 0; i < r; i++) { elms_2d[i].resize(aC); }
	}
	void setSize(uint32_t aR, uint32_t aC) {
		setRow(aR);
		setColumn(aC);
	}
	void setElement(uint32_t i, uint32_t j, const T& v) {
		if (i >= getRow()) { setRow(i + 1); }
		if (j >= getColumn()) { setColumn(j + 1); }
		elms_2d[i][j] = v;
	}
	void setRowElements(uint32_t r, const vector<T>& v) {
		assert(v.size() == getColumn());
		if (r >= getRow()) { setRow(r + 1); }
		elms_2d[r] = v;
	}
	void setColumnElements(uint32_t c, const vector<T>& v) {
		uint32_t r = getRow();
		assert(v.size() == r);
		if (c >= getColumn()) { setColumn(c + 1); }
		for (uint32_t i = 0; i < r; i++) { elms_2d[i][c] = v[i]; }
	}
	T getElement(uint32_t i, uint32_t j) const {
		assert(0 <= i && i < getRow() && 0 <= j && j < getColumn());
		return elms_2d[i][j];
	}
	vector<T> getRowElements(uint32_t r) const {
		assert(0 <= r && r < getRow());
		return elms_2d[r];
	}
	vector<T> getColumnElements(uint32_t c) const {
		assert(0 <= c && c < getColumn());
		uint32_t r = getRow();
		vector<T> ret(r);
		for (uint32_t i = 0; i < r; i++) { ret[i] = elms_2d[i][c]; }
		return ret;
	}
	void swapRow(uint32_t i, uint32_t j) {
		assert(0 <= i && i < getRow() && 0 <= j && j < getRow());
		if (i == j) { return; }
		uint32_t c = getColumn();
		for (uint32_t k = 0; k < c; k++) {
			T t = elms_2d[i][k];
			elms_2d[i][k] = elms_2d[j][k];
			elms_2d[j][k] = t;
		}
	}
	void swapColumn(uint32_t i, uint32_t j) {
		assert(0 <= i && i < getColumn()() && 0 <= j && j < getColumn());
		if (i == j) { return; }
		uint32_t r = getRow();
		for (uint32_t k = 0; k < r; k++) {
			T t = elms_2d[k][i];
			elms_2d[k][i] = elms_2d[k][j];
			elms_2d[k][j] = t;
		}
	}
	Matrix& transpose() {
		uint32_t r = getRow();
		uint32_t c = getColumn();
		if (r < c) {
			setRow(c);
			for (uint32_t i = 0; i < r; i++) {
				for (uint32_t j = r; j < c; j++) { elms_2d[j][i] = elms_2d[i][j]; }
			}
		}
		else if (r > c) {
			setColumn(r);
			for (uint32_t j = 0; j < c; j++) {
				for (uint32_t i = c; i < r; i++) { elms_2d[j][i] = elms_2d[i][j]; }
			}
		}
		uint32_t m = min(r, c);
		for (uint32_t i = 0; i < m; i++) {
			for (uint32_t j = i + 1; j < m; j++) {
				T tmp = elms_2d[i][j];
				elms_2d[i][j] = elms_2d[j][i];
				elms_2d[j][i] = tmp;
			}
		}
		if (r < c) { setColumn(r); }
		else if (r > c) { setRow(c); }
		return *this;
	}
	Matrix getTransposed() const {
		uint32_t r = getRow();
		uint32_t c = getColumn();
		Matrix ret = Matrix(c, r);
		ret.state = state;
		for (uint32_t i = 0; i < r; i++) {
			for (uint32_t j = 0; j < c; j++) { ret.elms_2d[j][i] = elms_2d[i][j]; }
		}
		return ret;
	}
	Matrix& invert();
	Matrix getInverted() const;
	virtual LinearEquationSolution<T> solveLinearEquationsSystem(const vector<T>& b) const {
		LinearEquationSolution<T> ret;
		return ret;
	}
	static Matrix getIdentity(uint32_t s) {
		Matrix ret(s, s);
		for (uint32_t i = 0; i < s; i++) { ret.elms_2d[i][i] = 1; }
		return ret;
	}
	static T det(const Matrix& mat) {
		assert(mat.getRow() == mat.getColumn());
		uint32_t n = (uint32_t)mat.getRow();
		Matrix m = mat;
		T ret = 1;
		vector<int> rows(mat.getRow());
		for (uint32_t i = 0; i < n; i++) {
			if (m[i][i] == 0) {
				ret = -ret;
				uint32_t j = i;
				while (j < n && m[j][i] == 0) { j++; }
				if (j == n) { return 0; }
				m.swapRow(i, j);
				
			}
			T v = m[i][i];
			ret *= v;
			for (uint32_t j = i + 1; j < n; j++) {
				T pv = m[j][i];
				for (uint32_t k = i + 1; k < n; k++) {
					m[j][k] -= m[i][k] * pv / v;
				}
			}
		}
		return ret;
	}
	template <typename S>
	static Matrix pow(const Matrix& mat, S aExp) {
		assert(mat.row == mat.column);
		Matrix ret = getIdentity(mat.row);
		Matrix b = (aExp < 0) ? mat : mat.getInverted();
		while (aExp) {
			if (aExp & 1) { ret *= b; }
			b *= b;
			aExp >>= 1;
		}
		return ret;
	}
	static Matrix getKroneckerProduct(const Matrix& mat1, const Matrix& mat2) {
		uint32_t r1 = mat1.getRow();
		uint32_t r2 = mat2.getRow();
		uint32_t c1 = mat1.getColumn();
		uint32_t c2 = mat2.getColumn();
		Matrix ret(r1 * r2, c1 * c2);
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

template <int MOD>
class ModintMatrix : public Matrix<modint_for_matrix<MOD>> {
	using mint = modint_for_matrix<MOD>;
public:
	ModintMatrix() : Matrix<mint>() {}
	ModintMatrix(uint32_t aRow, uint32_t aColumn)
	: Matrix<mint>(aRow, aColumn) {}
	ModintMatrix(uint32_t aRow, uint32_t aColumn, const vector<mint>& aElms)
	: Matrix<mint>(aRow, aColumn, aElms) {}
	ModintMatrix(const vector<vector<mint>>& aMat) : Matrix<mint>(aMat) {}
	ModintMatrix(const ModintMatrix& other) : Matrix<mint>(other) {}
	static mint det(const ModintMatrix& mat) {
		assert(mat.getRow() == mat.getColumn());
		uint32_t n = (uint32_t)mat.getRow();
		Matrix m = mat;
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
			mint v = m[i][i];
			mint iv = v.inv();
			for (uint32_t j = i; j < n; j++) { m[i][j] *= iv; }
			ret *= v;
			for (uint32_t j = i + 1; j < n; j++) {
				mint pv = m[j][i];
				for (uint32_t k = i + 1; k < n; k++) {
					m[j][k] -= m[i][k] * pv;
				}
			}
		}
		return ret;
	}
	virtual LinearEquationSolution<mint> solveLinearEquationsSystem(const vector<mint>& b) const override {
		LinearEquationSolution<mint> ret;
		return ret;
	}
};

#endif /* Matrix_hpp */
