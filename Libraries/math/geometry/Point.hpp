//
//  Point.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef Point_hpp
#define Point_hpp

template <typename T, size_t DIM>
class Point {
private:
	T elms[DIM];
public:
	Point() { initialize(); }
	Point(const Point& other) { operator=(other); }
	void initialize() {
		for (int i = 0; i < DIM; i++) { elms[i] = (T)0; }
	}
	Point& operator=(const Point& other) {
		for (int i = 0; i < DIM; i++) { elms[i] = other.elms[i]; }
		return *this;
	}
	T& operator[](int i) { return elms[i]; }
	const T& operator[](int i) const { return elms[i]; }
	bool operator==(const Point& other) const {
		for (int i = 0; i < DIM; i++) {
			if (elms[i] != other.elms[i]) { return false; }
		}
		return true;
	}
	bool operator!=(const Point& other) const { return !(*this == other); }
	Point& operator+=(const Point& other) {
		for (int i = 0; i < DIM; i++) { elms[i] += other[i]; }
		return *this;
	}
	Point& operator-=(const Point& other) {
		for (int i = 0; i < DIM; i++) { elms[i] -= other[i]; }
		return *this;
	}
	Point& operator*=(const T& other) {
		for (int i = 0; i < DIM; i++) { elms[i] *= other; }
		return *this;
	}
	Point& operator/=(const T& other) {
		for (int i = 0; i < DIM; i++) { elms[i] /= other; }
		return *this;
	}
	Point operator+(const Point& other) const { return Point(*this) += other; }
	Point operator-(const Point& other) const { return Point(*this) -= other; }
	Point operator*(const T& other) const { return Point(*this) *= other; }
	friend Point operator*(const T& lhs, const Point& rhs) {
		Point ret = rhs;
		for (int i = 0; i < DIM; i++) { ret.elms[i] *= lhs; }
		return ret;
	}
	Point operator/(const T& other) const { return Point(*this) /= other; }
	T dot(const Point& other) const {
		T ret = (T)0;
		for (int i = 0; i < DIM; i++) { ret += elms[i] * other.elms[i]; }
		return ret;
	}
	static T dot(const Point& lhs, const Point& rhs) { return lhs.dot(rhs); }
};

template <typename T>
class Point2D : public Point<T, 2> {
public:
	Point2D(T x, T y) {
		(*this)[0] = x;
		(*this)[1] = y;
	}
	T cross(const Point2D& other) const {
		return (*this)[0] * other[1] - (*this)[1] * other[0];
	}
	static T cross(const Point2D& lhs, const Point2D& rhs) { return lhs.cross(rhs); }
};

template <typename T>
class Point3D : public Point<T, 3> {
public:
	Point3D(T x, T y, T z) {
		(*this)[0] = x;
		(*this)[1] = y;
		(*this)[2] = z;
	}
	Point3D cross(const Point3D& other) const {
		T x = (*this)[1] * other[2] - (*this)[2] * other[1];
		T y = (*this)[2] * other[0] - (*this)[0] * other[2];
		T z = (*this)[0] * other[1] - (*this)[1] * other[0];
		return Point3D(x, y, z);
	}
	static Point3D cross(const Point3D& lhs, const Point3D& rhs) { return lhs.cross(rhs); }
};

#endif /* Point_hpp */
