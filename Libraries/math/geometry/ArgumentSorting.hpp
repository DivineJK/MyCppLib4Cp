//
//  ArgumentSorting.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef ArgumentSorting_hpp
#define ArgumentSorting_hpp

namespace argsort_internal {
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
	Point2D() {}
	Point2D(T x, T y) {
		(*this)[0] = x;
		(*this)[1] = y;
	}
	T cross(const Point2D& other) const {
		return (*this)[0] * other[1] - (*this)[1] * other[0];
	}
	static T cross(const Point2D& lhs, const Point2D& rhs) { return lhs.cross(rhs); }
};
} // namespace argsort_internal

template <typename T>
class ArgumentSorting {
	using Point = argsort_internal::Point2D<T>;
private:
	Point start_vec = Point(1, 0);
	bool is_vec_start = true;
	bool is_ccw_pos = true;
	Point origin_arg = Point(1, 0);
private:
	int getZone(const Point& vec) const {
		Point rv = (vec[0] == 0 && vec[1] == 0) ? origin_arg : vec;
		T cr = start_vec.cross(rv);
		T dt = start_vec.dot(rv);
		assert(cr != 0 || dt != 0);
		if (cr == 0) { return (dt < 0) ? 4 : ((is_vec_start) ? 0 : 8); }
		bool c = cr > 0 && is_ccw_pos;
		if (dt > 0) { return (c) ? 1 : 7; }
		else if (dt == 0) { return (c) ? 2 : 6; }
		else { return (c) ? 3 : 5; }
	}
	bool compare(const Point& lhs, const Point& rhs, bool is_contain_eq = true) const {
		if (lhs == rhs) { return is_contain_eq; }
		int lz = getZone(lhs), rz = getZone(rhs);
		if (lz != rz) { return lz < rz; }
		T dl = lhs.dot(lhs), dr = rhs.dot(rhs);
		if (lz % 2 == 0) { return dl < dr; }
		T cr = lhs.cross(rhs);
		return (cr > 0 && is_ccw_pos) || (cr < 0 && !is_ccw_pos) || (cr == 0 && dl < dr);
		
	}
public:
	ArgumentSorting() {}
	ArgumentSorting(const Point& vec, const Point& org = Point(1, 0),
					bool is_start = true, bool is_ccw = true)
	: start_vec(vec), origin_arg(org), is_vec_start(is_start), is_ccw_pos(is_ccw) {}
	void sort(vector<Point>* points) const {
		function<bool(const Point&, const Point&)> cmp;
		cmp = [this](const Point& l, const Point& r) {
			return compare(l, r, false);
		};
		std::sort(points->begin(), points->end(), cmp);
	}
	static void sort(vector<Point>* points, const Point& vec,
					 const Point& org = Point(1, 0), bool is_start = true,
					 bool is_ccw = true) {
		ArgumentSorting(vec, org, is_start, is_ccw).sort(points);
	}
};

#endif /* ArgumentSorting_hpp */
