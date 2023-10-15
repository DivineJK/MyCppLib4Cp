//
//  DirichletSeries.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/10/09.
//

#ifndef DirichletSeries_hpp
#define DirichletSeries_hpp

template <typename T>
class DirichletSeries {
private:
	vector<T> coeffs;
public:
	DirichletSeries() {}
	DirichletSeries(uint32_t n) { coeffs.resize(n); }
	DirichletSeries(const vector<T>& vec) { operator=(vec); }
	DirichletSeries(const DirichletSeries& other) { operator=(other); }
	DirichletSeries& operator=(const vector<T>& vec) { coeffs = vec; }
	DirichletSeries& operator=(const DirichletSeries& other) { coeffs = other.coeffs; }
	bool operator==(const DirichletSeries& other) { return coeffs == other.coeffs; }
	bool operator!=(const DirichletSeries& other) { return !(*this == other); }
	static DirichletSeries getZeta(uint32_t n) {
		DirichletSeries ret(n);
		for (uint32_t i = 0; i < n; i++) { ret[i] = (T)1; }
	}
	static DirichletSeries getMobius(uint32_t n) {
		DirichletSeries ret(n);
		return ret;
	}
};

#endif /* DirichletSeries_hpp */
