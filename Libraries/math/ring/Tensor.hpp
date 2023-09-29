//
//  Tensor.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/28.
//

#ifndef Tensor_hpp
#define Tensor_hpp

template <typename T>
class Tensor {
private:
	static constexpr uint64_t TOTAL_SIZE_LIMIT = UINT32_MAX;
private:
	vector<uint32_t> elms_size;
	vector<bool> isUpper;
	vector<T> elms;
public:
	Tensor() {}
	Tensor(const vector<int>& sizeArr) {
		int n = (int)sizeArr.size();
		uint64_t totSize = 1;
		for (int i = 0; i < n; i++) { totSize *= sizeArr[i]; }
		assert(totSize <= TOTAL_SIZE_LIMIT);
		elms_size = sizeArr;
		elms.resize(totSize);
		isUpper.resize(n);
		for (int i = 0; i < n; i++) { isUpper[i] = true; }
	}
	Tensor(const vector<int>& sizeArr, const vector<bool>& aIsUpper) : Tensor(sizeArr) {
		isUpper = aIsUpper;
	}
	Tensor(const Tensor& other) { operator=(other); }
	Tensor& operator=(const Tensor& other) {
		elms_size = other.elms_size;
		isUpper = other.isUpper;
		elms = other.elms;
		return *this;
	}
	Tensor& operator=(const vector<T>& aArr) {
		elms_size = {aArr.size()};
		isUpper = {true};
		elms = aArr;
		return *this;
	}
	bool operator==(const Tensor& other) const {
		return elms_size == other.elms_size && isUpper == other.isUpper && elms == other.elms;
	}
	bool operator!=(const Tensor& other) const { return !(*this == other); }
	Tensor& operator+=(const Tensor& other);
	Tensor& operator-=(const Tensor& other);
	Tensor& operator*=(const Tensor& other);
	Tensor& operator/=(const T& other);
	uint32_t getRank() const { return elms_size.size(); }
	T getElement(const vector<uint32_t>& idxs) const {
		uint32_t n = getRank();
		assert(idxs.size() == n);
		uint32_t place = 0;
		for (int i = 0; i < n; i++) {
			assert(idxs[i] < elms_size[i]);
			place *= elms_size[i];
			place += idxs[i];
		}
		return elms[place];
	}
	void setElement(const vector<uint32_t>& idxs, const T& val) {
		assert(idxs.size() == elms_size.size());
		uint32_t n = getRank();
		assert(idxs.size() == n);
		uint32_t place = 0;
		for (int i = 0; i < n; i++) {
			assert(idxs[i] < elms_size[i]);
			place *= elms_size[i];
			place += idxs[i];
		}
		elms[place] = val;
	}
	Tensor& contract(const vector<uint32_t>& u, const vector<uint32_t>& l);
	Tensor getContracted(const vector<uint32_t>& u, const vector<uint32_t>& l);
};

#endif /* Tensor_hpp */
