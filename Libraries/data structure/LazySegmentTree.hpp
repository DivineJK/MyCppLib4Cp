//
//  LazySegmentTree.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef LazySegmentTree_hpp
#define LazySegmentTree_hpp

template <typename MAIN_TYPE, typename LAZY_TYPE>
class LazySegmentTree {
	using MAIN_FUNC_TYPE = function<MAIN_TYPE(const MAIN_TYPE&, const MAIN_TYPE&)>;
	using UPDATE_FUNC_TYPE = function<MAIN_TYPE(uint32_t, uint32_t,
												const LAZY_TYPE&, const MAIN_TYPE&)>;
	using LAZY_FUNC_TYPE = function<LAZY_TYPE(const LAZY_TYPE&, const LAZY_TYPE&)>;
private:
	MAIN_FUNC_TYPE main_func;
	UPDATE_FUNC_TYPE update_func;
	LAZY_FUNC_TYPE lazy_func;
	MAIN_TYPE main_unit;
	LAZY_TYPE lazy_unit;
	mutable vector<MAIN_TYPE> main_tree;
	mutable vector<LAZY_TYPE> lazy_tree;
	static constexpr uint32_t getBinMin(uint32_t a) {
		uint32_t ret = 1;
		while (ret < a) { ret <<= 1; }
		return ret;
	}
	void setArray(const vector<MAIN_TYPE>& arr) {
		uint32_t n = (uint32_t)arr.size();
		uint32_t b = getBinMin(n);
		main_tree = vector<MAIN_TYPE>(b << 1, main_unit);
		for (uint32_t i = 0; i < n; i++) { main_tree[i + b] = arr[i]; }
		for (uint32_t i = b; i > 1; i--) {
			main_tree[i - 1] = main_func(main_tree[(i - 1) << 1],
										 main_tree[((i - 1) << 1) | 1]);
		}
	}
	void extendBottom(uint32_t s) {
		uint32_t b = bottom_size();
		if (s <= b) { return; }
		getSegment(0, b);
		uint32_t nb = getBinMin(s);
		main_tree.resize(nb << 1);
		for (uint32_t i = 0; i < b; i++) {
			main_tree[i + nb] = update_func(i, i + 1, lazy_tree[i + b], main_tree[i + b]);
			lazy_tree[i + nb] = lazy_unit;
		}
		for (uint32_t i = b; i < nb; i++) {
			main_tree[i + nb] = main_unit;
			lazy_tree[i + nb] = lazy_unit;
		}
		for (uint32_t i = nb; i > 0; i--) {
			main_tree[i - 1] = main_func(main_tree[(i - 1) << 1],
										 main_tree[((i - 1) << 1) | 1]);
		}
	}
	void get_update_interval(uint32_t l, uint32_t r,
							 vector<uint32_t>* upd, vector<uint32_t>* asc) const {
		uint32_t n = bottom_size();
		uint32_t regl = (l >= n) ? n : l;
		uint32_t regr = (r >= n) ? n : r;
		uint32_t lid = n + regl;
		uint32_t rid = n + regr;
		vector<uint32_t> ls, rs;
		while (lid < rid) {
			if (lid & 1) {
				ls.push_back(lid);
				lid++;
			}
			if (rid & 1) {
				rs.push_back(rid - 1);
				rid--;
			}
			lid >>= 1;
			rid >>= 1;
		}
		upd->clear();
		upd->reserve(ls.size() + rs.size());
		for (vector<uint32_t>::const_iterator it = ls.cbegin(); it != ls.cend(); ++it) {
			upd->push_back(*it);
		}
		for (vector<uint32_t>::const_reverse_iterator it = rs.crbegin(); it != rs.crend(); ++it) {
			upd->push_back(*it);
		}
		asc->clear();
		lid = regl + n;
		rid = regr + n;
		bool lf = false, rf = false;
		while (lid && rid) {
			if (rf) { asc->push_back(rid); }
			if (lf && lid != rid) { asc->push_back(lid); }
			lf |= lid & 1;
			rf |= rid & 1;
			lid >>= 1;
			rid >>= 1;
		}
		reverse(asc->begin(), asc->end());
	}
	void propagate(uint32_t l, uint32_t r,
				   vector<uint32_t>* upd, vector<uint32_t>* asc) const {
		get_update_interval(l, r, upd, asc);
		for (vector<uint32_t>::const_iterator it = asc->cbegin(); it != asc->cend(); ++it) {
			uint32_t li, ri;
			uint32_t idx = *it;
			get_interval(*it, &li, &ri);
			main_tree[idx] = update_func(li, ri, lazy_tree[idx], main_tree[idx]);
			lazy_tree[idx << 1] = lazy_func(lazy_tree[idx << 1], lazy_tree[idx]);
			lazy_tree[(idx << 1) | 1] = lazy_func(lazy_tree[(idx << 1) | 1], lazy_tree[idx]);
			lazy_tree[idx] = lazy_unit;
		}
	}
	void get_interval(uint32_t idx, uint32_t* l, uint32_t* r) const {
		uint32_t n = bottom_size();
		*l = idx;
		*r = idx + 1;
		while (*l < n) {
			*l <<= 1;
			*r <<= 1;
		}
	}
public:
	LazySegmentTree(uint32_t n, const MAIN_TYPE& mu, const MAIN_FUNC_TYPE& mf,
					const UPDATE_FUNC_TYPE& uf,
					const LAZY_TYPE& lu, const LAZY_FUNC_TYPE& lf) {
		main_func = mf;
		update_func = uf;
		lazy_func = lf;
		main_unit = mu;
		lazy_unit = lu;
		main_tree = vector<MAIN_TYPE>(getBinMin(n) << 1, mu);
		lazy_tree = vector<LAZY_TYPE>(getBinMin(n) << 1, lu);
		
	}
	LazySegmentTree(const vector<MAIN_TYPE>& arr, const MAIN_TYPE& mu, const MAIN_FUNC_TYPE& mf,
					const UPDATE_FUNC_TYPE& uf,
					const LAZY_TYPE& lu, const LAZY_FUNC_TYPE& lf) {
		main_func = mf;
		update_func = uf;
		lazy_func = lf;
		main_unit = mu;
		lazy_unit = lu;
		lazy_tree = vector<LAZY_TYPE>(getBinMin((uint32_t)arr.size()) << 1, lu);
		setArray(arr);
	}
	LazySegmentTree(const LazySegmentTree& other) { operator=(other); }
	LazySegmentTree& operator=(const LazySegmentTree& other) {
		main_func = other.main_func;
		update_func = other.update_func;
		lazy_func = other.lazy_func;
		main_unit = other.main_unit;
		lazy_unit = other.lazy_unit;
		main_tree = other.main_tree;
		lazy_tree = other.lazy_tree;
	}
	uint32_t size() const { return (uint32_t)main_tree.size(); }
	uint32_t bottom_size() const { return (uint32_t)main_tree.size() >> 1; }
	void assign(uint32_t idx, const MAIN_TYPE& val) {
		uint32_t n = bottom_size();
		if (idx >= n) { extendBottom(idx + 1); }
		n = bottom_size();
		uint32_t l = 0, r = n;
		uint32_t len = n;
		uint32_t c = 1;
		while (len > 1) {
			len >>= 1;
			main_tree[c] = update_func(l, r, lazy_func[c], main_func[c]);
			lazy_func[c << 1] = lazy_func(lazy_tree[c << 1], lazy_tree[c]);
			lazy_func[(c << 1) | 1] = lazy_func(lazy_tree[(c << 1) | 1], lazy_tree[c]);
			lazy_tree[c] = lazy_unit;
			c <<= 1;
			if (idx < l + len) {
				r -= len;
			} else {
				l += len;
				c |= 1;
			}
		}
		c = n + idx;
		main_tree[c] = val;
		lazy_tree[c] = lazy_unit;
		l = idx;
		r = idx;
		len = 1;
		c >>= 1;
		while (c) {
			if (c & 1) { l -= len; } else { r += len; }
			MAIN_TYPE lv, rv;
			lv = update_func(l, l + len, lazy_tree[c << 1], main_tree[c << 1]);
			rv = update_func(l + len, r, lazy_tree[(c << 1) | 1], main_tree[(c << 1) | 1]);
			main_tree[c] = main_func(lv, rv);
			lazy_tree[c] = lazy_unit;
			c >>= 1;
		}
	}
	void update(uint32_t l, uint32_t r, const LAZY_TYPE& lv) {
		if (l >= r) { return; }
		uint32_t n = bottom_size();
		if (l >= n) { return; }
		vector<uint32_t> upd, asc;
		propagate(l, r, &upd, &asc);
		for (vector<uint32_t>::const_iterator it = upd.cbegin(); it != upd.cend(); ++it) {
			lazy_tree[*it] = lazy_func(lazy_tree[*it], lv);
		}
		for (vector<uint32_t>::const_reverse_iterator it = asc.crbegin(); it != asc.crend(); ++it) {
			uint32_t li, mi, ri;
			get_interval((*it) << 1, &li, &mi);
			ri = (mi << 1) - li;
			MAIN_TYPE lv, rv;
			lv = update_func(li, mi, lazy_tree[(*it) << 1], main_tree[(*it) << 1]);
			rv = update_func(mi, ri, lazy_tree[((*it) << 1) | 1], main_tree[((*it) << 1) | 1]);
			main_tree[*it] = main_func(lv, rv);
		}
	}
	MAIN_TYPE getSegment(uint32_t l, uint32_t r) const {
		if (l >= r) { return main_unit; }
		uint32_t n = bottom_size();
		if (l >= n) { return main_unit; }
		vector<uint32_t> upd, asc;
		propagate(l, r, &upd, &asc);
		MAIN_TYPE ret = main_unit;
		for (vector<uint32_t>::const_iterator it = upd.cbegin(); it != upd.cend(); ++it) {
			uint32_t rl, rr;
			get_interval(*it, &rl, &rr);
			MAIN_TYPE val = update_func(rl, rr, lazy_tree[*it], main_tree[*it]);
			ret = main_func(ret, val);
		}
		return ret;
	}
};

#endif /* LazySegmentTree_hpp */
