//
//  LazySegmentTree.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef LazySegmentTree_hpp
#define LazySegmentTree_hpp

template <typename MAIN_TYPE, typename LAZY_TYPE,
const MAIN_TYPE& main_unit, const LAZY_TYPE& lazy_unit,
const function<MAIN_TYPE(const MAIN_TYPE&, const MAIN_TYPE&)>& main_func,
const function<MAIN_TYPE(uint32_t, uint32_t,
						 const LAZY_TYPE&, const MAIN_TYPE&)>& map_func,
const function<LAZY_TYPE(const LAZY_TYPE&, const LAZY_TYPE&)>& lazy_func>
class LazySegmentTree {
	struct MainStruct {
		uint32_t l = 0;
		uint32_t r = 0;
		MAIN_TYPE val = main_unit;
		MainStruct() {}
		MainStruct(uint32_t aL, uint32_t aR, const MAIN_TYPE& aVal) {
			l = aL;
			r = aR;
			val = aVal;
		}
	};
private:
	mutable vector<MainStruct> main_tree;
	mutable vector<LAZY_TYPE> lazy_tree;
	static constexpr uint32_t getBinMin(uint32_t a) {
		uint32_t ret = 1;
		while (ret < a) { ret <<= 1; }
		return ret;
	}
	void set_default(uint32_t s) {
		uint32_t b = getBinMin(s);
		main_tree = vector<MainStruct>(b << 1);
		for (uint32_t i = 0; i < b; i++) { main_tree[i + b] = MainStruct(i, i + 1, main_unit); }
		for (uint32_t i = b; i > 1; i--) {
			uint32_t l = main_tree[(i - 1) << 1].l;
			uint32_t r = main_tree[((i - 1) << 1) | 1].r;
			main_tree[i - 1] = MainStruct(l, r, main_unit);
		}
	}
	void setArray(const vector<MAIN_TYPE>& arr) {
		uint32_t n = (uint32_t)arr.size();
		uint32_t b = getBinMin(n);
		main_tree = vector<MainStruct>(b << 1);
		for (uint32_t i = 0; i < n; i++) { main_tree[i + b] = MainStruct(i, i + 1, arr[i]); }
		for (uint32_t i = n; i < b; i++) { main_tree[i + b] = MainStruct(i, i + 1, main_unit); }
		for (uint32_t i = b; i > 1; i--) {
			main_tree[i - 1].l = main_tree[(i - 1) << 1].l;
			main_tree[i - 1].r = main_tree[((i - 1) << 1) | 1].r;
			main_tree[i - 1].val = main_func(main_tree[(i - 1) << 1].val,
											 main_tree[((i - 1) << 1) | 1].val);
		}
	}
	void extendBottom(uint32_t s) {
		uint32_t b = bottom_size();
		if (s <= b) { return; }
		getSegment(0, b);
		uint32_t nb = getBinMin(s);
		main_tree.resize(nb << 1);
		for (uint32_t i = 0; i < b; i++) {
			main_tree[i + nb].l = i;
			main_tree[i + nb].r = i + 1;
			main_tree[i + nb].val = map_func(i, i + 1, lazy_tree[i + b], main_tree[i + b]);
			lazy_tree[i + nb] = lazy_unit;
		}
		for (uint32_t i = b; i < nb; i++) {
			main_tree[i + nb] = MainStruct(i, i + 1, main_unit);
			lazy_tree[i + nb] = lazy_unit;
		}
		for (uint32_t i = nb; i > 0; i--) {
			main_tree[i - 1].l = main_tree[(i - 1) << 1].l;
			main_tree[i - 1].r = main_tree[((i - 1) << 1) | 1].r;
			main_tree[i - 1].val = main_func(main_tree[(i - 1) << 1].val,
											 main_tree[((i - 1) << 1) | 1].val);
		}
	}
	void propagate(uint32_t l, uint32_t r) const {
		uint32_t n = bottom_size();
		uint32_t lg = n;
		uint32_t idxl = 1, idxr = 1;
		uint32_t ll = 0, rr = n;
		bool ul = ll != l || ll + lg > r, ur = rr != r || l + lg > rr;
		while (ul || ur) {
			if (ul) { descend(idxl); }
			if ((idxl != idxr || !ul) && ur) { descend(idxr); }
			lg >>= 1;
			if (ul) {
				idxl <<= 1;
				if (l >= ll + lg) {
					idxl++;
					ll += lg;
				}
			}
			if (ur) {
				idxr <<= 1;
				if (rr >= r + lg) { rr -= lg; }
				else { idxr++; }
			}
			ul &= ll != l || ll + lg > r;
			ur &= rr != r || l + lg > rr;
		}
	}
	MAIN_TYPE apply(uint32_t idx) const {
		return map_func(main_tree[idx].l, main_tree[idx].r, lazy_tree[idx], main_tree[idx].val);
	}
	void fetch(uint32_t idx) const {
		MAIN_TYPE lv = apply(idx << 1);
		MAIN_TYPE rv = apply((idx << 1) | 1);
		main_tree[idx].val = main_func(lv, rv);
	}
	void descend(uint32_t idx) const {
		main_tree[idx].val = apply(idx);
		lazy_tree[idx << 1] = lazy_func(lazy_tree[idx << 1], lazy_tree[idx]);
		lazy_tree[(idx << 1) | 1] = lazy_func(lazy_tree[(idx << 1) | 1], lazy_tree[idx]);
		lazy_tree[idx] = lazy_unit;
	}
public:
	LazySegmentTree(uint32_t n) {
		lazy_tree = vector<LAZY_TYPE>(getBinMin(n) << 1, lazy_unit);
		set_default(n);
	}
	LazySegmentTree(const vector<MAIN_TYPE>& arr) {
		lazy_tree = vector<LAZY_TYPE>(getBinMin((uint32_t)arr.size()) << 1, lazy_unit);
		setArray(arr);
	}
	LazySegmentTree(const LazySegmentTree& other) { operator=(other); }
	LazySegmentTree& operator=(const LazySegmentTree& other) {
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
			descend(c);
			c <<= 1;
			if (idx < l + len) {
				r -= len;
			} else {
				l += len;
				c |= 1;
			}
		}
		c = n + idx;
		main_tree[c].val = val;
		lazy_tree[c] = lazy_unit;
		l = idx;
		r = idx;
		len = 1;
		c >>= 1;
		while (c) {
			if (c & 1) { l -= len; } else { r += len; }
			MAIN_TYPE lv, rv;
			lv = map_func(l, l + len, lazy_tree[c << 1], main_tree[c << 1].val);
			rv = map_func(l + len, r, lazy_tree[(c << 1) | 1], main_tree[(c << 1) | 1].val);
			main_tree[c].val = main_func(lv, rv);
			lazy_tree[c] = lazy_unit;
			c >>= 1;
		}
	}
	void update(uint32_t l, uint32_t r, const LAZY_TYPE& lv) {
		if (l >= r) { return; }
		uint32_t n = bottom_size();
		uint32_t regl = (l >= n) ? n : l;
		uint32_t regr = (r >= n) ? n : r;
		if (regl >= regr) { return; }
		propagate(regl, regr);
		uint32_t rl = regl + n, rr = regr + n;
		while (rl < rr) {
			if (rl & 1) {
				lazy_tree[rl] = lazy_func(lazy_tree[rl], lv);
				rl++;
			}
			if (rr & 1) {
				lazy_tree[rr - 1] = lazy_func(lazy_tree[rr - 1], lv);
				rr--;
			}
			rl >>= 1;
			rr >>= 1;
		}
		uint32_t lml = regl + n, rmr = regr + n;
		uint32_t lg = 2;
		while (!(lml & 1) && regl + lg < regr) {
			lml >>= 1;
			lg <<= 1;
		}
		lg = 2;
		while (!(rmr & 1) && regl + lg < regr) {
			rmr >>= 1;
			lg <<= 1;
		}
		rmr--;
		lml >>= 1;
		rmr >>= 1;
		while (lml) {
			fetch(lml);
			if (lml != rmr) {
				fetch(rmr);
			}
			if (lml > rmr) {
				lml >>= 1;
			} else if (lml <= (rmr >> 1)) {
				rmr >>= 1;
			} else {
				lml >>= 1;
				rmr >>= 1;
			}
		}
	}
	MAIN_TYPE getSegment(uint32_t l, uint32_t r) const {
		if (l >= r) { return main_unit; }
		uint32_t n = bottom_size();
		uint32_t regl = (l >= n) ? n : l;
		uint32_t regr = (r >= n) ? n : r;
		if (regl >= regr) { return main_unit; }
		propagate(regl, regr);
		MAIN_TYPE retl = main_unit, retr = main_unit;
		uint32_t rl = regl + n, rr = regr + n;
		while (rl < rr) {
			if (rl & 1) { retl = main_func(retl, apply(rl++)); }
			if (rr & 1) { retr = main_func(apply(--rr), retr); }
			rl >>= 1;
			rr >>= 1;
		}
		return main_func(retl, retr);
	}
};

#endif /* LazySegmentTree_hpp */
