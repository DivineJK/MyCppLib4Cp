//
//  LiChaoTree.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/09/10.
//

#ifndef LiChaoTree_hpp
#define LiChaoTree_hpp

template <typename T, bool MIN_MODE, T X_MIN_INF, T X_MAX_INF, T VAL_INF>
class LiChaoTree {
	static_assert(X_MIN_INF < X_MAX_INF, "X_MIN_INF must be less than X_MAX_INF");
	struct Line {
		T a;
		T b;
		Line() {}
		Line(T aA, T aB) : a(aA), b(aB) {}
		inline T evaluate(T x) const { return a * x + b; }
	};
	struct Node {
		Line line;
		unique_ptr<Node> left;
		unique_ptr<Node> right;
		Node() : line(0, (MIN_MODE) ? VAL_INF : -VAL_INF), left(nullptr), right(nullptr) {}
		Node(const Line& l) : line(l), left(nullptr), right(nullptr) {}
	};
private:
	unique_ptr<Node> root;
private:
	void addLine_internal(const T& il, const T& ir, Line* l, unique_ptr<Node>* nd) {
		if (!*nd) {
			*nd = make_unique<Node>(*l);
			return;
		}
		T vl = l->evaluate(il), vr = l->evaluate(ir);
		T nl = (*nd)->line.evaluate(il), nr = (*nd)->line.evaluate(ir);
		if (nl <= vl && nr <= vr) { return; }
		else if (nl >= vl && nr >= vr) {
			(*nd)->line = *l;
			return;
		} else {
			T m = il + (ir - il) / 2;
			T nm = (*nd)->line.evaluate(m), vm = l->evaluate(m);
			if (m == il) {
				if (vm < nm || (vm == nm && vr < nr)) { swap((*nd)->line, *l); }
				return;
			}
			if (vm < nm) {
				swap((*nd)->line, *l);
				if (vl >= nl) { addLine_internal(il, m, l, &(*nd)->left); }
				else { addLine_internal(m, ir, l, &(*nd)->right); }
			} else {
				if (vl <= nl) { addLine_internal(il, m, l, &(*nd)->left); }
				else { addLine_internal(m, ir, l, &(*nd)->right); }
			}
		}
	}
	void addSegment_internal(const T& l, const T& r, const T& il, const T& ir,
							 Line* s, unique_ptr<Node>* nd) {
		if (r <= il || ir <= l) { return; }
		if (l <= il && ir <= r) {
			Line ns(*s);
			addLine_internal(il, ir, &ns, nd);
			return;
		}
		if (!*nd) { *nd = make_unique<Node>(); }
		T m = il + (ir - il) / 2;
		addSegment_internal(l, r, il, m, s, &(*nd)->left);
		addSegment_internal(l, r, m, ir, s, &(*nd)->right);
	}
	T getValue_internal(const T& x, const T& il, const T& ir, const Node* nd) const {
		if (!nd || il >= ir || x < il || x >= ir) { return VAL_INF; }
		if (ir - il == 1) { return nd->line.evaluate(x); }
		T m = il + (ir - il) / 2;
		if (x < m) {
			return min(nd->line.evaluate(x), getValue_internal(x, il, m, nd->left.get()));
		} else {
			return min(nd->line.evaluate(x), getValue_internal(x, m, ir, nd->right.get()));
		}
	}
public:
	LiChaoTree() {}
	void addLine(T a, T b) {
		Line l = (MIN_MODE) ? Line(a, b) : Line(-a, -b);
		addLine_internal(X_MIN_INF, X_MAX_INF, &l, &root);
	}
	void addSegment(T a, T b, T l, T r) {
		Line s = (MIN_MODE) ? Line(a, b) : Line(-a, -b);
		addSegment_internal(l, r, X_MIN_INF, X_MAX_INF, &s, &root);
	}
	T getValue(const T& x) const {
		T val = getValue_internal(x, X_MIN_INF, X_MAX_INF, root.get());
		return (MIN_MODE) ? val : -val;
	}
};

#endif /* LiChaoTree_hpp */
