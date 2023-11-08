//
//  DynamicSegmentTree.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2023/11/05.
//

#ifndef DynamicSegmentTree_hpp
#define DynamicSegmentTree_hpp

template <typename T>
class DynamicSegmentTree {
	struct Node {
		size_t index;
		T value;
		T product;
		unique_ptr<Node> left;
		unique_ptr<Node> right;
		Node(size_t i, const T& v) : index(i), value(v), product(v) {}
	};
private:
	function<T(const T&, const T&)> func;
	T unit;
	size_t n;
	unique_ptr<Node> root;
private:
	void set_internal(size_t p, size_t il, size_t ir, const T& v, unique_ptr<Node>* nd) {
		if (!*nd) {
			*nd = make_unique<Node>(p, v);
			return;
		}
		if (p == (*nd)->index) {
			(*nd)->value = v;
			update_node(nd->get());
			return;
		}
		size_t mid = il + ((ir - il) >> 1);
		size_t np = p;
		T nv = v;
		if (p < mid) {
			if ((*nd)->index < p) {
				np = (*nd)->index;
				nv = (*nd)->value;
				(*nd)->index = p;
				(*nd)->value = v;
			}
			set_internal(np, il, mid, nv, &(*nd)->left);
		} else {
			if (p < (*nd)->index) {
				np = (*nd)->index;
				nv = (*nd)->value;
				(*nd)->index = p;
				(*nd)->value = v;
			}
			set_internal(np, mid, ir, nv, &(*nd)->right);
		}
		update_node(nd->get());
	}
	T get_internal(size_t p, size_t il, size_t ir, const Node* nd) const {
		if (!nd || p < il || ir <= p) { return unit; }
		if (p == nd->index) { return nd->value; }
		size_t mid = il + ((ir - il) >> 1);
		if (p < mid) { return get_internal(p, il, mid, nd->left.get()); }
		else { return get_internal(p, mid, ir, nd->right.get()); }
	}
	T get_internal(size_t l, size_t r, size_t il, size_t ir, const Node* nd) const {
		if (!nd || l >= r || r <= il || ir <= l) { return unit; }
		if (l <= il && ir <= r) { return nd->product; }
		size_t mid = il + ((ir - il) >> 1);
		T ret = get_internal(l, r, il, mid, nd->left.get());
		if (l <= nd->index && nd->index < r) { ret = func(ret, nd->value); }
		return func(ret, get_internal(l, r, mid, ir, nd->right.get()));
	}
	void update_node(Node* nd) const {
		assert(nd);
		T lv = func((nd->left) ? nd->left->product : unit, nd->value);
		T rv = (nd->right) ? nd->right->product : unit;
		nd->product = func(lv, rv);
	}
public:
	DynamicSegmentTree(size_t aN, const T& aUnit,
					   const function<T(const T&, const T&)>& aFunc)
	: n(aN), func(aFunc), unit(aUnit) {}
	void setValue(size_t p, const T& v) { set_internal(p, 0, n, v, &root); }
	T getValue(size_t p) const { return get_internal(p, 0, n, root.get()); }
	T getSegment(size_t l, size_t r) const { return get_internal(l, r, 0, n, root.get()); }
};

#endif /* DynamicSegmentTree_hpp */
