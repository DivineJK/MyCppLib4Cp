//
//  Dijkstra.hpp
//  MyCppLib4Cp
//
//  Created by DivineJK on 2024/01/22.
//

#ifndef Dijkstra_h
#define Dijkstra_h

class Dijkstra {
private:
	size_t start = 0;
	vector<size_t> prev;
	vector<uint64_t> cost;
	vector<unordered_map<size_t, uint64_t>> graph;
public:
	Dijkstra(size_t n) {
		graph.resize(n);
		prev.resize(n);
		cost.resize(n);
		for (size_t i = 0; i < n; i++) { cost[i] = UINT64_MAX; }
	}
	Dijkstra(const vector<unordered_map<size_t, uint64_t>>& g) {
		size_t n = g.size();
		prev.resize(n);
		cost.resize(n);
		for (size_t i = 0; i < n; i++) { cost[i] = UINT64_MAX; }
		graph = g;
	}
	void addEdge(size_t a, size_t b, uint64_t c) {
		graph[a][b] = (graph[a].count(b)) ? min(c, graph[a][b]) : c;
	}
	void calculateCost(size_t st) {
		for (size_t i = 0; i < cost.size(); i++) { cost[i] = UINT64_MAX; }
		start = st;
		cost[st] = 0;
		using pair_ci = pair<uint64_t, size_t>;
		priority_queue<pair_ci, vector<pair_ci>, greater<pair_ci>> pq;
		pq.push(make_pair(0, st));
		while (!pq.empty()) {
			pair_ci b = pq.top();
			pq.pop();
			size_t idx = b.second;
			uint64_t r = b.first;
			if (cost[idx] < r) { continue; }
			for (const pair<const size_t, uint64_t>& p : graph[idx]) {
				size_t i = p.first;
				uint64_t v = p.second;
				if (cost[i] > r + v) {
					cost[i] = r + v;
					pq.push(make_pair(cost[i], i));
					prev[i] = idx;
				}
			}
		}
	}
	uint64_t getCost(size_t g) const { return cost[g]; }
	vector<size_t> getShortestPath(size_t g) const {
		if (cost[g] == UINT64_MAX) { return {}; }
		vector<size_t> ret;
		size_t c = g;
		while (c != start) {
			ret.push_back(c);
			c = prev[c];
		}
		ret.push_back(start);
		reverse(ret.begin(), ret.end());
		return ret;
	}
};

#endif /* Dijkstra_h */
