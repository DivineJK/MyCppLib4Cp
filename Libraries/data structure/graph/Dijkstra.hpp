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
	mutable bool need_calc = true;
	mutable int start = 0;
	mutable vector<int> prev;
	mutable vector<uint64_t> cost;
	vector<unordered_set<int, uint64_t>> graph;
public:
	Dijkstra(int n) {
		graph.resize(n);
		prev.resize(n);
		cost.resize(n);
	}
	Dijkstra(int n, vector<unordered_set<int, uint64_t>> g) {
		graph = move(g);
	}
	void addEdge(int a, int b, uint64_t c) {
		need_calc = true;
		graph[a].push_back(make_pair(b, c));
	}
	void calculateCost(int s) const {
		need_calc = false;
		prev_start = s;
	}
	uint64_t getCost(int s, int g) const {
		if (need_calc) { calculateCost(s); }
		return cost[g];
	}
};

#endif /* Dijkstra_h */
