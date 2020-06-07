#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>
#include <cstddef>
#include <string>
#include <utility>
#include <algorithm>
#include <limits>
#include <optional>
#include <exception>
#include <stdexcept>

#include "directed_graph.hpp"

using namespace std;

/*
 * Computes the shortest distance from u to v in graph g.
 * The shortest path corresponds to a sequence of vertices starting from u and ends at v,
 * which has the smallest total weight of edges among all possible paths from u to v.
 */

template<typename T>
int min_distance(directed_graph<T>& g, const std::map<T, T>& dijkstras, const std::unordered_set<T>& spt_set) {
	int min_vertex;
	// Set to max int value
	int min = std::numeric_limits<int>::max();
	// For each vertex
	for (vertex<T> g_it : g.get_vertices()) {
		// If it is not apart of the spt_set yet, and it is less than the current distance
		if (spt_set.count(g_it.id) == 0 && dijkstras.at(g_it.id) <= min) {
			// Set the vertex as the current minimum
			min = dijkstras.at(g_it.id); 
			min_vertex = g_it.id;
		}
	}
	return min_vertex;
}

template <typename T>
vector<vertex<T>> shortest_path(directed_graph<T> g, int u_id, int v_id) {

  map<T, T> dijkstras; //holds the distance values from source (Type key,Type value) : VertexIds, VertexDistanceFromSource
  unordered_set<T> spt_set; //spanning tree set, keeps tracks of vertices included in spt, i.e., whose minimum distance from source is calculated and finalized. Initially, this set is empty.
  vector<vertex<T>> result; //result of vertices in sequential shortest path visiting order from source(u_id) to target(v_id)
  int current = u_id; //index of which vertex currently on
  int numVertices = 0; //counter for number of reachable vertices

  // Initialise all vertex distances as infinite
	for (vertex<T> v : g.get_vertices()) {
		dijkstras.insert(std::pair<int, int>(v.id, std::numeric_limits<int>::max()));	
	}

  // Set source distance value to 0
  dijkstras[u_id] = 0;

  // Calculate how much vertices are reachable from u_id
  for(vertex<T> j : g.get_vertices()){
    if(g.reachable(u_id, j.id)){
      numVertices += 1;
    }
  }

  // Find shortest path for all vertices
  while(spt_set.size() != numVertices){ //while spt doesn't include all the vertices that are reachable
    // Pick the minimum distance vertex from the set of vertices not yet processed (minimum distance value vertex not in spt_set)
    int current = min_distance(g, dijkstras, spt_set);
    /// Mark the picked vertex as processed
    spt_set.insert(current);
    // update the distance value of the adjacent vertices of the picked vertex
    for(vertex<T> n : g.get_neighbours(current)){
       if(dijkstras[current]+g.get_edgeWeight(current, n.id) < dijkstras[n.id]){
         dijkstras[n.id] = dijkstras[current]+g.get_edgeWeight(current, n.id);
       }
    }
  }
  
  //test function
  for(auto elem: dijkstras){
    cout << elem.first << " " << elem.second << endl;
  }

  return result;
}

/*
 * Computes the strongly connected components of the graph.
 * A strongly connected component is a subset of the vertices
 * such that for every pair u, v of vertices in the subset,
 * v is reachable from u and u is reachable from v.
 */

template <typename T>
vector<vector<vertex<T>>> strongly_connected_components(directed_graph<T> g) {

  return vector<vector<vertex<T>>>();

}

/*
 * Computes a topological ordering of the vertices.
 * For every vertex u in the order, and any of its
 * neighbours v, v appears later in the order than u.
 * You will be given a DAG as the argument.
 */
template <typename T>
vector<vertex<T>> topological_sort(directed_graph<T> g) {

  return vector<vertex<T>>();

}

/*
 * Computes the lowest cost-per-person for delivery over the graph.
 * u is the source vertex, which send deliveries to all other vertices.
 * vertices denote cities; vertex weights denote cities' population;
 * edge weights denote the fixed delivery cost between cities, which is irrelevant to 
 * the amount of goods being delivered. 
 */
template <typename T>
T low_cost_delivery(directed_graph<T> g, int u_id) {

  return 0;

}

