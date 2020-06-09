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

/*
 * Computes the shortest distance from u to v in graph g.
 * The shortest path corresponds to a sequence of vertices starting from u and ends at v,
 * which has the smallest total weight of edges among all possible paths from u to v.
 */

template <typename T>
vector<vertex<T>> shortest_path(directed_graph<T> g, int u_id, int v_id) {

  map<T, T> dijkstras; //holds the distance values from source (Type key,Type value) : VertexIds, VertexDistanceFromSource
  unordered_set<T> spt_set; //spanning tree set, keeps tracks of vertices included in spt, i.e., whose minimum distance from source is calculated and finalized. Initially, this set is empty.
  vector<vertex<T>> path; //path of vertices in sequential shortest path visiting order from source(u_id) to target(v_id)
  int current = u_id; //index of which vertex currently on
  int numVertices = 0; //counter for number of reachable vertices
  map<T, T> prev; //holds the previous vertex from vertexid
  stack<int> reverse;
  int n = v_id;

  // Initialise all vertex distances as infinite
	for (vertex<T> v : g.get_vertices()) {
		dijkstras.insert(std::pair<int, int>(v.id, std::numeric_limits<int>::max()));	
    prev.insert(std::pair<int,int>(v.id, v.id));
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
       if(spt_set.count(n.id) == 0 && dijkstras[current]+g.get_edgeWeight(current, n.id) < dijkstras[n.id]){
         dijkstras[n.id] = dijkstras[current]+g.get_edgeWeight(current, n.id);
         prev[n.id] = current; //grab previous element of n.id
       }
    }
  }
  
  /*
  //TEST FUNCTION FOR DIJKSTRAS
  for(auto elem: dijkstras){
    cout << elem.first << " " << elem.second << endl;
  }
  cout << "----------------------" << endl;
  
  //TEST FUNCTION FOR PREV
  for(auto elem: prev){
    cout << elem.first << " " << elem.second << endl;
  }
  */

  //trace back path, let's say from vertex 4 (source would be 0)
  reverse.push(n);
  //cout << "path test" << endl;
  if(g.reachable(u_id,v_id)){
    while(n!= u_id){
      n = prev.at(n);
      reverse.push(n);
    }
  }

  //push vertex order path in vector
  while(!reverse.empty()){
    path.push_back(vertex<T>(reverse.top(),g.get_vertexWeight(reverse.top())));
    reverse.pop();
  }

  return path;
}

/*
 * Computes the strongly connected components of the graph.
 * A strongly connected component is a subset of the vertices
 * such that for every pair u, v of vertices in the subset,
 * v is reachable from u and u is reachable from v.
 */

template<typename T>
vector<vector<vertex<T>>> tarjan(int u, int dfn[], int low[], stack<int> *s, bool onStack[], directed_graph<T>& g) {

  static int dfn_cnt = 0;
  int w = 0;
  low[u] = dfn[u] = ++dfn_cnt;  //initialise low and order (dfn) value
  s->push(u); //add u to stack
  onStack[u] = true;
  set<set<int>> all_sccs; // algorithm’s result: all SCCs
  static vector<vector<vertex<T>>> result;
  vector<vertex<T>> av;

  for (vertex<T> v: g.get_neighbours(u)) { //get neighbours of u
    if (dfn[v.id] == -1){ // case1: v unvisited: : continue to search for v in depth. If not visited yet, recur for it
      tarjan(v.id, dfn, low, s, onStack, g); // recur for v
      low[u] = std::min(low[u], low[v.id]); // retrospective
    } 
    else if (onStack[v.id] == true) { // case2: v visited & in stack
      low[u] = std::min(low[u], dfn[v.id]); // update low[u]
    }
  }

  if (low[u]==dfn[u]){  // this only happens when u is the root of a SCC, indicating start of SCC, or if neighbour of current is equal to 
    set<int> scc; // create a new scc
    while (s->top() != u){ // add all vertices above u in stack to scc
      //for group sccs
      w = s->top();
      scc.insert(w);
      av.push_back(vertex<T>(w, g.get_vertexWeight(w)));
      onStack[w] = false;
      s->pop();
    }
    //for single sccs
    w = s->top();
    scc.insert(w);
    av.push_back(vertex<T>(w, g.get_vertexWeight(w)));
    onStack[w] = false;
    s->pop(); 
    all_sccs.insert(scc); // put scc to the result
    result.push_back(av);
  }

     /* set<set<int>>::iterator it_ex; // iterator for the "outer" structure
      set<int>::iterator it_in; // iterator for the "inner" structure
      for(it_ex = all_sccs.begin(); it_ex != all_sccs.end(); it_ex++) {
        for(it_in = it_ex->begin(); it_in != it_ex->end(); it_in++)   
            cout << *it_in << ", ";
        cout << endl;
    }*/
    
    return result;
}

template <typename T>
vector<vector<vertex<T>>> strongly_connected_components(directed_graph<T> g){

  vector<vector<vertex<T>>> result;
  int *dfn = new int[g.get_adjMatrixSize()]; 
  int *low = new int[g.get_adjMatrixSize()]; 
  bool *onStack = new bool[g.get_adjMatrixSize()];
  stack<int> *s = new stack<int>(); 

    for (int i = 0; i < g.get_adjMatrixSize(); i++) { 
        dfn[i] = -1;  
        low[i] = -1; 
        onStack[i] = false; 
    }

    for (vertex<T> a: g.get_vertices()){
      if (dfn[a.id] == -1){// u unvisited
      result = tarjan(a.id, dfn, low, s, onStack, g);
      }
    }

  // tarjan executes multiple times so result changes each time
   return result;

}

/*
 * Computes a topological ordering of the vertices.
 * For every vertex u in the order, and any of its
 * neighbours v, v appears later in the order than u.
 * You will be given a DAG as the argument.
 */

template <typename T>
void topologicalSort(directed_graph<T> g, int u, bool visited[], stack<int> &Stack){
  // Mark the current node as visited. 
    visited[u] = true; 
  
    // Recur for all the vertices adjacent to this vertex  
    for (auto i : g.get_neighbours(u)){
        if (!visited[i.id]){
            topologicalSort(g, i.id, visited, Stack); 
        }
    }
  
    // Push current vertex to stack which stores result 
    Stack.push(u); 
}

template <typename T>
vector<vertex<T>> topological_sort(directed_graph<T> g) { //Kahn's Algorithm or DFS

  stack<int> Stack;
  bool *visited = new bool[g.get_adjMatrixSize()];
  vector<vertex<T>> result;

  for (int i = 0; i < g.get_adjMatrixSize(); i++){
    visited[i] = false;
  }

  for (auto v : g.get_vertices()){
    if (visited[v.id] == false)
      topologicalSort(g, v.id, visited, Stack);
  }

  while (Stack.empty() == false) {
    result.push_back(vertex<T>(Stack.top(), g.get_vertexWeight(Stack.top())));
    Stack.pop(); } 

    return result;

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

  static vector<vector<vertex<T>>> result;
  stack<int> hold;
  T value;
  vector<int> cost;
  vector<int> population;
  int sumCost = 0;
  int sumPopulation = 0;

  for(auto a : g.get_vertices()){
    if(a.id !=u_id){
     result.push_back(shortest_path(g, u_id, a.id)); //return vector of shortest path from root to other nodes
      //if vertices return greater than 2, only return the last 2 
    }
  }

  for(auto vector: result){
    while(vector.size() > 2){ //if contains more than 2 vertices
      vector.erase(vector.begin());
    }
    for(auto v : vector){
        //cout << v.id << ",";
        hold.push(v.id);
    }
    //cout << " Edge Weight = ";
    value = hold.top();
    hold.pop();
    //cout << g.get_edgeWeight(hold.top(), value) << " , ";
    cost.push_back(g.get_edgeWeight(hold.top(), value));
    //cout << "Vertex Weight: " << g.get_vertexWeight(value);
    population.push_back(g.get_vertexWeight(value));
    //cout << endl;
  }

    for(int n : cost){
      sumCost += n;
    }

    for(int n : population){
      sumPopulation += n;
    }

  return sumCost/sumPopulation;

}
