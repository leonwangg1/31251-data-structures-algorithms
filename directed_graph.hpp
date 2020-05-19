 #ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <queue>
#include <stack>
#include <limits>
#include <utility> 
#include <algorithm>
#include <string>

using namespace std; 


template <typename T>
class vertex {

public:
	int id;
	T weight;

	vertex(int v_id, T v_weight) : id(v_id), weight(v_weight) {
	}
};


template <typename T>
class directed_graph {

private:
	vector<vector<T>> adj_matrix; // dj_matrix[u_id][v_id] = the weight for edge (u_id, u_id).
	vector<T> vertex_weights; // vertex_weights[u_id] stores the weight of vertex u_id.
	int edgeCount = 0;

public:

	directed_graph(); //A constructor for directed_graph. The graph should start empty.
	~directed_graph(); //A destructor. Depending on how you do things, this may not be necessary.

	void increaseCapacity();

	bool contains(const int&) const; //Returns true if the graph contains the given vertex_id, false otherwise.
	bool adjacent(const int&, const int&) const; //Returns true if the first vertex is adjacent to the second, false otherwise.

	void add_vertex(const vertex<T>&); //Adds the passed in vertex to the graph (with no edges).
	void add_edge(const int&, const int&, const T&); //Adds a weighted edge from the first vertex to the second.

	void remove_vertex(const int&); //Removes the given vertex. Should also clear any incident edges.
	void remove_edge(const int&, const int&); //Removes the edge between the two vertices, if it exists.

	size_t in_degree(const int&) const; //Returns number of edges coming in to a vertex.
	size_t out_degree(const int&) const; //Returns the number of edges leaving a vertex.
	size_t degree(const int&) const; //Returns the degree of the vertex (both in edges and out edges).

	size_t num_vertices() const; //Returns the total number of vertices in the graph.
	size_t num_edges() const; //Returns the total number of edges in the graph.

	vector<vertex<T>> get_vertices(); //Returns a vector containing all the vertices.
	vector<vertex<T>> get_neighbours(const int&); //Returns a vector containing all the vertices reachable from the given vertex. The vertex is not considered a neighbour of itself.
	vector<vertex<T>> get_second_order_neighbours(const int&); // Returns a vector containing all the second_order_neighbours (i.e., neighbours of neighbours) of the given vertex.
															  // A vector cannot be considered a second_order_neighbour of itself.
	bool reachable(const int&, const int&) const; //Returns true if the second vertex is reachable from the first (can you follow a path of out-edges to get from the first to the second?). Returns false otherwise.
	bool contain_cycles() const; // Return true if the graph contains cycles (there is a path from any vertices directly/indirectly to itself), false otherwise.

	vector<vertex<T>> depth_first(const int&); //Returns the vertices of the graph in the order they are visited in by a depth-first traversal starting at the given vertex.
	vector<vertex<T>> breadth_first(const int&); //Returns the vertices of the graph in the order they are visisted in by a breadth-first traversal starting at the given vertex.

	directed_graph<T> out_tree(const int&); //Returns a spanning tree of the graph starting at the given vertex using the out-edges. This means every vertex in the tree is reachable from the root.

	vector<vertex<T>> pre_order_traversal(const int&, directed_graph<T>&); // returns the vertices in the visiting order of a pre-order traversal of the minimum spanning tree starting at the given vertex.
	vector<vertex<T>> in_order_traversal(const int&, directed_graph<T>&); // returns the vertices in the visiting order of an in-order traversal of the minimum spanning tree starting at the given vertex.
	vector<vertex<T>> post_order_traversal(const int&, directed_graph<T>&); // returns the vertices in ther visitig order of a post-order traversal of the minimum spanning tree starting at the given vertex.

	vector<vertex<T>> significance_sorting(); // Return a vector containing a sorted list of the vertices in descending order of their significance.
};

template<typename T>
directed_graph<T>::directed_graph(){

	int initial_capacity = 100;

	// initialise a matrix of [initial_capacity*initial_capacity], which stores vertex ids from 0, 1, ..., initial_capacity-1.
	adj_matrix.resize(initial_capacity);
	for(int i=0; i<adj_matrix.size(); i++){
		adj_matrix[i].resize(initial_capacity);
		for(int j=0; j<adj_matrix[i].size(); j++){
			adj_matrix[i][j] = 0; // 0 indicates there is no edge from i to j in the graph
		}
	}

	// initialise the matrix to contain no vertex
	vertex_weights.resize(initial_capacity);
	for(int i=0; i<vertex_weights.size(); i++){
		vertex_weights[i] = 0;  // 0 indicate the vertex with the id of i is not in the graph 
	}
}

template <typename T>
directed_graph<T>::~directed_graph() {}

template <typename T>
void directed_graph<T>::increaseCapacity(){

	int old_capacity = vertex_weights.size();
	int new_capacity = 2 * old_capacity;

	// Step 1a: expand adj_matrix from [old_capacity*old_capacity] to [new_capacity*new_capacity]
	adj_matrix.resize(new_capacity);
	for(int i=0; i<adj_matrix.size(); i++){
		adj_matrix[i].resize(new_capacity);
	}

	// Step 1b: initialise the new space to contain no edge
	for(int i=0; i<old_capacity; i++){
		for(int j=old_capacity; j<new_capacity; j++){
			adj_matrix[i][j] = 0; // 0 indicates there is no edge from i to j in the graph 
		}
	}
	for(int i=old_capacity; i<new_capacity; i++){
		for(int j=0; j<new_capacity; j++){
			adj_matrix[i][j] = 0; /// 0 indicates there is no edge from i to j in the graph
		}
	}

	// Step 2a: expand the size of vertex_weights from capacity to new_capacity
	vertex_weights.resize(new_capacity);
	// Step 2b: initialise the new space to contain no vertex
	for(int i=old_capacity; i<new_capacity; i++){
		vertex_weights[i] = 0; // 0 indicate the vertex with the id of i is not in the graph 
	}
}

template <typename T>
bool directed_graph<T>::contains(const int& u_id) const {
	if (vertex_weights[u_id] > 0){ // 0 means the vertex is not in the graph
		return true;
	}
	return false;
}

template <typename T>
bool directed_graph<T>::adjacent(const int& u_id, const int& v_id) const { 
	//if vertices are valid, return whether or not u_id contains an edge to v_id, else return false.
	return (contains(u_id) && contains(v_id)) ? adj_matrix[u_id][v_id] > 0
		: false;
}

template <typename T>
void directed_graph<T>::add_vertex(const vertex<T>& u) {
	while(u.id > vertex_weights.size()-1){
		increaseCapacity();
	}
	vertex_weights[u.id] = u.weight;
}

template <typename T>
void directed_graph<T>::add_edge(const int& u_id, const int& v_id, const T& edge_weight) {
	if(contains(u_id) && contains(v_id)){ // check if the vertices are in the graph
		adj_matrix[u_id][v_id] = edge_weight; // this demo requires edge_weight != 0
		edgeCount += 1;
	}
}

template <typename T>
void directed_graph<T>::remove_vertex(const int& u_id) {
	vertex_weights[u_id] = 0;
}

template <typename T>
void directed_graph<T>::remove_edge(const int& u_id, const int& v_id) {
	if(contains(u_id) && contains(v_id)){ // check if the vertices are in the graph
		adj_matrix[u_id][v_id] = 0; // this demo requires edge_weight != 0
		edgeCount -= 1;
	}
}

template <typename T>
size_t directed_graph<T>::in_degree(const int& u_id) const { 
	int total_in=0;
	if(contains(u_id)){
		for(int i=0;i<adj_matrix[u_id].size();i++){ //how much columns does row u_id contain?
			//cout << adj_matrix[u_id][i] << endl;
			if(adj_matrix[i][u_id] > 0){ //adj_matrix at row u_id and i column
				total_in = total_in + 1;
			}
		}
	return total_in;
	}
	return 0;
 }

template <typename T>
size_t directed_graph<T>::out_degree(const int& u_id) const { //referencing get neighbors concept
	int total_out=0;
	if(contains(u_id)){
		for(int i=0;i<adj_matrix[u_id].size();i++){ //how much columns does row u_id contain?
			//cout << adj_matrix[u_id][i] << endl;
			if(adj_matrix[u_id][i] > 0){ //adj_matrix at row u_id and i column
				total_out = total_out + 1;
			}
		}
	return total_out;
	}
	return 0;
 }

template <typename T>
size_t directed_graph<T>::degree(const int& u_id) const { 
	return in_degree(u_id)+out_degree(u_id); 
}

template <typename T>
size_t directed_graph<T>::num_vertices() const { 
	int total = 0;
	for(int i=0; i < vertex_weights.size(); i++){
		if (vertex_weights[i] > 0){
		total = total + 1;
		}
	}
	return total;	
} 

template <typename T>
size_t directed_graph<T>::num_edges() const { 
	return edgeCount;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::get_vertices() {

	vector<vertex<T>> vertice_list;

	for(int i=0; i<vertex_weights.size(); i++){
		if (vertex_weights[i] > 0){
			vertice_list.push_back(vertex<T>(i, vertex_weights[i])); // construct vertex<T> from vertex_id
		}
	}

	return vertice_list;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::get_neighbours(const int& u_id) {
	vector<vertex<T>> result;
	vector<vertex<T>> dup;
	if(contains(u_id)){ // first make sure the vertex is in the graph
		for (int i=0; i<adj_matrix[u_id].size(); i++){
			if(adj_matrix[u_id][i] > 0){ // check if there is an edge in that row every column
					dup.push_back(vertex<T>(i, vertex_weights[i]));
			}
		}
		//check for duplicates
	  int pos = 0;
		for (vertex<T> a : dup){
			if(a.id == u_id){
				dup.erase (dup.begin()+pos);
			}
			pos +=1;
			result = dup;
		} 
	}
	result = dup;
	return result;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::get_second_order_neighbours(const int& u_id) { 

	vector<vertex<T>> first_orders = get_neighbours(u_id);
	vector<vertex<T>> second_orders;

	for (auto x: first_orders) {
		for (int i=0; i<adj_matrix[x.id].size();i++){
			if(adj_matrix[x.id][i] > 0){
				if(i != u_id){
					bool is_duplicate = false;
					for (auto v: second_orders){
						if (v.id == i){
							is_duplicate = true;
							break;
						}
					}
					if(!is_duplicate){
					second_orders.push_back(vertex<T>(i, vertex_weights[i]));
					}
				}
			}
		}
	}
	return second_orders;
}

template <typename T>
bool directed_graph<T>::reachable(const int& u_id, const int& v_id) const { return false; }

template <typename T>
bool directed_graph<T>::contain_cycles() const { return false; }

template <typename T>
vector<vertex<T>> directed_graph<T>::depth_first(const int& u_id) 
{ 	
	stack<T>unprocessed;
	vector<vertex<T>> visited_vertex;
	bool visited[adj_matrix.size()];

    for (unsigned i = 0; i<adj_matrix.size(); i++) // Intialise all vertices in the matrix to not visited or un visited.
    	{
        visited[i] = false;
    	} 

	unprocessed.push(u_id);

    while (!unprocessed.empty())
    {
        int top_stack = unprocessed.top();
        unprocessed.pop();
    
        if (!visited[top_stack])
        {
            visited[top_stack]=true;
            visited_vertex.push_back(vertex<T>(top_stack, vertex_weights[top_stack]));
			//check for exsiting neighbours
            for (unsigned i =num_vertices() ; i!= 0; i--)
            {
                if(adj_matrix[top_stack][i-1] != 0)
                {
                    unprocessed.push(i-1);
                }
            }
        }
    }
	return visited_vertex;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::breadth_first(const int& u_id) { return vector<vertex<T>>(); }

template <typename T>
directed_graph<T> directed_graph<T>::out_tree(const int& u_id) { return directed_graph<T>(); }

template <typename T>
vector<vertex<T>> directed_graph<T>::pre_order_traversal(const int& u_id, directed_graph<T>& mst) { return vector<vertex<T>>(); }

template <typename T>
vector<vertex<T>> directed_graph<T>::in_order_traversal(const int& u_id, directed_graph<T>& mst) { return vector<vertex<T>>(); }

template <typename T>
vector<vertex<T>> directed_graph<T>::post_order_traversal(const int& u_id, directed_graph<T>& mst) { return vector<vertex<T>>(); }

template <typename T>
vector<vertex<T>> directed_graph<T>::significance_sorting() { return vector<vertex<T>>(); }

#endif