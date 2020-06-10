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

	int get_edgeWeight(const int&, const int&);
	int get_vertexWeight(const int&);
	int get_adjMatrixSize();

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
	bool reachable(const int&,const int&); //Returns true if the second vertex is reachable from the first (can you follow a path of out-edges to get from the first to the second?). Returns false otherwise.
	bool contain_cycles(); // Return true if the graph contains cycles (there is a path from any vertices directly/indirectly to itself), false otherwise.

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
	return (contains(u_id) && contains(v_id)) ? adj_matrix[u_id][v_id] > 0 : false;
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
	//check if vertex exist in graph
	if(contains(u_id)){
		for(int i=0;i<adj_matrix[u_id].size();i++){ //how much columns does row u_id contain?
			if(adj_matrix[i][u_id] > 0){ //adj_matrix at row u_id and i column
				total_in = total_in + 1;
			}
		}
	}
	return total_in;
 }

template <typename T>
size_t directed_graph<T>::out_degree(const int& u_id) const { //referencing get neighbors concept
	int total_out=0;
	if(contains(u_id)){
		for(int i=0;i<adj_matrix[u_id].size();i++){ //how much columns does row u_id contain?
			if(adj_matrix[u_id][i] > 0){ //adj_matrix at row u_id and i column
				total_out = total_out + 1;
			}
		}
	}
	return total_out;
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
bool directed_graph<T>::reachable(const int& u_id,const int& v_id) { 

	bool visited[num_vertices()]; //initialize a boolean array of num of vertices
	stack<T> unprocessed;
	vector<vertex<T>> ordered;

	//if index of start_vertex is valid
	if (contains(u_id && v_id)){
		//set all index values to represent that they have not been visited yet
		for(unsigned i=0; i < num_vertices(); i++){
			visited[i] = false;
		}
		//push the start vertex id 'u' to the unprocessed stack
		unprocessed.push(u_id);
		//while there is still values in the unprocssed stack
		while (!unprocessed.empty()){
			//get the index of the top vertex and pop it from stack
			int index = unprocessed.top();
			unprocessed.pop();
			//if it hasn't been visited yet
			if(!visited[index]){
				//set it to visited
				visited[index] = true;
				//add the vertex to the ordered list
				ordered.push_back(vertex<T>(index, vertex_weights[index]));
					//if the vertex contains a neighbour(s), push it's unvisited neighbours to the stack
					if (out_degree(index) > 0){
						//add it's unvisited neighbour(s) to the unprocessed stack
						vector<vertex<T>> neighbour_list = get_neighbours(index);
						for(vertex<T> nb : neighbour_list){
						unprocessed.push(nb.id);
						}
					}
				}
			if(unprocessed.empty()){ //modified where all dfs wont traverse all nodes just the given neighbours of u_id.
				for(vertex<T> v : ordered){
					if(v.id == v_id){
						return true;
					}
				}
			}
		}
	}
	return false;
}

template <typename T>
bool directed_graph<T>::contain_cycles()
{	
	vector<vertex<T>> a = get_vertices();
	vector<vertex<T>> b = get_vertices();
	for (vertex<T> v : a){
		for (vertex<T> z : b){
			if ((v.id != z.id)){
				if(reachable(v.id,z.id) && reachable(z.id,v.id)){
					return true;
				}
			}
		}
	}
	return false;
}
    
template <typename T>
vector<vertex<T>> directed_graph<T>::depth_first(const int& u_id) 
{ 	
	stack<T>unprocessed;
	vector<vertex<T>> visited_vertex;
	bool visited[adj_matrix.size()];

    // Intialise all vertices in the matrix to not visited or un visited.
    for (unsigned i = 0; i<adj_matrix.size(); i++){
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
                if(adj_matrix[top_stack][i] != 0)
                {
                    unprocessed.push(i);
                }
            }
        }
		if(unprocessed.empty()){
			for (int i=1;i<adj_matrix.size();i++){
				for(int j=1;j<adj_matrix.size();j++){
					if(adj_matrix[i][j] > 0 && visited[i] == false){
						unprocessed.push(i);
					}
				}
			}
		}
    }
	return visited_vertex;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::breadth_first(const int& u_id) 
{ 
		bool visited[adj_matrix.size()];
		queue<T> unprocessed;
		vector<vertex<T>> ordered;

		// if the index of the start_vertex is valid
		if (contains(u_id)) {
			// set all index values to represent that they have not been visited yet 
			    for (unsigned i = 0; i<adj_matrix.size(); i++) {
       				 visited[i] = false;
    			} 
				// add the start_vertex to the unprocessed queue
				unprocessed.push(u_id);
				// while there is still values in the unprocessed stack
				while (!unprocessed.empty()){
				// get the index of the vertex at the front of the queue and remove it
				int index = unprocessed.front();
				unprocessed.pop();
				// if it hasn't been visited yet
				if (!visited[index]){
					// set it to visited
					visited[index] = true;
					// add the vertex to the ordered list
					ordered.push_back(vertex<T>(index, vertex_weights[index]));
					for (unsigned i = num_vertices(); i != 0; i--){
						// if the vertex contains a neighbour
						if (adj_matrix[index][i] != 0){
							// add the neighbour to the end of the unprocessed queue.
							unprocessed.push(i);
						}
					}
				}
				if(unprocessed.empty()){
					for (int i=1;i<adj_matrix.size();i++){
						for(int j=1;j<adj_matrix.size();j++){
							if(adj_matrix[i][j] > 0 && visited[i] == false){
								unprocessed.push(i);
							}
							if(adj_matrix[i][j] > 0 && visited[j] == false){
								unprocessed.push(j);
							}
						}
					}
				}
			}
		}
		return ordered;
}

template <typename T>
directed_graph<T> directed_graph<T>::out_tree(const int& u_id) { 

	/* 1. Contains no cycle 
	   2. Contains all vertices and all connect 
	   3. Edges = n - 1; where n is num of vertices 
	   4. All vertices reachable from root (doesn't need to contain all vertices)
	   5. Each vertex in tree will have 0, 1 2 out-edges only */
	   
	stack<T>unprocessed;
	vector<vertex<T>> visited_vertex;
	bool visited[adj_matrix.size()];

	directed_graph<T> stree;
	int iteration = 0;
	vector<T> arr_v;
	int n = 0;
	int limit = 0;
	int start = 0;

	// Intialise all vertices in the matrix to not visited or un visited.
    for (unsigned i = 0; i<adj_matrix.size(); i++) {
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
			stree.add_vertex(vertex<T>(top_stack, vertex_weights[top_stack]));
			arr_v.push_back(top_stack);
			iteration += 1;

			if(iteration > 1){
				//everytime we add a vertex, we want to add a edge between the two pair of vertices, this ensures edges = n - 1.
				T & prev = arr_v[n];
				if(reachable(prev,top_stack)==true && (adj_matrix[prev][top_stack] > 0)){
					stree.add_edge(prev, top_stack, adj_matrix[prev][top_stack]);
					n+=1;
				}
				if(stree.num_vertices()-1 != stree.num_edges()){
					for(vertex<T> v : stree.get_vertices()){
						if(stree.reachable(u_id,v.id) == false && reachable(u_id,v.id)){ //if vertex on tree not reachable from root but reachable in graph root, add edge to reach
								//fix below
								limit += 1;
								for(T w : arr_v){
									if(adjacent(w,v.id) && (start < limit) && (stree.out_degree(w)< 3)){
										stree.add_edge(w,v.id, adj_matrix[w][v.id]);
										start += 1;
									}
								}
							}
						}
					}
				}

			//check for exsiting neighbours
            for (unsigned i =num_vertices() ; i!= 0; i--)
            {
                if(adj_matrix[top_stack][i] != 0)
                {
                    unprocessed.push(i);
                }
            }
        }
	}
	return stree;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::pre_order_traversal(const int& u_id, directed_graph<T>& mst) { 
	
	stack<T>unprocessed;
	vector<vertex<T>> visited_vertex;
	bool visited[adj_matrix.size()];

	// Intialise all vertices in the matrix to not visited or un visited.
    for (unsigned i = 0; i<adj_matrix.size(); i++) {
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
                if(adj_matrix[top_stack][i] != 0)
                {
                    unprocessed.push(i);
                }
            }
        }
		if(unprocessed.empty()){
			for (int i=1;i<adj_matrix.size();i++){
				for(int j=1;j<adj_matrix.size();j++){
					if(adj_matrix[i][j] > 0 && visited[i] == false){
						unprocessed.push(i);
					}
				}
			}
		}
    }
	return visited_vertex;
}

struct Node 
{ 
    int data; 
    struct Node* left;
	struct Node* right; 
    Node(int data) 
    { 
        this->data = data; 
        left = right = NULL; 
    } 
}; 

template <typename T>
vector<vertex<T>> directed_graph<T>::in_order_traversal(const int& u_id, directed_graph<T>& mst) { 
//note that tutor has specified to ensure that out_tree will produce a binary tree 
//and also we can consider any (1 child as left or right) and (branch as left and other as right).

//1. if node == null then return
//2. then inorder(node.left)
//3. next visit(node)
//4. finally inorder(node.right)

	vector<vertex<T>> in_order;
	vector<vertex<T>> mst_neighbour = mst.get_neighbours(u_id); //get mst neighbour
	vector<vertex<T>> mst_vertices = mst.get_vertices(); //get mst vertices
	vertex<T> *root = 0;//pointer to vertex<T> that's set to 0

	for(int i=0; i<mst_vertices.size(); i++){ 
		if(mst_vertices[i].id == u_id){
			root = &mst_vertices[i]; 

			if(mst_neighbour.size()>0){
				vector<vertex<T>> left = in_order_traversal(mst_neighbour[0].id, mst); //adds nodes from left
				in_order.instert(in_order.end(), left.begin(), left.end());
			}

			in_order.push_back(*root); //adds the root

			if(mst_neighbour.size()>1){
				vector<vertex<T>> right = in_order_traversal(mst_neighbour[1].id,mst);
				in_order.insert(in_order.end(), right.begin(), right.end()); //adds nodes from right
			}
		}
	}
	return in_order;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::post_order_traversal(const int& u_id, directed_graph<T>& mst) { 

	vector<vertex<T>> post_order;
	vector<vertex<T>> mst_neighbour = mst.get_neighbours(u_id); //get mst neighbour
	vector<vertex<T>> mst_vertices = mst.get_vertices(); //get mst vertices
	vertex<T> *root = 0; //pointer to vertex<T> that's set to 0

	for(int i=0; i<mst_vertices.size(); i++){ 
		if(mst_vertices[i].id == u_id){ 
			root = &mst_vertices[i];//vertex id is passed to pointer(root) only when it matches u_id index 

			if(mst_neighbour.size()>0){
				vector<vertex<T>> left = post_order_traversal(mst_neighbour[0].id, mst); //adds nodes from left
				post_order.instert(post_order.end(), left.begin(), left.end());
			}

			if(mst_neighbour.size()>1){
				vector<vertex<T>> right = post_order_traversal(mst_neighbour[1].id,mst);
				post_order.insert(post_order.end(), right.begin(), right.end()); //adds nodes from right
			}

			post_order.push_back(*root); //adds the root
		}
	}
	return post_order;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::significance_sorting() { 
	//Using bubble sort
	vector<T> arrv;
	int i, j;
	int n = arrv.size();
	int temp;
	vector<vertex<T>> v = get_vertices();
	vector<vertex<T>> sorted;

	for(vertex<T> a : v){
		arrv.push_back(a.id);
	}
	for(i=0;i < n-1;i++){
		for(j=0; j < n-i-1; j++){
			if(arrv[j] > arrv[j+1]){
				temp = arrv[j];
				arrv[j]=arrv[j+1];
				arrv[j+1] = temp;
			}
		}
	}
	for(int i=0;i<arrv.size();i++){
		sorted.push_back(vertex<T>(arrv[i],vertex_weights[arrv[i]]));
	}
	return sorted;
}  

template <typename T>
int directed_graph<T>::get_edgeWeight(const int& u, const int& v) {
	return adj_matrix[u][v];
}

template <typename T>
int directed_graph<T>::get_vertexWeight(const int& u) {
	return vertex_weights[u];
}

template <typename T>
int directed_graph<T>::get_adjMatrixSize() {
	return adj_matrix.size();
}

#endif