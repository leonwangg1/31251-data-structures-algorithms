#include "directed_graph.hpp"

int main() {

	directed_graph<double> g1;

	vertex<double> v0(1, 0.25); 
	vertex<double> v1(2, 1.41);
	vertex<double> v2(3, 2.32);
	vertex<double> v3(4, 3.66);
	vertex<double> v4(5, 4.12);

	g1.add_vertex(v0); // add vertex class
	g1.add_vertex(v1);
	g1.add_vertex(v2);
	g1.add_vertex(v3);
	g1.add_vertex(v4);

	g1.remove_vertex(1); // delete vertex by id

	vector<vertex<double>> vertex_list = g1.get_vertices();
	cout << "all vertices: ";
	for (vertex<double> vt : vertex_list) {
	 	cout << "(" << vt.id << ", " << vt.weight << ") ";
	}
	cout << endl;

	g1.add_edge(2, 3, 10); //(u, v, edge weight)
	g1.add_edge(2, 4, 20);
	g1.add_edge(3, 4, 30);
	g1.add_edge(4, 2, 40);
	g1.add_edge(4, 3, 50);
	g1.add_edge(5, 3, 60);

	cout << "all neighbours of 5: ";
	vector<vertex<double>> neighbour_list = g1.get_neighbours(5);
	for (vertex<double> nb : neighbour_list) {
	 	cout << "(" << nb.id << ", " << nb.weight << ") ";
	}
	cout << endl;

	if (g1.adjacent(3,2)){
		cout << "yay";
	}

	g1.out_degree(2);
	
	cout << "all second-order neighbours of 5: ";
    vector<vertex<double>> two_neighbour_list = g1.get_second_order_neighbours(5);
	for (vertex<double> nb : two_neighbour_list){
		cout <<  "(" << nb.id << ", " << nb.weight << ") ";
	}
	cout << endl;
	
}