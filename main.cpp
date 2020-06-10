#include <iostream>
#include "directed_graph.hpp"
#include "directed_graph_algorithms.cpp"

int main() {

	directed_graph<double> g1;

	vertex<double> A(1, 0.25); 
	vertex<double> B(2, 1.41);
	vertex<double> C(3, 2.32);
	vertex<double> D(4, 3.66);
	vertex<double> E(5, 4.12);
	vertex<double> F(6, 2.12);
	vertex<double> G(7, 2.12);

	g1.add_vertex(A); // add vertex class
	g1.add_vertex(B);
	g1.add_vertex(C);
	g1.add_vertex(D);
	g1.add_vertex(E);
	g1.add_vertex(F);
	g1.add_vertex(G);

	//g1.remove_vertex(1); // delete vertex by id
	cout << "num of vertices: " << g1.num_vertices() << endl;
	vector<vertex<double>> vertex_list = g1.get_vertices();
	cout << "all vertices: ";
	for (vertex<double> vt : vertex_list) {
	 	cout << "(" << vt.id << ", " << vt.weight << ") ";
	}
	cout << endl;

	g1.add_edge(1, 2, 10); //(u, v, edge weight)
	g1.add_edge(1, 3, 10); //(u, v, edge weight)
	g1.add_edge(1, 4, 10); 
	g1.add_edge(2, 1, 20);
	g1.add_edge(2, 5, 30);
	g1.add_edge(4, 2, 10);
	g1.add_edge(4, 7, 10);
	g1.add_edge(5, 4, 40);
	g1.add_edge(5, 6, 40);
	g1.add_edge(6, 5, 30);
	g1.add_edge(6, 7, 30);
	g1.add_edge(7, 1, 30);

	cout << "num of edges: " << g1.num_edges() << endl;

	cout << "all neighbours of 2: ";
	vector<vertex<double>> neighbour_list = g1.get_neighbours(2);
	for (vertex<double> nb : neighbour_list) {
	 	cout << "(" << nb.id << ", " << nb.weight << ") ";
	}
	cout << endl;

	if (g1.adjacent(2,3)){
		cout << "Edges 2, 3 are adjacent" << endl;
	}

	g1.out_degree(2);
	cout << "In-degree of vertex 3: " << g1.in_degree(3) << endl;
	
	for (int i = 2; i < 6; i++){
		cout << "all second-order neighbours of " << i << ": ";
		vector<vertex<double>> two_neighbour_list = g1.get_second_order_neighbours(i);
		//cout << "size: " << two_neighbour_list.size() << endl;
			for (vertex<double> nb : two_neighbour_list){
				cout <<  "(" << nb.id << ", " << nb.weight << ") ";
			}
			cout << endl;
	}

	cout << "Num of vertices " << g1.num_vertices() << endl;

	cout << "depth first search below:" << endl;
	vector<vertex<double>> dfs = g1.depth_first(1);
	for (vertex<double> ver : dfs){
		cout <<  "(" << ver.id << ", " << ver.weight << ") ";
	}
	cout << endl;

	cout << "breadth first search below:" << endl;
	vector<vertex<double>> bfs = g1.breadth_first(2);
	for (vertex<double> ver : bfs){
		cout <<  "(" << ver.id << ", " << ver.weight << ") ";
	}

	g1.contain_cycles();
	g1.reachable(3,4);

	directed_graph<double> g2;
	g2 = g1.out_tree(1);
	cout << endl << "trees num of vertices: " << g2.num_vertices();
	cout << endl << "trees num of edges: " << g2.num_edges() << endl;
	vector<vertex<double>> vertex_list2 = g2.get_vertices();
	cout << "all vertices of tree: ";
	for (vertex<double> vc : vertex_list2) {
	 	cout << "(" << vc.id << ", " << vc.weight << ") ";
	}
	cout << endl;
	cout << "Is 3 reachable from 1?: " << g2.reachable(1,3) << endl;

	cout << "Shortest Path" << endl;
    vector<vertex<int>> short_path = shortest_path(g1, 1, 4);
    for (auto vertex : short_path) {
            cout << "(" << vertex.id << ", " << vertex.weight << ")";
    }
    cout << endl;

    vector<vector<vertex<int>>> scc = strongly_connected_components(g1);

    cout << "Topological Sort" << endl;
    vector<vertex<int>> top_sort = topological_sort(g1);
    for (auto vert : top_sort) {
        cout << "(" << vert.id << ", " << vert.weight << ")";
    }
    cout << endl;

    cout << "Low Cost Delivery: " << low_cost_delivery(g1,1);
}