#include <iostream>
#include "directed_graph.hpp"
#include "directed_graph_algorithms.cpp"

int main() {

    directed_graph<int> g1;

    vertex<int> va = vertex<int>(0, 800);        //A
    vertex<int> vb = vertex<int>(1, 300);        //B
    vertex<int> vc = vertex<int>(2, 400);        //C
    vertex<int> vd = vertex<int>(3, 710);        //D
    vertex<int> ve = vertex<int>(4, 221);        //E;

    g1.add_vertex(va);
    g1.add_vertex(vb);
    g1.add_vertex(vc);
    g1.add_vertex(vd);
    g1.add_vertex(ve);

	//g1.remove_vertex(1); // delete vertex by id
	cout << "num of vertices: " << g1.num_vertices() << endl;
	vector<vertex<int>> vertex_list = g1.get_vertices();
	cout << "all vertices: ";
	for (vertex<int> vt : vertex_list) {
		cout << "(" << vt.id << ", " << vt.weight << ") ";
	}
	cout << endl;

    g1.add_edge(va.id, vb.id, 600); //A-B 6
    g1.add_edge(va.id, vc.id, 900); //A-C 9
    g1.add_edge(vb.id, ve.id, 3000); //B-E 3
    g1.add_edge(vc.id, vd.id, 4000); //C-D 4
    g1.add_edge(vd.id, va.id, 1); //D-A 1
    g1.add_edge(vd.id, vc.id, 700); //D-C 7
    g1.add_edge(vd.id, ve.id, 500); //D-E 5

	cout << "num of edges: " << g1.num_edges() << endl;

	cout << "all neighbours of 2: ";
	vector<vertex<int>> neighbour_list = g1.get_neighbours(2);
	for (vertex<int> nb : neighbour_list) {
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
		vector<vertex<int>> two_neighbour_list = g1.get_second_order_neighbours(i);
		//cout << "size: " << two_neighbour_list.size() << endl;
			for (vertex<int> nb : two_neighbour_list){
				cout <<  "(" << nb.id << ", " << nb.weight << ") ";
			}
		cout << endl;
	}

	cout << "Num of vertices " << g1.num_vertices() << endl;

	cout << "depth first search below:" << endl;
	vector<vertex<int>> dfs = g1.depth_first(1);
	for (vertex<int> ver : dfs){
		cout <<  "(" << ver.id << ", " << ver.weight << ") ";
	}
	cout << endl;

	cout << "breadth first search below:" << endl;
	vector<vertex<int>> bfs = g1.breadth_first(2);
	for (vertex<int> ver : bfs){
		cout <<  "(" << ver.id << ", " << ver.weight << ") ";
	}

	g1.contain_cycles();
	g1.reachable(3,4);

	directed_graph<int> g2;
	g2 = g1.out_tree(1);
	cout << endl << "trees num of vertices: " << g2.num_vertices();
	cout << endl << "trees num of edges: " << g2.num_edges() << endl;
	vector<vertex<int>> vertex_list2 = g2.get_vertices();
	cout << "all vertices of tree: ";
	for (vertex<int> vc : vertex_list2) {
	 	cout << "(" << vc.id << ", " << vc.weight << ") ";
	}
	cout << endl;
	cout << "Is 3 reachable from 1?: " << g2.reachable(1,3) << endl;

	cout << "Shortest Path" << endl;
    vector<vertex<int>> short_path = shortest_path(g1, 0, 4);
	for (auto vertex : short_path) {
		cout << "(" << vertex.id << ", " << vertex.weight << ")";
	}
	cout << endl;

    cout << "Strongly Connected Components" << endl;
    vector<vector<vertex<int>>> scc = strongly_connected_components(g1);
    for (auto vector : scc) {
      for (auto vertex : vector) {
        cout << "(" << vertex.id << ", " << vertex.weight << ")";
      }
    cout << endl;
    }

    cout << "Topological Sort" << endl;
    vector<vertex<int>> top_sort = topological_sort(g1);
    for (auto vert : top_sort) {
        cout << "(" << vert.id << ", " << vert.weight << ")";
    }
    cout << endl;

    cout << "Low Cost Delivery: " << low_cost_delivery(g1,0); //D 
}
