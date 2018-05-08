/* Testing Boost Graph Lib.
 * Author: Raphaël Dang-Nhu.
 * Date: 08/05 */

#include <iostream>
#include <utility>
#include <algorithm>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace std;
using namespace boost;

void test0(){

    // create a typedef for the graph type
    typedef adjacency_list<vecS,vecS, bidirectionalS> Graph;

    // Create evrtices labels
    enum { A,B,C,D,E,N};
    const int num_vertices = N;
    const char* name = "ABCDE";

    // writing out edges in the graph
    typedef std::pair<int, int> Edge;
    Edge edge_array[] = {Edge(A,B),Edge(A,D),Edge(B,C)};
    const int num_edges = sizeof(edge_array)/sizeof(edge_array[0]);

    // Declare graph object
    Graph g(num_vertices);

    // add the edges to the graph object
    for(int i = 0; i < num_edges; i++){
        add_edge(edge_array[i].first,edge_array[i].second, g);
    }
}

int main(){
    test0();
    cout << endl << "Finished test" << endl;
    return 0;
}
