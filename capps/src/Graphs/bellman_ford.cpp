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
    
    // Create a struct to hold the weights
    struct edge_info{
        string name;
        int weight;
    };
    
    // create a typedef for the graph type
    typedef adjacency_list<vecS,vecS, bidirectionalS,no_property,edge_info> Graph;
    typedef graph_traits<Graph>::edge_iterator edge_iterator; 
    
    // Declare a graph with 10 vertices
    Graph g(10);
    
    // Add some edges
    edge_info i1;
    i1.name = "First edge";
    i1.weight = 10;

    edge_info i2;
    i2.name = "Second edge";
    i2.weight = -5;

    add_edge(0,1,g);
    add_edge(1,2,g);

    // Print edges
    pair<edge_iterator,edge_iterator> ei = edges(g);
    cout << "Number of edges: " << num_edges(g) << endl;
    cout << "Edge list: " << endl;

    for(;ei.first != ei.second;++ei.first){
        cout << "Name: " << (*ei.first.name) << endl;
        cout << "Weight: " << (*ei.first.weight) << endl;
        cout << endl << "(" << source(*ei.first,g) << ",";
        cout << target(*ei.first,g) << ")" << endl;
    }
}

int main(){
    test0();
    cout << endl << "Finished test" << endl;
    return 0;
}
