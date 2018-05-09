/* Testing Boost Graph Lib.
 * Author: Raphaël Dang-Nhu.
 * Date: 08/05 */

#include <iostream>
#include <utility>
#include <algorithm>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/random.hpp>

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

    add_edge(0,1,i1,g);
    add_edge(1,2,i2,g);

    // Print edges
    pair<edge_iterator,edge_iterator> ei = edges(g);
    cout << "Number of edges: " << num_edges(g) << endl;
    cout << "Edge list: " << endl;

    for(;ei.first != ei.second;++ei.first){
        cout << endl << "Name: " << g[*ei.first].name << endl;
        cout << "Weight: " << g[*ei.first].weight << endl;
        cout <<  "(" << source(*ei.first,g) << ",";
        cout << target(*ei.first,g) << ")" << endl;
    }
}

// Function to test random graph generation */
void test1(){
    
    struct edge_info{
        edge_info(int x) : weight(x) {}
        int weight;
    };

    typedef adjacency_list<vecS,vecS,bidirectionalS,no_property,edge_info> Graph;
    typedef graph_traits<Graph>::edge_iterator edge_iterator; 
    typedef sorted_erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;
   
    // Generate Erdos-Renyi random graph
    boost::minstd_rand gen(0);
    /*Graph g(ERGen(gen,10,0.5),ERGen(),10);

    // Give random weights to the edges
    //randomize_property<edge_bundle_t>(g,gen);

    // Print edges
    pair <edge_iterator,edge_iterator> ei = edges(g);
    cout << "Number of edges: " << num_edges(g) << endl;
    cout << "Edge list: " << endl;

    for(;ei.first != ei.second;++ei.first){
        cout << endl << "(" << source(*ei.first,g) << ",";
        cout << target(*ei.first,g) << ")" << endl;
        cout << "Weight: " << g[*ei.first].weight << endl;
    }*/
    
}

int main(){
    test0();
    cout << endl << "Finished test" << endl;
    return 0;
}
