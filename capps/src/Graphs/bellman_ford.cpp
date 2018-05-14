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

#include <random>

using namespace std;
using namespace boost;

double randDouble(int emin, int emax, int neg_ratio) {
    // Uniform mantissa
    double x = double(rand()) / double(RAND_MAX * .99) + 1.;
    // Uniform exponent
    int e = (rand() % (emax - emin)) + emin;
    // Sign
    if(neg_ratio > 1 && rand() % neg_ratio == 0) {
        x = -x;
    }
    return ldexp(x, e);
}

void test0(){
    
    // Create property
    typedef property<edge_weight_t,double> EdgeProperty; 

    // create a typedef for the graph type
    typedef adjacency_list<vecS,vecS, bidirectionalS,no_property,EdgeProperty> Graph;
    typedef graph_traits<Graph>::edge_iterator edge_iterator; 
    
    // Declare a graph with 10 vertices
    Graph g(10);
   
    // Obtain property map
    property_map<Graph, edge_weight_t>::type w = get(edge_weight,g);

    // Add some edges
    add_edge(0,1,10,g);
    add_edge(1,2,20,g);

    // Print edges
    pair<edge_iterator,edge_iterator> ei = edges(g);
    cout << "Number of edges: " << num_edges(g) << endl;
    cout << "Edge list: " << endl;

    for(;ei.first != ei.second;++ei.first){
        cout << endl << "Vertice 1: " << source(*ei.first,g)<< endl;
        cout << endl << "Vertice 2: " << target(*ei.first,g)<< endl;
        cout << "Weight: " << w[*ei.first] << endl;
    }
}

// Function to test random graph generation */
void test1(){
    
    // Create property
    typedef boost::property<boost::edge_weight_t, double> WeightProperty;

    typedef adjacency_list<vecS,vecS,bidirectionalS,no_property,WeightProperty> Graph;
    typedef graph_traits<Graph>::edge_iterator edge_iterator; 
    typedef sorted_erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;
   
    
    // Generate Erdos-Renyi random graph
    boost::minstd_rand gen(0);
    Graph g(ERGen(gen,10,0.5),ERGen(),10);

    // Properties
    property_map<Graph, edge_weight_t>::type w = get(edge_weight,g);
    // To store parents

      IndexMap indexMap = boost::get(boost::vertex_index, g);
        PredecessorMap predecessorMap(&predecessors[0], indexMap);
          DistanceMap distanceMap(&distances[0], indexMap);

    // Give random weights to the edges
    pair <edge_iterator,edge_iterator> eiw = edges(g);
    for(;eiw.first != eiw.second;++eiw.first){
        w[*eiw.first] = randDouble(-100,100,2);      
    }

    // Print edges
    pair <edge_iterator,edge_iterator> ei = edges(g);
    cout << "Number of edges: " << num_edges(g) << endl;
    cout << "Edge list: " << endl;

    for(;ei.first != ei.second;++ei.first){
        cout << endl << "(" << source(*ei.first,g) << ",";
        cout << target(*ei.first,g) << ")" << endl;
        cout << "Weight: " << w[*ei.first] << endl;
    }

    // Bellman-Ford algorithm
    bellman_ford_shortest_paths(g,0,boost::distance_map(
}

int main(){
    test1();
    cout << endl << "Finished test" << endl;
    return 0;
}
