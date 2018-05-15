/* Implmentation of Bellman Ford algorithm.
 * Author: Raphaël Dang-Nhu.
 * Date: 15/05 */

#include <vector>

#include "interval_arithmetic.hpp"

using namespace std;

class edge {
    
    public:
    edge(int source_, double weight_) :
        source(source_),
        weight(weight_)
    {}
    void print();

    int source;
    double weight;
};

class node {

    public:
    node(int index_) :
        index(index_)
    {}

    int index;
    vector<edge*> adjacentEdges;
};

class Graph {

    public:
    // Constructor to generate a Erdos-Renyi random graph with random edge weights
    Graph(int nVertices, double edgeProba, int emin, int emax, int negratio);
    // Function to print a graph
    void printGraph();
    // Function to perform Bellman-Ford algorithm on the graph. Takes as a parameter the index of the origin, an array to put distances, and an array to put predecessors
    bool bellmanFord(int origin, vector<double> &distances, vector<int> &predecessors);
    // Same function with interval arithmetic
    boolean intervalBellmanFord(int origin, vector<__m128d> &distance, vector<int> &predecessors, vector<vector<int>>& totalc, vector<vector<int>>& undefinedc);
   
    int nVertices;
    vector<node*> nodes;
};

