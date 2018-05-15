/* Implmentation of Bellman Ford algorithm.
 * Author: Raphaël Dang-Nhu.
 * Date: 15/05 */

#include <vector>

using namespace std;

class edge {
    
    public:
    edge(int target_, double weight_) :
        target(target_),
        weight(weight_)
    {}
    void print();

    int target;
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
    void print();
    // Function to perform Bellman-Ford algorithm on the graph. Takes as a parameter the index of the origin, an array to put distances, and an array to put predecessors
    void bellmanFord(int origin, vector<double> &distances, vector<int> &predecessors);
   
    int nVertices;
    vector<node*> nodes;
};

