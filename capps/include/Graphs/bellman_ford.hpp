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

// Enum to represent a comparison status
enum CompResult{
    newOptimum,
    noChange,
    undefined,
    useless
};

// typedef to memorize information
typedef vector<vector<vector<CompResult>>> compType;

// Struct returned by the bellman Ford method
struct bellmanFordResult{
    bool b;
    vector<double> distances;
    vector<int> pred;
};

// Struct returned by the interval bellman Ford method
struct intervalBellmanFordResult{
    boolean b;
    vector<__m128d> distances;
    vector<int> pred;
    vector<vector<int>> totalc;
    vector<vector<int>> undefinedc;
    compType memo;
};

// Struct returned by the superacc bellman Ford result
struct mpfrBellmanFordResult{
    bool b;
    vector<double> distances;
    vector<int> pred;
};
class Graph {

    public:

    // Constructor to generate a Erdos-Renyi random graph with random edge weights
    Graph(int nVertices, double edgeProba, int emin, int emax, int negratio, double offset);

    // Function to print a graph
    void printGraph();

    // Function to perform Bellman-Ford algorithm on the graph. Takes as a parameter the index of the origin, an array to put distances, and an array to put predecessors
    bellmanFordResult bellmanFord(int origin);

    // Same function with interval arithmetic
    intervalBellmanFordResult intervalBellmanFord(int origin);
   
    // Reverse processing
    compType reverseBellmanFord(int origin,compType c);

    // Same function with mpfr
    mpfrBellmanFordResult mpfrBellmanFord(int origin);
    
    // Same function with lazy
    mpfrBellmanFordResult lazyMpfrBellmanFord(int origin, compType c);
    
    int nVertices;
    vector<node*> nodes;
};

// Function to print comparison results 
void printComps(vector<vector<int>> undefc, vector<vector<int>> totalc);
