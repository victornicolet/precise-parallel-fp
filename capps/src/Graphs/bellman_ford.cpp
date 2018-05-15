/* Class to implement Bellam Ford algorithm.
 * Author: Raphaël Dang-Nhu.
 * Date: 08/05 */

#include <iostream>
#include <random>
#include<math.h>

#include "common.hpp"
#include "bellman_ford.hpp"

using namespace std;

template<class t> void  printVector(vector<t> v){
    for(typename vector<t>::iterator it = v.begin(); it != v.end(); it++){
        cout << *it << ",";
    }
    cout << endl;
}

void edge::print(){
    cout << "Target: " << target << endl;
    cout << "Weight: " << weight << endl;
}

Graph::Graph(int nVertices, double edgeProba, int emin, int emax, int negratio) : 
    nVertices(nVertices),
    nodes(nVertices)    
{
    // Generate all nodes
    for(int i = 0; i != nVertices; i++){
        nodes[i] = new node(i);
    }

    // For all pair of nodes, generate random edges
    srand(time(NULL));
    for(int i = 0; i !=nVertices; i++){
        for(int j = i+1; j!= nVertices; j++){

            double test = rand() / (double) (RAND_MAX);
            if(test <= edgeProba){
                // Generate weight
                double w = randDouble(emin,emax,negratio);
                // Create edges
                edge* ei = new edge(j,w);
                nodes[i]->adjacentEdges.push_back(ei);
                edge* ej = new edge(i,w);
                nodes[j]->adjacentEdges.push_back(ej);
            }
        }
    }
}

void Graph::print(){
    for(vector<node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        cout << endl << "Vertice: " << (*it)->index << endl;
        for(vector<edge*>::iterator it0 = (*it)->adjacentEdges.begin(); it0 != (*it)->adjacentEdges.end(); it0++){
            (*it0)->print();
        }
    }
}

void Graph::bellmanFord(int origin, vector<double> &distances, vector<int> &predecessors){
    
    // Create structure to memoize distances
    vector<vector<double>> distancesAux;
    for(int i = 0; i != nVertices; i++){
        vector<double> distancesLine(nVertices,INFINITY);
        distancesAux.push_back(distancesLine);
    }
    distancesAux[0][origin] = 0;
    
    // Main loop
    for(int i = 1; i != nVertices; i++){
        // For each node
        for(vector<node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
            int nodeindex = (*it)->index;
            distancesAux[i][nodeindex] = distancesAux[i-1][nodeindex];
            // For each edge adjacent to this node
            for(vector<edge*>::iterator it0 = (*it)->adjacentEdges.begin(); it0 != (*it)->adjacentEdges.end(); it0++){

                // If using this edge makes the distance shorter...
                int sourceIndex = (*it0)->target;
                double aux = distancesAux[i-1][sourceIndex] + (*it0)-> weight;
                if(aux < distancesAux[i][nodeindex]){
                    distancesAux[i][nodeindex] = aux;
                    predecessors[nodeindex] = sourceIndex;
                }
            }
        }
    }
    distances = distancesAux[nVertices-1];
}
                
int main(){
    int n = 5;

    Graph g(n,0.01,-100,100,1);    
    g.print();

    vector<double> distances;
    vector<int> pred(n,-1);
    pred[0] = 0;
    g.bellmanFord(0,distances,pred);
   
    cout << endl;
    cout << "Distances: " << endl;
    printVector(distances);
    cout << "Predecessors: " << endl;
    printVector(pred);
    
}



