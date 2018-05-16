/* Class to implement Bellam Ford algorithm.
 * Author: Raphaël Dang-Nhu.
 * Date: 08/05 */

#include "debug.hpp"

#include <iostream>
#include <random>
#include<math.h>
#include <emmintrin.h>

#include <mpreal.h>
#include <mpfr.h>
#include <gmp.h>

#include "pfpdefs.hpp"
#include "common.hpp"
#include "bellman_ford.hpp"
#include "interval_arithmetic.hpp"

using namespace std;
using mpfr::mpreal;

struct mpfr_wrapper{
    mpfr_t val;
};

template<class t> void  printVector(vector<t> v){
    for(typename vector<t>::iterator it = v.begin(); it != v.end(); it++){
        cout << *it << ",";
    }
    cout << endl;
}

void  printVectorInt(vector<__m128d> v){
    for(vector<__m128d>::iterator it = v.begin(); it != v.end(); it++){
        print(*it);
        cout << endl;
    }
}

void printComps(vector<vector<int>> undefc, vector<vector<int>> totalc){
    for(unsigned int i = 0; i!=undefc.size(); i++){

        for(unsigned int j = 0; j!= undefc[i].size(); j++){
            cout << undefc[i][j] << "/" << totalc[i][j] << " ";
        }
    cout << endl;
    }
    cout << endl;
}

void edge::print(){
    cout << "source: " << source << endl;
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
    for(int i = 0; i !=nVertices; i++){
        for(int j = i+1; j!= nVertices; j++){

            double test = rand() / (double) (RAND_MAX);
            if(test <= edgeProba){
                // Generate weight
                double w = randDouble(emin,emax,negratio);
                // Create edges
                edge* ei = new edge(i,w);
                nodes[j]->adjacentEdges.push_back(ei);
            }
            test = rand() / (double) (RAND_MAX);
            if(test <= edgeProba){
                double w = randDouble(emin,emax,negratio);
                edge* ej = new edge(j,w);
                nodes[i]->adjacentEdges.push_back(ej);
            }
        }
    }
}

void Graph::printGraph(){
    for(vector<node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        cout << endl << "Vertice: " << (*it)->index << endl;
        for(vector<edge*>::iterator it0 = (*it)->adjacentEdges.begin(); it0 != (*it)->adjacentEdges.end(); it0++){
            (*it0)->print();
        }
    }
}

bool Graph::bellmanFord(int origin, vector<double> &distances, vector<int> &predecessors){
    
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
                int sourceIndex = (*it0)->source;
                double aux = distancesAux[i-1][sourceIndex] + (*it0)-> weight;
                if(aux < distancesAux[i][nodeindex]){
                    distancesAux[i][nodeindex] = aux;
                    predecessors[nodeindex] = sourceIndex;
                }
            }
        }
    }
    distances = distancesAux[nVertices-1];

    // Check that no negative cycles
    for(vector<node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        // For each edge adjacent to this node
        for(vector<edge*>::iterator it0 = (*it)->adjacentEdges.begin(); it0 != (*it)->adjacentEdges.end(); it0++){

            // If using this edge makes the distance shorter...
            int sourceIndex = (*it0)->source;
            double aux = distancesAux[nVertices-1][sourceIndex] + (*it0)-> weight;
            if(aux < distancesAux[nVertices-1][(*it)->index]){
                return false;
            }
        }
    }
    return true;
}

boolean Graph::intervalBellmanFord(int origin, vector<__m128d> &distances, vector<int> &predecessors, vector<vector<int>> & totalComps, vector<vector<int>> & undefComps){
    
    // Create structure to memoize distances
    vector<vector<__m128d>> distancesAux;
    for(int i = 0; i != nVertices; i++){
        vector<__m128d> distancesLine(nVertices,in2_create(INFINITY,INFINITY));
        distancesAux.push_back(distancesLine);
    }
    distancesAux[0][origin] = in2_create(0,0);

    // Create structure to memoize undefined comparisons
    for(int i = 0; i != nVertices+1; i++){
        vector<int> undefCompsLine(nVertices,0);
        vector<int> totalCompsLine(nVertices,0);
        undefComps.push_back(undefCompsLine);
        totalComps.push_back(totalCompsLine);
    }
    
    // Main loop
    for(int i = 1; i != nVertices; i++){
        // For each node
        for(vector<node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
            // Counting undefined comparisons
            int undefinedc = 0;
            int totalc = 0;

            int nodeindex = (*it)->index;
            distancesAux[i][nodeindex] = distancesAux[i-1][nodeindex];

            // For each edge adjacent to this node
            for(vector<edge*>::iterator it0 = (*it)->adjacentEdges.begin(); it0 != (*it)->adjacentEdges.end(); it0++){
    
                // If using this edge makes the distance shorter...
                int sourceIndex = (*it0)->source;
                __m128d aux = in2_add_double(distancesAux[i-1][sourceIndex],(*it0)-> weight);
                /*cout << "node: " << nodeindex<< endl;
                cout << "source: " << sourceIndex << endl;
                print(aux);
                cout << endl;
                print(distancesAux[i][nodeindex]);
                cout << endl;*/
                boolean b = inferior(distancesAux[i][nodeindex],aux);
                if(b == False){
                    distancesAux[i][nodeindex] = aux;
                    predecessors[nodeindex] = sourceIndex;
                }else if (b == Undefined){
                    /*cout<< endl;
                    cout << nodeindex;
                    cout<< endl;
                    cout << predecessors[nodeindex];
                    cout << endl;
                    cout << sourceIndex;
                    cout << endl;
                    print(distancesAux[i][nodeindex]);
                    cout << endl;
                    print(aux);
                    cout << endl;
                */
                    distancesAux[i][nodeindex] = in2_min(aux,distancesAux[i][nodeindex]);
                    undefinedc += 1;
                    // What do we do for predecessors ? 
                }
                totalc += 1;
            }

            totalComps[i][nodeindex] = totalc;
            undefComps[i][nodeindex] = undefinedc;
        }
    }
    distances = distancesAux[nVertices-1];
    
    boolean result = True;
    // Check that no negative cycles
    for(vector<node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        int undefinedc = 0;
        int totalc = 0;
        int nodeindex = (*it)->index;
        // For each edge adjacent to this node
        for(vector<edge*>::iterator it0 = (*it)->adjacentEdges.begin(); it0 != (*it)->adjacentEdges.end(); it0++){

            // If using this edge makes the distance shorter...
            int sourceIndex = (*it0)->source;
            __m128d aux = in2_add_double(distancesAux[nVertices-1][sourceIndex],(*it0)-> weight);

            boolean b = inferior(distancesAux[nVertices-1][nodeindex],aux);
            if(b == False) result = False;
            else if(b == Undefined && result == True){
                /*cout << endl;
                print(aux);
                cout << endl;
                print(distancesAux[nVertices-1][nodeindex]);
                cout << endl;
                */
                undefinedc += 1;
                result = Undefined;
            }
            totalc += 1;
        }
        totalComps[nVertices][nodeindex] = totalc;
        undefComps[nVertices][nodeindex] = undefinedc;
    }
    return result;
}

bool Graph::mpfrBellmanFord(int origin, vector<double> &distances, vector<int> &predecessors){

    // Create structure to memoize distances
    vector<vector<mpreal>> distancesAux(nVertices);
    for(int i = 0; i != nVertices; i++){
        
        distancesAux[i] = vector<mpreal>(nVertices);
        for(int j = 0; j!= nVertices; j++){
            distancesAux[i][j] = mpfr::const_infinity(1,mpreal::get_default_prec());
        }
        
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

                int sourceIndex = (*it0)->source;
                // If using this edge makes the distance shorter...
                
                mpreal aux = distancesAux[i-1][sourceIndex]+(*it0)-> weight;
                if(distancesAux[i][nodeindex]>aux){
                    distancesAux[i][nodeindex] = aux;
                    predecessors[nodeindex] = sourceIndex;
                }
            }
        }
    }
    
    for(unsigned int i = 0; i != distancesAux[nVertices-1].size(); i++){
        distances[i] = (double)distancesAux[nVertices-1][i];
    }

    // Check that no negative cycles
    for(vector<node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        // For each edge adjacent to this node
        for(vector<edge*>::iterator it0 = (*it)->adjacentEdges.begin(); it0 != (*it)->adjacentEdges.end(); it0++){

            int sourceIndex = (*it0)->source;
            // If using this edge makes the distance shorter...
            mpreal aux = distancesAux[nVertices-1][sourceIndex]+(*it0)-> weight;
            if(distancesAux[nVertices-1][(*it)->index]>aux){
                return false;
            }
        }
    }
    return true;
    
}
                
void test(){
    int n = 200;

    Graph g(n,1,-1000,1000,2);    
    g.printGraph();
    
    // Test with doubles
    vector<double> distances;
    vector<int> pred(n,-1);
    pred[0] = 0;
  
    cout << endl;
    bool b = g.bellmanFord(0,distances,pred);
    if(b){
        cout << "No negative cycles" << endl;
    }
    else{
        cout << "Negative cycle" << endl;
    }
    
    //cout << "Distances: " << endl;
    //printVector(distances);
    cout << "Predecessors: " << endl;
    printVector(pred);

    // Test with interval arithmetic
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    vector<__m128d> distancesI;
    vector<int> predI(n,-1);
    vector<vector<int>> totalc;
    vector<vector<int>> undefc;
    predI[0] = 0;
  
    cout << endl;
    boolean b1 = g.intervalBellmanFord(0,distancesI,predI,totalc,undefc); 
    if(b1 == True){
        cout << "No negative cycles" << endl;
    }
    else if ( b1 == False) {
        cout << "Negative cycle" << endl;
    }
    else if (b1 == Undefined){
        cout << "Uncertain presence of negative cycles" << endl;
    }
    //cout << "Distances: " << endl;
    //printVectorInt(distancesI);
    cout << "Predecessors: " << endl;
    printVector(predI);
    //cout << "Undefined comparisons" << endl;
    //printComps(undefc,totalc);
    _MM_SET_ROUNDING_MODE(0);

    // Test with mpfr
    vector<double> distancesM(n);
    vector<int> predM(n,-1);
    predM[0] = 0;
  
    cout << endl;
    if(g.mpfrBellmanFord(0,distancesM,predM)){
        cout << "No negative cycles" << endl;
    }
    else{
        cout << "Negative cycle" << endl;
    }
    //cout << "Distances: " << endl;
    //printVector(distancesM);
    cout << "Predecessors: " << endl;
    printVector(predM);

    // Check if different result
    int d = 0;
    for(unsigned int i = 0; i != pred.size(); i++){
        if(pred[i] != predM[i]){
            d++;
        }
    }
    cout << endl << "Number of different results: "<< d <<  endl;
}

void runtimeTest(){
   
    // Generate a graph
    int n = 100;
    Graph g(n,0.5,-1000,1000,2);    

    // Variables
    vector<double> distances;
    vector<__m128d> distancesI;
    vector<double> distancesM(n);
    vector<int> predecessors(n,-1);
    vector<int> predecessorsI(n,-1);
    vector<int> predecessorsM(n,-1);
    vector<vector<int>> undefc;
    vector<vector<int>> totalc;

    // Tests
    double start = 0.0, time_double = 0.0, time_interval = 0.0, time_mpfr = 0.0;

    PFP_TIME(g.bellmanFord(0,distances,predecessors),start,time_double);
    PFP_TIME(g.intervalBellmanFord(0,distancesI,predecessorsI,undefc,totalc),start,time_interval);
    PFP_TIME(g.mpfrBellmanFord(0,distancesM,predecessorsM),start,time_mpfr);

    cout << endl << time_interval/time_double << endl;
    cout << endl << time_mpfr/time_double << endl;

}

int main(){
    srand(time(NULL));
    for(int i = 0; i != 1; i++){
        test();
    }
}

