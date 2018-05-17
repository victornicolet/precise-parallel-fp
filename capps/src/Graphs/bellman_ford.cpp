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

#include "common.hpp"
#include "bellman_ford.hpp"
#include "interval_arithmetic.hpp"

using namespace std;
using mpfr::mpreal;

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

Graph::Graph(int nVertices, double edgeProba, int emin, int emax, int negratio, double offset) : 
    nVertices(nVertices),
    nodes(nVertices)    
{
    // Generate all nodes
    for(int i = 0; i != nVertices; i++){
        nodes[i] = new node(i);
    }

    // For all pair of nodes, generate random edges
    for(int j = 1; j!= nVertices; j++){

        double test = rand() / (double) (RAND_MAX);
        if(test <= edgeProba){
            // Generate weight
            double w = offset;
            // Create edges
            edge* ei = new edge(0,w);
            nodes[j]->adjacentEdges.push_back(ei);
        }
        test = rand() / (double) (RAND_MAX);
        if(test <= edgeProba){
            double w = offset;
            edge* ej = new edge(j,w);
            nodes[0]->adjacentEdges.push_back(ej);
        }
    }

    for(int i = 1; i !=nVertices; i++){
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

bellmanFordResult Graph::bellmanFord(int origin){

    // Store predecessors
    vector<int> pred(nVertices,-1);
    pred[0] = 0;

    // Store distances
    vector<vector<double>> distancesAux(nVertices);
    for(int i = 0; i != nVertices; i++){
        vector<double> distancesLine(nVertices,INFINITY);
        distancesAux[i] = distancesLine;
    }
    distancesAux[0][origin] = 0;
    
    // Main loop
    for(int i = 1; i != nVertices; i++){

        // For each node
        for(vector<node*>::iterator it = nodes.begin(); it != nodes.end(); it++){

            int nodeindex = (*it)->index;

            // Set initial distance 
            distancesAux[i][nodeindex] = distancesAux[i-1][nodeindex];

            // For each edge adjacent to this node
            for(vector<edge*>::iterator it0 = (*it)->adjacentEdges.begin(); it0 != (*it)->adjacentEdges.end(); it0++){

                // If using this edge makes the distance shorter...
                int sourceIndex = (*it0)->source;
                double aux = distancesAux[i-1][sourceIndex] + (*it0)-> weight;
                if(aux < distancesAux[i][nodeindex]){
                    distancesAux[i][nodeindex] = aux;
                    pred[nodeindex] = sourceIndex;
                }
            }
        }
    }

    // Check that no negative cycles
    bool b;

    for(vector<node*>::iterator it = nodes.begin(); it != nodes.end(); it++){

        // For each edge adjacent to this node
        for(vector<edge*>::iterator it0 = (*it)->adjacentEdges.begin(); it0 != (*it)->adjacentEdges.end(); it0++){

            // If using this edge makes the distance shorter...
            int sourceIndex = (*it0)->source;
            double aux = distancesAux[nVertices-1][sourceIndex] + (*it0)-> weight;
            if(aux < distancesAux[nVertices-1][(*it)->index]){
                b = false;
            }
        }
    }

    // Final result
    bellmanFordResult R;
    R.pred = pred;
    R.distances = distancesAux[nVertices-1];
    R.b = b;
    return R;
}

intervalBellmanFordResult Graph::intervalBellmanFord(int origin){

    // Memorize pred
    vector<int> predI(nVertices,-1);
    predI[0] = 0;
    
    // Create structure to memoize distances
    vector<vector<__m128d>> distancesAux(nVertices);
    for(int i = 0; i != nVertices; i++){
        vector<__m128d> distancesLine(nVertices,in2_create(INFINITY,INFINITY));
        distancesAux[i] = distancesLine;
    }
    distancesAux[0][origin] = in2_create(0,0);

    // Create structure to memoize undefined comparisons
    vector<vector<int>> totalc(nVertices+1);
    vector<vector<int>> undefc(nVertices+1);
    for(int i = 0; i != nVertices+1; i++){
        vector<int> undefCompsLine(nVertices,0);
        vector<int> totalCompsLine(nVertices,0);
        undefc[i] = undefCompsLine;
        totalc[i] = totalCompsLine;
    }

    // Generate comparison memorization
    vector<vector<vector<CompResult>>> c(nVertices);
    for(int i = 1; i != nVertices; i++){
        vector<vector<CompResult>> c0(nVertices);
        for(int j = 0; j!= nVertices; j++){
            vector<CompResult> c1(nodes[j]->adjacentEdges.size(),noChange);
            c0[j] = c1;
        }
        c[i] = c0;
    }
    
    // Main loop
    for(int i = 1; i != nVertices; i++){
        // For each node
        for(vector<node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
            // Counting undefined comparisons
            int uc = 0;
            int tc = 0;

            int nodeindex = (*it)->index;
            distancesAux[i][nodeindex] = distancesAux[i-1][nodeindex];

            // For each edge adjacent to this node
            for(int k = 0; k != (*it)->adjacentEdges.size(); k++){
    
                // If using this edge makes the distance shorter...
                int sourceIndex = (*it)->adjacentEdges[k]->source;
                __m128d aux = in2_add_double(distancesAux[i-1][sourceIndex],(*it)->adjacentEdges[k]->weight);
                /*cout << "node: " << nodeindex<< endl;
                cout << "source: " << sourceIndex << endl;
                print(aux);
                cout << endl;
                print(distancesAux[i][nodeindex]);
                cout << endl;*/
                boolean b = inferior(distancesAux[i][nodeindex],aux);
                if(b == False){
                    distancesAux[i][nodeindex] = aux;
                    predI[nodeindex] = sourceIndex;
                    c[i][nodeindex][k] = newOptimum;
                }else if (b == Undefined){
                    c[i][nodeindex][k] = undefined;
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
                    uc += 1;
                    // What do we do for predecessors ? 
                }
                tc += 1;
            }

            totalc[i][nodeindex] = tc;
            undefc[i][nodeindex] = uc;
        }
    }
    
    boolean result = True;
    // Check that no negative cycles
    for(vector<node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        int uc = 0;
        int tc = 0;
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
                uc += 1;
                result = Undefined;
            }
            tc += 1;
        }
        totalc[nVertices][nodeindex] = tc;
        undefc[nVertices][nodeindex] = uc;
    }

    // Final result
    intervalBellmanFordResult R;
    R.b = result;
    R.distances =distancesAux[nVertices-1];
    R.pred = predI;
    R.totalc = totalc;
    R.undefinedc = undefc;
    R.memo = c;
    return R;
}

mpfrBellmanFordResult Graph::lazyMpfrBellmanFord(int origin, vector<vector<vector<CompResult>>> c){
   
    // Memorize pred
    vector<int> predM(nVertices,-1);
    predM[0] = 0;

    // Create structure to memoize distances
    vector<vector<mpreal>> distancesAux(nVertices);
    for(int i = 0; i != nVertices; i++){
        
        distancesAux[i] = vector<mpreal>(nVertices);
        for(int j = 0; j!= nVertices; j++){
            distancesAux[i][j] = mpfr::const_infinity(1,1000);
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
            for(int k = 0; k != (*it)->adjacentEdges.size(); k++){
    
                // Take the result of the lazy comparison
                CompResult cr = c[i][nodeindex][k];
    
                if(cr == undefined){
                    // If using this edge makes the distance shorter...
                    int sourceIndex = (*it)->adjacentEdges[k]->source;
                    mpreal aux = distancesAux[i-1][sourceIndex]+(*it)->adjacentEdges[k]->weight;
                    
                    

                    if(distancesAux[i][nodeindex]>aux){
                        distancesAux[i][nodeindex] = aux;
                        predM[nodeindex] = sourceIndex;
                    }
                }
                else if(cr == newOptimum){
                    int sourceIndex = (*it)->adjacentEdges[k]->source;
                    mpreal aux = distancesAux[i-1][sourceIndex]+(*it)->adjacentEdges[k]->weight;
                    distancesAux[i][nodeindex] = aux;
                    predM[nodeindex] = sourceIndex;
                }
            }
        }
    }

    // Check that no negative cycles
    bool b = true;
    for(vector<node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        // For each edge adjacent to this node
        for(vector<edge*>::iterator it0 = (*it)->adjacentEdges.begin(); it0 != (*it)->adjacentEdges.end(); it0++){

            int sourceIndex = (*it0)->source;
            // If using this edge makes the distance shorter...
            mpreal aux = distancesAux[nVertices-1][sourceIndex]+(*it0)-> weight;
            if(distancesAux[nVertices-1][(*it)->index]>aux){
                b = false;
            }
        }
    }
    
    // Construct result
    mpfrBellmanFordResult R;
    R.pred = predM;
    R.distances = vector<double>(nVertices);
    for(unsigned int i = 0; i != distancesAux[nVertices-1].size(); i++){
        R.distances[i] = (double)distancesAux[nVertices-1][i];
    }
    return R;
    
}

mpfrBellmanFordResult Graph::mpfrBellmanFord(int origin){
   
    // Memorize pred
    vector<int> predM(nVertices,-1);
    predM[0] = 0;

    // Create structure to memoize distances
    vector<vector<mpreal>> distancesAux(nVertices);
    for(int i = 0; i != nVertices; i++){
        
        distancesAux[i] = vector<mpreal>(nVertices);
        for(int j = 0; j!= nVertices; j++){
            distancesAux[i][j] = mpfr::const_infinity(1,1000);
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
            for(int k = 0; k != (*it)->adjacentEdges.size(); k++){
    
    
                // If using this edge makes the distance shorter...
                int sourceIndex = (*it)->adjacentEdges[k]->source;
                mpreal aux = distancesAux[i-1][sourceIndex]+(*it)->adjacentEdges[k]->weight;
                
                

                if(distancesAux[i][nodeindex]>aux){
                    distancesAux[i][nodeindex] = aux;
                    predM[nodeindex] = sourceIndex;
                }
            }
        }
    }

    // Check that no negative cycles
    bool b = true;
    for(vector<node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        // For each edge adjacent to this node
        for(vector<edge*>::iterator it0 = (*it)->adjacentEdges.begin(); it0 != (*it)->adjacentEdges.end(); it0++){

            int sourceIndex = (*it0)->source;
            // If using this edge makes the distance shorter...
            mpreal aux = distancesAux[nVertices-1][sourceIndex]+(*it0)-> weight;
            if(distancesAux[nVertices-1][(*it)->index]>aux){
                b = false;
            }
        }
    }
    
    // Construct result
    mpfrBellmanFordResult R;
    R.pred = predM;
    R.distances = vector<double>(nVertices);
    for(unsigned int i = 0; i != distancesAux[nVertices-1].size(); i++){
        R.distances[i] = (double)distancesAux[nVertices-1][i];
    }
    return R;
    
}
