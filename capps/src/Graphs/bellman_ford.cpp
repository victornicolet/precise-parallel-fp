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
            double w = randDouble(emin,emax,negratio);
            // Create edges
            edge* ei = new edge(0,w);
            nodes[j]->adjacentEdges.push_back(ei);
        }
        test = rand() / (double) (RAND_MAX);
        if(test <= edgeProba){
            double w = randDouble(emin,emax,negratio);
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
        for(int j = 0; j!= nVertices; j++){

            // Set initial distance 
            distancesAux[i][j] = distancesAux[i-1][j];

            // For each edge adjacent to this node
            for(int k = 0; k != nodes[j]->adjacentEdges.size(); k++){

                int sourceIndex = nodes[j]->adjacentEdges[k]->source;

                double aux = distancesAux[i-1][sourceIndex] + nodes[j]->adjacentEdges[k]-> weight;
                /*
                cout << "Node Index: " << j << endl;
                cout << "Pred: " << pred[j] << endl;
                cout << "Edge source: " << sourceIndex << endl;
                cout << "Temp distance: " << distancesAux[i][j] << endl;
                cout << "Aux: " << aux << endl << endl;
                */
                if(aux < distancesAux[i][j]){
                    distancesAux[i][j] = aux;
                    pred[j] = sourceIndex;
                }
            }
        }
    }

    // Final result
    bellmanFordResult R;
    R.pred = pred;
    R.distances = distancesAux[nVertices-1];
    R.b = True;
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

    // Generate comparison memorization
    compType c(nVertices);
    for(int i = 1; i != nVertices; i++){
        vector<vector<memo>> c0(nVertices);
        for(int j = 0; j!= nVertices; j++){
            vector<memo> c1(nodes[j]->adjacentEdges.size());
            c0[j] = c1;
        }
        c[i] = c0;
    }
    
    // Main loop
    for(int i = 1; i != nVertices; i++){
        // For each node
        for(int j = 0; j != nodes.size(); j++){

            distancesAux[i][j] = distancesAux[i-1][j];

            // For each edge adjacent to this node
            for(int k = 0; k != nodes[j]->adjacentEdges.size(); k++){
    
                // If using this edge makes the distance shorter...
                int sourceIndex = nodes[j]->adjacentEdges[k]->source;
                __m128d aux = in2_add_double(distancesAux[i-1][sourceIndex],nodes[j]->adjacentEdges[k]->weight);
                /*
                cout << "Node Index: " << j << endl;
                cout << "Pred: " << predI[j] << endl;
                cout << "Edge source: " << sourceIndex << endl;
                cout << "Temp distance: ";
                print(distancesAux[i][j]);
                cout << endl << "Aux: ";
                print(aux);
                cout << endl << endl;
                */
                boolean b = inferior(distancesAux[i][j],aux);
                if(b == False){
                    distancesAux[i][j] = aux;
                    predI[j] = sourceIndex;
                }else if (b == Undefined){
                    distancesAux[i][j] = in2_min(aux,distancesAux[i][j]);
                }
                c[i][j][k].useful2 = b;
            }
        }
    }
   

    // Final result
    intervalBellmanFordResult R;
    R.b = False;
    R.distances =distancesAux[nVertices-1];
    R.pred = predI;
    R.mem = c;
    return R;
}

compType Graph::reverseBellmanFord(int origin, compType T){

    // Memorize pred
    vector<bool> predI(nVertices,false);
    predI[nVertices-1] = true;
    
    // Create structure to memoize distances
    vector<vector<bool>> distancesAux(nVertices);
    for(int i = 0; i != nVertices; i++){
        vector<bool> distancesLine(nVertices,false);
        distancesAux[i] = distancesLine;
    }

    bool aux = false;

    // Main loop
    for(int i = nVertices-1; i > 0; i--){
        // For each node
        for(int j = nodes.size()-1; j >= 0; j--){
            
            // For each edge adjacent to this node
            for(int k = nodes[j]->adjacentEdges.size() - 1; k >= 0; k--){

                int sourceIndex = nodes[j]->adjacentEdges.size();
    
                memo b = T[i][j][k];

                // Comparison
                if(b.useful2 == False){
                    if(!(predI[j] || distancesAux[i][j])){
                        T[i][j][k].useful2 = Useless;
                    }
                    if(predI[j]){
                        predI[j] = false;
                        predI[sourceIndex] = true;
                    }
                    if(distancesAux[i][j]){
                        aux = true;
                        distancesAux[i][j] = false;
                    }
                }else if (b.useful2 == True){
                    T[i][j][k].useful2 = Useless;
                }else if (b.useful2 == Undefined){
                    bool t = predI[j] || distancesAux[i][j];
                    if(!t){
                        T[i][j][k].useful2 = Useless;
                    }
                    if(predI[j]){
                        predI[sourceIndex] = true;
                    }
                    if(distancesAux[i][j]){
                        aux = true;
                    }
                    if(t){
                        aux = true;
                    }
                }

                // Computation of aux 
                if(aux){
                    aux = false;
                    distancesAux[i-1][sourceIndex] = true;
                    b.useful1 = true;
                }else{
                    b.useful1 = false;
                }

            }
        // Assignment of previous value
        if(distancesAux[i][j]){
            distancesAux[i-1][j] = true;
            distancesAux[i][j] = false;
        }
        }
    }
    return T;
}

mpfrBellmanFordResult Graph::lazyMpfrBellmanFord(int origin, compType c){
   
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
    
    mpreal aux;
    int sourceIndex;

    // Main loop
    for(int i = 1; i != nVertices; i++){
        // For each node
        for(int j = 0; j!= nVertices; j++){

            distancesAux[i][j] = distancesAux[i-1][j];

            // For each edge adjacent to this node
            for(int k = 0; k != nodes[j]->adjacentEdges.size(); k++){
    
                memo cr = c[i][j][k];
                
                if(cr.useful1){
                    sourceIndex = nodes[j]->adjacentEdges[k]->source;
                    aux = distancesAux[i-1][sourceIndex]+nodes[j]->adjacentEdges[k]->weight;
                }
    
                if(cr.useful2 == False){
                    distancesAux[i][j] = aux;
                    predM[j] = sourceIndex;
                }
                else if(cr.useful2 == Undefined){
                    if(aux < distancesAux[i][j]){
                        distancesAux[i][j] = aux;
                        predM[j] = sourceIndex;
                    }
                }
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
        for(int j = 0; j!= nVertices; j++){

            distancesAux[i][j] = distancesAux[i-1][j];

            // For each edge adjacent to this node
            for(int k = 0; k != nodes[j]->adjacentEdges.size(); k++){
    
                // If using this edge makes the distance shorter...
                int sourceIndex = nodes[j]->adjacentEdges[k]->source;
                mpreal aux = distancesAux[i-1][sourceIndex]+nodes[j]->adjacentEdges[k]->weight;
                
                if(distancesAux[i][j]>aux){
                    distancesAux[i][j] = aux;
                    predM[j] = sourceIndex;
                }
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
