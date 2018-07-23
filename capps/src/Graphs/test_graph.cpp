/* Main function for the Bellman-Ford tests.
 * Author: Raphaël Dang-Nhu.
 * Date: 15/05 */

#include <fstream>
#include <iostream>
#include <cmath>

#include "debug.hpp"
#include "pfpdefs.hpp"
#include "bellman_ford.hpp"
#include "printFunctions.hpp"

using namespace std;
                
void test(){
    int n = 30;
    Graph g(n,0.5,-1000,1000,1,1);    

    if(PRINT){
        g.printGraph();
    }
    
    // Test with doubles
    bellmanFordResult R = g.bellmanFord(0);

    if(PRINT){
        if(R.b){
            cout << "No negative cycles" << endl;
        }
        else{
            cout << "Negative cycle" << endl;
        }
    
        cout << "Distances: " << endl;
        printVector(R.distances);
        cout << "Predecessors: " << endl;
        printVector(R.pred);
    }

    // Test with interval arithmetic
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);

    intervalBellmanFordResult R0 = g.intervalBellmanFord(0); 
    if(PRINT){
        if(R0.b == True){
                cout << "No negative cycles" << endl;
        }
        else if ( R0.b == False) {
            cout << "Negative cycle" << endl;
        }
        else if (R0.b == Undefined){
            cout << "Uncertain presence of negative cycles" << endl;
        }
        cout << "Distances: " << endl;
        printVectorInt(R0.distances);
        cout << "Predecessors: " << endl;
        printVector(R0.pred);
        cout << "Undefined comparisons" << endl;
        printComps(R0.undefinedc,R0.totalc);
    }
    _MM_SET_ROUNDING_MODE(0);

    // Test with mpfr
    mpfrBellmanFordResult mR= g.mpfrBellmanFord(0);
    if(PRINT){
        if(mR.b){
            cout << "No negative cycles" << endl;
        }
        else{
            cout << "Negative cycle" << endl;
        }
        cout << "Distances: " << endl;
        printVector(mR.distances);
        cout << "Predecessors: " << endl;
        printVector(mR.pred);
    }

    // Check if different result
    int d = 0;
    for(unsigned int i = 0; i != R.pred.size(); i++){
        if(R.pred[i] != mR.pred[i]){
            d++;
        }
    }
    cout << endl << "Number of different results: "<< d <<  endl;
}

void runtimeTest(){
   
    // for each dynamic range
    vector<int> graphSizes  {1200};
    unsigned long s = graphSizes.size();
    
    // Store results to plot
    fstream results;
    results.open("Plots/bellman_ford.csv", ofstream::out | ofstream::trunc);

    for(unsigned long r = 0; r != s; r++){
        
        // Generate a graph
        int nV = graphSizes[r];
        Graph g(nV,0.5,-1000,1000,1,0);    

        // Tests
        double start = 0.0, time_double = 0.0, time_interval = 0.0, time_reverse = 0.0, time_mpfr = 0.0, time_lazy_mpfr = 0.0;

        PFP_TIME(bellmanFordResult R = g.bellmanFord(0),start,time_double);
        PFP_TIME(intervalBellmanFordResult R0 = g.intervalBellmanFord(0),start,time_interval);
        PFP_TIME(g.reverseBellmanFord(0,R0.mem,R0.memsize),start,time_reverse);
        //PFP_TIME(mpfrBellmanFordResult R1 = g.mpfrBellmanFord(0),start,time_mpfr);
        PFP_TIME(mpfrBellmanFordResult R2 = g.lazyMpfrBellmanFord(0,R0.mem),start,time_lazy_mpfr);

        double time_double_aux = time_double / time_double;
        time_interval = time_interval / time_double;
        time_mpfr = time_mpfr / time_double;
        time_lazy_mpfr = time_lazy_mpfr / time_double;
        time_reverse = time_reverse / time_double;

        // Compute mean discrepancy between result
        /*
        double mean1 = 0.0, mean2 = 0.0;
        for(int i = 0; i != nV; i++){
            cout << "Floating-Point pred: " << R.pred[i] << endl;
            cout << "Interval pred: " << R0.pred[i] << endl;
            //cout << "Precise pred: " << R1.pred[i] << endl;
            cout << "Lazy pred: " << R2.pred[i] << endl << endl;
            //mean1 += abs(R1.distances[i] - R.distances[i]);
            //mean2 += abs(R2.distances[i] - R1.distances[i]);
        }
        cout << "Mpfr result - double result: " << mean1 << endl;
        cout << "Lazy mpfr result - Mpfr result: " << mean2 << endl;
        */
        results << to_string(nV) << "," << to_string(time_double_aux) << "," << to_string(time_mpfr) << "," << to_string(time_lazy_mpfr) << "," << to_string(time_interval) << "," << to_string(time_reverse) << "," << to_string(time_interval+time_reverse+time_lazy_mpfr) << endl;
    }
    results.close();

}

int main(int argc, char** argv){
    srand(time(NULL));
    if(argc >= 1){
        int long a = atoi(argv[1]);
        switch (a){
            case 0 :
                test();
                break;
            case 1 :
                runtimeTest();
                break;
        }
    }
}

