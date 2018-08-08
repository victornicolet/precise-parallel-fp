/* Main file for comparison of steep methods.
 * Author: Raphael Dang-Nhu.
 * Date: 03/08/2018 */


#include <iostream>
#include <fstream>
#include <time.h>

#include "superaccumulator.hpp"
#include "common.hpp"
#include "pfpdefs.hpp"

#include "tbb/tbb.h"

#include "parallel_steep.hpp"
#include "summation.hpp"

using namespace tbb;
using namespace std;

// Runtime comparison for the new steep hybrid lazy computation
void runtime_comparison_parallel_steep_hybrid(){
    // Parameters
    double start;
    int size = pow(10,6);
    int N = 100;
    vector<int> depths = {1,3,5,7,10,11,12,13,14,15};
    vector<double> hybrid(depths.size());
    vector<double> lazy(depths.size());
    
    // Store results to plot
    fstream results;
    results.open("Plots/csv/steep_hybrid.csv", ofstream::out | ofstream::trunc);

    // Generating array
    srand(time(NULL));
    double* drray = new double[size];
    init_fpuniform(size, drray, 200, 100);
    for(int j = 0; j < size ; j++){
         drray[j] = (rand() % 2) ? drray[j] : -drray[j];
    }
     
    // Hybrid Reduction
    for(int it = 0; it != depths.size(); it++){

        double mean_hybrid = 0.;
        for(int i = 0; i < N; i++){
            
            // Declare result variables
            double time_hybrid = 0.0;
            PFP_TIME(parallel_steep_hybrid(drray,size,depths[it]),start,time_hybrid);
            mean_hybrid += time_hybrid;

        }
        mean_hybrid = mean_hybrid / N;
        hybrid[it] = mean_hybrid;
    }

    // Standard reduction
    double mean_tbb = 0.;

    for(int i = 0; i < N; i++){
        
        // Declare result variables
        double time_tbb = 0.0;
        PFP_TIME(parallel_steep_double(drray,size),start,time_tbb);
   
        mean_tbb += time_tbb;
        
    }
    mean_tbb = mean_tbb/N ;

    
    // Lazy hybrid reduction
    for(int it = 0; it != depths.size(); it++){

        double mean_hybrid_interval = 0.;
        for(int i = 0; i < N; i++){
            
            // Declare result variables
            double time_hybrid_interval = 0.0;
            PFP_TIME(parallel_steep_hybrid_interval(drray,size,depths[it]),start,time_hybrid_interval);
            mean_hybrid_interval += time_hybrid_interval;
            
        }
        mean_hybrid_interval = mean_hybrid_interval / N;
        lazy[it] = mean_hybrid_interval;
    }

    // Writing results to a file
    results << to_string(mean_tbb) << endl;
    for(int it = 0; it != depths.size(); it++){
        results << to_string(depths[it]) << "," <<
        to_string(lazy[it]) << "," <<
        to_string(hybrid[it]) << endl;
    }

    results.close();
}

void optimal_depth(){
    // Parameters
    double start;
    int N = 10;
    vector<int> sizes = {3*pow(10,2),pow(10,3),3*pow(10,3),pow(10,4),3*pow(10,4),pow(10,5),3*pow(10,5),pow(10,6),3*pow(10,6),pow(10,7),3*pow(10,7),pow(10,8),3*pow(10,8),pow(10,9)};
    vector<int> ntrials = {100000,30000,10000,3000,1000,1000,1000,1000,1000,300,100,30,10,3};
    
    // Store results to plot
    fstream results;
    results.open("Plots/csv/optimal_depth.csv", ofstream::out | ofstream::trunc);

    for(int a = 0; a != sizes.size(); a++){
        int size = sizes[a];
        cout << "Size: " << size << endl;

        // Generating array
        srand(time(NULL));
        double* drray = new double[size];
        init_fpuniform(size, drray, 200, 100);
        for(int j = 0; j < size ; j++){
             drray[j] = (rand() % 2) ? drray[j] : -drray[j];
        }

        // Standard reduction
        double mean_tbb = 0.;

        for(int i = 0; i < ntrials[a]; i++){
            
            // Declare result variables
            double time_tbb = 0.0;
            PFP_TIME(parallel_steep_double(drray,size),start,time_tbb);
       
            mean_tbb += time_tbb;
            
        }
        mean_tbb = mean_tbb/ntrials[a] ;

        // Hybrid Reduction
        int d = 1;
        double mean_hybrid = 0;
        while(mean_hybrid <= 1.2){

            mean_hybrid = 0.;
            for(int i = 0; i < ntrials[a]; i++){
                
                // Declare result variables
                double time_hybrid = 0.0;
                PFP_TIME(parallel_steep_hybrid(drray,size,d),start,time_hybrid);
                mean_hybrid += time_hybrid;

            }
            mean_hybrid = mean_hybrid / ntrials[a];
            mean_hybrid /= mean_tbb;
            d++;
        }
        d--;

        // Writing results to a file
        results <<
        to_string(size) << "," <<
        to_string(d) << endl;

    }
    results.close();
}

int main(int argc, char** argv){
    if(argc >= 1){
        int long a = atoi(argv[1]);
        switch (a){
            case 0 :
                runtime_comparison_parallel_steep_hybrid();
                break;
            case 1 :
                optimal_depth();
        }
    }
}
