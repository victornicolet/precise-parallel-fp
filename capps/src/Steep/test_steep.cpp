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

void debug(){
    // Parameters
    double start;
        int size = 100;

        // Generating array
        srand(time(NULL));
        double* drray = new double[size];
        init_fpuniform(size, drray, 1, 1);
        for(int j = 0; j < size ; j++){
             drray[j] = (rand() % 2) ? drray[j] : -drray[j];
        }
         
        cout << "Starting standard: " << size << endl;
        // Standard reduction

            // Declare result variables
            parallel_steep_double(drray,size);
       
        cout << "Starting interval 2: " << size << endl;
            double time_hybrid_interval;
            double time_hybrid_total;
            parallel_steep_hybrid_lazy(drray,size,3,time_hybrid_interval,time_hybrid_total);
            
}

// Runtime comparison for the new steep hybrid lazy computation
void runtime_comparison_parallel_steep_hybrid(){
    // Parameters
    double start;
    vector<int> sizes = {3*pow(10,2),pow(10,3),3*pow(10,3),pow(10,4),3*pow(10,4),pow(10,5),3*pow(10,5),pow(10,6),3*pow(10,6),pow(10,7),3*pow(10,7),pow(10,8),3*pow(10,8),pow(10,9)};
    vector<int> depths = {5,5,5,5,5,5,6,7,8,10,11,13,14,14};
    vector<int> ntrials = {10000,30000,10000,3000,3000,3000,1000,1000,300,100,30,10,3,1};
    int n = depths.size();
    vector<double> interval(n);
    vector<double> total(n);
    vector<double> intervalref(n);
    
    // Store results to plot
    fstream results;
    results.open("Plots/csv/steep_hybrid.csv", ofstream::out | ofstream::trunc);

    for(int i = 0; i != n; i++){
        int size = sizes[i];
        cout << "Size: " << size << endl;
        int N = ntrials[i];

        // Generating array
        srand(time(NULL));
        double* drray = new double[size];
        init_fpuniform(size, drray, 200, 100);
        for(int j = 0; j < size ; j++){
             drray[j] = (rand() % 2) ? drray[j] : -drray[j];
        }
         
        cout << "Starting hybrid interval: " << size << endl;
        // Hybrid Reduction
        double mean_hybrid = 0.;
        for(int k = 0; k < N; k++){
            
            // Declare result variables
            double time_hybrid = 0.0;
            PFP_TIME(parallel_steep_hybrid_interval(drray,size,1),start,time_hybrid);
            mean_hybrid += time_hybrid;

        }
        mean_hybrid = mean_hybrid / N;
        mean_hybrid = mean_hybrid / size;
        intervalref[i] = mean_hybrid;

        cout << "Starting standard: " << size << endl;
        // Standard reduction
        double mean_tbb = 0.;

        for(int k = 0; k < N; k++){
            
            // Declare result variables
            double time_tbb = 0.0;
            PFP_TIME(parallel_steep_double(drray,size),start,time_tbb);
       
            mean_tbb += time_tbb;
            
        }
        mean_tbb = mean_tbb/N ;
        mean_tbb = mean_tbb/size ;

        cout << "Starting interval 2: " << size << endl;
        double mean_hybrid_interval = 0.;
        double mean_hybrid_total = 0.;
        for(int k = 0; k < N; k++){
            
            // Declare result variables
            double time_hybrid_interval;
            double time_hybrid_total;
            parallel_steep_hybrid_lazy(drray,size,depths[i],time_hybrid_interval,time_hybrid_total);
            mean_hybrid_interval += time_hybrid_interval;
            mean_hybrid_total += time_hybrid_total;
            
        }
        mean_hybrid_interval = mean_hybrid_interval / size;
        mean_hybrid_interval = mean_hybrid_interval / N;
        mean_hybrid_total = mean_hybrid_total / size;
        mean_hybrid_total = mean_hybrid_total / N;
        interval[i] = mean_hybrid_interval;
        total[i] = mean_hybrid_total;

        cout << "Starting sequential: " << size << endl;
        // Lazy hybrid reduction
        double mean_sequential = 0.;
        // Declare result variables
        double time_seq = 0.0;
        double s,c;
        bool b;
        PFP_TIME(sequential_steep(drray,size,s,c,b),start,time_seq);
        mean_sequential += time_seq;
        mean_sequential /= size;

        interval[i] = mean_hybrid_interval;
        results << to_string(sizes[i]) << "," <<
        to_string(mean_tbb) << "," <<
        to_string(interval[i]) << "," <<
        to_string(intervalref[i]) << "," <<
        to_string(mean_sequential) << "," <<
        to_string(total[i]) << endl;
    }
    results.close();
}

void optimal_depth(){
    // Parameters
    double start;
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
                break;
            case 2 :
                debug();
        }
    }
}
