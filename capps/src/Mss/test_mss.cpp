/* Main file for comparison of mps methods.
 * Author: Raphael Dang-Nhu.
 * Date: 27/04/2018 */


#include <iostream>
#include <fstream>
#include <time.h>

#include "superaccumulator.hpp"
#include "common.hpp"
#include "pfpdefs.hpp"

#include "tbb/tbb.h"

#include "sequential_mss.hpp"
#include "parallel_mss.hpp"
#include "summation.hpp"

using namespace tbb;
using namespace std;

/* This function compares the runtime of the different sequential mss implementations */
void runtime_comparison_sequential_mss(){
    
    // Variables declaration and initialisation 
    double start;
    long size = pow(10,2);
    long N = 1;

    // for each dynamic range
    vector<long> dynRanges  {100,400,700,1000,1300,1600,1900};
    unsigned long s = dynRanges.size();
    
    // Store results to plot
    fstream results;
    results.open("Plots/sequential_mss.csv", ofstream::out | ofstream::trunc);
    vector<double> x(s),r0(s),r1(s),r2(s),r3(s),r4(s);

    // Random seed
    srand(time(NULL));
    
    for(unsigned long r = 0; r < dynRanges.size(); r++){

        // initialization of means
        double mean_double = 0.,mean_sum_superacc = 0., mean_superacc = 0., mean_lazy = 0.; 
        

        for(long i = 0; i < N; i++){
            
            // Generating array
            double* drray = new double[size];
            init_fpuniform(size, drray, dynRanges[r], dynRanges[r]/2);

            // Randomly change signs
            for(long j = 0; j < size ; j++){
                 drray[j] = (rand() % 2) ? drray[j] : -drray[j];
            }
            
            // Declare result variables
            double sum;
            double mss;
            long pos,pos0;
            
            double time_double = 0.0;
            PFP_TIME(sequential_mss_double(drray,size),start,time_double);
            double time_superacc = 0.0;
            PFP_TIME(sequential_mss_superacc(drray,size,&mss,&pos,&pos0),start,time_superacc);
            double time_lazy = 0.0;
            PFP_TIME(sequential_mss_lazy(drray,size,&mss,&pos,&pos0),start,time_lazy);
            double time_sum_superacc = 0.0;
            PFP_TIME(sequential_summation_superacc(drray,size,&sum),start,time_sum_superacc);
        
            mean_double += time_double;
            mean_sum_superacc += time_sum_superacc;
            mean_superacc += time_superacc;
            mean_lazy += time_lazy;
            
            delete[] drray;
        }
        // Finalize mean computation
        mean_double = mean_double / N;
        mean_sum_superacc = mean_sum_superacc /N;
        mean_superacc = mean_superacc /N;
        mean_lazy = mean_lazy /N;
        
        x[r]= dynRanges[r];
        r0[r]= mean_double/mean_double;
        r1[r]= mean_sum_superacc/mean_double;
        r2[r]= mean_superacc/mean_double;
        r3[r]= mean_lazy/mean_double;

        // Writing results to a file
        results << to_string(x[r]) << "," << to_string(r0[r]) << "," << to_string(r1[r]) << "," << to_string(r2[r]) << "," << to_string(r3[r])<< endl;
    }
    results.close();
}

/* This function compares the runtime of the different parallel mss implementations */
void runtime_comparison_parallel_mss(){
    
    // Variables declaration and initialisation 
    double start;
    long size = pow(10,5);
    long N = 1;

    // for each dynamic range
    vector<long> dynRanges  {1000,1900};
    unsigned long s = dynRanges.size();
    
    // Store results to plot
    fstream results;
    results.open("Plots/parallel_mss.csv", ofstream::out | ofstream::trunc);
    vector<double> x(s),r0(s),r1(s),r2(s),r3(s),r4(s);

    // Random seed
    srand(time(NULL));
    
    for(unsigned long r = 0; r < dynRanges.size(); r++){

        // initialization of means
        double mean_double = 0., mean_parallel_double = 0., mean_dynamic_lazy = 0.; 
        

        for(long i = 0; i < N; i++){
            
            // Generating array
            double* drray = new double[size];
            init_fpuniform(size, drray, dynRanges[r], dynRanges[r]/2);

            // Randomly change signs
            for(long j = 0; j < size ; j++){
                 drray[j] = (rand() % 2) ? drray[j] : -drray[j];
            }
            
            // Declare result variables
            double  mss;
            long pos, pos0;
           
            // Check result
            //sequential_mss_superacc(drray,size,&mss,&pos,&pos0);
            double time_seq_double = 0.0;
            PFP_TIME(sequential_mss_double(drray,size),start,time_seq_double);
            double time_double = 0.0;
            PFP_TIME(parallel_mss_double(drray,size),start,time_double);
            double time_interval = 0.0;
            PFP_TIME(parallel_mss_interval(drray,size),start,time_interval);
        
            mean_double += time_seq_double;
            mean_parallel_double += time_double;
            mean_dynamic_lazy += time_interval;
            
            delete[] drray;
        }
        // Finalize mean computation
        mean_double = mean_double / N;
        mean_parallel_double = mean_parallel_double / N;
        mean_dynamic_lazy = mean_dynamic_lazy / N;
        
        x[r]= dynRanges[r];
        r0[r]= mean_double/mean_double;
        r1[r] = mean_parallel_double / mean_double;
        r2[r] = mean_dynamic_lazy / mean_double;

        // Writing results to a file
        results << to_string(x[r]) << "," << to_string(r0[r]) << "," << to_string(r1[r]) << "," << to_string(r2[r]) << endl;
    }
    results.close();
}

// Runtime comparison for the new mss hybrid lazy computation
void runtime_comparison_parallel_mss_hybrid(){
    // Parameters
    double start;
    int size = pow(10,6);
    int N = 1;
    vector<int> depths = {1,5,11,12,13,14,15};
    vector<double> hybrid(depths.size());
    vector<double> lazy(depths.size());
    
    // Store results to plot
    fstream results;
    results.open("Plots/csv/mss_hybrid.csv", ofstream::out | ofstream::trunc);

    // Generating array
    srand(time(NULL));
    double* drray = new double[size];
    init_fpuniform(size, drray, 100, 50);
    for(int j = 0; j < size ; j++){
         drray[j] = (rand() % 2) ? drray[j] : -drray[j];
    }
     
    // Hybrid Reduction
    for(int it = 0; it != depths.size(); it++){

        double mean_hybrid = 0.;
        for(int i = 0; i < N; i++){
            
            // Declare result variables
            double time_hybrid = 0.0;
            PFP_TIME(parallel_mss_hybrid(drray,size,depths[it]),start,time_hybrid);
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
        PFP_TIME(parallel_mss_double(drray,size),start,time_tbb);
   
        mean_tbb += time_tbb;
        
    }
    mean_tbb = mean_tbb / N;

    
    // Lazy hybrid reduction
    for(int it = 0; it != depths.size(); it++){

        double mean_hybrid_interval = 0.;
        for(int i = 0; i < N; i++){
            
            // Declare result variables
            double time_hybrid_interval = 0.0;
            PFP_TIME(parallel_mss_hybrid_interval(drray,size,depths[it]),start,time_hybrid_interval);
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

void runtime_comparison_parallel_mss_hybrid_final(){
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
    results.open("Plots/csv/mss_hybrid_final_with_recomp.csv", ofstream::out | ofstream::trunc);

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
            PFP_TIME(parallel_mss_hybrid_interval(drray,size,1),start,time_hybrid);
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
            PFP_TIME(parallel_mss_double(drray,size),start,time_tbb);
       
            mean_tbb += time_tbb;
            
        }
        mean_tbb = mean_tbb/N ;
        mean_tbb = mean_tbb/size ;


        cout << "Starting interval 2: " << size << endl;
        // Lazy hybrid reduction
        double mean_hybrid_interval = 0.;
        double mean_hybrid_total = 0.;
        for(int k = 0; k < N; k++){
            
            // Declare result variables
            double time_hybrid_interval;
            double time_hybrid_total;
            parallel_mss_hybrid_lazy(drray,size,depths[i],time_hybrid_interval,time_hybrid_total);
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
        PFP_TIME(sequential_mss_double(drray,size),start,time_seq);
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

int main(int argc, char** argv){
    if(argc >= 1){
        int long a = atoi(argv[1]);
        switch (a){
            case 0 :
                runtime_comparison_sequential_mss();
                break;
            case 1 :
                runtime_comparison_parallel_mss();
                break;
            case 2 :
                runtime_comparison_parallel_mss_hybrid();
                break;
            case 3 :
                runtime_comparison_parallel_mss_hybrid_final();
                break;
        }
    }
}
