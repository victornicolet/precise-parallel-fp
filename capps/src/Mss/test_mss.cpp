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
    long size = pow(10,7);
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
    // Variables declaration and initialisation 
    double start;
    int size = pow(10,8);
    int N = 5;

    
    // Store results to plot
    fstream results;
    results.open("Plots/mss_hybrid.csv", ofstream::out | ofstream::trunc);

    // Random seed
    srand(time(NULL));
    // Generating array
    double* drray = new double[size];
    init_fpuniform(size, drray, 100, 50);

    // Randomly change signs
    for(int j = 0; j < size ; j++){
         drray[j] = (rand() % 2) ? drray[j] : -drray[j];
    }
    
    // initialization of means
    double mean_hybrid = 0.;

    for(int i = 0; i < N; i++){
        
        // Declare result variables
        double time_hybrid = 0.0;
        PFP_TIME(parallel_mss_hybrid(drray,size),start,time_hybrid);
   
        mean_hybrid += time_hybrid;

    }
    mean_hybrid = mean_hybrid / N;

    // initialization of means
    double mean_tbb = 0.;

    for(int i = 0; i < N; i++){
        
        // Declare result variables
        double time_tbb = 0.0;
        PFP_TIME(parallel_mss_double(drray,size),start,time_tbb);
   
        mean_tbb += time_tbb;
        
    }
    mean_tbb = mean_tbb / N;

    
        cout << endl << "Performance ratio: " << mean_hybrid / mean_tbb << endl;
        // Writing results to a file
        results << to_string(mean_tbb) << "," << to_string(mean_hybrid) << "," << endl;

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
        }
    }
}
