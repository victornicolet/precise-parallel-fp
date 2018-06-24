/* Main file for comparison of mps methods.
 * Author: Raphael Dang-Nhu.
 * Date: 27/04/2018 */


#include <iostream>
#include <fstream>
#include <time.h>

#include "superaccumulator.hpp"
#include "summation.hpp"
#include "lazy_mps_implementations.hpp"
#include "parallel_mps.hpp"
#include "2_lazy_mps_implementations.hpp"
#include "common.hpp"
#include "pfpdefs.hpp"
#include "tbb/tbb.h"
#include "sequential_mps.hpp"
#include "sequential_mps_alt.hpp"
#include "sequential_mts.hpp"
#include "debug.hpp"

using namespace tbb;
using namespace std;


// This function helps debug a mps implementation 
void debug_test(){
    
    // Test interval arithmetic
    double brray[] = {1.75,-10.,20.,-3.,1.};
    task_scheduler_init init(1);
    parallel_mps_mpfr_lazy_2(brray,5,2);

    double brray4[] = {1.76,-10.,20.,-3.,1.};
    parallel_mps_mpfr_lazy_2(brray4,5,2);
    
    double brray2[] = {1.,pow(2,-50),-1.};
    parallel_mps_mpfr_lazy_2(brray2,3,2);

    double brray3[] = {1.,pow(2,-53),-1.};
    parallel_mps_mpfr_lazy_2(brray3,3,2);

    // Test what happens in case of uncertain comparison
    double crray[] = {1.,pow(2,-53),-1,1.,pow(2,-54)};
    parallel_mps_mpfr_lazy_2(crray,5,2);

    // Test big arrays
    srand(time(NULL));
    for(int i = 2; i < 5; i++){
        int size = pow(10,i);
        double* drray = new double[size];
        init_fpuniform(size, drray, 2000, 1000);
        // Randomly change signs
        for(int j = 0; j < size ; j++){
             drray[j] = (rand() % 2) ? drray[j] : -drray[j];
        }
        for(int k = 0; k < 5; k++){ 
            parallel_mps_mpfr_lazy_2(drray,size,100);
        }
        delete[] drray;
    }
}

/* This function compares the runtime of mps with superaccumulators or mpfr, and its lazy implementation with interval arithmetic */
void runtime_comparison_parallel(){
    
    // Variables declaration and initialisation 
    double start;
    int size = pow(10,7);
    int N = 2;

    // for each dynamic range
    vector<int> dynRanges  {6,10,20,30,60,100,300,600,1000,2000};
    int s = dynRanges.size();
    
    // Store results to plot
    fstream results;
    results.open("Plots/plot1.csv", ofstream::out | ofstream::trunc);
    vector<double> x(s),r0(s),r1(s),r2(s),r3(s),r4(s),r5(s),r6(s);

    // Random seed
    srand(time(NULL));
    
    for(int r = 0; r < dynRanges.size(); r++){

        // initialization of means
        double mean_float = 0, mean_parallel_float = 0, mean_superacc = 0, mean_superacc_lazy = 0, mean_mpfr = 0, mean_lazy_mpfr = 0, mean_lazy_mpfr_2 = 0;
        
        // 1000 trials with N = 10^5
        for(int i = 0; i < N; i++){
            
            // Generating array
            double* drray = new double[size];
            init_fpuniform(size, drray, dynRanges[r], dynRanges[r]/2);

            // Randomly change signs
            for(int j = 0; j < size ; j++){
                 drray[j] = (rand() % 2) ? drray[j] : -drray[j];
            }

            // Variables for result
            double sum, mps;
            int pos;
            
            double time_float = 0.0;
            PFP_TIME(sequential_mps_double(drray,size,&sum,&mps,&pos),start,time_float);
            double time_parallel_float = 0.0;
            PFP_TIME(parallel_mps_float(drray,size),start,time_parallel_float);
            double time_superacc = 0.0;
            //PFP_TIME(parallel_mps_superacc(drray,size),start,time_superacc);
            double time_superacc_lazy = 0.0;
            //PFP_TIME(parallel_mps_superacc_lazy(drray,size),start,time_superacc_lazy);
            double time_mpfr = 0.0;
            //PFP_TIME(parallel_mps_mpfr(drray,size),start,time_mpfr);
            double time_mpfr_lazy = 0.0;
            PFP_TIME(parallel_mps_mpfr_lazy(drray,size),start,time_mpfr_lazy);
            double time_mpfr_lazy_2 = 0.0;
            PFP_TIME(parallel_mps_mpfr_lazy_2(drray,size,10000),start,time_mpfr_lazy_2);
        
            mean_float += time_float;
            mean_parallel_float += time_parallel_float;
            mean_superacc += time_superacc;
            mean_superacc_lazy += time_superacc_lazy;
            mean_mpfr += time_mpfr;
            mean_lazy_mpfr += time_mpfr_lazy;
            mean_lazy_mpfr_2 += time_mpfr_lazy_2;
            
            delete[] drray;
        }
        // Finalize mean computation
        mean_float = mean_float/N;
        mean_parallel_float = mean_parallel_float/N;
        mean_superacc = mean_superacc/N;
        mean_superacc_lazy = mean_superacc_lazy/N;
        mean_mpfr = mean_mpfr/N;
        mean_lazy_mpfr = mean_lazy_mpfr/N;
        mean_lazy_mpfr_2 = mean_lazy_mpfr_2/N;
        
        x[r]= dynRanges[r];
        r0[r]= mean_float/mean_float;
        r1[r] = mean_parallel_float/mean_float;
        r2[r]= mean_superacc/mean_float;
        r3[r]= mean_superacc_lazy/mean_float;
        r4[r]= mean_mpfr/mean_float;
        r5[r]= mean_lazy_mpfr/mean_float;
        r6[r]= mean_lazy_mpfr_2/mean_float;

        // Writing results to a file
        results << to_string(x[r]) << "," << to_string(r0[r]) << "," << to_string(r1[r]) << "," << to_string(r2[r]) << "," << to_string(r3[r])<< "," << to_string(r4[r]) << "," << to_string(r5[r]) << "," << to_string(r6[r]) << endl;;
        
        
    }

    results.close();
}

/* This function compares the runtime of the different sequential mps implementations */
void runtime_comparison_sequential(){
    
    // Variables declaration and initialisation 
    double start;
    int size = pow(10,6);
    int N = 10;

    // for each dynamic range
    vector<int> dynRanges  {300,600,900,1200,1500,1800};
    int s = dynRanges.size();
    
    // Store results to plot
    fstream results;
    results.open("Plots/lazympsseq.csv", ofstream::out | ofstream::trunc);
    vector<double> x(s),r0(s),r1(s),r2(s),r3(s),r4(s);

    // Random seed
    srand(time(NULL));
   
    for(int r = 0; r < dynRanges.size(); r++){

        // initialization of means
        double mean_double = 0.,mean_sum_superacc = 0., mean_superacc = 0., mean_interval = 0.,mean_reverse_mps = 0., mean_reverse_pos = 0., mean_lazy_superacc_mps = 0., mean_lazy_superacc_pos = 0., mean_lazy_mpfr_mps = 0., mean_lazy_mpfr_pos = 0.; 
        

        for(int i = 0; i < N; i++){

            if(PRINT){
                cout << endl << endl <<"***********************************" << endl << "New Input Array" << endl;
            }
            
            // Generating array
            double* drray = new double[size];
            init_fpuniform(size, drray, dynRanges[r], dynRanges[r]/2);

            // Randomly change signs
            for(int j = 0; j < size ; j++){
                 drray[j] = (rand() % 2) ? drray[j] : -drray[j];
            }
            
            // Declare result variables
            double sum;
            double mps;
            int pos;
            memo* da;
            
            double time_double = 0.0;
            PFP_TIME(sequential_mps_double(drray,size,&sum,&mps,&pos),start,time_double);
            double time_superacc = 0.0;
            PFP_TIME(sequential_mps_superacc(drray,size,&sum,&mps,&pos),start,time_superacc);

            double time_sum_superacc = 0.0;
            PFP_TIME(sequential_summation_superacc(drray,size,&sum),start,time_sum_superacc);

            /* Lazy computation, mps superacc */
            double time1 = 0.0;
            PFP_TIME(sequential_mps_interval_memorized(drray,size,&sum,&mps,&pos,&da),start,time1);
            double time2 = 0.0;
            PFP_TIME(sequential_mps_iterate_reverse_mps(drray,size,&sum,&mps,&pos,&da),start,time2);
            double time3 = 0.0;
            PFP_TIME(sequential_mps_lazy_superacc(drray,size,&sum,&mps,&pos,&da),start,time3);
            /* Lazy computation, pos superacc */
            double time4 = 0.0;
            PFP_TIME(sequential_mps_interval_memorized(drray,size,&sum,&mps,&pos,&da),start,time4);
            double time5 = 0.0;
            PFP_TIME(sequential_mps_iterate_reverse_pos(drray,size,&sum,&mps,&pos,&da),start,time5);
            double time6 = 0.0;
            PFP_TIME(sequential_mps_lazy_superacc(drray,size,&sum,&mps,&pos,&da),start,time6);
            
            /* Lazy computation, mps mpfr */
            double time7 = 0.0;
            PFP_TIME(sequential_mps_interval_memorized(drray,size,&sum,&mps,&pos,&da),start,time7);
            double time8 = 0.0;
            PFP_TIME(sequential_mps_iterate_reverse_mps(drray,size,&sum,&mps,&pos,&da),start,time8);
            double time9 = 0.0;
            PFP_TIME(sequential_mps_lazy_mpfr(drray,size,&sum,&mps,&pos,&da),start,time9);
            /* Lazy computation, pos mpfr */
            double time10 = 0.0;
            PFP_TIME(sequential_mps_interval_memorized(drray,size,&sum,&mps,&pos,&da),start,time10);
            double time11 = 0.0;
            PFP_TIME(sequential_mps_iterate_reverse_pos(drray,size,&sum,&mps,&pos,&da),start,time11);
            double time12 = 0.0;
            PFP_TIME(sequential_mps_lazy_mpfr(drray,size,&sum,&mps,&pos,&da),start,time12);

            mean_double += time_double;
            mean_superacc += time_superacc;
            mean_sum_superacc += time_sum_superacc;
            mean_interval += (time1+time4+time7+time10)/4;
            mean_reverse_mps += (time2+time8)/2;
            mean_reverse_pos += (time5+time11)/2;
            mean_lazy_superacc_mps += time3;
            mean_lazy_superacc_pos += time6;
            mean_lazy_mpfr_mps += time9;
            mean_lazy_mpfr_pos += time12;
            
            delete[] drray;
        }
        // Finalize mean computation
        mean_double = mean_double / N;
        mean_sum_superacc = mean_sum_superacc / N;
        mean_superacc = mean_superacc / N;
        mean_interval = mean_interval/N;
        mean_reverse_mps = mean_reverse_mps / N;
        mean_reverse_pos = mean_reverse_pos / N;
        mean_lazy_superacc_mps = mean_lazy_superacc_mps / N;
        mean_lazy_superacc_pos = mean_lazy_superacc_pos / N;
        mean_lazy_mpfr_mps = mean_lazy_mpfr_mps / N;
        mean_lazy_mpfr_pos = mean_lazy_mpfr_pos / N;
        
        x[r]= dynRanges[r];
        r0[r]= mean_double/mean_double;
        r1[r]= mean_sum_superacc/mean_double;
        r2[r]= mean_superacc/mean_double;
        r3[r]= mean_interval/mean_double;
        r4[r]= mean_reverse_mps/mean_double;

        // Writing results to a file
        results << to_string(x[r]) << "," << to_string(r0[r]) << "," << to_string(r1[r]) << "," << to_string(r2[r]) << "," << to_string(r3[r])<< "," << to_string(r4[r]) << endl;
        
        
    }

    results.close();
}

/* This function compares the runtime of the different sequential mts implementations */
void runtime_comparison_sequential_mts(){
    
    // Variables declaration and initialisation 
    double start;
    long size = pow(10,6);
    int N = 5;

    // for each dynamic range
    vector<int> dynRanges  {300,600,900,1200,1500,1800};
    int s = dynRanges.size();
    
    // Store results to plot
    fstream results;
    results.open("Plots/plot3.csv", ofstream::out | ofstream::trunc);
    vector<double> x(s),r0(s),r1(s),r2(s),r3(s),r4(s);

    // Random seed
    srand(time(NULL));
    
    for(int r = 0; r < dynRanges.size(); r++){

        // initialization of means
        double mean_double = 0.,mean_sum_superacc = 0., mean_superacc = 0., mean_lazy = 0.,mean_lazy_opt = 0.; 
        

        for(int i = 0; i < N; i++){
            
            // Generating array
            double* drray = new double[size];
            init_fpuniform(size, drray, dynRanges[r], dynRanges[r]/2);

            // Randomly change signs
            for(int j = 0; j < size ; j++){
                 drray[j] = (rand() % 2) ? drray[j] : -drray[j];
            }
            
            // Declare result variables
            double sum;
            double mts;
            long pos;
            
            double time_double = 0.0;
            PFP_TIME(sequential_mts_double(drray,size,&mts,&pos),start,time_double);
            double time_superacc = 0.0;
            PFP_TIME(sequential_mts_superacc(drray,size,&mts,&pos),start,time_superacc);
            double time_lazy = 0.0;
            PFP_TIME(sequential_mts_lazy(drray,size,&mts,&pos,0),start,time_lazy);
            double time_lazy_opt = 0.0;
            PFP_TIME(sequential_mts_lazy(drray,size,&mts,&pos,1),start,time_lazy_opt);
            double time_sum_superacc = 0.0;
            PFP_TIME(sequential_summation_superacc(drray,size,&sum),start,time_sum_superacc);
        
            mean_double += time_double;
            mean_sum_superacc += time_sum_superacc;
            mean_superacc += time_superacc;
            mean_lazy += time_lazy;
            mean_lazy_opt += time_lazy_opt;
            
            delete[] drray;
        }
        // Finalize mean computation
        mean_double = mean_double / N;
        mean_sum_superacc = mean_sum_superacc /N;
        mean_superacc = mean_superacc /N;
        mean_lazy = mean_lazy /N;
        mean_lazy_opt = mean_lazy_opt /N;
        
        x[r]= dynRanges[r];
        r0[r]= mean_double/mean_double;
        r1[r]= mean_sum_superacc/mean_double;
        r2[r]= mean_superacc/mean_double;
        r3[r]= mean_lazy/mean_double;
        r4[r]= mean_lazy_opt/mean_double;

        // Writing results to a file
        results << to_string(x[r]) << "," << to_string(r0[r]) << "," << to_string(r1[r]) << "," << to_string(r2[r]) << "," << to_string(r3[r])<< "," << to_string(r4[r]) << endl;
        
        
    }

    results.close();
}
int main(int argc, char** argv){
    if(argc >= 1){
        int a = atoi(argv[1]);
        switch (a){
            case 0 :
                debug_test();
                break;
            case 1 : 
                runtime_comparison_parallel();
                break;
            case 2 :
                runtime_comparison_sequential();
                break;
            case 3 :
                runtime_comparison_sequential_mts();
                break;
        }
    }
}
