/* Main file for comparison of mps methods.
 * Author: Raphael Dang-Nhu.
 * Date: 27/04/2018 */

#include <iostream>
#include <fstream>

#include "matplotlibcpp.h"
#include "superaccumulator.hpp"
#include "lazy_mps_implementations.hpp"
#include "2_lazy_mps_implementations.hpp"
#include "common.hpp"
#include "pfpdefs.hpp"

using namespace std;

namespace plt = matplotlibcpp;

// This function tests the lazy computation of mps with interval arithmetic and superaccumulators
void debug_test(){
    // Small tests
    /* Control experimets first
    double array[] = {1.1,-10.,20.,-3.,1.};
    parallel_superacc_mps(array,5);
    
    double array2[] = {1.,pow(2,-50),-1.};
    parallel_superacc_mps(array2,3);

    double array3[] = {1.,pow(2,-53),-1.};
    parallel_mps_superacc(array3,3);
    */

    // Test regular addition rounding
    //cout << endl << (1 + pow(2,-60)) - 1 << endl;
    
    // Test how superacc deals with negative number
    Superaccumulator test = Superaccumulator();
    test.Accumulate(1);
    test.Accumulate(pow(2,-53));
    test.Accumulate(-1);

    cout << endl << test.Round() << endl;

    // Test interval arithmetic
    double brray[] = {1.75,-10.,20.,-3.,1.};
    parallel_mps_mpfr_lazy_2(brray,5);

    double brray4[] = {1.76,-10.,20.,-3.,1.};
    parallel_mps_mpfr_lazy_2(brray4,5);
    
    double brray2[] = {1.,pow(2,-50),-1.};
    parallel_mps_mpfr_lazy_2(brray2,3);

    double brray3[] = {1.,pow(2,-53),-1.};
    parallel_mps_mpfr_lazy_2(brray3,3);

    // Test what happens in case of uncertain comparison
    double crray[] = {1.,pow(2,-53),-1,1.,pow(2,-54)};
    parallel_mps_mpfr_lazy_2(crray,5);

    // Test big arrays
    srand(time(NULL));
    for(int i = 2; i < 2; i++){
        int size = pow(10,i);
        double* drray = new double[size];
        init_fpuniform(size, drray, 30, 15);
        // Randomly change signs
        for(int j = 0; j < size ; j++){
             drray[j] = (rand() % 2) ? drray[j] : -drray[j];
        }
        for(int k = 0; k < 5; k++){ 
            parallel_mps_superacc(drray,size);
        }
        delete[] drray;
    }
}

/* This function compares the runtime of mps with superaccumulators or mpfr, and its lazy implementation with interval arithmetic */
void runtime_comparison(){
    
    // Variables declaration and initialisation 
    double start;
    int size = pow(10,5);
    int N = 1;

    // for each dynamic range
    //vector<int> dynRanges {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
    vector<int> dynRanges  {10,20,40,60,80,100,200,400,600,800};
    int s = dynRanges.size();
    
    // Store results to plot
    fstream results;
    results.open("Plots/plot1.csv", ofstream::out | ofstream::trunc);
    vector<double> x(s),r1(s),r2(s),r3(s),r4(s),r5(s),r6(s);

    // Random seed
    srand(time(NULL));
    
    for(int r = 0; r < dynRanges.size(); r++){

        // initialization of means
        double mean_float = 0, mean_superacc = 0, mean_superacc_lazy = 0, mean_mpfr = 0, mean_lazy_mpfr = 0, mean_lazy_mpfr_2 = 0;
        
        // 1000 trials with N = 10^5
        for(int i = 0; i < N; i++){
            
            // Generating array
            double* drray = new double[size];
            init_fpuniform(size, drray, dynRanges[r], dynRanges[r]/2);

            // Randomly change signs
            for(int j = 0; j < size ; j++){
                 drray[j] = (rand() % 2) ? drray[j] : -drray[j];
            }
            
            double time_float = 0.0;
            PFP_TIME(parallel_mps_float(drray,size),start,time_float);
            double time_superacc = 0.0;
            PFP_TIME(parallel_mps_superacc(drray,size),start,time_superacc);
            double time_superacc_lazy = 0.0;
            PFP_TIME(parallel_mps_superacc_lazy(drray,size),start,time_superacc_lazy);
            double time_mpfr = 0.0;
            PFP_TIME(parallel_mps_mpfr(drray,size),start,time_mpfr);
            double time_mpfr_lazy = 0.0;
            PFP_TIME(parallel_mps_mpfr_lazy(drray,size),start,time_mpfr_lazy);
            double time_mpfr_lazy_2 = 0.0;
            PFP_TIME(parallel_mps_mpfr_lazy_2(drray,size),start,time_mpfr_lazy_2);
        
            /* Print times
            cout << endl << "Float time: " << time_float << " / Superacc time: " << time_superacc << " / Lazy superacc time: " << time_superacc_lazy << endl;
            cout << "Ratio superacc/float " << time_superacc/time_float << " / Ratio lazy superacc/superacc: " << time_superacc_lazy/time_superacc << endl;  
            cout << "Ratio mpfr/float " << time_mpfr/time_float << " / Ratio lazy mpfr/mpfr: " << time_mpfr_lazy/time_mpfr << endl;  
            */

            mean_float += time_float;
            mean_superacc += time_superacc;
            mean_superacc_lazy += time_superacc_lazy;
            mean_mpfr += time_mpfr;
            mean_lazy_mpfr += time_mpfr_lazy;
            mean_lazy_mpfr_2 += time_mpfr_lazy_2;
            
            delete[] drray;
        }
        // Finalize mean computation
        mean_float = mean_float/N;
        mean_superacc = mean_superacc/N;
        mean_superacc_lazy = mean_superacc_lazy/N;
        mean_mpfr = mean_mpfr/N;
        mean_lazy_mpfr = mean_lazy_mpfr/N;
        mean_lazy_mpfr_2 = mean_lazy_mpfr_2/N;
        
        x[r]= dynRanges[r];
        r1[r]= mean_float/mean_float;
        r2[r]= mean_superacc/mean_float;
        r3[r]= mean_superacc_lazy/mean_float;
        r4[r]= mean_mpfr/mean_float;
        r5[r]= mean_lazy_mpfr/mean_float;
        r6[r]= mean_lazy_mpfr_2/mean_float;

        // Writing results to a file
        results << to_string(x[r]) << "," << to_string(r1[r]) << "," << to_string(r2[r]) << "," << to_string(r3[r])<< "," << to_string(r4[r]) << "," << to_string(r5[r]) << "," << to_string(r6[r]) << endl;;
        
        
    }

    results.close();
}

int main(int argc, char** argv){
    if(argc >= 1){
        int a = atoi(argv[1]);
        cout <<endl<<a<<endl;
        switch (a){
            case 0 :
                debug_test();
                break;
            case 1 : 
                runtime_comparison();
                break;
        }
    }
}
