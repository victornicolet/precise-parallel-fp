/* Main file for comparison of mps methods.
 * Author: Raphael Dang-Nhu.
 * Date: 27/04/2018 */

#include <iostream>

#include "superaccumulator.hpp"
#include "lazy_mps_implementations.hpp"
#include "common.hpp"
#include "pfpdefs.hpp"

using namespace std;

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
    parallel_mps_mpfr(brray,5);

    double brray4[] = {1.76,-10.,20.,-3.,1.};
    parallel_mps_mpfr(brray4,5);
    
    double brray2[] = {1.,pow(2,-50),-1.};
    parallel_mps_mpfr(brray2,3);

    double brray3[] = {1.,pow(2,-53),-1.};
    parallel_mps_mpfr(brray3,3);

    // Test what happens in case of uncertain comparison
    double crray[] = {1.,pow(2,-53),-1,1.,pow(2,-54)};
    parallel_mps_mpfr(crray,5);

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

/* This function compares the runtime of mps with superaccumulators, and its lazy implementation with interval arithmetic */
void runtime_comparison(){
    
    double start;
    
    srand(time(NULL));
    for(int i = 6; i < 7; i++){
        // Generating arrays
        int size = pow(10,i);
        double* drray = new double[size];
        init_fpuniform(size, drray, 10, 0);
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
        PFP_TIME(parallel_mps_superacc_lazy(drray,size),start,time_mpfr_lazy);
    
        // Print times
        cout << endl << "Float time: " << time_float << " / Superacc time: " << time_superacc << " / Lazy superacc time: " << time_superacc_lazy << endl;
        cout << "Ratio superacc/float " << time_superacc/time_float << " / Ratio lazy superacc/superacc: " << time_superacc_lazy/time_superacc << endl;  
        cout << "Ratio mpfr/float " << time_mpfr/time_float << " / Ratio lazy mpfr/mpfr: " << time_mpfr_lazy/time_mpfr << endl;  
        delete[] drray;
    }
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
