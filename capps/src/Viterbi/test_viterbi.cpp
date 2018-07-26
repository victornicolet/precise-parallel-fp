/* Main file for comparison of mps methods.
 * Author: Raphael Dang-Nhu.
 * Date: 27/04/2018 */


#include <iostream>
#include <fstream>
#include <time.h>

#include "summation.hpp"
#include "common.hpp"
#include "pfpdefs.hpp"
#include "debug.hpp"
#include "printFunctions.hpp"
#include "sequential_viterbi.hpp"

using namespace tbb;
using namespace std;

void runtime_comparison(){
    
    // Variables declaration and initialisation 
    double start;
    int size = pow(10,5);
    int N = 5;

    // for each dynamic range
    vector<int> dynRanges  {2000};
    int s = dynRanges.size();
    
    // Store results to plot
    fstream results;
    results.open("Plots/csv/lazyviterbi.csv", ofstream::out | ofstream::trunc);
    vector<double> x(s),r0(s),r1(s),r2(s),r3(s),r4(s);

    // Random seed
    srand(time(NULL));
   
    for(int r = 0; r < dynRanges.size(); r++){
        if(PRINT){
            cout << endl << endl <<"***********************************" << endl << "New Input Array, Mps Alt" << endl;
        }
        // Generating array
        double* drray = new double[size];
        init_fpuniform(size, drray, dynRanges[r], dynRanges[r]/2);

        // Randomly change signs
        for(int j = 0; j < size ; j++){
             drray[j] = (rand() % 2) ? drray[j] : -drray[j];
        }

        // initialization of means
        double mean_double = 0.;
        double mean_mpfr = 0.;
        double mean_interval = 0.;
        double mean_reverse = 0.;
        double mean_lazy_mpfr = 0.;
        

        for(int i = 0; i < N; i++){

            // Declare result variables
            double sum;
            bool b;
            boolean* da;
            
            double time_double = 0.0;
            PFP_TIME(viterbi_double(drray,size,&sum,&b),start,time_double);
            double time_mpfr = 0.0;
            PFP_TIME(viterbi_mpfr(drray,size,&sum,&b),start,time_mpfr);

            /* Lazy computation, mpfr */
            double time1 = 0.0;
            PFP_TIME(viterbi_interval(drray,size,&sum,&b,&da),start,time7);
            double time2 = 0.0;
            PFP_TIME(viterbi_reverse(drray,size,&sum,&b,&da),start,time8);
            double time3 = 0.0;
            PFP_TIME(viterbi_lazy(drray,size,&sum,&b,&da),start,time9);

            mean_double += time_double;
            mean_mpfr += time_mpfr;
            mean_interval += time1;
            mean_reverse += time2;
            mean_lazy += time13;
            
        }
        // Finalize mean computation
        mean_double = mean_double / N;
        mean_mpfr = mean_mpfr / N;
        mean_interval = mean_interval/N;
        mean_reverse = mean_reverse / N;
        mean_lazy = mean_lazy / N;
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
                runtime_comparison();
                break;
        }
    }
}
