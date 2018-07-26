/* File for sequential viterbi implementations.
 * Author: Raphaël Dang-Nhu
 * Date: July 26th */

#include <iostream>
#include <emmintrin.h>
#include <mpreal.h>
#include <mpfr.h>
#include <gmp.h>

#include "interval_arithmetic.hpp"
#include "debug.hpp"

using mpfr::mpreal;


void sequential_viterbi(double*array, long size, double* sum, double* mps, int* pos){
    mpreal dsumA(0.,1000);
    mpreal mpsA(0.,1000);
    int t = 0;

    for(int i = 0; i != size; i++){
        dsumA += array[i];
        if(dsumA >= 0){
            mpsA += dsumA;
            dsumA = 0;
            t = i+1;
        }
    }
    *deltasum = dsumA.toDouble();
    *mps = mpsA.toDouble();
    *pos = t;
    
    if(PRINT){
        cout << endl << "Mps with mpfr" << endl;
        cout << "Delta Sum: " << *deltasum << endl;
        cout << "Mps: " << *mps << endl;
        cout << "Pos: " << *pos << endl;
    }

}

void sequential_mps_double(double* array, int size, double* deltasum, double* mps, int* pos){
    int t = 0; 
    double sumt = 0.;
    double mpst = 0.;
    for(int i = 0; i != size; i++){
        sumt += array[i];
        if(sumt >= 0){
            mpst += sumt;
            sumt = 0;
            t = i+1;
        }
    }
    *deltasum = sumt;
    *mps = mpst;
    *pos = t;

    if(PRINT){
        cout << endl << "Mps with doubles" << endl;
        cout << "Delta Sum: " << *deltasum << endl;
        cout << "Mps: " << *mps << endl;
        cout << "Pos: " << *pos << endl;
    }
}


void sequential_mps_interval_memorized(double* array, int size, double* deltasum, double* mps, int* pos, memo** da){
    // Internal variables and memorization of decisions
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    __m128d sumI = in2_create(0.);     
    __m128d mpsI = in2_create(0.);     
    memo* d = new memo[size](); 
    int t = 0;
    
    /* First iteration with interval arithmetic */
    for(int i = 0; i != size; i++){
        sumI = in2_add(sumI,array[i]);
        boolean b = in2_le(0,sumI);

        if(b == True){
            mpsI = in2_add(mpsI,sumI);
            sumI = in2_create(0.);
            t = i+1;
        }
        else if(b == Undefined){
            sumI = in2_merge(sumI,in2_create(0.));
            mpsI = in2_add(mpsI,sumI);
            t = i+1;
        }
        d[i].useful2 = b;
    }

    *da = d;

    if(PRINT){
        cout << endl << "Interval mps"  << endl;
        cout << "Delta Sum:";
        print(sumI);
        cout << endl << "Mps:";
        print(mpsI);
        cout << endl << "Pos:" << t << endl;
        /*for(int i = 0;  i != size; i++){
            Dec d = da[i]; 
            if(d == True) cout << "X";
            else if(d == False) cout << "_";
            else cout << "|";
        }*/
        cout << endl;
    }
    _MM_SET_ROUNDING_MODE(0);
}

void sequential_viterbi_iterate_reverse_mps(double* array, int size, double* sum, double* mps, int* pos, memo** da){
    /* Iterate in da in reverse order */
    bool sum_useful = false, mps_useful = true;
    memo* dap = *da;
    for(int i = size - 1; i >= 0 ; i--){
        memo d = dap[i];

        
        // Handling step2
        if(d.useful2 == True){
            if(sum_useful){
                sum_useful = false;
            }
            if(mps_useful){
                sum_useful = true;
            }
            else{
                dap[i].useful2 = Useless;
            }
        }
        else if(d.useful2 == False){
            dap[i].useful2 = Useless;
        }
        else{
            bool t = mps_useful || sum_useful;
            if(mps_useful){
                sum_useful = true;
            }
            if(!t){
              sum_useful = true;
            }
            else{ 
                dap[i].useful2 = Useless;
            }
        }

        // Handling step1
        if(sum_useful){
            dap[i].useful1 = true;
        }
    }
    if(PRINT){
        cout << "Finished reverse mps" << endl;
    }
}

                
void sequential_viterbi_lazy_mpfr(double* array, int size, double* sum, double* mps, int* pos, memo** da){

    /* Second iteration with superaccumulators */
    mpreal sumA(0.,1000);
    mpreal mpsA(0.,1000);
    double post = 0;
    memo* dap = *da;

    for(int i = 0; i != size; i++){
        memo d = dap[i];
        if(d.useful1){
            sumA += array[i];
        }
     
        if(d.useful2 == True){
            mpsA += sumA;
            sumA = 0.;
            post = i+1;
        }
        else if (d.useful2 == False){
        }
        else if (d.useful2 == Undefined){
            // Redo the comparison
            if(sumA >= 0){
                mpsA += sumA;
                sumA = 0.;
                post = i+1;
            }
        }
    }
    delete[] dap;
    *sum = sumA.toDouble();
    *mps = mpsA.toDouble();
    *pos = post;

    if(PRINT){
        cout << endl << "Exact computation with mpfr" << endl;
        cout << "Delta Sum: " << *sum << endl;
        cout << "Mps: " << *mps << endl;
        cout << "Pos: " << *pos << endl;
    }
}

