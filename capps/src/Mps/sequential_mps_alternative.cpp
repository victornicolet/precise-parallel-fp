/* File for sequential mps implementations.
 * Author: Raphaël Dang-Nhu
 * Date: May 7th */

#include <iostream>
#include <emmintrin.h>

#include "interval_arithmetic.hpp"
#include "superaccumulator.hpp"
#include "debug.hpp"

void sequential_mps_superacc_alt(double*array, int size, double* sum, double* mps, int* pos){
    Superaccumulator sumA = Superaccumulator();
    Superaccumulator mpsA = Superaccumulator();
    int t = 0;
    for(int i = 0; i != size; i++){
        sumA.Accumulate(array[i]);
        if(!sumA.comp(mpsA)){
            mpsA = Superaccumulator(sumA.get_accumulator());
            t = i+1;
        }
    }
    *sum = sumA.Round();
    *mps = mpsA.Round();
    *pos = t;
    
    if(PRINT){
        cout << endl << "Mps with superaccumulators" << endl;
        cout << "Sum: " << *sum << endl;
        cout << "Mps: " << *mps << endl;
        cout << "Pos: " << *pos << endl;
    }

}

void sequential_mps_double_alt(double* array, int size, double* sum, double* mps, int* pos){
    int t = 0; 
    double sumt = 0.;
    double mpst = 0.;
    for(int i = 0; i != size; i++){
        sumt += array[i];
        if(sumt >= mpst){
            mpst = sumt;
            t = i+1;
        }
    }
    *sum = sumt;
    *mps = mpst;
    *pos = t;

    if(PRINT){
        cout << endl << "Mps with doubles" << endl;
        cout << "Sum: " << *sum << endl;
        cout << "Mps: " << *mps << endl;
        cout << "Pos: " << *pos << endl;
    }
}


void sequential_mps_interval_memorized_alt(double* array, int size, double* sum, double* mps, int* pos, memo* da){
    // Internal variables and memorization of decisions
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    __m128d sumI = in2_create(0.,0.);     
    __m128d mpsI = in2_create(0.,0.);     
    da = new memo[size]; 
    int t = 0;
    
    /* First iteration with interval arithmetic */
    for(int i = 0; i != size; i++){
        sumI = in2_add_double(sumI,array[i]);
        boolean b = inferior(mpsI,sumI);

        if(b == True){
            da[i].useful2 = True;
            mpsI = sumI;
            t = i+1;
        }
        else if(b == False){
            da[i].useful2 = False;
        }
        else {
            da[i].useful2 = Undefined;
            mpsI = in2_max(sumI,mpsI);
            t = i+1;
        }
    }

    if(PRINT){
        cout << endl << "Lazy mps, first results"  << endl;
        cout << "Sum: ";
        print(sumI);
        cout << endl << "Mps: ";
        print(mpsI) ;
        cout << endl << "Pos: " << t << endl;
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

void sequential_mps_iterate_reverse_mps_alt(double* array, int size, double* sum, double* mps, int* pos, memo* da){
    /* Iterate in da in reverse order */
    bool sum_useful = false, mps_useful = true;
    for(int i = size - 1; i >= 0 ; i--){
        memo d = da[i];

        // Handling step1
        if(sum_useful){
            d.useful1 = true;
        }
        
        // Handling step2
        if(d.useful2 == True){
            if(mps_useful){
                sum_useful = true;
                mps_useful = false;
            }
            else{
                da[i].useful2 = Useless;
            }
        }
        else if(d.useful2 == False){
            da[i].useful2 = Useless;
        }
        else{
            if(mps_useful){
                sum_useful = true;
                mps_useful = true;
            }
            else{
                da[i].useful2 = Useless;
            }
        }
    }
}

void sequential_mps_iterate_reverse_pos_alt(double* array, int size, double* sum, double* mps, int* pos, memo* da){
    /* Iterate in da in reverse order */
    bool sum_useful = false, mps_useful = false, pos_useful = true;
    for(int i = size - 1; i >= 0 ; i--){
        memo d = da[i];

        // Handling step1
        if(sum_useful){
            d.useful1 = true;
        }
        
        // Handling step2
        if(d.useful2 == True){
            if(mps_useful){
                sum_useful = true;
                mps_useful = false;
            }
            if(pos_useful){
                pos_useful = false;
            }
            if(!(mps_useful || pos_useful)){
                da[i].useful2 = Useless;
            }
        }
        else if(d.useful2 == False){
            da[i].useful2 = Useless;
        }
        else{
            if(mps_useful){
                sum_useful = true;
                mps_useful = true;
            }
            if(pos_useful){
                sum_useful = true;
                mps_useful = true;
            }
            if(!(mps_useful || pos_useful)){
                da[i].useful2 = Useless;
            }
        }
    }
}

void sequential_mps_lazy_superacc_alt(double* array, int size, double* sum, double* mps, int* pos, memo* da){

    /* Second iteration with superaccumulators */
    Superaccumulator sumA = Superaccumulator();
    Superaccumulator mpsA = Superaccumulator();
    double post = 0;

    for(int i = 0; i != size; i++){
        memo d = da[i];
        if(d.useful1){
            sumA.Accumulate(array[i]);
        }
     
        if(d.useful2 == True){
            mpsA = Superaccumulator(sumA.get_accumulator());
            post = i+1;
        }
        else if (d.useful2 == False){
        }
        else if (d.useful2 == Undefined){
            // Redo the comparison
            if(!sumA.comp(mpsA)){
                mpsA = Superaccumulator(sumA.get_accumulator());
                post = i+1;
            }
        }
    }
    delete[] da;
    *sum = sumA.Round();
    *mps = mpsA.Round();
    *pos = post;

    if(PRINT){
        cout << endl << "Lazy mps superacc, precise results" << endl;
        cout << "Sum: " << *sum << endl;
        cout << "Mps: " << *mps << endl;
        cout << "Pos: " << *pos << endl;
    }
}
                
void sequential_mps_lazy_mpfr_alt(double* array, int size, double* sum, double* mps, int* pos, memo* da){
    /* Second iteration with superaccumulators */
    Superaccumulator sumA = Superaccumulator();
    Superaccumulator mpsA = Superaccumulator();
    double post = 0;

    for(int i = 0; i != size; i++){
        memo d = da[i];
        if(d.useful1){
            sumA.Accumulate(array[i]);
        }
     
        if(d.useful2 == True){
            mpsA = Superaccumulator(sumA.get_accumulator());
            post = i+1;
        }
        else if (d.useful2 == False){
        }
        else if (d.useful2 == Undefined){
            // Redo the comparison
            if(!sumA.comp(mpsA)){
                mpsA = Superaccumulator(sumA.get_accumulator());
                post = i+1;
            }
        }
    }
    delete[] da;
    *sum = sumA.Round();
    *mps = mpsA.Round();
    *pos = post;

    if(PRINT){
        cout << endl << "Lazy mps superacc, precise results" << endl;
        cout << "Sum: " << *sum << endl;
        cout << "Mps: " << *mps << endl;
        cout << "Pos: " << *pos << endl;
    }
}


    
        


