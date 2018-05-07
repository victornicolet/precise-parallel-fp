/* File for sequential mps implementations.
 * Author: Raphaël Dang-Nhu
 * Date: May 7th */

#include <iostream>
#include <emmintrin.h>

#include "interval_arithmetic.hpp"
#include "superaccumulator.hpp"

#define VERBOSE 1


void sequential_summation_superacc(double* array, int size, double* sum){
    Superaccumulator sumA = Superaccumulator();
    for(int i = 0; i != size; i++){
        sumA.Accumulate(array[i]);
    }
    *sum = sumA.Round();

    if(VERBOSE){
        cout << endl << "Summation with superaccumulators" << endl;
        cout << "Sum: " << *sum << endl;
    }

}

void sequential_mps_superacc(double*array, int size, double* sum, double* mps, int* pos){
    Superaccumulator sumA = Superaccumulator();
    Superaccumulator mpsA = Superaccumulator();
    int t = 0;
    for(int i = 0; i != size; i++){
        sumA.Accumulate(array[i]);
        if(sumA.comp(mpsA)){
            mpsA = Superaccumulator(sumA.get_accumulator());
            t = i+1;
        }
    }
    *sum = sumA.Round();
    *mps = mpsA.Round();
    *pos = t;
    
    if(VERBOSE){
        cout << endl << "Mps with superaccumulators" << endl;
        cout << "Sum: " << *sum << endl;
        cout << "Mps: " << *mps << endl;
        cout << "Pos: " << *pos << endl;
    }

}

void sequential_mps_double(double* array, int size, double* sum, double* mps, int* pos){
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

    if(VERBOSE){
        cout << endl << "Mps with doubles" << endl;
        cout << "Sum: " << *sum << endl;
        cout << "Mps: " << *mps << endl;
        cout << "Pos: " << *pos << endl;
    }
}

// Enum to represent decisions
enum Dec{
    Dsum,
    Dmps,
    Dundef,
};

void sequential_mps_lazy(double* array, int size, double* sum, double* mps, int* pos){
    // Internal variables and memorization of decisions
    __m128d sumI = in2_create(0.,0.);     
    __m128d mpsI = in2_create(0.,0.);     
    Dec* da = new Dec[size]; 
    int t = 0;
    
    /* First iteration with interval arithmetic */
    for(int i = 0; i != size; i++){
        sumI = in2_add_double(sumI,array[i]);
        boolean b = inferior(mpsI,sumI);

        if(b == True){
            da[i] = Dsum;
            mpsI = sumI;
            t = i+1;
        }
        else if(b == False){
            da[i] = Dmps;
        }
        else {
            da[i] = Dundef;
            mpsI = in2_merge(sumI,mpsI);
            t = i+1;
        }
    }

    if(VERBOSE){
        cout << endl << "Lazy mps, first results"  << endl;
        cout << "Sum: ";
        print(sumI);
        cout << endl << "Mps: ";
        print(mpsI) ;
        cout << endl << "Pos: " << t << endl;
        for(int i = 0; i++; i != size){
            Dec d = da[i]; 
            if(d == Dsum) cout << "X";
            else if(d = Dmps) cout << "_";
            else cout << "|";
        }
        cout << endl;
    }
    
    /* Second iteration with superaccumulators */
    Superaccumulator sumA = Superaccumulator();
    Superaccumulator mpsA = Superaccumulator();

    for(int i = 0; i != size; i++){
        sumA.Accumulate(array[i]);
    
        Dec d = da[i];
        if(d == Dsum){
            mpsA = Superaccumulator(sumA.get_accumulator());
            t = i+1;
        }
        else if (d == Dundef){
            // Redo the comparison
            if(sumA.comp(mpsA)){
                mpsA = Superaccumulator(sumA.get_accumulator());
                t = i+1;
            }
        }
    }
    delete[] da;
    *sum = sumA.Round();
    *mps = mpsA.Round();
    *pos = t;

    if(VERBOSE){
        cout << endl << "Lazy mps, precise results" << endl;
        cout << "Sum: " << *sum << endl;
        cout << "Mps: " << *mps << endl;
        cout << "Pos: " << *pos << endl;
    }
}
                


    
        



