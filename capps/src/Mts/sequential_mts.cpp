/* File for sequential mts implementations.
 * Author: Raphaël Dang-Nhu
 * Date: May 7th */

#include <iostream>
#include <emmintrin.h>

#include "interval_arithmetic.hpp"
#include "superaccumulator.hpp"

#define VERBOSE 0


void sequential_mts_superacc(double*array, long size, double* mts, long* pos){
    Superaccumulator mtsA = Superaccumulator();
    long t = 0;
    for(long i = 0; i != size; i++){
        mtsA.Accumulate(array[i]);
        if(mtsA.Normalize()){
            mtsA = Superaccumulator();
            t = i+1;
        }
    }
    *mts = mtsA.Round();
    *pos = t;
    
    if(VERBOSE){
        cout << endl << "mts with superaccumulators" << endl;
        cout << "mts: " << *mts << endl;
        cout << "Pos: " << *pos << endl;
    }

}

void sequential_mts_double(double* array, long size, double* mts, long* pos){
    long t = 0; 
    double mtst = 0.;
    for(long i = 0; i != size; i++){
        mtst += array[i];
        if(mtst < 0){
            mtst = 0;
            t = i+1;
        }
    }
    *mts = mtst;
    *pos = t;

    if(VERBOSE){
        cout << endl << "mts with doubles" << endl;
        cout << "mts: " << *mts << endl;
        cout << "Pos: " << *pos << endl;
    }
}

// Enum to represent decisions
enum Dec{
    D0,
    Dmts,
    Dundef,
    Dirrelevant
};

void sequential_mts_lazy(double* array, long size,  double* mts, long* pos,long opt){
    // longernal variables and memorization of decisions
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    __m128d mtsI = in2_create(0.,0.);     
    Dec* da = new Dec[size]; 
    long t = 0;
    
    /* First iteration with longerval arithmetic */
    for(long i = 0; i != size; i++){
        mtsI = in2_add(mtsI,array[i]);
        boolean b = inferior_double(0.,mtsI);

        if(b == False){
            da[i] = D0;
            mtsI = in2_create(0.,0.);
            t = i+1;
        }
        else if(b == True){
            da[i] = Dmts;
        }
        else {
            da[i] = Dundef;
            t = i+1;
        }
    }

    if(VERBOSE){
        cout << endl << "Lazy mts, first results"  << endl;
        cout <<  "mts: ";
        print(mtsI) ;
        cout << endl << "Pos: " << t << endl;
        for(long i = 0;  i != size; i++){
            Dec d = da[i]; 
            if(d == D0) cout << "X";
            else if(d == Dmts) cout << "_";
            else cout << "|";
        }
        cout << endl;
    }

    /* Iterate in da in reverse order */
    if(opt){
        bool seenD0 = 0;
        for(long i = size - 1; i >= 0 ; i--){
            Dec d = da[i];
            if(!seenD0 && d == D0){
                seenD0 = 1;
            }
            else if(seenD0){
                da[i] = Dirrelevant;
            }
        }
    }

    /* Second iteration with superaccumulators */
    _MM_SET_ROUNDING_MODE(0);
    Superaccumulator mtsA = Superaccumulator();

    for(long i = 0; i != size; i++){
    
        Dec d = da[i];
        if(d == D0){
            mtsA = Superaccumulator();
            t = i+1;
        }
        else if (d == Dmts){
            mtsA.Accumulate(array[i]);
        }
        else if (d == Dundef){
            // Redo the comparison
            if(mtsA.Normalize()){
                mtsA = Superaccumulator();
                t = i+1;
            }
            else {
                mtsA.Accumulate(array[i]);
            }
        }
    }
    delete[] da;
    *mts = mtsA.Round();
    *pos = t;

    if(VERBOSE){
        cout << endl << "Lazy mts, precise results" << endl;
        cout << "mts: " << *mts << endl;
        cout << "Pos: " << *pos << endl;
    }
}
                


    
        



