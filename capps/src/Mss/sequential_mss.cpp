/* File for sequential mss implementations.
 * Author: Raphaël Dang-Nhu
 * Date: May 7th */

#include <iostream>
#include <emmintrin.h>

#include "interval_arithmetic.hpp"
#include "superaccumulator.hpp"

#define VERBOSE 1


void sequential_mss_superacc(double*array, long size, double* mss, long* posl, long* posr){
    // Positions for mss
    long poslt = 0; 
    long posrt = 0;
    
    // Position for mts
    long post = 0; 

    Superaccumulator mssA = Superaccumulator();
    Superaccumulator mtsA = Superaccumulator();

    for(long i = 0; i != size; i++){

        // Update msst
        mtsA.Accumulate(array[i]);
        if(!mtsA.comp(mssA)){
            mssA = Superaccumulator(mtsA.get_accumulator());
            poslt = post;
            posrt = i+1;
        }

        // Update mtst
        if(mtsA.Normalize()){
            mtsA = Superaccumulator();
            post = i+1;
        }
    }
    *mss = mssA.Round();
    *posl = poslt;
    *posr = posrt;
    
    if(VERBOSE){
        cout << endl << "mss with superaccumulators" << endl;
        cout << "mss: " << *mss << endl;
        cout << "Posl: " << *posl << endl;
        cout << "Posr: " << *posr << endl;
    }

}

void sequential_mss_double(double* array, long size, double* mss, long* posl, long*posr){
    // Positions for mss
    long poslt = 0; 
    long posrt = 0;
    
    // Position for mts
    long post = 0; 

    double msst = 0.;
    double mtst = 0.;
    
    for(long i = 0; i != size; i++){

        // Update msst
        mtst += array[i];
        if(msst <= mtst){
            msst = mtst;
            poslt = post;
            posrt = i+1;
        }

        // Update mtst
        if(mtst < 0){
            mtst = 0.;
            post = i+1;
        }
    }
    *mss = msst;
    *posl = poslt;
    *posr = posrt;

    if(VERBOSE){
        cout << endl << "mss with doubles" << endl;
        cout << "mss: " << *mss << endl;
        cout << "Posl: " << *posl << endl;
        cout << "Posr: " << *posr << endl;
    }
}

// Enum to represent decisions for mss
enum Dec{
    Dmss,
    Dmts,
    Dundef
};

// Enum to represent decisions for mts
enum Decc{
    Dcmts,
    Dc0,
    Dcundef
};

// Enum to indicate which variables to compute
enum Var{
    Vmss,
    Vmts,
    Vnull
};

void sequential_mss_lazy(double* array, long size,  double* mss, long* posl, long* posr){
    // Rounding mode
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);

    __m128d mssI = in2_create(0.,0.);     
    __m128d mtsI = in2_create(0.,0.);     

    // Positions for mss
    long poslt = 0; 
    long posrt = 0;
    
    // Position for mts
    long post = 0; 

    Dec* da = new Dec[size]; 
    Decc* dca = new Decc[size]; 
    
    /* First iteration with longerval arithmetic */
    for(long i = 0; i != size; i++){
        mtsI = in2_add_double(mtsI,array[i]);

        // Update mssI
        boolean b = inferior(mssI,mtsI);

        if(b == True){
            da[i] = Dmts;
            mssI = mtsI;
            poslt = post;
            posrt = i+1;
        }
        else if(b == False){
            da[i] = Dmss;
        }
        else {
            da[i] = Dundef;
            mssI = in2_merge(mssI,mtsI);
            poslt = post;
            posrt = i+1;
        }

        // Update mtsI
        boolean b0 = inferior_double(0,mtsI);

        if(b0 == False){
            dca[i] = Dc0;
            mtsI = in2_create(0.,0.);
            post = i+1;
        }
        else if (b0 == True){
            dca[i] = Dcmts;
        }
        else{
            dca[i] = Dcundef;
            post = i+1;
        }
    }


    /* Iterate in reverse order */
    Var* va = new Var[size];
    Var v = Vmss;
    for(long i = size - 1; i >= 0 ; i--){
        va[i] = v;
        if(v == Vmss && da[i] == Dmts){
            v = Vmts;
        }
        else if(v == Vmts && dca[i] == Dc0){
            v = Vnull;
        }
    }

    if(VERBOSE){
        cout << endl << "Lazy mss, first results"  << endl;
        cout <<  "mss: ";
        print(mssI) ;
        cout << endl << "Posl: " << poslt ;
        cout << endl << "Posr: " << posrt << endl << endl;
        for(long i = 0;  i != size; i++){
            Dec d = da[i]; 
            if(d == Dmts) cout << "X";
            else if(d == Dmss) cout << "_";
            else cout << "|";
        }
        cout << endl << endl << endl << endl;
        for(long i = 0;  i != size; i++){
            Decc dc = dca[i]; 
            if(dc == Dc0) cout << "X";
            else if(dc == Dcmts) cout << "_";
            else cout << "|";
        }
        cout << endl << endl << endl << endl;
        for(long i = 0;  i != size; i++){
            Var v = va[i]; 
            if(v == Vmss) cout << "2";
            else if(v == Vmts) cout << "1";
            else cout << "0";
        }
        cout << endl;
    }

    /* Second iteration with superaccumulators */
    // Positions for mss
    poslt = 0; 
    posrt = 0;
    
    // Position for mts
    post = 0; 

    _MM_SET_ROUNDING_MODE(0);
    Superaccumulator mssA = Superaccumulator();
    Superaccumulator mtsA = Superaccumulator();

    for(long i = 0; i != size; i++){
        
        // Check which variables need to be computed
        Var v = va[i];
        
        if(v == Vmss){
            mtsA.Accumulate(array[i]);
            // Check decision
            Dec d = da[i];
            if(d == Dmts){
                mssA = Superaccumulator(mtsA.get_accumulator());
                poslt = post;
                posrt = i+1;
            }
            else if (d == Dundef){
                // Update msst
                if(!mtsA.comp(mssA)){
                    mssA = Superaccumulator(mtsA.get_accumulator());
                    poslt = post;
                    posrt = i+1;
                }
            }
            // Check decision
            Decc dc = dca[i];

            if(dc == Dc0){
                mtsA = Superaccumulator();
                post = i+1;
            }
            else if (dc == Dcundef){
                // Redo the comparison
                if(mtsA.Normalize()){
                    mtsA = Superaccumulator();
                    post = i+1;
                }
                else {
                    mssA.Accumulate(array[i]);
                }
            }
        }
        else if(v == Vmts){
            // Check decision
            Decc dc = dca[i];

            if(dc == Dc0){
                mtsA = Superaccumulator();
                post = i+1;
            }
            else if (dc == Dcmts){
                mtsA.Accumulate(array[i]);
            }
            else if (dc == Dcundef){
                // Redo the comparison
                if(mtsA.Normalize()){
                    mtsA = Superaccumulator();
                    post = i+1;
                }
                else {
                    mssA.Accumulate(array[i]);
                }
            }
        }
    }
    delete[] da;
    delete[] dca;
    delete[] va;
    *mss = mssA.Round();
    *posl = poslt;
    *posr = posrt;

    if(VERBOSE){
        cout << endl << "Lazy mss, precise results" << endl;
        cout << "mss: " << *mss << endl;
        cout << "Posl: " << *posl << endl;
        cout << "Posr: " << *posr << endl;
    }
}
                


    
        



