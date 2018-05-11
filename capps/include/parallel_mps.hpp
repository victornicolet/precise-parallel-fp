/* Interface for several precise mps implementations.
 * Author: Raphael Dang-Nhu
 * Date: 26/04/18 */

#ifndef PRECISE_MPS_IMPLEMENTATIONS_H_
#define PRECISE_MPS_IMPLEMENTATIONS_H_

#include <mpfr.h>
#include "tbb/blocked_range.h"
#include "superaccumulator.hpp"

using namespace tbb;
using namespace std;


// Naive mps structure
struct __mps_naive{
    // pointer to the array
    double* array;
    // Superaccumulators for sum and mps
    double sum;
    double mps;
    // Position of the maximum prefix sum
    int position;
    // Constructor
    __mps_naive(double* a);
    // Splitting constructor
    __mps_naive(__mps_naive&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<int>&);
    // Join operation for the reduction
    void join(__mps_naive& rightMps); 
    // Printing function
    void print_mps();
};

// Same structure with Collange implementation. 
template<int N> struct __mps_precise{
    // pointer to the array
    double* array;
    // Arrays for FPE
    double sum[N];
    double mps[N];
    // Superpreciseumulators for sum and mps
    Superaccumulator sumA;
    Superaccumulator mpsA;
    // Position of the maximum prefix sum
    int position;
    // Constructor
    __mps_precise(double* a);
    // Splitting constructor
    __mps_precise(__mps_precise&,split);
    // preciseumulate result for subrange
    void operator()(const blocked_range<int>&);
    // Join operation for the reduction
    void join(__mps_precise& rightMps); 
    // Printing function
    void print_mps();
    // Rounding functions
    double getSumDoubleDown();
    double getMpsDoubleDown();
    double getSumDoubleUp();
    double getMpsDoubleUp();
};

// Same structure with super accumulators only. 
struct __mps_acc{
    // pointer to the array
    double* array;
    // Superaccumulators for sum and mps
    Superaccumulator sum;
    Superaccumulator mps;
    // Position of the maximum prefix sum
    int position;
    // Constructor
    __mps_acc(double* a);
    // Splitting constructor
    __mps_acc(__mps_acc&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<int>&);
    // Join operation for the reduction
    void join(__mps_acc& rightMps); 
    // Printing function
    void print_mps();
    // Rounding functions
    double getSumDoubleDown();
    double getMpsDoubleDown();
    double getSumDoubleUp();
    double getMpsDoubleUp();
};

// Same structure with mpfr library.
struct __mps_mpfr{
    double* array;
    mpfr_t sum;
    mpfr_t mps;
    int position;
    __mps_mpfr(double* a);
    __mps_mpfr(__mps_mpfr&,split);
    void operator()(const blocked_range<int>&);
    void join(__mps_mpfr& right);
    void print_mps();
    // Rounding functions
    double getSumDoubleDown();
    double getMpsDoubleDown();
    double getSumDoubleUp();
    double getMpsDoubleUp();
};


inline static double twosum(double a, double b, double & s)
{
    double r = a + b;
    double z = r - a;
    s = (a - (r - z)) + (b - z);
    return r;
}


 template<int N>  void __mps_precise<N>::operator()(const blocked_range<int>& r){

    _MM_SET_ROUNDING_MODE(0);

    if(position == -1){
        position = r.begin();
    }
    for(int i = r.begin(); i!= r.end(); i++){
       // Accumulate array[i]
       double x = array[i];
       double s1;
       double sumtest = 0, mpstest = 0;
       for(int j = 0; j != N; j++){
            sum[j] = twosum(sum[j],x,s1);
            sumtest += sum[j];
            x = s1;
            mpstest += mps[j];
            cout << sum[j];
       }
       if(mpstest<=sumtest){
           for(unsigned int j = 0; j != N; ++j) {
               mps[j] = sum[j];
           }
           position = i;
       }
       if(x != 0){
           sumA.Accumulate(x);
           mpsA.Accumulate(x);
       }
    }

    // Flushing the results into the superaccumulator
    for(int j = 0; j!= N; j++){
        sumA.Accumulate(sum[j]);
        mpsA.Accumulate(mps[j]);
    }

    // To be completed, not ready yet

    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);

}

 template<int N> void __mps_precise<N>::join(__mps_precise& right){
    _MM_SET_ROUNDING_MODE(0);

   // computing sum-l + mps-r
   right.sumA.Accumulate(sumA);
   right.mpsA.Accumulate(sumA);
   // comparison of mpsCandidate and mps-l
   if(right.mpsA.Round() >= mpsA.Round()){
       mpsA = right.mpsA;
       position = right.position;
   }
        
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
}


 template<int N> __mps_precise<N>::__mps_precise(double* a) :
    array(a),
    position(-1)
{
    sumA = Superaccumulator();
    mpsA = Superaccumulator();
}

template<int N> __mps_precise<N>::__mps_precise(__mps_precise& parent,split) :
    array(parent.array),
    position(-1)
{
    sumA = Superaccumulator();
    mpsA = Superaccumulator();
}

 template<int N> void __mps_precise<N>::print_mps(){
    cout << "sum: " << sumA.Round();
    cout << endl << "mps: " << mpsA.Round();
    cout << endl << "position: " << position << endl;
}
#endif
