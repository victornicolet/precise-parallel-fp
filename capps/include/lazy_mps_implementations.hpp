/* Implementation of lasy parallel maximum prefix sum with interval arithmetic and superaccumulators.
 * Author: Raphael Dang-Nhu 
 * Date: 04/19/2018 */

#ifndef LAZY_MPS_H_
#define LAZY_MPS_H_

#include <emmintrin.h>
#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"
#include "interval_arithmetic.hpp"

using namespace tbb;
using namespace std;

// Same structure with lazy computation, i.a. and superaccumulators. 
template <typename __mps_high_precision> struct __mps{
    // pointer to the array
    double* array;
    // Both sum and maximum prefix sum are represented as intervals  
    __m128d sum_interval;
    __m128d mps_interval;
    // Position of the maximum prefix sum
    int position;
    // size
    int size;
    int left;
    int right;
    // Constructor
    __mps(double* a);
    // Splitting constructor
    __mps(__mps&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<int>&);
    // Join operation for the reduction
    void join(__mps& rightMps); 
    // Printing function
    void print_mps();
};

/* MPS computation with floats */
void parallel_mps_float(double*,int);

/* Function computing in parallel the mps with superaccumulators, same arguments as previous function */
void parallel_mps_superacc(double*,int);

/* Function computing mps with mpfr */
void parallel_mps_mpfr(double*,int);

/* Function computimg mps with Collange mixed FPE-superacc solution (not working yet)*/
void parallel_mps_Collange(double*,int);

/* Lazy computation of mps, with superaccs for precise computations */
void parallel_mps_superacc_lazy(double*,int);

/* Lazy computations of mps, with mpfr for precise computations */
void parallel_mps_mpfr_lazy(double*,int);

template <typename __mps_high_precision> __mps<__mps_high_precision>::__mps(double* a) :
    array(a),
    position(-1),
    size(0),
    left(-1),
    right(-1)
{
    sum_interval = in2_create(0.,0.);
    mps_interval = in2_create(0.,0.);
}

template <typename __mps_high_precision> __mps<__mps_high_precision>::__mps(__mps& x, split) :
    array(x.array),
    position(-1),
    size(0),
    left(-1),
    right(-1)
{
    sum_interval = in2_create(0.,0.);
    mps_interval = in2_create(0.,0.);
}

template <typename __mps_high_precision> void __mps<__mps_high_precision>::print_mps(){
    cout << "sum: ";
    print(sum_interval);
    cout << endl << "mps: ";
    print(mps_interval);
    cout << endl << "size: " << size << endl << "position: " << position << endl;
}

template <typename __mps_high_precision> void __mps<__mps_high_precision>::join(__mps& rightMps){
   size += rightMps.size;
   right = rightMps.right;
   // computing sum-l + mps-r
   __m128d mpsCandidate = in2_add(sum_interval,rightMps.mps_interval);
   // adding two sums
   sum_interval = in2_add(sum_interval,rightMps.sum_interval);
   // comparison of mpsCandidate and mps-l
   boolean b = inferior(mps_interval,mpsCandidate);
   if (b == True){
       mps_interval = mpsCandidate;
       position = rightMps.position;
   }
   // In case of undefined comparison:
   else if(b == Undefined){
       __mps_high_precision precise_result(array);
       parallel_reduce(blocked_range<int>(left,right),precise_result);
       double mpsaux = precise_result.getMpsDouble();
       mps_interval = in2_create(mpsaux,mpsaux);
       double sumaux = precise_result.getSumDouble();
       sum_interval = in2_create(sumaux,sumaux);
       position = precise_result.position;
   } 
}

template <typename __mps_high_precision> void __mps<__mps_high_precision>::operator()(const blocked_range<int>& r){
    if(position == -1){
        position = r.begin();
    }
    if(left == -1){
        left = r.begin();
    }
    right = r.end();
    size = right - left;
    // iterating over the subrange
    for(int i = r.begin(); i != r.end(); i++){
        sum_interval = in2_add_double(sum_interval,array[i]);
        boolean b = inferior(mps_interval,sum_interval);
        if (b == True){
            mps_interval = sum_interval;
            position = i; 
        }
        else if(b == Undefined){
           __mps_high_precision precise_result(array);
           parallel_reduce(blocked_range<int>(left,right),precise_result);
           double mpsaux = precise_result.getMpsDouble();
           mps_interval = in2_create(mpsaux,mpsaux);
           double sumaux = precise_result.getSumDouble();
           sum_interval = in2_create(sumaux,sumaux);
           position = precise_result.position;
           break;
        }
    }
}

#endif
