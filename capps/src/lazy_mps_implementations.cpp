/* Implementation of lazy parallel maximum prefix sum with interval arithmetic and superaccumulators.
 * Author: Raphael Dang-Nhu 
 * Date: 04/19/2018 */

#include <iostream>
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"

#include "precise_mps_implementations.hpp"
#include "lazy_mps_implementations.hpp"
#include "interval_arithmetic.hpp"

using namespace std;
using namespace tbb;


__mps::__mps(double* a) :
    array(a),
    position(-1),
    size(0),
    left(-1),
    right(-1)
{
    sum_interval = in2_create(0.,0.);
    mps_interval = in2_create(0.,0.);
}

__mps::__mps(__mps& x, split) :
    array(x.array),
    position(-1),
    size(0),
    left(-1),
    right(-1)
{
    sum_interval = in2_create(0.,0.);
    mps_interval = in2_create(0.,0.);
}

void __mps::print_mps(){
    cout << "sum: ";
    print(sum_interval);
    cout << endl << "mps: ";
    print(mps_interval);
    cout << endl << "size: " << size << endl << "position: " << position << endl;
}

void __mps::join(__mps& rightMps){
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
       __mps_precise<4> precise_result(array);
       parallel_reduce(blocked_range<int>(left,right),precise_result);
       double mpsaux = precise_result.mpsA.Round();
       mps_interval = in2_create(mpsaux,mpsaux);
       double sumaux = precise_result.sumA.Round();
       sum_interval = in2_create(sumaux,sumaux);
       position = precise_result.position;
   } 
}

void __mps::operator()(const blocked_range<int>& r){
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
           __mps_precise<4> precise_result(array);
           parallel_reduce(blocked_range<int>(left,right),precise_result);
           double mpsaux = precise_result.mpsA.Round();
           mps_interval = in2_create(mpsaux,mpsaux);
           double sumaux = precise_result.sumA.Round();
           sum_interval = in2_create(sumaux,sumaux);
           position = precise_result.position;
           break;
        }
    }
}
void printA(double* array, int size){
    // printing the array
    cout << endl << "{";
    if(size <= 10){
        for(int i = 0; i < size; i++){
            cout << array[i] << ",";
        }
        cout << "}" << endl;
    }
    else{
        for(int i = 0; i < 10; i++){
            cout << array[i] << ",";
        }
        cout << "... (size: " << size << ")" << endl;
    }

}

void parallel_mps(double* array, int size){
    printA(array,size);
    // Computing sum and mp__mps_acc result(array);
    __mps_naive result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    result.print_mps();
}
void parallel_superacc_mps(double* array, int size){
    printA(array,size);
    // Computing sum and mps
    __mps_precise<4> result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    result.print_mps();
}

void parallel_lazy_superacc_mps(double* array, int size){
    printA(array,size);
    // Computing sum and mps
    __mps result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    result.print_mps();
}
