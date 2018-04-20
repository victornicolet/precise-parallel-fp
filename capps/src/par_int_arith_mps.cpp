/* Implementation of parallel maximum prefix sum with interval arithmetic.
 * Author: Raphael Dang-Nhu 
 * Date: 04/19/2018 */

#include <iostream>
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"

#include "par_int_arith_mps.hpp"
#include "interval_arithmetic.hpp"

using namespace std;
using namespace tbb;

__mps::__mps(double* a) {
    array = a;
    sum_interval = in2_create(0.,0.);
    mps_interval = in2_create(0.,0.);
    position = 0;
    size = 0;
}

__mps::__mps(__mps& x, split){
    array = x.array;
    sum_interval = in2_create(0.,0.);
    mps_interval = in2_create(0.,0.);
    position = 0;
    size = 0;
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
   // computing sum-l + mps-r
   __m128d mpsCandidate = in2_add(sum_interval,rightMps.mps_interval);
   // comparison of mpsCandidate and mps-l
   boolean b = inferior(mps_interval,mpsCandidate);
   if (b == True){
       mps_interval = mpsCandidate;
       position = rightMps.position;
   }
   else if(b == Undefined){
      cout << "Undefined comparison result in join operation, size: " << size << endl;
   }
   // adding two sums
   sum_interval = in2_add(sum_interval,rightMps.sum_interval);
}

void __mps::operator()(const blocked_range<int>& r){
    // iterating over the subrange
    for(int i = r.begin(); i != r.end(); i++){
        size++;
        sum_interval = in2_add_double(sum_interval,array[i]);
        boolean b = inferior(mps_interval,sum_interval);
        if (b == True){
            mps_interval = sum_interval;
            position = i; 
        }
        else if(b == Undefined){
            int l = 0;
        }
    }
    // Print debugging information
    //cout << "Processed subrange of length: " << size << endl;
}
        

void parallel_ia_mps(double* array, int size){
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
    // Computing sum and mps
    __mps result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    result.print_mps();
}
