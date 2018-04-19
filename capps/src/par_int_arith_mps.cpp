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

__mps::__mps() {
    sum_interval = in2_create(0.,0.);
    mps_interval = in2_create(0.,0.);
    position = 0;
    size = 0;
}

__mps::__mps(__mps& r, split){
    sum_interval = in2_create(0.,0.);
    mps_interval = in2_create(0.,0.);
    position = 0;
    size = 0;
}


void __mps::print_mps(){
    cout << endl << "sum: ";
    print(sum_interval);
    cout << endl << "mps: ";
    print(mps_interval);
    cout << endl << "position: " << position << endl;
}


void __mps::join(__mps& rightMps){
   // computing sum-l + mps-r
   __m128d mpsCandidate = in2_add(sum_interval,rightMps.mps_interval);
   // comparison of mpsCandidate and mps-l
   boolean b = inferior(mps_interval,mpsCandidate);
   if (b == True){
       mps_interval = mpsCandidate;
       position = size + rightMps.position;
   }
   else if(b == Undefined){
      cout << "Undefined comparison result" << endl;
   }
   // adding two sums and sizes 
   size += rightMps.size;
   sum_interval = in2_add(sum_interval,rightMps.sum_interval);

}

void __mps::operator()(const blocked_range<double*>& r){
    size = r.end() - r.begin();
    // iterating over the subrange
    for(double* a = r.begin(); a != r.end(); a++){
        sum_interval = in2_add_double(sum_interval,*a);
        boolean b = inferior(mps_interval,sum_interval);
        if (b == True){
            mps_interval = sum_interval;
            position = a - r.begin(); 
        }
        else if(b == Undefined){
           cout << "Undefined comparison result" << endl;
        }
    }
}
        

__mps parallel_ia_mps(double* array, int size){
    __mps result;
    parallel_reduce(blocked_range<double*>(array,array+size),result);
    result.print_mps();
}
