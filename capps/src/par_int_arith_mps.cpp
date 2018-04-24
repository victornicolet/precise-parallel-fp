/* Implementation of lazy parallel maximum prefix sum with interval arithmetic and superaccumulators.
 * Author: Raphael Dang-Nhu 
 * Date: 04/19/2018 */

#include <iostream>
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"

#include "par_int_arith_mps.hpp"
#include "interval_arithmetic.hpp"

using namespace std;
using namespace tbb;


__mps_naive::__mps_naive(double* a):
    array(a),
    position(-1),
    sum(0),
    mps(0)

{}

__mps_naive::__mps_naive(__mps_naive& x, split) :
    array(x.array),
    position(-1),
    sum(0),
    mps(0)
{}
void __mps_naive::print_mps(){
    cout << "sum: " << sum;
    cout << endl << "mps: " << mps;
    cout << endl << "position: " << position << endl;
}

void __mps_naive::operator()(const blocked_range<int>& r){
    if(position == -1){
        position = r.begin();
    }
    for(int i = r.begin(); i != r.end(); i++){
        sum += array[i];
        if(sum> mps){
            mps = sum;
            position = i; 
        }
    }
}

void __mps_naive::join(__mps_naive& rightMps){
   // computing sum-l + mps-r
   rightMps.mps += sum;
   // adding two sums
   sum += rightMps.sum;
   // comparison of mpsCandidate and mps-l
   if(rightMps.mps > mps){
       mps = rightMps.mps;
       position = rightMps.position;
   }
}
__mps::__mps(double* a) :
    array(a),
    size(0),
    position(-1),
    left(-1),
    right(-1)
{
    sum_interval = in2_create(0.,0.);
    mps_interval = in2_create(0.,0.);
}

__mps::__mps(__mps& x, split) :
    array(x.array),
    size(0),
    position(-1),
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
       __mps_acc precise_result(array);
       parallel_reduce(blocked_range<int>(left,right),precise_result);
       double mpsaux = precise_result.mps.Round();
       mps_interval = in2_create(mpsaux,mpsaux);
       double sumaux = precise_result.sum.Round();
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
           __mps_acc precise_result(array);
           parallel_reduce(blocked_range<int>(left,right),precise_result);
           double mpsaux = precise_result.mps.Round();
           mps_interval = in2_create(mpsaux,mpsaux);
           double sumaux = precise_result.sum.Round();
           sum_interval = in2_create(sumaux,sumaux);
           position = precise_result.position;
           break;
        }
    }
}

__mps_acc::__mps_acc(double* a):
    array(a),
    position(-1)
{
    sum = Superaccumulator();
    neg_sum = Superaccumulator();
    mps = Superaccumulator();
    neg_mps = Superaccumulator();
}

__mps_acc::__mps_acc(__mps_acc& x, split) :
    array(x.array),
    position(-1)
{
    sum = Superaccumulator();
    neg_sum = Superaccumulator();
    mps = Superaccumulator();
    neg_mps = Superaccumulator();
}

void __mps_acc::print_mps(){
    cout << "sum: " << sum.Round();
    cout << endl << "mps: " << mps.Round();
    cout << endl << "position: " << position << endl;
}

void __mps_acc::operator()(const blocked_range<int>& r){
    if(position == -1){
        position = r.begin();
    }
    for(int i = r.begin(); i != r.end(); i++){
        //cout << endl << array[i] << endl;
        sum.Accumulate(array[i]);
        neg_sum.Accumulate(-array[i]);
        // Trick to do the comparison
        Superaccumulator test = Superaccumulator();
        test.Accumulate(sum);
        test.Accumulate(neg_mps);
        test.Normalize();
        double testAux = test.Round();
        //cout << endl << testAux << endl;
        if(testAux >= 0){
            mps = Superaccumulator(sum.get_accumulator());
            neg_mps = Superaccumulator(neg_sum.get_accumulator());
            position = i; 
        }
    }
}

void __mps_acc::join(__mps_acc& rightMps){
   // computing sum-l + mps-r
   rightMps.mps.Accumulate(sum);
   rightMps.neg_mps.Accumulate(neg_sum);
   // adding two sums
   sum.Accumulate(rightMps.sum);
   neg_sum.Accumulate(rightMps.neg_sum);
   // comparison of mpsCandidate and mps-l
   Superaccumulator test = Superaccumulator();
   test.Accumulate(rightMps.mps);
   test.Accumulate(neg_mps);
   if(test.Round() >= 0){
       mps = rightMps.mps;
       neg_mps = rightMps.neg_mps;
       position = rightMps.position;
   }
}

void parallel_mps(double* array, int size){
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
    // Computing sum and mp__mps_acc result(array);
    __mps_naive result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    result.print_mps();
}
void parallel_superacc_mps(double* array, int size){
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
    __mps_acc result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    result.print_mps();
}

void parallel_lazy_superacc_mps(double* array, int size){
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
