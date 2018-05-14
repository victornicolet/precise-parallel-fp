/* Author: Raphael Dang-Nhu */

#include <iostream>
#include <gmp.h>
#include <mpfr.h>

#include "tbb/blocked_range.h"
#include "superaccumulator.hpp"
#include "parallel_mps.hpp"

using namespace std;

__mps_naive::__mps_naive(double* a):
    array(a),
    sum(0),
    mps(0),
    position(0)
{}

__mps_naive::__mps_naive(__mps_naive& x, split) :
    array(x.array),
    sum(0),
    mps(0),
    position(0)
{}

void __mps_naive::print_mps(){
    cout << "sum: " << sum;
    cout << endl << "mps: " << mps;
    cout << endl << "position: " << position << endl;
}

void __mps_naive::operator()(const blocked_range<int>& r){
    if(position == 0){
        position = r.begin();
    }
    for(int i = r.begin(); i != r.end(); i++){
        sum += array[i];
        if(sum >= mps){
            mps = sum;
            position = i+1; 
        }
    }
}

void __mps_naive::join(__mps_naive& rightMps){
   // computing sum-l + mps-r
   rightMps.mps += sum;
   // adding two sums
   sum += rightMps.sum;
   // comparison of mpsCandidate and mps-l
   if(rightMps.mps >= mps){
       mps = rightMps.mps;
       position = rightMps.position;
   }
}

__mps_acc::__mps_acc(double* a):
    array(a),
    position(0)
{
    sum = Superaccumulator();
    mps = Superaccumulator();
}

__mps_acc::__mps_acc(__mps_acc& x, split) :
    array(x.array),
    position(0)
{
    sum = Superaccumulator();
    mps = Superaccumulator();
}

void __mps_acc::print_mps(){
    cout << "sum: " << sum.Round();
    cout << endl << "mps: " << mps.Round();
    cout << endl << "position: " << position << endl;
}

void __mps_acc::operator()(const blocked_range<int>& r){
    _MM_SET_ROUNDING_MODE(0);
    if(position == 0){
        position = r.begin();
    }
    for(int i = r.begin(); i != r.end(); i++){
        sum.Accumulate(array[i]);
        if(!sum.comp(mps)){
            mps = Superaccumulator(sum.get_accumulator());
            position = i+1; 
        }
    }
}

void __mps_acc::join(__mps_acc& rightMps){
   //Set rounding mode
   _MM_SET_ROUNDING_MODE(0);
   // computing sum-l + mps-r
   rightMps.mps.Accumulate(sum);
   // adding two sums
   sum.Accumulate(rightMps.sum);
   // comparison of mpsCandidate and mps-l
   if(!rightMps.mps.comp(mps)){
       mps = rightMps.mps;
       position = rightMps.position;
   }
}

double __mps_acc::getSumDoubleUp(){
    return sum.Round();
}

double __mps_acc::getMpsDoubleUp(){
    return mps.Round();
}


double __mps_acc::getSumDoubleDown(){
    return sum.Round();
}

double __mps_acc::getMpsDoubleDown(){
    return mps.Round();
}
__mps_mpfr::__mps_mpfr(double* array) :
    array(array),
    position(0)
{
    mpfr_init2(sum,30000);
    mpfr_init2(mps,30000);
    mpfr_set_d(sum,0.,MPFR_RNDN);
    mpfr_set_d(mps,0.,MPFR_RNDN);
}

__mps_mpfr::__mps_mpfr(__mps_mpfr& x,split) :
    array(x.array),
    position(0)
{
    mpfr_init2(sum,30000);
    mpfr_init2(mps,30000);
    mpfr_set_d(sum,0.,MPFR_RNDN);
    mpfr_set_d(mps,0.,MPFR_RNDN);
}

void __mps_mpfr::operator()(const blocked_range<long>& range){
    if(position == 0){
        position = range.begin();
    }
    for(long i = range.begin(); i != range.end(); i++){
        mpfr_add_d(sum,sum,array[i],MPFR_RNDN);
        if(mpfr_cmp(sum,mps)>=0){
            mpfr_set(mps,sum,MPFR_RNDN);
            position = i+1;
        }
    }
}

void __mps_mpfr::join(__mps_mpfr& right){
    mpfr_add(right.mps,right.mps,sum,MPFR_RNDN);
    mpfr_add(sum,right.sum,sum,MPFR_RNDN);
    if(mpfr_cmp(right.mps,mps) >= 0){
            mpfr_set(mps,right.mps,MPFR_RNDN);
            position = right.position;
    }
    mpfr_clear(right.mps);
    mpfr_clear(right.sum);
}

void __mps_mpfr::print_mps(){
    cout << "sum: " << getSumDoubleUp();
    cout << endl << "mps: " << getMpsDoubleUp();
    cout << endl << "position: " << position << endl;
}

double __mps_mpfr::getSumDoubleDown(){
    return mpfr_get_d(sum,MPFR_RNDD);
}

double __mps_mpfr::getMpsDoubleDown(){
    return mpfr_get_d(mps,MPFR_RNDD);
}

double __mps_mpfr::getSumDoubleUp(){
    return mpfr_get_d(sum,MPFR_RNDU);
}

double __mps_mpfr::getMpsDoubleUp(){
    return mpfr_get_d(mps,MPFR_RNDU);
}


