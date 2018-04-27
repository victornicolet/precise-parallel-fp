/* Author: Raphael Dang-Nhu */

#include <iostream>
#include <gmp.h>
#include <mpfr.h>

#include "tbb/blocked_range.h"
#include "superaccumulator.hpp"
#include "precise_mps_implementations.hpp"

using namespace std;

__mps_naive::__mps_naive(double* a):
    array(a),
    sum(0),
    mps(0),
    position(-1)
{}

__mps_naive::__mps_naive(__mps_naive& x, split) :
    array(x.array),
    sum(0),
    mps(0),
    position(-1)
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
        //Set rounding mode
        _MM_SET_ROUNDING_MODE(0);

        sum.Accumulate(array[i]);
        neg_sum.Accumulate(-array[i]);
        // Trick to do the comparison
        Superaccumulator test = Superaccumulator();
        test.Accumulate(sum);
        test.Accumulate(neg_mps);
        double testAux = test.Round();
        //cout << endl << testAux << endl;
        if(testAux >= 0){
            mps = Superaccumulator(sum.get_accumulator());
            neg_mps = Superaccumulator(neg_sum.get_accumulator());
            position = i; 
        }
        //Set back rounding mode
        _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    }
}

void __mps_acc::join(__mps_acc& rightMps){
    //Set rounding mode
    _MM_SET_ROUNDING_MODE(0);
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
    //Set back rounding mode
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
}

double __mps_acc::getSumDouble(){
    return sum.Round();
}

double __mps_acc::getMpsDouble(){
    return sum.Round();
}

__mps_mpfr::__mps_mpfr(double* array) :
    array(array),
    position(-1)
{
    mpfr_init2(sum,MPFR_PREC_MAX);
    mpfr_init2(mps,MPFR_PREC_MAX);
    mpfr_set_d(sum,0.,MPFR_RNDN);
    mpfr_set_d(sum,0.,MPFR_RNDN);
}

__mps_mpfr::__mps_mpfr(__mps_mpfr& x,split) :
    array(x.array),
    position(-1)
{
    mpfr_init2(sum,MPFR_PREC_MAX);
    mpfr_init2(mps,MPFR_PREC_MAX);
    mpfr_set_d(sum,0.,MPFR_RNDN);
    mpfr_set_d(sum,0.,MPFR_RNDN);
}

void __mps_mpfr::operator()(const blocked_range<int>& range){
    if(position == -1){
        position = range.begin();
    }
    for(int i = range.begin(); i != range.end(); i++){
        mpfr_add_d(sum,sum,array[i],MPFR_RNDN);
        if(mpfr_cmp(sum,mps)>=0){
            mpfr_set(mps,sum,MPFR_RNDN);
            position = i;
        }
    }
}

void __mps_mpfr::join(__mps_mpfr& right){
    mpfr_add(right.sum,right.sum,sum,MPFR_RNDN);
    mpfr_add(right.mps,right.mps,sum,MPFR_RNDN);
    if(mpfr_cmp(right.mps,mps) >= 0){
            mpfr_set(mps,right.mps,MPFR_RNDN);
            position = right.position;
    }
    mpfr_clear(right.mps);
    mpfr_clear(right.sum);
}

void __mps_mpfr::print_mps(){
    cout << "sum: ";
    mpfr_out_str(NULL,10,0,sum,MPFR_RNDN);
    cout << endl << "mps: " ;
    mpfr_out_str(NULL,10,0,mps,MPFR_RNDN);
    cout << endl << "position: " << position << endl;
}

double __mps_mpfr::getSumDouble(){
    return mpfr_get_d(sum,MPFR_RNDN);
}

double __mps_mpfr::getMpsDouble(){
    return mpfr_get_d(mps,MPFR_RNDN);
}
