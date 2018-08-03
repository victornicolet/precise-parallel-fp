
#include <iostream>
#include <algorithm>

#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/task_group.h"
#include "tbb/task_scheduler_init.h"

#include "parallel_steep.hpp"
#include "parallel_mps.hpp"
#include "sequential_mts.hpp"

#include "debug.hpp"

using namespace tbb;
using namespace std;


void printA(double* array, int size){
    // printing the array
    cout << "{";
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

void parallel_steep_double(double* array, long size){
    task_scheduler_init init;
    __steep_naive result(array);
    parallel_reduce(blocked_range<long>(0,size),result);
    if(PRINT){
        cout << endl << "Parallel double" << endl;
        printA(array,size);
        result.print_steep();
    }
    init.terminate();
}

void parallel_steep_interval(double* array, long size){
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    task_scheduler_init init;
    __steep_interval result(array);
    parallel_reduce(blocked_range<long>(0,size),result);
    init.terminate();
    _MM_SET_ROUNDING_MODE(0);

    // If non precise result redo computation
    long posl = -1;
    long posr = -1;
    if(result.possteepl1 != result.possteepl2 && result.possteepr1 != result.possteepr2){

       // First case : possteepl2 <= possteepr1
       if(result.possteepl2 <= result.possteepr1){
           task_scheduler_init init2(2);
           __mps_mpfr resultmps(array);
                       
           parallel_reduce(blocked_range<long>(result.possteepr1,result.possteepr2),resultmps);

           posr = resultmps.position;
           
           // Creation of parameters
           long length = result.possteepl2 - result.possteepl1;
           double* a = new double[length];
           for(long i = 0; i != length; i++){
               a[i] = array[i+result.possteepl1];
           }
           double mts;
           sequential_mts_superacc(a,length,&mts,&posl);
           posl += result.possteepl1;
        }
       // Second case : possteepl2 > possteepr1
       else {
           cout << endl <<  "Recomputation not yet implemented, because super rare" << endl;
        }
    }
    else if(result.possteepr1 != result.possteepr2){
           task_scheduler_init init2(2);
           __mps_mpfr resultmps(array);
                       
           parallel_reduce(blocked_range<long>(result.possteepr1,result.possteepr2),resultmps);

           posr = resultmps.position;
    }
    else if(result.possteepl1 != result.possteepl2){
           // Creation of parameters
           long length = result.possteepl2 - result.possteepl1;
           double* a = new double[length];
           for(long i = 0; i != length; i++){
               a[i] = array[i+result.possteepl1];
           }
           double mts;
           sequential_mts_superacc(a,length,&mts,&posr);
           posl += result.possteepl1;
    }
    if(PRINT){
        cout << endl << "Dynamic lazy computation first results" << endl;
        printA(array,size);
        result.print_steep();
        cout << endl << "Precise results" << endl;
        cout << "Left pos: " << posl << endl;
        cout << "Right pos: " << posr << endl;
    }
}

__steep_naive::__steep_naive(double* a) :
    array(a),
    sum(0.),
    capacity(0.),
    b(true)
{}

__steep_naive::__steep_naive(__steep_naive& x, split) :
    array(x.array),
    sum(0.),
    capacity(0.),
    b(true)
{}

void __steep_naive::operator()(const blocked_range<long>& r){

    for(int i = r.begin(); i != r.end(); i++){
        double aux = a[i]-sum;
        if(aux < capacity){
            capacity = aux;
        }
        sum += a[i];
        b &= capacity >= 0;
    }
}

void __steep_naive::join(__steep_naive& right){

    capacity = min(capacity,right.capacity - sum);
    sum += right.sum;
    b &= capacity >= 0;

}

void __steep_naive::print_steep(){
    cout << "Sum: " << sum << endl;
    cout << "Capacity: " << capacity << endl;
    cout << "Boolean: " << b;
}

__steep_mpfr::__steep_mpfr(double* a) :
    array(a),
    sum(0.,1000),
    capacity(0.,1000),
    b(true)
{}

__steep_mpfr::__steep_mpfr(__steep_mpfr& x, split) :
    array(x.a),
    sum(0.,1000),
    capacity(0.,1000),
    b(true)
{}

void __steep_mpfr::operator()(const blocked_range<long>& r){

    for(int i = r.begin(); i != r.end(); i++){
        mpreal aux = a[i]-sum;
        if(aux < capacity){
            capacity = aux;
        }
        sum += a[i];
        b &= capacity >= 0;
    }

}

void __steep_mpfr::join(__steep_mpfr& right){

    mpreal aux = right.capacity - sum;
    capacity = min(capacity,aux);
    sum += right.sum;
    b &= aux >= 0;
}

void __steep_mpfr::print_steep(){
    cout << "Sum: " << sum << endl;
    cout << "Capacity: " << capacity << endl;
    cout << "Boolean: " << b << endl;
}

__steep_interval::__steep_interval(double* array_) :
    array(a),
    b(True)
{
    sum = in2_create(0.);
    capacity = in2_create(0.);
}

__steep_interval::__steep_interval(__steep_interval& x ,split s) :
    array(x.array),
    b(True)
{
    sum = in2_create(0.);
    capacity = in2_create(0.);
}

void __steep_interval::operator()(const blocked_range<long>& r){

    boolean dec;
    for(int i = r.begin(); i != r.end(); i++){
        __m128d aux = in2_sub(a[i],sum);

        dec = in2_ge(aux,0.);
        b = booleanAnd(b,dec); 

        dec = in2_lt(aux,capacity);
        switch(dec){
            case True:
                capacity = aux;
                break;
            case Undefined:
                in2_merge(aux,capacity);
                break;
        }


        sum = in2_add(sum,a[i]);
    }

}

void __steep_interval::join(__steep_interval& right){

    __m128d aux = in2_sum(right.capacity,sum);
    capacity = in2_min(capacity,aux);
    sum = in2_add(sum,right.sum);
    
    b = booleanAnd(b,in2_ge(aux,0.));
}

void __steep_interval::print_steep(){
    cout << "Sum: ";
    print(sum);
    cout << endl << "Capacity: ";
    print(steep);
    cout << endl << "Boolean: ";
    print(b);
}

