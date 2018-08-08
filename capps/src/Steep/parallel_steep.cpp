
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

void sequential_steep(double* a, int size, double& s, double& c, bool& b){

    double sum = 0.;
    double aux;
    double capacity = 0.;
    bool bo = true;
    for(int i = 0; i != size; i++){
        aux = a[i]-sum;
        if(aux < capacity){
            capacity = aux;
        }
        sum += a[i];
        bo &= aux >= 0;
    }
    s = sum;
    c = capacity;
    b = bo;
}

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
        cout << endl;
    }
    init.terminate();
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
        double aux = array[i]-sum;
        if(aux < capacity){
            capacity = aux;
        }
        sum += array[i];
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
    array(x.array),
    sum(0.,1000),
    capacity(0.,1000),
    b(true)
{}

void __steep_mpfr::operator()(const blocked_range<long>& r){

    for(int i = r.begin(); i != r.end(); i++){
        mpreal aux = array[i]-sum;
        if(aux < capacity){
            capacity = aux;
        }
        sum += array[i];
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
    array(array_),
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
    __m128d aux;
    for(int i = r.begin(); i != r.end(); i++){

        aux = in2_sub(array[i],sum);
        dec = in2_ge(aux,0.);
        b = booleanAnd(b,dec); 
        
        dec = in2_lt(aux,capacity);
        switch(dec){
            case False:
                break;
            case True:
                capacity = aux;
                break;
            case Undefined:
                capacity = in2_merge(aux,capacity);
                break;
        }
        sum = in2_add(sum,array[i]);
    }

}

void __steep_interval::join(__steep_interval& right){

    __m128d aux = in2_add(right.capacity,sum);
    capacity = in2_min(capacity,aux);
    sum = in2_add(sum,right.sum);
    b = booleanAnd(b,in2_ge(aux,0.));
}

void __steep_interval::print_steep(){
    cout << "Sum: ";
    print(sum);
    cout << endl << "Capacity: ";
    print(capacity);
    cout << endl << "Boolean: ";
    print(b);
}

