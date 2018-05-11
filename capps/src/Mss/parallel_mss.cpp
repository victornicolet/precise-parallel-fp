
#include <iostream>

#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"

#include "parallel_mss.hpp"

#define PRINT 1

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

void parallel_mss_double(double* array, long size){
    __mss_naive result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    if(PRINT){
        cout << endl << "Parallel double" << endl;
        printA(array,size);
        result.print_mss();
    }
}

__mss_naive::__mss_naive(double* a) :
    array(a),
    sum(0.),
    mps(0.),
    mts(0.),
    mss(0.),
    posmps(0),
    posmts(0),
    posmssl(0),
    posmssr(0)
{}

__mss_naive::__mss_naive(__mss_naive& x, split) :
    array(x.array),
    sum(0.),
    mps(0.),
    mts(0.),
    mss(0.),
    posmps(0),
    posmts(0),
    posmssl(0),
    posmssr(0)
{}

void __mss_naive::operator()(const blocked_range<int>& r){
    if(posmps == 0){
        posmps = r.begin();
    }
    if(posmts == 0){
        posmts = r.begin();
    }
    if(posmssl == 0){
        posmssl = r.begin();
    }
    if(posmssr == 0){
        posmssr = r.begin();
    }

    for(int i = r.begin(); i != r.end(); i++){

        // Update sum
        sum += array[i];
        
        // Update mts
        mts += array[i];
        if(mts < 0){
            mts = 0.;
            posmts = i+1;
        }
        
        // Update mps
        if(sum >= mps){
            mps = sum;
            posmps = i+1;
        }

        // Update mss
        if(mts >= mss){
            mss = mts;
            posmssl = posmts;
            posmssr = i+1;
        }
    }
}

void __mss_naive::join(__mss_naive& right){

    // Join mss
    double aux = mts + right.mps;
    bool test1 = aux >= mss;
    bool test2 = right.mss >= aux;
    bool test3 = right.mss >= mss;

    if(test2 && test3){
        mss = right.mss;
        posmssl = right.posmssl;
        posmssr = right.posmssr;
    }
    else if (test1 && !test2){
        mss = aux;
        posmssl = posmts;
        posmssr = right.posmps;
    }
    
    // Join mps
    right.mps += sum;
    if(right.mps >= mps){
        mps = right.mps;
        posmps = right.posmps;
    }

    // Join mts
    mts += right.sum;
    if(right.mts >= mts){
        mts = right.mts;
        posmts = right.posmts;
    }
}


