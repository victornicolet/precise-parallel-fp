/* File for small parallel reduce tests.
 * Author: Raphaël Dang-Nhu.
 * Date: May 9th */

#include <iostream>

#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/atomic.h"
#include "tbb/mutex.h"

using namespace std;
using namespace tbb;

/* Lock for print function */
mutex countMutex;

// Naive mps structure
struct __mps_naive{
    // pointer to the array
    double* array;
    // Storing which body processed the element
    int* memo;
    // sum
    double sum;
    // Pointer to global counter
    tbb::atomic<int>* counter; 
    // Index
    int index;
    // Constructor
    __mps_naive(double* a,int* memo,tbb::atomic<int>* c);
    // Splitting constructor
    __mps_naive(__mps_naive&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<int>&);
    // Join operation for the reduction
    void join(__mps_naive& rightMps); 
    // Printing function
    void print_mps();
};

__mps_naive::__mps_naive(double* a,int* m,tbb::atomic<int>* c):
    array(a),
    memo(m),
    sum(0),
    counter(c)
{
    index = counter->fetch_and_add(1);
}

__mps_naive::__mps_naive(__mps_naive& x, split) :
    array(x.array),
    memo(x.memo),
    counter(x.counter),
    sum(0)
{
    index = counter->fetch_and_add(1);
}

void __mps_naive::print_mps(){
    cout << "sum: " << sum;
}

void __mps_naive::operator()(const blocked_range<int>& r){
    for(int i = r.begin(); i != r.end(); i++){

        //cout << endl << "Body " << index << " processed element of index " << i;
        sum += array[i];
        
        // Useless operation to make things longer
        int aux = 0;
        for(int i = 0; i != 1000; i++){
            aux = aux + sum;
        }
        array[i] = aux;

        // Memorize the body 
        memo[i] = index;
    }
    // Special memo for the first index
    memo[r.begin()] = -index;
}

void __mps_naive::join(__mps_naive& rightMps){
    countMutex.lock();
    cout << endl << "Body " << index << " joined Body " << rightMps.index ;
    countMutex.unlock();
   sum += rightMps.sum;
}

void test(){
    
    // Initialization
    const int size = 1000;
    double* array = new double[size];
    
    // Put random elements in the array
    for(int i = 0; i != size; i++){
        array[i] = (double)(rand()%100);
    }

    int* memo = new int[size];
    tbb::atomic<int> counter = 1;
    __mps_naive result(array,memo,&counter);

    parallel_reduce(blocked_range<int>(0,size),result);

    // Printing memo
    cout << endl << endl;
    for(int i = 0; i != size ; i++){
        cout << memo[i] << ",";
    }
    cout << endl;
}
/*
int main(){
    // Init task scheduler
    task_scheduler_init init(3);

    for(int i = 0; i != 10; i++) 
    test();

    return 0;

}*/
