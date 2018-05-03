/* Implementation of lazy parallel maximum prefix sum with interval arithmetic and superaccumulators.
 * Author: Raphael Dang-Nhu 
 * Date: 04/19/2018 */

#include <iostream>
#include <math.h>
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/task_group.h"
#include "tbb/task_scheduler_init.h"

#include "precise_mps_implementations.hpp"
#include "lazy_mps_implementations.hpp"
#include "2_lazy_mps_implementations.hpp"
#include "interval_arithmetic.hpp"

using namespace std;
using namespace tbb;

#define PRINT 0

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

void sequential_mps(double* array, int size){
    double sum = 0., mps = 0.;
    int position = 0;
    if(PRINT){
        cout << endl << "Sequential mps" << endl;
        printA(array,size);
        cout <<  "sum: " << sum;
        cout << endl << "mps: " << mps;
        cout << endl << "pos: " << position << endl;
    }
    sequentialMps(array,size,&sum,&mps,&position);
}

void parallel_mps_float(double* array, int size){
    __mps_naive result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    if(PRINT){
        cout << endl << "Parallel double" << endl;
        printA(array,size);
        result.print_mps();
    }
}

void parallel_mps_Collange(double* array, int size){
    __mps_precise<4> result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    if(PRINT){
        cout << endl << "Parallel Collange" << endl;
        printA(array,size);
        result.print_mps();
    }
}

void parallel_mps_superacc(double* array, int size){
    __mps_acc result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    if(PRINT){
        cout << endl << "Parallel superacc" << endl;
        printA(array,size);
        result.print_mps();
    }
}

void parallel_mps_superacc_lazy(double* array, int size){
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    __mps<__mps_acc> result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    if(PRINT){
        cout << endl << "Parallel superacc lazy" << endl;
        printA(array,size);
        result.print_mps();
    }
}

void parallel_mps_mpfr(double* array, int size){
    __mps_mpfr result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    if(PRINT){
        cout << endl << "Parallel mpfr" << endl;
        printA(array,size);
        result.print_mps();
    }
} 

void parallel_mps_mpfr_lazy(double* array, int size){
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    __mps<__mps_mpfr> result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    _MM_SET_ROUNDING_MODE(0);
    if(PRINT){
        cout << endl << "Parallel mpfr lazy" << endl;
        printA(array,size);
        result.print_mps();
    }
}

void parallel_mps_mpfr_lazy_2(double* array, int size, int grainsize){

    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    
    // Parameters
    int Cutoff = grainsize;
    
    // Create variables for result
    __m128d sum_interval = in2_create(0.,0.);
    __m128d mps_interval = in2_create(0.,0.);
    int position = 0;

    // Create matrix to store the information
    double ratio = (double) (size)/(double)Cutoff;
    int ceilRatio = ceil(log2(ratio));
    const int maxDepth = ceilRatio + 1;
    
    int maxIndex = 1;
    Status** memo = new Status*[maxDepth]; 
    for(int i = 0; i < maxDepth;i++){
        memo[i] = new Status[maxIndex];
        maxIndex = 2*maxIndex;
    }

    int validity = 0;

    MpsTask1& root = *new(task::allocate_root()) MpsTask1(Cutoff,array,size,&validity,&sum_interval,&mps_interval,&position,memo);

    task::spawn_root_and_wait(root);
    
    if(PRINT){
        // Printing result
        cout << endl << "Parallel mpfr lazy 2"<< endl;
        printA(array,size);
            if(validity == 0){
                cout <<"Valid position" << endl;
            cout << "Sum: ";
            print(sum_interval) ;
            cout << endl << "Mps: ";
            print(mps_interval) ;
            cout << endl << "Position: " << position << endl;
        }else{
            cout << "Unvalid position" << endl;
            __mps_mpfr result(array);
            parallel_reduce(blocked_range<int>(0,size),result);
            result.print_mps();
        }

        // Printing memo
        cout << endl;
        maxIndex= 1;
        for(int i = 0; i!=maxDepth; i++){
            for(int j = 0; j!= maxIndex; j++){
                switch(memo[i][j]){
                    case undefinedComparison :
                        cout << "0";
                        break;
                    case cutoffPrecise : 
                        cout << "0";
                        break;
                    case leftChild :
                        cout << "|";
                        break;
                    case rightChild :
                        cout << "|";
                        break;
                    case cutoff :
                        cout << "|";
                        break;
                
                }
            }
            maxIndex = 2*maxIndex;
            cout << endl;
        }
    }
}

