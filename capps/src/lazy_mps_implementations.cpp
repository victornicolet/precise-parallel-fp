/* Implementation of lazy parallel maximum prefix sum with interval arithmetic and superaccumulators.
 * Author: Raphael Dang-Nhu 
 * Date: 04/19/2018 */

#include <iostream>
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/task_group.h"

#include "precise_mps_implementations.hpp"
#include "lazy_mps_implementations.hpp"
#include "2_lazy_mps_implementations.hpp"
#include "interval_arithmetic.hpp"

using namespace std;
using namespace tbb;

void printA(double* array, int size){
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

}

void parallel_mps_float(double* array, int size){
    //printA(array,size);
    __mps_naive result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    //result.print_mps();
}

void parallel_mps_Collange(double* array, int size){
    //printA(array,size);
    __mps_precise<4> result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    //result.print_mps();
}

void parallel_mps_superacc(double* array, int size){
    //printA(array,size);
    __mps_acc result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    //result.print_mps();
}

void parallel_mps_superacc_lazy(double* array, int size){
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    //printA(array,size);
    __mps<__mps_acc> result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    //result.print_mps();
    _MM_SET_ROUNDING_MODE(0);
}

void parallel_mps_mpfr(double* array, int size){
    //printA(array,size);
    __mps_mpfr result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    //result.print_mps();
} 

void parallel_mps_mpfr_lazy(double* array, int size){
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    //printA(array,size);
    __mps<__mps_mpfr> result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    //result.print_mps();
    _MM_SET_ROUNDING_MODE(0);
}

void parallel_mps_mpfr_lazy_2(double* array, int size){
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    //printA(array,size);
    // Create containers for result
    __m128d sum_interval = in2_create(0.,0.);
    __m128d mps_interval = in2_create(0.,0.);
    int position = 0;

    MpsTask1& root = *new(task::allocate_root()) MpsTask1(10000,array,size,&sum_interval,&mps_interval,&position);

    task::spawn_root_and_wait(root);

    
    // Printing result
    /*cout << "Sum: ";
    print(sum_interval) ;
    cout << endl << "Mps: ";
    print(mps_interval) ;
    cout << endl << "Position: " << position << endl;
    */

    
}

