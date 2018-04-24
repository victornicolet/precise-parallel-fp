/* Implementation of lasy parallel maximum prefix sum with interval arithmetic and superaccumulators.
 * Author: Raphael Dang-Nhu 
 * Date: 04/19/2018 */

#include <emmintrin.h>
#include "tbb/blocked_range.h"
#include "superaccumulator.hpp"

using namespace tbb;

// Naive mps structure. 
struct __mps_naive{
    // pointer to the array
    double* array;
    // Superaccumulators for sum and mps
    double sum;
    double mps;
    // Position of the maximum prefix sum
    int position;
    // Constructor
    __mps_naive(double* a);
    // Splitting constructor
    __mps_naive(__mps_naive&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<int>&);
    // Join operation for the reduction
    void join(__mps_naive& rightMps); 
    // Printing function
    void print_mps();
};

// Same structure with lazy computation, i.a. and superaccumulators. 
struct __mps{
    // pointer to the array
    double* array;
    // Both sum and maximum prefix sum are represented as intervals  
    __m128d sum_interval;
    __m128d mps_interval;
    // Position of the maximum prefix sum
    int position;
    // size
    int size;
    int left;
    int right;
    // Constructor
    __mps(double* a);
    // Splitting constructor
    __mps(__mps&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<int>&);
    // Join operation for the reduction
    void join(__mps& rightMps); 
    // Printing function
    void print_mps();
};

// Same structure with super accumulators. 
struct __mps_acc{
    // pointer to the array
    double* array;
    // Superaccumulators for sum and mps
    Superaccumulator sum;
    Superaccumulator mps;
    Superaccumulator neg_sum;
    Superaccumulator neg_mps;
    // Position of the maximum prefix sum
    int position;
    // Constructor
    __mps_acc(double* a);
    // Splitting constructor
    __mps_acc(__mps_acc&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<int>&);
    // Join operation for the reduction
    void join(__mps_acc& rightMps); 
    // Printing function
    void print_mps();
};

/* Naive parallel mps function */
void parallel_mps(double*,int);

/* Function computing in parallel the maximum prefix sum in a lazy way, with interval arithmetic and superaccumulators.
 * Arguments: input array, size of the array
 */
void parallel_lazy_superacc_mps(double*,int);

/* Function computing in parallel the mps with superaccumulators, same arguments as previous function */
void parallel_superacc_mps(double*,int);
