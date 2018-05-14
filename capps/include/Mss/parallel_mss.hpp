/* Parallel mss implementations.
 * Author: Raphaël Dang-Nhu
 * Date: 11/05/18 */

#ifndef PARALLEL_MPS_H
#define PARALLEL_MPS_H
#include "interval_arithmetic.hpp"

#include "tbb/blocked_range.h"

using namespace tbb;
using namespace std;

// Main function for naive parallel mss
void parallel_mss_double(double* array, long size);

// Main function for parallel mss with intervals
void parallel_mss_interval(double* array, long size);

// Naive mss structure
struct __mss_naive{
    // pointer to the array
    double* array;
    // Superaccumulators for sum and mss
    double sum;
    double mps;
    double mts;
    double mss;
    // Position of the maximum prefix sum
    long posmps;
    long posmts;
    long posmssl;
    long posmssr;
    
    // Constructor
    __mss_naive(double* a);
    // Splitting constructor
    __mss_naive(__mss_naive&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<long>&);
    // Join operation for the reduction
    void join(__mss_naive& rightmss); 
    // Printing function
    void print_mss();
};

// Mss structure with interval arithmetic
struct __mss_interval{
    double* array;

    __m128d sum;
    __m128d mss;
    __m128d mps;
    __m128d mts;

    long posmssl1;
    long posmssl2;
    long posmssr1;
    long posmssr2;
    long posmts1;
    long posmts2;
    long posmps1;
    long posmps2;

    __mss_interval(double* a);

    __mss_interval(__mss_interval& x, split s);

    void operator()(const blocked_range<long>&);
    
    void join(__mss_interval& right);

    void print_mss();
};

#endif
