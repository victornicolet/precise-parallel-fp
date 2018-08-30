/* Parallel steep implementations.
 * Author: Raphaël Dang-Nhu
 * Date: 11/05/18 */

#ifndef PARALLEL_MPS_H
#define PARALLEL_MPS_H

#include "interval_arithmetic.hpp"

#include "tbb/blocked_range.h"
#include <mpreal.h>
#include <mpfr.h>
#include <gmp.h>

using mpfr::mpreal;
using namespace tbb;
using namespace std;

void sequential_steep(double* a, int size,double&,double&,bool&);

// Main function for naive parallel steep
void parallel_steep_double(double* array, long size);

// Main function for parallel steep with intervals
void parallel_steep_interval(double* array, long size);

// Main function for parallel steep with intervals
void parallel_steep_hybrid(double* array, long size, int maxDepth);

// Main function for parallel steep with intervals
boolean*** parallel_steep_hybrid_interval(double* array, long size, int maxDepth);

void parallel_steep_hybrid_exact(double* array, long size, int maxDepth, boolean***);

void parallel_steep_hybrid_lazy(double* array, long size, int maxDepth, double&,double&);

// Naive steep structure
struct __steep_naive{
    double* array;
    double sum;
    double capacity;
    bool b;
    
    // Constructor
    __steep_naive(double* a);
    // Splitting constructor
    __steep_naive(__steep_naive&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<long>&);
    // Join operation for the reduction
    void join(__steep_naive& rightsteep); 
    // Printing function
    void print_steep();
};

// Mpfr steep structure
struct __steep_mpfr{
    double* array;
    mpreal sum;
    mpreal capacity;
    bool b;
    
    // Constructor
    __steep_mpfr(double* a);
    // Splitting constructor
    __steep_mpfr(__steep_mpfr&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<long>&);
    // Join operation for the reduction
    void join(__steep_mpfr& rightsteep); 
    // Printing function
    void print_steep();
};

// steep structure with interval arithmetic
struct __steep_interval{
    double* array;
    __m128d sum;
    __m128d capacity;
    boolean b;

    __steep_interval(double* a);

    __steep_interval(__steep_interval& x, split s);

    void operator()(const blocked_range<long>&);
    
    void join(__steep_interval& right);

    void print_steep();
};

#endif
