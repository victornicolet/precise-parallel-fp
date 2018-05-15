/* Optimized functions for interval arithmetic 
 * Source: https://locklessinc.com/articles/interval_arithmetic/ 
 * */

#ifndef INTERVAL_ARITHMETIC_H
#define INTERVAL_ARITHMETIC_H

#include <emmintrin.h>

// Function to get the lower bound of an interval
double in2_min(__m128d);

// Function to get the upper bound of an interval
double in2_max(__m128d);

// Function to get the real value stored as the lower bound (opposite of the lower bound)
double in2_rmin(__m128d);

// Function to create an interval, first argument is lower bound, second argument is upper bound
__m128d in2_create(double,double);

// Same function, except that the first argument must be the opposite of the lower bound
__m128d in2_raw_create(double,double);

// Function to add to intervals. CAUTION: THE ROUNDING MODE MUST BE SET TOWARDS INFINITY, with  _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
__m128d in2_add(__m128d,__m128d);

// Functions below have been implemented by Raphael Dang-Nhu 

// Function to print a __m128d
void print(__m128d);

// datatype to return the result of a partial comparison
enum boolean{
    True,
    False,
    Undefined
};

// Function to add an interval and a double
__m128d in2_add_double(__m128d, double);

// Function to compare two intervals. Return True if a is inferior or equal to b, False if a is strictly superior to b, Undefined if the two intervals intersect in more than one point
boolean inferior(__m128d a,__m128d b);

// Function to compare an interval and a double. Return True if a is inferior or equal to b, False if a is strictly superior to b, undefined is a belongs to the interval b
boolean inferior_double(double a, __m128d b);

// Function to merge two intervals in case of undefined comparison
__m128d in2_max(__m128d,__m128d);

// Function to merge two intervals in case of undefined comparison
__m128d in2_min(__m128d,__m128d);

#endif
