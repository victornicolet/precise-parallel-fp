/* Optimized functions for interval arithmetic 
 * Source: https://locklessinc.com/articles/interval_arithmetic/ 
 * */

#include <iostream>
#include <emmintrin.h>
#include "interval_arithmetic.hpp"


double in2_min(__m128d x)
{
        return -x[0];
}

double in2_max(__m128d x)
{
        return x[1];
}

/* Get the real value stored in the SSE register */
double in2_rmin(__m128d x)
{
        return x[0];
}

__m128d in2_raw_create(double nmin, double max)
{
        return (__m128d){nmin, max};
}


void print(__m128d x){
    cout << "[" << in2_min(x) << "," << in2_max(x) << "]" << ", length: " << in2_max(x) - in2_min(x);
}



