/* Optimized functions for interval arithmetic 
 * Source: https://locklessinc.com/articles/interval_arithmetic/ 
 * */

#include <iostream>
#include <math.h>
#include <emmintrin.h>
#include "interval_arithmetic.hpp"

using namespace std;

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
__m128d in2_create(double min, double max)
{
        return  (__m128d){-min, max};
}

__m128d in2_raw_create(double nmin, double max)
{
        return (__m128d){nmin, max};
}

__m128d in2_add(__m128d x, __m128d y)
{
        /* Component-wise addition */
        return  x + y;
}

void print(__m128d x){
    cout << "[" << in2_min(x) << "," << in2_max(x) << "]" << endl;
}


__m128d in2_add_double(__m128d x, double y){
    __m128d y0 = in2_create(y,y);
    return in2_add(x,y0);
}

boolean inferior(__m128d a, __m128d b){
    if(in2_max(a) <= in2_min(b)) return True;
    else if (in2_min(a) > in2_max(b)) return False;
    else return Undefined;
}
