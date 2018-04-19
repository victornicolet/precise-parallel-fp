/* Optimized functions for interval arithmetic 
 * Source: https://locklessinc.com/articles/interval_arithmetic/ 
 * */

#include <math.h>
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
