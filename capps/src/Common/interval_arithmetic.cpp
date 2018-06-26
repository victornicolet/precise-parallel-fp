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
        //if(_MM_GET_ROUNDING_MODE() == 0) cout << endl <<"false rm for ia" <<endl; 
        /* Component-wise addition */
        return  x + y;
}

void print(__m128d x){
    cout << "[" << in2_min(x) << "," << in2_max(x) << "]" << ", length: " << in2_max(x) - in2_min(x);
}


__m128d in2_add_double(__m128d x, double y){
    //__m128d y0 = in2_create(y,y);
    //return in2_add(x,y0);
    return x + (__m128d){-y,y};
}

boolean inferior(__m128d a, __m128d b){
    if(a[1] <= -b[0]) return True;
    // Should it be equal or inferior or equal there ?
    else if (-a[0] > b[1]) return False;
    else return Undefined;
}

boolean inferior_double(double a, __m128d b){
    if(a <= -b[0]) return True;
    else if( a > b[1]) return False;
    else return Undefined;
}

__m128d in2_max(__m128d a,__m128d b){
    return (__m128d){min(a[0],b[0]),max(a[1],b[1])};
}

__m128d in2_min(__m128d a,__m128d b){
    return (__m128d){max(a[0],b[0]),min(a[1],b[1])};
}
