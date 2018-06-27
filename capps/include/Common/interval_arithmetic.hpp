/* Optimized functions for interval arithmetic 
 * Source: https://locklessinc.com/articles/interval_arithmetic/ 
 * */

#ifndef INTERVAL_ARITHMETIC_H
#define INTERVAL_ARITHMETIC_H

#include <emmintrin.h>
#include <math.h>

using namespace std;

// Function to get the lower bound of an interval
double in2_min(__m128d);

// Function to get the upper bound of an interval
double in2_max(__m128d);

// Function to get the real value stored as the lower bound (opposite of the lower bound)
double in2_rmin(__m128d);

// Function to create an interval, first argument is lower bound, second argument is upper bound
inline __m128d in2_create(double,double);

inline __m128d in2_create(double min, double max)
{
        return  (__m128d){-min, max};
}

// Same function, except that the first argument must be the opposite of the lower bound
__m128d in2_raw_create(double,double);

// Function to add to intervals. CAUTION: THE ROUNDING MODE MUST BE SET TOWARDS INFINITY, with  _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
inline __m128d in2_add(__m128d,__m128d);

inline __m128d in2_add(__m128d x, __m128d y)
{
        //if(_MM_GET_ROUNDING_MODE() == 0) cout << endl <<"false rm for ia" <<endl; 
        /* Component-wise addition */
        return  x + y;
}

// Functions below have been implemented by Raphael Dang-Nhu 

// Function to print a __m128d
void print(__m128d);

// datatype to return the result of a partial comparison
enum boolean{
    True,
    False,
    Undefined,
    Useless
};

struct memo{
    bool useful1;
    boolean useful2;
};

// Function to add an interval and a double
inline __m128d in2_add_double(__m128d, double);


inline __m128d in2_add_double(__m128d x, double y){
    //__m128d y0 = in2_create(y,y);
    //return in2_add(x,y0);
    return x + (__m128d){-y,y};
}

// Function to compare two intervals. Return True if a is inferior or equal to b, False if a is strictly superior to b, Undefined if the two intervals intersect in more than one point
inline boolean inferior(__m128d a,__m128d b);

inline boolean inferior(__m128d a, __m128d b){
    if(a[1] <= -b[0]) return True;
    // Should it be equal or inferior or equal there ?
    else if (-a[0] > b[1]) return False;
    else return Undefined;
}

// Function to compare an interval and a double. Return True if a is inferior or equal to b, False if a is strictly superior to b, undefined is a belongs to the interval b
inline boolean inferior_double(double a, __m128d b);

inline boolean inferior_double(double a, __m128d b){
    if(a <= -b[0]) return True;
    else if( a > b[1]) return False;
    else return Undefined;
}

// Function to merge two intervals in case of undefined comparison
inline __m128d in2_max(__m128d,__m128d);

inline __m128d in2_max(__m128d a,__m128d b){
    return (__m128d){min(a[0],b[0]),max(a[1],b[1])};
}

// Function to merge two intervals in case of undefined comparison
inline __m128d in2_min(__m128d,__m128d);

inline __m128d in2_min(__m128d a,__m128d b){
    return (__m128d){max(a[0],b[0]),min(a[1],b[1])};
}

#endif
