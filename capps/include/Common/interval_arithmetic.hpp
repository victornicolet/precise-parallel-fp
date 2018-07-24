/* Optimized functions for interval arithmetic 
 * Source: https://locklessinc.com/articles/interval_arithmetic/ 
 * */

#ifndef INTERVAL_ARITHMETIC_H
#define INTERVAL_ARITHMETIC_H

#include <emmintrin.h>
#include <math.h>
#include "struct.hpp"

using namespace std;

// Function to get the lower bound of an interval
inline double in2_min(__m128d x)
{
        return -x[0];
}

// Function to get the upper bound of an interval
inline double in2_max(__m128d x)
{
        return x[1];
}

// Function to get the real value stored as the lower bound (opposite of the lower bound)
inline double in2_rmin(__m128d x)
{
        return x[0];
}

// Function to create an interval, first argument is lower bound, second argument is upper bound
inline __m128d in2_create(double min, double max)
{
        return  (__m128d){-min, max};
}

inline __m128d in2_create(double x)
{
        return  (__m128d){-x, x};
}

inline __m128d in2_create(int x)
{
        return  (__m128d){-(double)x,(double) x};
}

inline __m128i in2i_create(int x)
{
        return  (__m128i){-x,x};
}

inline __m128d in2_create(__m128d x)
{
        return  x;
}
// Same function, except that the first argument must be the opposite of the lower bound
inline __m128d in2_raw_create(double nmin, double max)
{
        return (__m128d){nmin, max};
}

// Function to multiply 
inline __m128d c_swap(__m128d x, __m128d cond)
{
	__m128d t = _mm_xor_pd((__m128d)_mm_shuffle_epi32((__m128i) x, 0x4e), x);
	
	cond = _mm_and_pd(cond, t);
	
	return _mm_xor_pd(cond, x);
}

inline __m128d in2_mul(__m128d x, __m128d y)
{
	__m128d t1 = (__m128d)_mm_shuffle_epi32((__m128i) x, 0xee);
	__m128d t2 = (__m128d)_mm_shuffle_epi32((__m128i) y, 0xee);
	
	__m128d t3 = _mm_xor_pd(x, t1);
	__m128d t4 = _mm_xor_pd(y, t2);
	
	if (_mm_movemask_pd(_mm_and_pd(t3, t4)))
	{
		__m128d c = {0.0, 0.0};
		__m128d c1 = _mm_cmple_pd(t2, c);
		__m128d c2 = _mm_cmple_pd(t1, c);
		__m128d c3 = (__m128d) {-0.0, 0.0};
		
		x = c_swap(_mm_xor_pd(x, c3), c1);
		y = c_swap(_mm_xor_pd(y, c3), c2);
		
		return x * _mm_xor_pd(y, c3);
	}
	
	/* There is a zero overlap */
	t1 = (__m128d)_mm_shuffle_epi32((__m128i) x, 0x4e) * _mm_unpacklo_pd(y, y);
	t2 *= x;

	return _mm_max_pd(t1, t2);
}
inline __m128d in2_mul(__m128d x, double y){
	
	return in2_mul(x,in2_create(y));

}

inline __m128d in2_mul(double x, __m128d y){
	
	return in2_mul(in2_create(x),y);

}
// Function to add to intervals. CAUTION: THE ROUNDING MODE MUST BE SET TOWARDS INFINITY, with  _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
inline __m128d in2_add(__m128d x, __m128d y)
{
        //if(_MM_GET_ROUNDING_MODE() == 0) cout << endl <<"false rm for ia" <<endl; 
        /* Component-wise addition */
        return  x + y;
}

inline __m128d in2_add(__m128d x, double y)
{
        //if(_MM_GET_ROUNDING_MODE() == 0) cout << endl <<"false rm for ia" <<endl; 
        /* Component-wise addition */
	return x + (__m128d){-y,y};
}

inline __m128d in2_add(double x, __m128d y)
{
        //if(_MM_GET_ROUNDING_MODE() == 0) cout << endl <<"false rm for ia" <<endl; 
        /* Component-wise addition */
	return y + (__m128d){-x,x};
}

inline __m128d in2_add(double x, double y)
{
        //if(_MM_GET_ROUNDING_MODE() == 0) cout << endl <<"false rm for ia" <<endl; 
        /* Component-wise addition */
	double aux = x+y;
        return  (__m128d){-aux,aux};
}

// Function to print a __m128d
inline void print(__m128d x){
    cout << "[" << in2_min(x) << "," << in2_max(x) << "]" << ", length: " << in2_max(x) - in2_min(x);
}



// Function to compare two intervals. Return True if a is inferior or equal to b, False if a is strictly superior to b, Undefined if the two intervals intersect in more than one point
inline boolean in2_le(__m128d a, __m128d b){
    if(a[1] <= -b[0]) return True;
    // Should it be equal or inferior or equal there ?
    else if (-a[0] > b[1]) return False;
    else return Undefined;
}

inline boolean in2_le(double a, __m128d b){
    if(a <= -b[0]) return True;
    // Should it be equal or inferior or equal there ?
    else if (a > b[1]) return False;
    else return Undefined;
}

inline boolean in2_ge(__m128d a, __m128d b){
	return in2_le(b,a);
}

inline boolean in2_ge(double b, __m128d a){
	return in2_le(in2_create(b),a);
}

inline boolean in2_ge(__m128d b, double a){
	return in2_le(b,in2_create(a));
}

inline boolean in2_lt(__m128d a, __m128d b){
    if(a[1] < -b[0]) return True;
    // Should it be equal or inferior or equal there ?
    else if (-a[0] >= b[1]) return False;
    else return Undefined;
}

inline boolean in2_lt(double a, __m128d b){
    if(a < -b[0]) return True;
    // Should it be equal or inferior or equal there ?
    else if (a >= b[1]) return False;
    else return Undefined;
}

inline boolean in2_lt(__m128d a, double b){
    if(a[1] < b) return True;
    // Should it be equal or inferior or equal there ?
    else if (-a[0] >= b) return False;
    else return Undefined;
}

inline boolean in2_lt(double a, double b){
    if(a < b) return True;
    // Should it be equal or inferior or equal there ?
    else if (a >= b) return False;
}

inline boolean in2_gt(__m128d a, __m128d b){
	return in2_lt(b,a);
}

// Function to compare an interval and a double. Return True if a is inferior or equal to b, False if a is strictly superior to b, undefined is a belongs to the interval b
inline boolean inferior_double(double a, __m128d b);

inline boolean inferior_double(double a, __m128d b){
    if(a <= -b[0]) return True;
    else if( a > b[1]) return False;
    else return Undefined;
}

// Function to merge two intervals in case of undefined comparison
inline __m128d in2_merge(__m128d a,__m128d b){
    return (__m128d){max(a[0],b[0]),max(a[1],b[1])};
}

inline __m128d in2_merge(double a,__m128d b){
    return (__m128d){max(a,b[0]),max(a,b[1])};
}

inline __m128i in2_merge(__m128i a,__m128i b){
    return (__m128i){max(a[0],b[0]),max(a[1],b[1])};
}
#endif
